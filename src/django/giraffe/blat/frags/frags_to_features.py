##############################################################################:
## frags_to_features.py
##
## frags_to_features --- Convert a list of feature fragments found on a plasmid
##                       into a list of features by joining fragments of the
##                       same feature into "trains" that later become features
##
## MYW 1/12/2011

##############################################################################
## Imports
from __future__ import division
from copy import deepcopy
from operator import attrgetter

from giraffe.blat.models import Feature_Type
from giraffe.blat.models import Feature_DB_Index
from giraffe.blat.models import Sequence_Feature

##############################################################################
## Global Constants
PCT_IDENTITY_ERROR_THRESHOLD = 0.25
WT_THRESHOLD = 0.05

##############################################################################
## Flags
_debug = False
_debug_features_to_observe = (98, 132)

##############################################################################
## Classes
class FeatureData(object):
    """ Caches feature data to prevent db access. """
    def __init__(self, fdb):
        self.__fdb = fdb
        self.__feature = fdb.feature
        self.__type = self.__feature.type_id
        self.__name = self.__feature.name
        self.__length = len(self.__feature.sequence)
        self.__clockwise = not fdb.antisense

    ## Accessors
    @property
    def type(self):
        return self.__type
    @property
    def length(self):
        return self.__length
    @property
    def clockwise(self):
        return self.__clockwise
    @property
    def name(self):
        return self.__name


class Frag(object):
    """
    Fragment of a feature recognized in a sequence.

    NOTE: Interpreting a fragment.
    
    Example: 56        3         1234           3
             feat_idx  frag_idx  seq_start_pos shift
    This is a fragment of a feature with feature index 56 (feature_index)
    It is the third fragment of that feature (fragment_index)  
    It was found at position 1234 in the sequence (seq_start_position)
    It has a shift of 3 (shift)
    """

    ## Constants
    SIZE = 12 # The size of an individual fragment, in bases

    def __init__(self, string = None, feature_index = -1, fragment_index = -1, 
            seq_start_position = -1, shift = 0):
        """ Create a fragment, either from a string or from the data itself. """
        if string == None:
            # Create from given arguments
            self.__feature_index  = feature_index
            self.__fragment_index = fragment_index
            self.__seq_start_position = seq_start_position
            self.__shift = shift
        else:
            # Parse a line of text (faster witout regexp)
            # to make the frag object
            frag_list = string.split()
            if len(frag_list) != 4: raise ValueError
            (self.__feature_index, self.__fragment_index, 
             self.__seq_start_position, self.__shift) = \
                  [ int(num) for num in frag_list ]

        # Incorporate the shift into the position at the beginning
        if self.__shift > 0:
            self.__seq_start_position += self.__shift

    ## Read-only basic fragment data
    @property
    def feature_index(self):
        return self.__feature_index
    @property
    def fragment_index(self):
        return self.__fragment_index
    @property
    def seq_start_position(self):
        return self.__seq_start_position
    @property
    def shift(self):
        return self.__shift

    ## Overloaded operators
    def __repr__(self):
        return 'Frag("%d %d %d %d")' % (self.feature_index,
                                        self.fragment_index,
                                        self.seq_start_position,
                                        self.shift)
    def __str__(self):
        return "[%d %d %d %d]" % (self.feature_index,
                                  self.fragment_index,
                                  self.seq_start_position,
                                  self.shift)


class FragTrain(object):
    """
    A train of fragments of a single feature to consider for upgrade to 
    "Feature" status.

    Contains of a list of Frags, counts of the hits, mutations, intserts, and
    deletes encountered by those frags, and methods to manipulate the data. Also
    contains information about the feature (such as its length).
    """

    ## Constants
    HIGH_FIDELITY_CUTOFF = 0.2 # Need at least 20% of the feature length
                               # matched to be high-fidelity
    MAX_INSERT_FRACTION = 0.75


    def __init__(self, feature_data, 
            frags = None, hits = 0, short = False, mutations = 0,
            inserts = 0, deletes = 0):
        """ Create a train with the supplied data, and add any frags given. """

        # Modifiable counts
        self.hits = hits
        self.short = short
        self.mutations = mutations
        self.inserts = inserts
        self.deletes = deletes

        # Read-only feature info
        self.__feature_data = feature_data

        # Keep the internals hidden
        self.__train = []
        if frags != None:
            self.__make_train(frags)

    ## Accessors
    # For read-only feature info
    @property
    def feature_index(self):
        return self.head.feature_index

    @property
    def feature(self):
        return self.__feature_data

    # For read-only train access
    @property
    def head(self):
        return self.__train[0]

    @property
    def tail(self):
        return self.__train[-1]

    ## Calculated properties
    @property
    def left_position(self):
        """
        The earliest possible place in the sequence that this feature
        starts, going clockwise.

        If it's clockwise, it's just the start. If not, it's the 
        startpos - length =  2*start - end
        """
        try:
            return self.__left_position
        except AttributeError:
            self.__left_position = self.start_position if self.feature.clockwise \
                                   else 2 * self.start_position - self.stop_position
            return self.__left_position

    @property
    def start_position(self):
        """ The starting position of the first fragment, regardless of direction. """
        return self.head.seq_start_position

    @property
    def stop_position(self):
        """ The stopping position of the last fragment, regardless of direction. """
        try:
            return self.__stop_position
        except AttributeError:
            self.__stop_position = self.tail.seq_start_position
            if (self.tail.fragment_index == \
                int((self.feature.length - 1)/Frag.SIZE)) and ( # If the last fragment
                self.feature.length % Frag.SIZE != 0):          # is not an even divisor
                                                                # of Frag.SIZE,
                                                                # add the remainder
                self.__stop_position += (self.feature.length % Frag.SIZE) - 1
            else:
                self.__stop_position +=  Frag.SIZE - 1          # Otherwise, add in 
                                                                # a whole Frag.SIZE
            return self.__stop_position

    @property
    def is_high_fidelity(self):
        """
        Determines if the FragTrain is high fidelity or not.

        A high fidelity FragTrain is a a FragTrain with perfect
        identity for more than 20% of the total size of the gene.
        """
        try:
            return self.__is_high_fidelity
        except AttributeError:
            self.__is_high_fidelity = \
               (self.inserts == 0 and self.deletes == 0 and 
                self.hits >= FragTrain.HIGH_FIDELITY_CUTOFF *
                self.feature.length)
            return self.__is_high_fidelity

    ## Action methods
    def extend(self, frag, add_hits = True):
        """ Adds one fragment to the end of the train. """
        self.__train.append(frag)

        if add_hits:
            self.hits += self.__calculate_hit_score(frag)

    def __deepcopy__(self, memo):
        # Don't want true deep copy: __feature_data stays shallow, as
        # do the Frags in the __train itself
        new_train = FragTrain(self.__feature_data,  hits=self.hits,
                              short=self.short,     mutations=self.mutations,
                              inserts=self.inserts, deletes=self.deletes,
                              frags=self.__train)

        return new_train

    ## Comparison methods
    def overlaps_with_preceeding(self, preceeding_train):
        """ 
        Checks for overlap with a train of the same feature that starts in the 
        sequence before this one.
        
        Checks for at least one fragment's worth of overlap.
        """
        return (preceeding_train.start_position <= self.start_position) and \
               (preceeding_train.stop_position  >= self.start_position + Frag.SIZE - 1)

    def overlaps_with_in_sequence(self, other_train):
        """ 
        Checks for overlap with another train in the final sequence.
        
        Checks for overlap with another train, but by using the left-most
        position, allows for the features to be different, and takes into 
        account whether the features are clockwise or not. Any overlap is
        enough.
        """
        return (self.left_position >= other_train.left_position and self.left_position <= other_train.stop_position) or \
                (other_train.left_position >= self.left_position and other_train.left_position <= self.stop_position)

    def has_similar_name(self, other_train):
        """ 
        Checks to see if two trains' names are similar.
        
        Checks for identical type and whether or not one name contains the other.
        """
        return self.feature.type == other_train.feature.type and \
                (self.feature.name.find(other_train.feature.name) >= 0 or
                 self.feature.name.find(other_train.feature.name) >= 0)

    ## Conversion methods
    def to_sequence_feature(self, seq_length):
        """
        Convert the train to a sequence feature (without a sequence).
        
        The feature is not saved.
        """
        seq_feat = Sequence_Feature()
        seq_feat.start = self.start_position
        seq_feat.end = self.stop_position
        seq_feat.clockwise = self.feature.clockwise

        # Features that cross the boundary should wrap their starting and ending
        # positions
        if seq_feat.start < 0:
            seq_feat.start += seq_length
        if seq_feat.end < 0:
            seq_feat.end += seq_length

        if seq_feat.start > seq_length:
            seq_feat.start %= seq_length

        if seq_feat.end > seq_length:
            seq_feat.end %= seq_length

        # Feature object itself is private, to prevent unnecessary db access
        seq_feat.feature_db_index = self.feature._FeatureData__fdb

        if globals()["_debug"]: 
            if self.feature_index in globals()["_debug_features_to_observe"]:
                print "local index %d is %s" % \
                      (self.feature_index, self.feature.name)

        # Variant labels for genes
        #  XXX FeatureType-Dependent Code
        if self.feature.type == Feature_Type.GENE:
            # A "high-fidelity" gene with a high deletion rate must just
            # be a fragment of a gene
            if not self.matches() and self.is_high_fidelity:
                s_start = self.head.fragment_index * Frag.SIZE
                s_end = s_start + self.stop_position - self.start_position
                seq_feat.subset_start = s_start
                seq_feat.subset_end = s_end
            elif self.score > WT_THRESHOLD or self.deletes > Frag.SIZE:
                seq_feat.is_variant = True
            elif self.inserts > 2 * Frag.SIZE:
                seq_feat.has_gaps = True
        #  XXX END FeatureType-Dependent Code
        return seq_feat

    ## Overloaded operators
    def __repr__(self):
        return '[' + ' '.join([repr(frag) for frag in self.__train]) + ']'

    def __str__(self):
        return '[' + ' '.join([str(frag) for frag in self.__train]) + ']'

    def __len__(self):
        return len(self.__train)
            
    ## Private utility methods
    def __make_train(self, frags):
        add_hits = (self.hits == 0) # only add in the hits if we have no
                                    # hit score yet
        for frag in frags:
            self.extend(frag, add_hits)

    def __calculate_hit_score(self, frag):
        """
        Calculates the hit score of a fragment of this feature type.
        
        The hit score is how many bases in that particualr fragment were
        actually matched
        """
        hit_score = self.feature.length - frag.fragment_index * Frag.SIZE
        if hit_score > Frag.SIZE: hit_score = Frag.SIZE
        return hit_score
            
    ##############################################################################
    ## Match Evaluation
    #  This method should be the only explicitly feature type-dependent
    #  code in this module
    #  TODO: Make this error/threshold function part of the FeatureType class
    #  XXX FeatureType-Dependent Code
    def matches(self):
        """ 
        Use a heuristic weighted sum to approximate how much of the feature's
        sequence is mutated.

        If the feature is exact, determine if the match is perfect and return either
        a perfect score or a perfect fail.
        """

        # Exact features must be exact: no scoring system
        if (self.feature.type == Feature_Type.EXACT_FEATURE) or \
           (self.feature.type == Feature_Type.ENZYME):
            if self.hits == self.feature.length and \
               self.inserts == 0 and self.deletes == 0 and self.mutations == 0:
                pct_identity_error = 0.0 # Perfect score
            else:
                pct_identity_error = 1.0 # Perfect fail
        else:
            # All other features: scoring system
            # scoring function:
            #   +1   for each exact match
            #   +0   for missing nucleotides at the start or end of sequence
            #   +0.3 for each nucleotide in a mutated fragment
            #        (i.e. assume roughly 70% of the fragment are mutated)
            #   -0.1 for each deleted nucleotide
            #   +0   for each inserted nucleotide in a gene and if total inserts
            #        is less than a threshold that's related to gene size
            #        (i.e. we don't penalize small inserts for genes)
            #   -0.1 for each inserted nucleotide not in a gene

            factor_match        =  1.0
           #factor_missing      =  0.0
            factor_mutations    =  0.3
            factor_deletes      = -0.1
           #factor_inserts_gene =  0.0
            factor_inserts      = -0.1

            net_matches = self.hits * factor_match

            # Missing nucleotides
            # net_matches += (self.feature.length - (net_matches + self.mutations)) * \
            #                 factor_missing

            # Mutations
            net_matches += self.mutations * factor_mutations

            # Deletions
            net_matches += self.deletes * factor_deletes

            # Insertions
            if self.feature.type != Feature_Type.GENE or \
               self.inserts > self.feature.length * FragTrain.MAX_INSERT_FRACTION:
                net_matches += self.inserts * factor_inserts
            else:
                pass
                # net_matches += self.inserts * factor_inserts_gene

            # Normalization
            pct_identity_error = 1 - net_matches / self.feature.length

        # Store the most recently calculated error score in the train
        self.score = pct_identity_error
        return pct_identity_error < PCT_IDENTITY_ERROR_THRESHOLD
    #  XXX End FeatureType-Dependent Code


class FragTrainLink(object):
    """ Links a Frag to a FragTrain, temporarily or permanently. """

    def __init__(self, frag, train):
        self.frag = frag
        self.train = train
        self.__calculate_diffs()

    ## Private utility methods
    def __calculate_diffs(self):
        """
        Calculate the sequence position and fragment differences between the
        fragment and the last fragment of the train.
        """
        self.frag_index_diff = self.frag.fragment_index - self.train.tail.fragment_index
        self.seq_pos_diff = self.frag.seq_start_position - (self.train.tail.seq_start_position + Frag.SIZE)

    ## Calculated properties
    @property
    def insert_size(self):
        """
        The number of bases that must have been inserted.
        
        The delete size is simply the negative of this number.
        """
        try:
            return self.__insert_size
        except AttributeError:
            self.__insert_size = self.seq_pos_diff - (self.frag_index_diff - 1) * Frag.SIZE 
            return self.__insert_size

    ## Relationship detection methods
    #  Detect whether or not a train-fragment link is consecutive, or if it has
    #  some mutations/inserts/deletes by examining the fragment and position
    #  differences
    def is_nonoverlapping(self):
        """ Basic sanity check. Should not have overlapping fragments. """
        return self.frag_index_diff > 0 and self.seq_pos_diff >= 0

    def is_consecutive(self):
        if globals()["_debug"]: 
            if self.frag.feature_index in globals()["_debug_features_to_observe"]:
                print "fdiff: %d, pdiff: %d" % (self.frag_index_diff,
                        self.seq_pos_diff)
            return self.frag_index_diff == 1 and self.seq_pos_diff == 0

    def has_mutation(self):
        """ Only called if not already consecutive. """
        return self.insert_size == 0

    def has_insert(self):
        """
        Heuristically determines if the observed gap is due to an insert.

        An insert is only valid if the disturbance it makes is sufficiently 
        small: if the disturbance (measured by seq_pos_diff) is too big, 
        it probably means there are two separate features, or the feature
        is really split
        """
        return self.insert_size >= 0 and self.seq_pos_diff < \
               int(self.train.feature.length * FragTrain.MAX_INSERT_FRACTION)
    
    ## Link solidification
    #  Extends the train and updates its mutation/deletion/insert counts
    def solidify(self):
        """ Just extend the train with the fragment. """
        self.train.extend(self.frag)

    def update_mutations(self):
        """
        Use the sequence position difference to estimate the number of
        mutations, which is actually just the number of nucleotides in the
        missing fragments.
        """
        self.train.mutations += self.seq_pos_diff

    def update_inserts(self):
        """
        Update the insert count heuristically.

        An insert should really only cause one fragment to be
        missing, but there may be some mutations around the
        splicing sites so we conservatively assume two
        fragments may be missing due to the insert. all other
        missing fragments must have mutations.
        """

        if self.frag_index_diff > 3:
            # (frag_index_diff-1)-2 is the number of fragments that may
            # have mutations, assuming we allow two fragments to
            # be mutated due to the insert
            self.train.mutations += (self.frag_index_diff - 3) * Frag.SIZE
        self.train.inserts += self.insert_size;

    def update_deletes(self):
        """
        Update the train's delete count with the insert_size.

        This is because "delete size" is just insert size when
        insert size is < 0: both measure how far off the actual position
        of the fragment is from its expected position.
        """
        self.train.deletes += abs(self.insert_size)

    def make_hypo_train(self):
        """
        Return a hypothetical empty train with only the counts: useful in
        calculating whether or not the missing fragments are due to a deletion
        or not.
        """
        hypo_train = FragTrain(self.train.feature,
                mutations = self.train.mutations,
                inserts = self.train.inserts)
        hypo_train.deletes = self.train.deletes + abs(self.insert_size)

        # Hypothetical hits is current hits plus all remaining
        # fragments
        hypo_train.hits = Frag.SIZE *(len(self.train) + 1) + \
                          self.train.feature.length - \
                          self.frag.fragment_index * Frag.SIZE

        return hypo_train



##############################################################################
## Global Functions
def _group_frags_by_feature_index(frags_strings):
    """ Parse the frag strings and group the fragments. """
    frags_by_feature = {}

    for fx, feature_frag_line in enumerate(frags_strings):
        # Make a feature from the line of text
        try:
            frag = Frag(string = feature_frag_line)
        except ValueError:
            continue

        # Hash the frags together by feature index
        if frag.feature_index in frags_by_feature:
            frags_by_feature[frag.feature_index].append(frag)
        else:
            frags_by_feature[frag.feature_index] = [frag]

    return frags_by_feature

def _frags_to_trains(frags, feature_data, seq_length):
    """
    Converts a list of fragments (of a single feature index) into a list of
    trains of fragments, which are candidates for feature identification.

    Arguments:
        frags         --- the list of fragments
    """
    if _debug: 
        if frags[0].feature_index in _debug_features_to_observe:
            print "--- Frags (%d)" % len(frags)
            for frag in frags:
                print frag

    trains = []

    # Main loop: iterate through the list of fragments
    for fx, frag in enumerate(frags):
       
        # Default to making a new train from the fragment, and if any links
        # are solidified (i.e., the fragment extends some other train),
        # make sure that the new train is marked as short.
        create_new_train = True
        new_train_is_short = False

        for train in trains:
            link = FragTrainLink(frag, train)

            # First: try to extend consecutive trains
            if link.is_consecutive():
                # Solidify the link by extending the train
                link.solidify()
                create_new_train = False

                if _debug: 
                    if train.feature_index in _debug_features_to_observe:
                        print "extended consecutive train"
                        print train


            # Second: try to extend trains with inserts/mutations/deletions
            elif not train.short and link.is_nonoverlapping():

                if link.has_mutation():
                    # Solidify the link by extending the train
                    # and updating the mutation count
                    link.solidify()
                    link.update_mutations()
                    create_new_train = False

                    if _debug: 
                        if train.feature_index in _debug_features_to_observe:
                            print "extended train with mutations"
                            print train


                # In case the sequence really represents two separate, 
                # features which we have misread as an insertion/deletion,
                # we hedge our bets. Whenever we detect an insertion or
                # a deletion in a train, we add a copy of that train
                # to the list of trains, and only extend one of the trains.
                elif link.has_insert():
                    if train.matches():
                        new_train = deepcopy(train)
                        new_train.short = True
                        trains.append(new_train)

                        if _debug: 
                            if train.feature_index in _debug_features_to_observe:
                                print "made dupliate train due to inserts"
                                print train


                    # Solidify the link by extending the train
                    # and updating the insert and/or mutation count
                    link.solidify()
                    link.update_inserts()
                    new_train_is_short = True
                    
                    if _debug: 
                        if train.feature_index in _debug_features_to_observe:
                            print "extended train with inserts"
                            print train


                # Link has deletion
                else:
                    # For deletions, use a "hypothetical" train, with the
                    # maximum possible number of hits
                    hypo_train = link.make_hypo_train()
                    # XXX I don't understand why this has to be the case

                    if hypo_train.matches():
                        # Append an identical copy of the train in its current
                        # state to the beginning of the train list.  We need to
                        # insert the new, identical copy before the iterator, so
                        # that it doesn't get iterated over again, or else an
                        # infinite loop results in certain cases
                        new_train = deepcopy(train)
                        trains.insert(fx, new_train)

                        # Solidify the link by extending the train
                        # and updating the delete count
                        link.solidify()
                        link.update_deletes()

                        if _debug: 
                            if train.feature_index in _debug_features_to_observe:
                                print "extended train with deletes"
                                print train


                    # The order of trains should therefore have the same
                    # properites as the order of trains in the perl code.

        # Make a new train, if necessary
        if create_new_train: 
            if frag.seq_start_position <= seq_length:
                if _debug: 
                    if frag.feature_index in _debug_features_to_observe:
                        if new_train_is_short:
                            print "new short train"
                        else:
                            print "new train"
                        print frag

                trains.append(FragTrain(feature_data, frags = [frag], \
                                        short = new_train_is_short))

    if _debug: 
        notprint = True
        for train in trains:
            if train.feature_index in _debug_features_to_observe:
                if notprint:
                    print "--- Trains (%d)" % len(trains)
                    notprint = False
                print train

    return trains

def _pick_good_trains(trains):
    """
    Pick only trains of sufficiently high quality that don't overlap.

    Quality is set by a call-back function that includes feature type-dependent
    code.
    """
    good_trains = []

    last_train = None

    # Remove overlapping preceding trains
    for train in trains:
        if last_train != None and train.overlaps_with_preceeding(last_train):
            continue

        # XXX Feature type-dependent name/variant modifications done in 
        #     check_features in old code would go here. Can be moved off
        #     to the JavaScript code
        # XXX FeatureType-Dependent Code
        if train.matches() or (train.feature.type == Feature_Type.GENE and
                train.is_high_fidelity):
            good_trains.append(train)
            last_train = train

    return good_trains

def _trains_to_features(trains, seq_length):
    """
    Prune interleaving trains and make features from the rest.
    """
    features = []

    trains.sort(key = attrgetter('left_position'))

    for outer_train in trains:
        add_feature = True

        # XXX FeatureType-Dependent Code
        #     Non-Enzyme features must be pruned.
        if outer_train.feature.type != Feature_Type.ENZYME:
            for inner_train in trains:
                if inner_train.left_position > outer_train.stop_position:
                    break

                if (inner_train.feature.type == Feature_Type.GENE or
                    inner_train.has_similar_name(outer_train)) and \
                    inner_train.overlaps_with_in_sequence(outer_train):

                    # If the conditions are met, check the scores, and if the 
                    # inner feature scores better, keep it instead of the outer
                    # one

                    if inner_train.score < outer_train.score:
                        if _debug: print "Pruning %s (%d - %d)" % \
                            (outer_train.feature.name, 
                             outer_train.start_position,
                             outer_train.stop_position)

                        add_feature = False
                        break
        # else: Enzymes always get added.
        # XXX End FeatureType-Dependent Code
           
        # If we found a better feature in the inner loop, don't add this outer
        # feature. When the outer loop gets to the better feature on its own, 
        # it will add it
        if add_feature:
            features.append(outer_train.to_sequence_feature(seq_length))

    return features


## Main entry point
def frags_to_features(frags_strings, db, seq_length):
    """
    Given a list of feature fragments, turns them into Sequence Features.
    """

    all_trains = [] 
    
    # Sort the fragments into groups of the same feature index
    frags_by_feature_index = _group_frags_by_feature_index(frags_strings)

    # Iterate over each group
    for (feature_index, single_index_frags) in frags_by_feature_index.items():

        # Pull the feature data from the DB
        feature_db = Feature_DB_Index.objects.select_related(
            'feature', 'feature__type'
        ).get(
            feature_index=feature_index,db=db
        )

        # Make the biggest possible trains you can
        single_index_trains  = _frags_to_trains(single_index_frags,
                                                FeatureData(feature_db),
                                                seq_length)

        # Prune them to keep only a relevant subset of sufficiently high
        # quality
        single_index_good_trains = _pick_good_trains(single_index_trains)

        # Add them to the global list
        all_trains.extend(single_index_good_trains)

    # Turn the global list into a set of non-overlapping SequenceFeature objects
    return _trains_to_features(all_trains, seq_length)

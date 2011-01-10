from __future__ import division
from copy import deepcopy
from giraffe.blat.models import Feature_DB_Index
from giraffe.blat.models import Sequence_Feature

##############################################################################
## Global Constants
PCT_IDENTITY_ERROR_THRESHOLD = 0.25
WT_THRESHOLD = 0.05

##############################################################################
## Flags
__debug = True


##############################################################################
## Classes

class AutoCalc(object):
    """
    Base class for FragTrains and FragTrainLinks, which calculate data members 
    on demand.
    """

    def __getattr__(self, name):
        """
        Automatically calculates sequence distances and overlaps on demand.

        Sees if there's a __calculate_<name> method for it. If so, calculates
        it and stores the result. If not, raises an exception. 
        """
        # XXX hasattr/getattr, as external functions, need to work with the
        # mangled name. Kind of silly, but so be it.
        calc_method_name = "_%s__calculate_%s" % (self.__class__.__name__, name)
        if hasattr(self, calc_method_name):
            setattr(self, name, getattr(self, calc_method_name)())
        else:
            raise AttributeError("Cannot calculate " + name)

        return getattr(self, name)

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

    # Read-only basic fragment data
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


class FragTrain(AutoCalc):
    """
    A train of fragments to consider for upgrade to "Feature" status.

    Contains of a list of Frags, and methods to manipulate it.
    """

    ## Constants
    HIGH_FIDELITY_CUTOFF = 0.2 # Need at least 20% of the feature length
                               # matched to be high-fidelity

    def __init__(self, feature_db, 
            frags = None, hits = 0, short = False, mutations = 0,
            inserts = 0, deletes = 0):

        # Modifiable counts
        self.hits = hits
        self.short = short
        self.mutations = mutations
        self.inserts = inserts
        self.deletes = deletes

        # Read-only feature info
        self.__load_feature_info(feature_db)

        # Keep the internals hidden
        self.__train = []
        if frags != None:
            self.__make_train(frags)

    def __load_feature_info(self, feature_db):
        # feature:        the feature object associated with the fragment
        # feature_length: the length of that feature, in bases
        # clockwise:      which way the feature has been found
        self.__feature = feature_db.feature
        self.__feature_length = len(self.__feature.sequence)
        self.__clockwise = not feature_db.antisense


    # Accessors for read-only feature info
    @property
    def feature(self):
        return self.__feature

    @property
    def feature_length(self):
        return self.__feature_length

    @property
    def clockwise(self):
        return self.__clockwise

    @property
    def start_position(self):
        return self.head().seq_start_position

    def __make_train(self, frags):
        for frag in frags:
            self.extend(frag)

    def __repr__(self):
        return '[' + ' '.join([repr(frag) for frag in self.__train]) + ']'

    def __str__(self):
        return '[' + ' '.join([str(frag) for frag in self.__train]) + ']'

    def __len__(self):
        return len(self.__train)
            
    # Not called by AutoCalc:
    def __frag_calculate_hit_score(self, frag):
        # hit_score: how many bases in that particualr fragment 
        #            were actually matched
        hit_score = self.feature_length - frag.fragment_index * Frag.SIZE
        if hit_score > Frag.SIZE: hit_score = Frag.SIZE
        return hit_score
            
    # Called on demand by AutoCalc: only one value needed per lifetime
    def __calculate_stop_position(self):
        stop_position = self.tail().seq_start_position
        if (self.tail().fragment_index == \
            int((self.feature_length - 1)/Frag.SIZE)) and ( # If the last fragment
            self.feature_length % Frag.SIZE != 0):          # is not an even divisor
                                                            # of Frag.SIZE,
                                                            # add the remainder
            stop_position += (self.feature_length % Frag.SIZE) - 1
        else:
            stop_position +=  Frag.SIZE - 1                 # Otherwise, add in 
                                                            # a whole Frag.SIZE
        return stop_position

    # A high fidelity FragTrain is a a FragTrain with perfect
    # identity for more than 20% of the total size of the gene
    # Called on demand by AutoCalc: only one value needed per lifetime
    def __calculate_is_high_fidelity(self):
        return (self.inserts == 0 and self.deletes == 0 and 
                self.hits >= FragTrain.HIGH_FIDELITY_CUTOFF *
                self.feature_length)

    def extend(self, frag):
        self.__train.append(frag)
        self.hits += self.__frag_calculate_hit_score(frag)

    def head(self):
        return self.__train[0]

    def tail(self):
        return self.__train[-1]

    def overlaps_with(self, other_train):
        return (other_train.head().seq_start_position <= \
                self.head().seq_start_position) and \
               (other_train.stop_position >= self.tail().seq_start_position + \
                       Frag.SIZE - 1)

    def to_feature(self):
        feat = Sequence_Feature()
        feat.feature = self.feature
        feat.start = self.start_position
        feat.end = self.stop_position
        feat.clockwise = self.clockwise
        return feat

class FragTrainLink(AutoCalc):
    """ Links a Frag to a FragTrain, temporarily or permanently. """

    ## Constants
    MAX_INSERT_FRACTION = 0.75

    def __init__(self, frag, train):
        self.frag = frag
        self.train = train


    """
    Calculate the sequence position and fragment differences between the 
    fragment and the last fragment of the train
    """
    def __calculate_frag_index_diff(self):
        self.frag_index_diff = self.frag.fragment_index - self.train.tail().fragment_index
    def __calculate_seq_pos_diff(self):
        self.seq_pos_diff = self.frag.seq_start_position - \
            (self.train.tail().seq_start_position + Frag.SIZE)

    def __calculate_insert_size(self):
        return self.seq_pos_diff - (self.frag_index_diff - 1) * Frag.SIZE 

    def __getattr__(self, name):
        """
        Automatically calculates sequence distances and overlaps on demand.

        Sees if there's a __calculate_<name> method for it. If so, calculates
        it and stores the result. If not, raises an exception. 
        """

        # XXX hasattr/getattr, as external functions, need to work with the
        # mangled name. Kind of silly, but so be it.
        calc_method_name = "_FragTrainLink__calculate_" + name
        if hasattr(self, calc_method_name):
            setattr(self, name, getattr(self, calc_method_name)())
        else:
            raise AttributeError(
                    "Cannot find %s and no method %s to calculate it." % 
                    (name, calc_method_name))

        return getattr(self, name)

    ## Detect consecutive fragments, inserts/deletions, or mutations
    def is_nonoverlapping(self):
        return self.frag_index_diff > 0 and self.seq_pos_diff >= 0

    def is_consecutive(self):
        return self.frag_index_diff == 1 and self.seq_pos_diff == 0

    def has_mutation(self):
        return self.insert_size == 0

    def has_insert(self):
        # An insert is only valid if the disturbance it makes is sufficiently 
        # small: if the disturbance (measured by seq_pos_diff) is too big, 
        # it probably means there are two separate features, or the feature
        # is really split
        return self.insert_size >= 0 and self.seq_pos_diff < \
               int(self.train.feature_length * FragTrainLink.MAX_INSERT_FRACTION)
    
    ## Solidify the link by extending the train
    def solidify(self):
        self.train.extend(self.frag)

    def update_mutations(self):
        self.train.mutations += self.seq_pos_diff

    def update_inserts(self):
        # an insert should really only cause one fragment to be
        # missing, but there may be some mutations around the
        # splicing sites so we conservatively assume two
        # fragments may be missing due to the insert. all other
        # missing fragments must have mutations.

        if self.frag_index_diff > 3:
            # (frag_index_diff-1)-2 is the number of fragments that may
            # have mutations, assuming we allow two fragments to
            # be mutated due to the insert
            self.train.mutations += (self.frag_index_diff - 3) * Frag.SIZE
        self.train.inserts += self.insert_size;

    def update_deletes(self):
        # This is because "delete size" is just insert size when
        # insert size is < 0: both measure how far off the actual position
        # of the fragment is from its expected position
        self.train.deletes += abs(self.insert_size)

    def make_hypo_train(self):
        hypo_train = deepcopy(self.train)
        hypo_train.deletes += abs(self.insert_size)

        # hypothetical hits is current hits plus all remaining
        # fragments
        hypo_train.hits = Frag.SIZE *(len(hypo_train) + 1) + \
                          self.train.feature_length - \
                          self.frag.fragment_index * Frag.SIZE

        return hypo_train

##############################################################################
## Global Functions

def _group_frags_by_feature_index(frags_strings):
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

# This callback function should be the only explicitly feature type-dependent
# code
# TODO: Make this threshold function part of the FeatureType class
def _percent_identity_threshold(train):

    # Exact features must be exact: no scoring system
    if (train.feature.type.type == "ExactFeature") or \
       (train.feature.type.type == "Cutter"):
        if train.hits == train.feature_length and \
           train.inserts == 0 and train.deletes == 0 and train.mutations == 0:
            return True
        else:
            return False
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
        #
        # $i: inserts
        # $d: deletes
        # $m: mutations (actually, # nucleotides from missing fragments)

        factor_match        =  1.0
       #factor_missing      =  0.0
        factor_mutations    =  0.3
        factor_deletes      = -0.1
       #factor_inserts_gene =  0.0
        factor_inserts      = -0.1

        net_matches = train.hits * factor_match

        # Missing nucleotides
        # net_matches += (train.feature_length - (net_matches + train.mutations)) * \
        #                 factor_missing

        # Mutations
        net_matches += train.mutations * factor_mutations

        # Deletions
        net_matches += train.deletes * factor_deletes

        # Insertions
        if train.feature.type.type != "Gene" or \
           train.inserts > train.feature_length * FragTrainLink.MAX_INSERT_FRACTION:
            net_matches += train.inserts * factor_inserts
        else:
            pass
            # net_matches += train.inserts * factor_inserts_gene

        # Normalization
        pct_identity_error = 1 - net_matches / train.feature_length

        return pct_identity_error < PCT_IDENTITY_ERROR_THRESHOLD

def _frags_to_trains(frags, feature_db, 
                     train_matches = _percent_identity_threshold):
    """
    Converts a list of fragments (of a single feature index) into a list of
    trains of fragments, which are candidates for feature identification.

    Arguments:
        frags         --- the list of fragments
        train_matches --- a call-back function to use when evaluating
                          whether or not a train is of high enough quality.
                          must take a (train) as its argument
    """
    trains = []

    # Main loop: iterate through the list of fragments
    for fx, frag in enumerate(frags):
       
        # Default to making a new train from the fragment, and if any links
        # are solidified (i.e., the fragment extends some other train),
        # make sure that the new train is marked as short.
        create_new_train = True
        new_train_is_short = False

        dx = 0;
        for train in trains:
            link = FragTrainLink(frag, train)

            # First: try to extend consecutive trains
            if link.is_consecutive():
                # Solidify the link by extending the train
                link.solidify()
                create_new_train = False

            # Second: try to extend trains with inserts/mutations/deletions
            elif not train.short and link.is_nonoverlapping():

                if link.has_mutation():
                    # Solidify the link by extending the train
                    # and updating the mutation count
                    link.solidify()
                    link.update_mutations()
                    create_new_train = False

                # In case the sequence really represents two separate, 
                # features which we have misread as an insertion/deletion,
                # we hedge our bets. Whenever we detect an insertion or
                # a deletion in a train, we add a copy of that train
                # to the list of trains, and only extend one of the trains.
                elif link.has_insert():
                    if train_matches(train):
                        new_train = deepcopy(train)
                        new_train.short = True
                        trains.append(new_train)

                    # Solidify the link by extending the train
                    # and updating the insert and/or mutation count
                    link.solidify()
                    link.update_inserts()
                    new_train_is_short = True
                    
                # Link has deletion
                else:
                    # For deletions, use a "hypothetical" train, with the
                    # maximum possible number of hits
                    hypo_train = link.make_hypo_train()
                    # XXX I don't understand why this has to be the case

                    if train_matches(hypo_train):
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
                    # The order of trains should therefore have the same
                    # properites as the order of trains in the perl code.

        # Make a new train, if necessary
        if create_new_train: 
            trains.append(FragTrain(feature_db, frags = [frag], \
                                    short = new_train_is_short))

    if __debug: 
            print '[' + ", ".join([str(train) for train in trains]) + ']'

    return trains

def _trains_to_features(trains,
                        train_matches = _percent_identity_threshold):
    features = []

    last_train = None

    # Remove overlapping trains, adjust details of particualr train
    for train in trains:
        if last_train != None and train.overlaps_with(last_train):
            continue

        # name/variant modifications done in check_features in old code
        # TODO rewrite them in a logical way consistent with the new
        #      object-oriented model
        if train_matches(train) or train.is_high_fidelity:
            features.append(train.to_feature())
            last_train = train

    return features

def _pick_good_features(features):
    # TODO: make this function a method of Seq_Feat objects
    def left_pos(feat):
        return feat.start if feat.clockwise else 2 * feat.start - feat.end

    # The returns left-most position of sequence. if it's
    # clockwise, it's just the start. if not, it's the startpos - length =
    # 2*start - end
    features.sort(key = left_pos)

    for feature in features:
        left_pos_outer = left_pos(feature)

### MAIN ENTRY POINT
def frags_to_features(frags_strings, db):
    """
    Given a list of feature fragments, turns them into features.
    """

    if __debug: print "Turning frags into features"

    all_features = [] 

    # Sort the fragments into groups of the same feature
    frags_by_feature_index = _group_frags_by_feature_index(frags_strings)

    # Iterate over each group
    for (feature_index, single_index_frags) in frags_by_feature_index.items():
        feature_db = Feature_DB_Index.objects.get(feature_index=feature_index,
                                                  db=db)
        single_index_trains   = _frags_to_trains(single_index_frags,
                                                 feature_db)
        single_index_features = _trains_to_features(single_index_trains)
        all_features.extend(single_index_features)

    _pick_good_features(all_features)

    return all_features

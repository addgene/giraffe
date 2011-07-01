import sys
sys.path.append('../../')
import re
import unittest
import giraffe.simple_test as st
st.setup()
import giraffe.blat.models as models
import django

class TestFeature(object):
    property_map = { 
      "name":       "feature",
      "id":         "feature_id",
      "start":      "start",
      "end":        "end",
      "cut":        "cut",
      "clockwise":  "clockwise"
    }

    def __init__(self, test, name = None, id = None, 
            start = None, end = None, cut = None, clockwise = None, sequence = ''):
            self.__test      = test
            self.__name      = name
            self.__id        = id
            self.__start     = start
            self.__end       = end
            self.__cut       = cut
            self.__clockwise = clockwise
            self.__sequence  = sequence

    @property
    def sequence(self):
        return self.__sequence


    def assertEqual(self, db_seq_feature):
        for property_name in TestFeature.property_map:
            property = getattr(self, "_TestFeature__" + property_name)

            if property != None:
                try:
                    self.__test.assertEqual(
                        property,
                        db_seq_feature[TestFeature.property_map[property_name]]
                    )
                except Exception as e:
                    modargs = list(e.args)
                    modargs[0] = (e.args[0] + " (%s, %s) " %
                          (db_seq_feature[TestFeature.property_map['name']],
                           TestFeature.property_map[property_name]))
                    e.args = tuple(modargs)
                    raise e

    def assertFound(self, db_seq_feature_list):
        featureFound= False

        for feature in db_seq_feature_list:
            matches = True

            for property_name in TestFeature.property_map:
                property = getattr(self, "_TestFeature__" + property_name)
                if property != None:
                    if property != feature[TestFeature.property_map[property_name]]:
                        matches = False
                        break
            
            if matches:
                featureFound = True
                break

        self.__test.assertTrue(featureFound)


class ItDetectsFeaturesInDNASequences(unittest.TestCase):

    # Setup code
    def setUp(self):
        """Turns off debugging to prevent spurious output in tests"""
        django.conf.settings.DEBUG = False

    # Utility methods
    def find_features(self, sequence, db_name = 'default'):
        """Returns a dictionary of feature data, given a sequence"""

        # Call the feature detection algorithm on the sequence and get out a
        # hash
        hash = models.Giraffe_Mappable_Model.detect_features(sequence, db_name)

        # Use the hash to query the DB to pull out the sequence model object
        db = models.Feature_Database.objects.get(name=db_name)
        sequence_obj = models.Sequence.objects.get(db=db,hash=hash)

        # Iterate through that sequence object's sequence_feature set,
        # i.e. the set of all sequence_feature objects that have this sequence
        # as their .sequence property
        features = []
        for f in sequence_obj          \
                 .sequence_feature_set \
                 .select_related(
                     'feature_db_index__feature',
                     'feature_db_index__feature__type',
                 ):

            features.append(f.to_dict())

        for f in sequence_obj          \
                 .sequence_feature_annotated_set \
                 .select_related('feature_type'):
            features.append(f.to_dict())

        return features


    def test_ItDetectsALoneShortFeature(self):
        """Tests that a very short feature-only sequence detects only itself."""
        # Detect the feature
        lone_feat = TestFeature(self,
            name = 'EK', id = 15, start = 1, end = 15, clockwise = True,
            sequence = 'GATGACGACGACAAG'
        ) 
        features = self.find_features(lone_feat.sequence)
        self.assertEqual(len(features), 1)
        lone_feat.assertEqual(features[0])


    def test_ItDetectsFeaturesInLowerCase(self):
        """Tests that lowercase sequences still detect features"""
        # Detect a lowercase feature
        lone_feat = TestFeature(self,
            name = 'EK', id = 15, start = 1, end = 15, clockwise = True,
            sequence = 'gatgacgacgacaag'
        ) 
        features = self.find_features(lone_feat.sequence)
        self.assertEqual(len(features), 1)
        lone_feat.assertEqual(features[0])

    def test_ItDetectsNothingInAFeatureLessSequence(self):
        """Tests that sequences which clearly have no features in them do not
        somehow generate features"""

        for seq in ['A', 'ATGC', 'T' * 4096 ]:
            features = self.find_features(seq)
            self.assertEqual(len(features), 0)

    def test_ItIgnoresInvalidTextInSequences(self):
        """Tests that FASTA comments and invalid characters are cut"""

        seqs = [
            'GA5TGA1CGAC2392GA2CAA1G',
            ' {+G( A%T[G]A1;CG  > >AC"2<3@&92,G~`A2C.A?A/G}', 
            """>EK | feature id 15    
01 GATG 04
; This is a comment!
;; So is this!
05 ACGA 08
;; this is too
; is this?
09 CGAC 10
11 AAG  13"""
        ]

        feat = TestFeature(self,
            name = 'EK', id = 15, start = 1, end = 15, clockwise = True,
        ) 

        for seq in seqs:
            features = self.find_features(seq)
            self.assertEqual(len(features), 1)
            feat.assertEqual(features[0])
    
    def test_ItRefusesMultipleSequences(self):
        self.assertRaises(models.BadSequenceError, self.find_features, """
>EK|15
GATGACGACGACAAG
>EK|15
GATGACGACGACAAG""")


    def test_ItDetectsFeatureOnlySequences(self):
        """Tests that sequences which consist of only one feature
        return at least that one feature and nothing else with the
        same length. """

        # Pick some random features from the default db
        feature_ids_to_test = frozenset([
            39,  184, 153, 54, 96,  183, 31,  22,  80,  71,  108, 73, 114, 147, 10,
            138, 21,  85,  41, 180, 150, 160, 167, 133, 168, 131, 47, 24,  96,  60
        ])
        # 23 and 82: HIV conflict. avoid for now
        # 37 and 38: T7 leader conflict. avoid for now


        with open('fixtures/features/generic.features') as features_file:
            for (line_num, line) in enumerate(features_file):

                # XXX not reliable, but works or for now
                feature_id = line_num + 1

                if feature_id in feature_ids_to_test:
                    line_items = line.split()

                    feat = TestFeature(self,
                            name = line_items[0].split(':')[1],
                            sequence = line_items[1],
                            id = feature_id,
                            start = 1,
                            end = len(line_items[1]))

                    # Detect the features
                    features = self.find_features(feat.sequence)

                    # How many of them have the same length?
                    same_length = [f for f in features
                         if abs(f['end'] - f['start']) + 1 == len(feat.sequence) and
                            not re.match('ORF', f['feature'])]
                            

                    # Only one feature has the same length
                    try:
                        self.assertEqual(len(same_length), 1)
                    except:
                        print same_length
                        raise

                    # That feature should have same ID and name
                    feat.assertEqual(same_length[0])

    def test_ItDetectsFeatureGroupSequences(self):
        """Tests that sequences which consist of only a triplet of features
        return all of those features. """

        # Pick some random features from the default db
        feature_ids_to_test = [
            39, 184, 153,  54,  96, 183, 31, 22, 80,  71,
           108,  73, 114, 147,  10, 138, 21, 85, 41, 180,
           150, 160, 167, 133, 168, 131, 47, 24, 96,  60
        ]
        feature_ids_to_test_set = frozenset(feature_ids_to_test)

        features_to_test = {}


        # Load the feature sequences to test
        with open('fixtures/features/generic.features') as features_file:
            for (line_num, line) in enumerate(features_file):

                # XXX not reliable, but works or for now
                feature_id = line_num + 1

                if feature_id in feature_ids_to_test_set:
                    line_items = line.split()
                    features_to_test[feature_id] = \
                        TestFeature(self,
                            id = feature_id,
                            sequence = line_items[1],
                            name = line_items[0].split(':')[1]
                        )

        # Make the combined feature sequences
        for group in zip(feature_ids_to_test[0:10],
                         feature_ids_to_test[10:20], 
                         feature_ids_to_test[20:30]):

            # Create the long sequence and groups of features to detect
            sequence = ''
            for id in group:
                sequence += features_to_test[id].sequence

            # Detect the features
            features = self.find_features(sequence)
            
            # Do we detect at least as many features as we put in?
            self.assertTrue(len(features) >= len(group))

            # Are all of the features there?
            for id in group:
                features_to_test[id].assertFound(features)

    def test_DetectCrossBoundaryORFs(self):
        sequence = "GAGCAGGACTGGGCGGCGGCCAAAGCGGTCGGACAGTGCTCCGAGAACGGGTGCGCATAGAAATTGCATCAACGCATATAGCGCTAGCAGCACGCCATAGTGACTGGCGATGCTGTCGGAATGGACGATATCCCGCAAGAGGCCCGGCAGTACCGGCATAACCAAGCCTATGCCTACAGCATCCAGGGTGACGGTGCCGAGGATGACGATGAGCGCATTGTTAGATTTCATAAAAAAAAAAAAAATCAGGTCGAGGTGGCCCGGCTCCATGCACCGCGACGCAACGCGGGGAGGCAGACAAGGTATAGGGCGGCGCCTACAATCCATGCCAACCCGTTCCATGTGCTCGCCGAGGCGGCATAAATCGCCGTGACGATCAGCGGTCCAATGATCGAAGTTAGGCTGGTAAGAGCCGCGAGCGATCCTTGAAGCTGTCCCTGATGGTCGTCATCTACCTGCCTGGACAGCATGGCCTGCAACGCGGGCATCCCGATGCCGCCGGAAGCGAGAAGAATCATAATGGGGAAGGCCATCCAGCCTCGCGTCGCGAACGCCAGCAAGACGTAGCCCAGCGCGTCGGCCGCCATGCCGGCGATAATGGCCTGCTTCTCGCCGAAACGTTTGGTGGCGGGACCAGTGACGAAGGCTTGAGCGAGGGCGTGCAAGATTCCGAATACCGCAAGCGACAGGCCGATCATCGTCGCGCTCCAGCGAAAGCGGTCCTCGCCGAAAATGACCCAGAGCGCTGCCGGCACCTGTCCTACGAGTTGCATGATAAAGAAGACAGTCATAAGTGCGGCGACGATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATGGGTCGGCGCTCTCCCTTATGCGACTCCTGCATTAGGAAGCAGCCCAGTAGTAGGTTGAGGCCGTTGAGCACCGCCGCCGCAAGGAATGGTGCATGTAAGGAGATGGCGCCCAACAGTCCCCCGGCCACGGGGCCTGCCACCATACCCACGCCGAAACAAGCGCTCATGAGCCCGAAGTGGCGAGCCCGATCTTCCCCATCGGTGATGTCGGCGATATAGGCGCCAGCAACCGCACCTGTGGCGCCGGTGATGCCGGCCACGATGCGTCCGGCGTAGAGAATCCACAGGACGGGTGTGGTCGCCATGATCGCGTAGTCGATAGTGGCTCCAAGTAGCGAAGC"
        features = self.find_features(sequence)
        orfs = [f for f in features if re.match('ORF', f['feature'])]
        found_cross_boundary_orf = False
        for f in orfs:
            if f['start'] == 246 and f['end'] == 231 and f['clockwise'] == False:
                found_cross_boundary_orf = True
                break
        self.assertTrue(found_cross_boundary_orf)

 
    def test_ItDetectsLoneORFs(self):
        """
        This test tries to make sure we only get one ORF, the entire sequence,
        from an ORF sequence. Note that there's actually an ORF in the sequence
        reverse complement, but not long enough to qualify.
        """
        lone_orf = TestFeature(self,
            name = 'ORF frame 1', start = 1, end = 840, clockwise = True,
            sequence = "atggcagcgcgccgaccgcgatgggctgtggccaatagcggctgctcagcagggcgcgccgagagcagcggccgggaaggggcggtgcgggaggcggggtgtggggcggtagtgtgggccctgttcctgcccgcgcggtgttccgcattctgcaagcctccggagcgcacgtcggcagtcggctccctcgttgaccgaatcaccgacctctctccccagggggatccaccggagcttaccatgaccgagtacaagcccacggtgcgcctcgccacccgcgacgacgtccccagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgacatcggcaaggtgtgggtcgcggacgacggcgccgcggtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggcgtctcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaggcggccgagcgcgccggggtgcccgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgcctga"
        ) 
        features = self.find_features(lone_orf.sequence)
        orfs = [f for f in features if re.match('ORF', f['feature'])]
        # check we only found one ORF, smaller orfs from reverse complement did not qualify
        self.assertTrue(len(orfs) == 1)
        # and the ORF we found, is the entire sequence
        lone_orf.assertEqual(orfs[0])

        
    def test_ItDetectsFeaturesThatCrossTheBoundary(self):
        cross_features = [
            # Enzyme
            TestFeature(self,
                name = 'DraI', start = 28, end = 3, cut = 30, clockwise = True,
                sequence = "aaatgaccctttgggatgaaagggcccttt"
            ),
            # Enzyme with wrapped cut site
            TestFeature(self,
                name = 'DraI', start = 29, end = 4, cut = 1, clockwise = True,
                sequence = "taaatgaccctttgggatgaaagggccctt"
            ),
            # Normal feature
            TestFeature(self,
                name = 'EK', id = 15, start = 4103, end = 6, clockwise = True,
                sequence = 'gacaag' + 't' * 4096  + 'gatgacgac'
            ),
            # ORF
            TestFeature(self,
                name = 'ORF frame 1', start = 430, end = 420,
                sequence = "atcggcaaggtgtgggtcgcggacgacggcgccgcggtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggcgtctcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaggcggccgagcgcgccggggtgcccgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgcctga" + \
                't' * 9 + \
"atggcagcgcgccgaccgcgatgggctgtggccaatagcggctgctcagcagggcgcgccgagagcagcggccgggaaggggcggtgcgggaggcggggtgtggggcggtagtgtgggccctgttcctgcccgcgcggtgttccgcattctgcaagcctccggagcgcacgtcggcagtcggctccctcgttgaccgaatcaccgacctctctccccagggggatccaccggagcttaccatgaccgagtacaagcccacggtgcgcctcgccacccgcgacgacgtccccagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgac"
            ) 
        ]

        for cf in cross_features:
            features = self.find_features(cf.sequence)
            cf.assertEqual(features[-1])



if __name__ == '__main__':
    unittest.main()


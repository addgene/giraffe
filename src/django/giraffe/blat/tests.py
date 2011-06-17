import sys
sys.path.append('../../')
import re
import unittest
import giraffe.simple_test as st
st.setup()
import giraffe.blat.models as models
import django

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
        features = self.find_features('GATGACGACGACAAG')
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]['feature'], 'EK')
        self.assertEqual(features[0]['feature_id'], 15)

    def test_ItDetectsFeaturesInLowerCase(self):
        """Tests that lowercase sequences still detect features"""
        # Detect the feature
        features = self.find_features('gatgacgacgacaag')
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]['feature'], 'EK')
        self.assertEqual(features[0]['feature_id'], 15)

    def test_ItDetectsNothingInAFeatureLessSequence(self):
        """Tests that sequences which clearly have no features in them do not
        somehow generate features"""

        features = self.find_features('A')
        self.assertEqual(len(features), 0)

        features = self.find_features('ATGC')
        self.assertEqual(len(features), 0)

        features = self.find_features('T'*4096)
        self.assertEqual(len(features), 0)

    def test_ItIgnoresInvalidTextInSequences(self):
        """Tests that FASTA comments and invalid characters are cut"""

        def asserts(features):
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0]['feature'], 'EK')
            self.assertEqual(features[0]['feature_id'], 15)

        features = self.find_features('GA5TGA1CGAC2392GA2CAA1G')
        asserts(features)

        features = self.find_features(' {+G( A%T[G]A1;CG  > >AC"2<3@&92,G~`A2C.A?A/G}')
        asserts(features)

        features = self.find_features("""
>EK | feature id 15    
01 GATG 04
; This is a comment!
;; So is this!
05 ACGA 08
;; this is too
; is this?
09 CGAC 10
11 AAG  13""")
        asserts(features)
    
        self.assertRaises(self.find_features, models.BadSequencError, """
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
                    feature_sequence = line_items[1]
                    feature_name = line_items[0].split(':')[1]

                    # Detect the features
                    features = self.find_features(feature_sequence)

                    # How many of them have the same length?
                    same_length = [f for f in features
                         if abs(f['end'] - f['start']) + 1 ==
                         len(feature_sequence)]

                    # Only one feature has the same length
                    try:
                        self.assertEqual(len(same_length), 1)
                    except:
                        print same_length
                        raise

                    # That feature should have same ID and name
                    self.assertEqual(same_length[0]['feature'], feature_name)
                    self.assertEqual(same_length[0]['feature_id'], feature_id)
                    self.assertEqual(same_length[0]['start'], 1)
                    self.assertEqual(same_length[0]['end'],
                            len(feature_sequence))

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
                    features_to_test[feature_id] = {
                        'feature_sequence': line_items[1],
                        'feature_name': line_items[0].split(':')[1]
                    }

        # Make the combined feature sequences
        for group in zip(feature_ids_to_test[0:10],
                         feature_ids_to_test[10:20], 
                         feature_ids_to_test[20:30]):

            # Create the long sequence and groups of features to detect
            sequence = ''
            names = []
            for id in group:
                sequence += features_to_test[id]['feature_sequence']
                names.append(features_to_test[id]['feature_name'])

            # Detect the features
            features = self.find_features(sequence)
            
            # Do we detect at least as many features as we put in?
            self.assertTrue(len(features) >= len(group))

            # Are all of the names there?
            for (nx, name) in enumerate(names):
                featureFound = False
                for feature in features:
                    if feature['feature']    == name and \
                       feature['feature_id'] == group[nx]:
                        featureFound = True
                        break

                self.assertTrue(featureFound)


    def test_ItDetectsLoneORFs(self):
        orf_seq = "atggcagcgcgccgaccgcgatgggctgtggccaatagcggctgctcagcagggcgcgccgagagcagcggccgggaaggggcggtgcgggaggcggggtgtggggcggtagtgtgggccctgttcctgcccgcgcggtgttccgcattctgcaagcctccggagcgcacgtcggcagtcggctccctcgttgaccgaatcaccgacctctctccccagggggatccaccggagcttaccatgaccgagtacaagcccacggtgcgcctcgccacccgcgacgacgtccccagggccgtacgcaccctcgccgccgcgttcgccgactaccccgccacgcgccacaccgtcgatccggaccgccacatcgagcgggtcaccgagctgcaagaactcttcctcacgcgcgtcgggctcgacatcggcaaggtgtgggtcgcggacgacggcgccgcggtggcggtctggaccacgccggagagcgtcgaagcgggggcggtgttcgccgagatcggcccgcgcatggccgagttgagcggttcccggctggccgcgcagcaacagatggaaggcctcctggcgccgcaccggcccaaggagcccgcgtggttcctggccaccgtcggcgtctcgcccgaccaccagggcaagggtctgggcagcgccgtcgtgctccccggagtggaggcggccgagcgcgccggggtgcccgccttcctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgcaagcccggtgcctga"
        features = self.find_features(orf_seq)

        orfs = [f for f in features if re.match('ORF', f['feature'])]

        self.assertTrue(len(orfs) == 1)
        self.assertEqual(orfs[0]['feature'], 'ORF frame 1')
        self.assertEqual(orfs[0]['start'], 1)
        self.assertEqual(orfs[0]['end'], len(orf_seq))
        self.assertTrue(orfs[0]['clockwise'])
        
    def test_ItDetectsFeaturesThatCrossTheBoundary(self):
        boundary_cross_seq = "aaatgaccctttgggatgaaagggcccttt"

        features = self.find_features(boundary_cross_seq)
        
        self.assertEqual(features[-1]['feature'], 'DraI')
        self.assertTrue (features[-1]['clockwise'])
        self.assertEqual(features[-1]['start'], 28)
        self.assertEqual(features[-1]['end'], 3)

        


if __name__ == '__main__':
    unittest.main()


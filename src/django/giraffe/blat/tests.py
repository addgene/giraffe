import sys
sys.path.append('../../')
import unittest
import giraffe.simple_test as st
st.setup()
import giraffe.blat.models as models

class ItDetectsFeaturesInDNASequences(unittest.TestCase):

    def find_features(self, sequence, db_name = 'default'):
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
                 .select_related( 'feature_db_index__feature'):

            features.append({
                'id': f.feature.id,
                'name': f.feature.name
            })

        return features


    def test_ItDetectsFeatureOnlySequences(self):
        """
        Tests detection of complete feature sequence
        """

        # EK feature in default db
        features = self.find_features('GATGACGACGACAAG');

        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]['id'],   15)
        self.assertEqual(features[0]['name'], 'EK')

if __name__ == '__main__':
    unittest.main()


import os
import tempfile

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature_DB_Index
from giraffe.blat.models import Feature
from giraffe.blat.frags.frags_to_features import frags_to_features

import time


def _get_frags(db_name,sequence):
    """
    Returns list of feature fragments for the sequence, after
    detecting features using the specified feature database. Each
    member of the list is a string with 4 numbers:

    feature index
    fragment index
    start position
    shift
    """

    BIN_PATH = '/usr/local/bin/'
    DATA_PATH = '/usr/local/share/giraffe/'

    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    # XXX should we write sequence twice to handle boundary?
    tmp_file.write(sequence)
    tmp_file.close()

    cmd = '%s/giraffe-frags %s/%s.data %s' % (
        BIN_PATH, DATA_PATH, db_name, tmp_file.name
    )
    f = os.popen(cmd)
    res = f.readlines()
    os.unlink(tmp_file.name)

    for i,line in enumerate(res):
        if i == 0:
            if not line.startswith("===="):
                raise Exception('Error: '+line)
        res[i] = line.strip()
    res.pop(0)
    # print str(res)

    return res

def _pick_good_features(seq_features, sequence):
    good_seq_features = []

    # TODO: make this function a method of Seq_Feat objects
    def left_pos(feat):
        return feat.start if feat.clockwise else 2 * feat.start - feat.end

    # The returns left-most position of sequence. if it's
    # clockwise, it's just the start. if not, it's the startpos - length =
    # 2*start - end
    seq_features.sort(key = left_pos)

    print len(sequence.sequence)

    for seq_feature_outer in seq_features:
        left_pos_outer = left_pos(seq_feature_outer)

        if left_pos_outer > len(sequence.sequence):
            continue

        # XXX FeatureType-Dependent Code!
        if seq_feature_outer.feature.type.type != "Enzyme":
            found = False

            for seq_feature_inner in seq_features:
                left_pos_inner = left_pos(seq_feature_inner)

                if left_pos_inner > seq_feature_outer.end:
                    break

                # This long test checks if it's a gene, or if it's not, if the
                # outer and inner features are sufficiently similar, and if that
                # is the case, looks to see if the features overlap for any part
                # of the sequence
                if (seq_feature_inner.feature.type.type == "Gene" or
                    (seq_feature_inner.feature.type.type == seq_feature_outer.feature.type.type and 
                     (seq_feature_inner.feature.name.find(seq_feature_outer.feature.name) >= 0 or
                      seq_feature_inner.feature.name.find(seq_feature_outer.feature.name) >= 0))) and \
                   ((left_pos_inner >= left_pos_outer and left_pos_inner <= seq_feature_outer.end) or 
                    (left_pos_inner <= left_pos_outer and left_pos_outer <= seq_feature_inner.end)):
                    # If the conditions are met, check the scores, and if the 
                    # inner feature scores better, keep it instead of the outer
                    # one

                    if seq_feature_inner.score < seq_feature_outer.score:
                        print "Pruning %s (%d - %d)" % (seq_feature_outer.feature.name, 
                                                        seq_feature_outer.start,
                                                        seq_feature_outer.end)

                        found = True
                        break

            # If we found a better feature in the inner loop, don't add this
            # feature. When the outer loop gets to it on its own, it will add it
            if found:
                continue
        # else: Enzymes always get added.

        good_seq_features.append(seq_feature_outer)

    return good_seq_features

def blat(db,sequence):
    sequence = Sequence.strip(sequence)

    # create sequence record
    s = Sequence()
    s.sequence = sequence
    s.db = db
    s.save()
    s.clear_features()

    frags = _get_frags(db.name,sequence)

    # Translate frags to features and store them into the database
    features = frags_to_features(frags, db)

    good_features = _pick_good_features(features, s)

    # Store features into database
    print "Storing %d features" % len(good_features)
    t0 = time.time()
    for fx, feature in enumerate(good_features):
        feature.sequence = s
        feature.save()
    t0 = time.time() - t0
    print "database add took %f seconds" % t0

    return s


import os
import tempfile

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature_DB_Index
from giraffe.blat.models import Feature
from giraffe.blat.frags.frags_to_features import frags_to_features

import time


def get_frags(db_name,sequence):
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

def blat(db,sequence):
    sequence = Sequence.strip(sequence)

    # create sequence record
    s = Sequence()
    s.sequence = sequence
    s.db = db
    s.save()
    s.clear_features()

    frags = get_frags(db.name,sequence)
    frag_names = set()

    for frag in frags:
        fdb = Feature_DB_Index.objects.get(
                       feature_index=int(frag.split()[0]),
                       db=db)
        frag_names.add("%s %s" % (fdb.feature.name, "antisense" if fdb.antisense else
            ''))
    for fname in frag_names:
        print fname

    # Translate frags to features and store them into the database
    features = frags_to_features(frags, db)
    print "exited frags_to_features with %d" % len(features)

    # Store features into database
    t0 = time.time()
    for fx, feature in enumerate(features):
        feature.sequence = s
        feature.save()
    t0 = time.time() - t0
    print "database add took %f seconds" % t0

    return s


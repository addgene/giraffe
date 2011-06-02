import os
import tempfile

from django.conf import settings

from giraffe.blat.models import Sequence
from giraffe.blat.frags.frags_to_features import frags_to_features

# For Debugging
import time
_debug = settings.DEBUG
def timer():
    t0 = 0
    while True:
        t = time.time()
        yield t - t0
        t0 = t


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

    DATA_PATH = BIN_PATH = os.path.dirname(__file__)

    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    tmp_file.write(sequence)
    tmp_file.close()

    cmd = '%s/bin/frags %s/data/%s.data %s' % (
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

def blat(db,sequence_obj):
    if _debug: 
        t = timer()
        t.next()
    frags = _get_frags(db.name,Sequence.convert_to_dna(sequence_obj.sequence))
    if _debug: print "get_frags took %f seconds" % t.next()

    # Translate frags to features and store them into the database
    seq_length = len(sequence_obj.sequence)
    features = frags_to_features(frags, db, seq_length)
    if _debug: print "frags_to_features took %f seconds" % t.next()

    # Store features into database
    if _debug: print "Storing %d features" % len(features)
    for feature in features:
        feature.sequence = sequence_obj
        feature.save()
    if _debug: print "Database adds took %f seconds" % t.next()


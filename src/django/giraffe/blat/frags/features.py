import os
import tempfile

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature_DB_Index
from giraffe.blat.models import Feature
from giraffe.blat.frags.frags_to_features import frags_to_features

from django.db import transaction

# For Debugging
import time
_debug = True
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

@transaction.commit_manually()
def blat(db,sequence):
    """ Analyze the features of a sequence. """
    sequence = Sequence.strip(sequence)

    # create sequence record
    try:
        s = Sequence()
        s.sequence = sequence
        s.db = db
        s.save()
        s.clear_features()
    except:
        transaction.rollback()
    else:
        transaction.commit()

    t = timer()
    t.next()
    frags = _get_frags(db.name,sequence)
    if _debug: print "get_frags took %f seconds" % t.next()

    # Translate frags to features and store them into the database
    seq_length = len(sequence)
    features = frags_to_features(frags, db, seq_length)
    if _debug: print "frags_to_features took %f seconds" % t.next()

    # Store features into database
    if _debug: print "Storing %d features" % len(features)
    try:
        for feature in features:
            feature.sequence = s
            feature.save()
    except Exception, e:
        transaction.rollback()
        raise e
    else:
        transaction.commit()

    if _debug: print "Database adds took %f seconds" % t.next()

    return s


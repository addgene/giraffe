import os
import tempfile

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature


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


def blat(db,sequence):
    sequence = Sequence.strip(sequence)

    frags = get_frags(db.name,sequence)

    print str(frags)
    # XXX translate frags to features

    # create sequence record
    s = Sequence()
    s.sequence = sequence
    s.db = db
    s.save()
    s.clear_features()

    f = Sequence_Feature()
    f.sequence = s
    f.feature = Feature.objects.all()[0:1][0]
    f.start = 0
    f.end = len(f.feature.sequence)
    f.clockwise = False
    f.save()
   
    return s


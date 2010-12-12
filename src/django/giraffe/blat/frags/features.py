import os
import tempfile
import hashlib

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature


def get_frags(db_name,sequence):

    # XXX how to find the binary?
    PATH = '/Users/benjie/git/giraffe/src/django/giraffe/blat/frags/'

    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    # XXX should we write sequence twice to handle boundary?
    tmp_file.write(sequence)
    tmp_file.close()

    cmd = '%s/frags %s/%s.data %s' % (
        PATH, PATH, db_name, tmp_file.name
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

    # XXX translate frags to features

    hash = hashlib.sha1(sequence).hexdigest()

    # create sequence record
    s = Sequence.objects.filter(hash=hash,db=db).delete()
    s = Sequence()
    s.sequence = sequence
    s.db = db
    s.save()

    f = Sequence_Feature()
    f.sequence = s
    f.feature = Feature.objects.all()[0:1][0]
    f.start = 0
    f.end = len(f.feature.sequence)
    f.clockwise = False
    f.save()
   
    return s


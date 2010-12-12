import os
import tempfile
import hashlib

from giraffe.blat.models import Sequence
from giraffe.blat.models import Sequence_Feature
from giraffe.blat.models import Feature


def get_frags(db_name,sequence):
    # XXX should we double up the sequence to detect features across
    # the boundary?

    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    print tmp_file.name
    tmp_file.write(sequence)
    tmp_file.close()

    f = os.popen('/Users/benjie/git/giraffe/src/django/giraffe/blat/frags/frags %s %s' % (db_name, tmp_file.name))
    res = f.readlines()

    print str(res)
    os.unlink(tmp_file.name)
    return res


def blat(db,sequence):
    frags = get_frags(db.name,sequence)
    # XXX translate frags to features

    seq = sequence.strip()
    hash = hashlib.sha1(seq).hexdigest()

    # create sequence record
    s = Sequence.objects.filter(hash=hash,db=db).delete()
    s = Sequence()
    s.sequence = seq
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


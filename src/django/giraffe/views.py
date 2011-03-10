import time
from django.shortcuts import render_to_response
from django.template import RequestContext
import blat.models


def test_draw(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.
    """
    db = blat.models.Feature_Database.objects.get(name=db_name)
    sequence = blat.models.Sequence.objects.get(db=db,hash=hash)
    ts = int(time.mktime(sequence.modified.timetuple()))

    return render_to_response(
        'test/draw.html', { "hash" : hash, "mtime" : ts },
        context_instance=RequestContext(request)
    )


def test_analyze(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.
    """
    db = blat.models.Feature_Database.objects.get(name=db_name)
    sequence = blat.models.Sequence.objects.get(db=db,hash=hash)
    ts = int(time.mktime(sequence.modified.timetuple()))

    return render_to_response(
        'test/analyze.html', { "hash" : hash, "mtime" : ts },
        context_instance=RequestContext(request)
    )




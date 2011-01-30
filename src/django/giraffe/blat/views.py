from django.http import HttpResponse
from django.core import serializers
import json
import httplib

import models
import frags.features

from django.shortcuts import redirect
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext


def post(request):
    if request.method == 'GET':
        return render_to_response(
            'blat/post.html', {}, 
            context_instance=RequestContext(request)
        )
    else:
        db_name = request.POST['db']
        sequence = request.POST['sequence']
        db = models.Feature_Database.objects.get(name=db_name)
        s = frags.features.blat(db,sequence)
        return redirect(reverse(get,args=[s.hash,db_name]))


def get(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.
    """
    db = models.Feature_Database.objects.get(name=db_name)
    sequence = models.Sequence.objects.get(db=db,hash=hash)

    res = [len(sequence.sequence),]
    for f in sequence.sequence_feature_set.order_by("start"):
        res.append(f.to_dict())

    j = json.JSONEncoder().encode(res)

    if 'jsonp' in request.GET:
        j = request.GET['jsonp']+'('+j+')'
        http_res = HttpResponse(
            j,mimetype="text/javascript",status=httplib.OK
        )
    else:
        http_res = HttpResponse(
            j,mimetype="application/json",status=httplib.OK
        )

    return http_res


def post_check(request):
    assert(request.method == 'POST')
    db_name = request.POST['db']
    sequence = request.POST['sequence']
    (sequence,hash) = models.Sequence.clean_and_hash(sequence)
    db = models.Feature_Database.objects.get(name=db_name)
    obj = models.Sequence.objects.filter(db=db,hash=hash)
    if len(obj) > 0:
        return get(request, hash, db_name)
    else:
        return post(request)


def draw(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.
    """
    db = models.Feature_Database.objects.get(name=db_name)
    sequence = models.Sequence.objects.get(db=db,hash=hash)

    return render_to_response(
        'blat/draw.html',
        { "sequence" : sequence }, 
        context_instance=RequestContext(request)
    )


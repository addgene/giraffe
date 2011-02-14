from django.http import HttpResponse
from django.core import serializers
import json
import httplib

import models

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
        hash = models.Giraffe_Mappable_Model.detect_features(sequence,db_name)
        return redirect(reverse(get,args=[hash,db_name]))


def get(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.

    If the 'sc' key appears in GET dictionary, then return single
    cutters and non-cutter features. Otherwise, return all cutters and
    non-cutter features.
    """
    db = models.Feature_Database.objects.get(name=db_name)
    sequence = models.Sequence.objects.get(db=db,hash=hash)

    res = [len(sequence.sequence),]

    if 'sc' in request.GET:
        features = []
        cutters = {}
        for f in sequence.sequence_feature_set.order_by("start").select_related(
            'feature', 'feature__type',
        ):
            features.append(f)
            if f.feature.type_id == models.Feature_Type.ENZYME:
                if f.feature.name in cutters:
                    cutters[f.feature.name] = cutters[f.feature.name]+1
                else:
                    cutters[f.feature.name] = 1

        for f in features:
            if f.feature.type_id == models.Feature_Type.ENZYME:
                if cutters[f.feature.name] == 1:
                    res.append(f.to_dict())
            else:
                res.append(f.to_dict())

    else:
        for f in sequence.sequence_feature_set.order_by("start").select_related(
            'feature', 'feature__type',
        ):
            res.append(f.to_dict())

    j = json.JSONEncoder().encode(res)

    if 'jsonp' in request.GET:
        j = request.GET['jsonp']+'('+j+')'
        http_res = HttpResponse(
            j,mimetype="text/javascript",status=httplib.OK
        )

    else:
        http_res = HttpResponse(
            #j,mimetype="application/json",status=httplib.OK
            j,status=httplib.OK
        )

    # we tell browser to cache this; if the sequence change, the hash would
    # change. the only danger is if we re-blat the sequence, in that case the
    # features list cached by browser will be out of date.
    http_res['Cache-Control'] = 'max-age=2592000'
    return http_res


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



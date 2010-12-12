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
        db = models.Feature_Database.objects.get(name=db_name)
        sequence = request.POST['sequence']
        s = frags.features.blat(db,sequence)
        return redirect(reverse(get,args=[s.hash,db_name]))


def get(request,hash,db_name):
    """
    Get features of a sequence, using the sequence's sha-1 hash as the
    identifier.
    """
    db = models.Feature_Database.objects.get(name=db_name)
    sequence = models.Sequence.objects.get(db=db,hash=hash)

    j = serializers.get_serializer("json")()
    data = j.serialize(sequence.sequence_feature_set.all(),
                       ensure_ascii=False,
                       use_natural_keys=True)
    http_res = HttpResponse(
        data,mimetype="application/json",status=httplib.OK
    )
    return http_res


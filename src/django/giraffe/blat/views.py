from django.http import HttpResponse
import json
import httplib

import models

from django.shortcuts import redirect
from django.core.urlresolvers import reverse


"""
post view: post a sequence and run the sequence through blat and orf
detection.

    expects: db and sequence
    response:

        1. if 'next' is specified as a CGI argument, redirect to that
        URL with '/hash/db_name/' appended, where 'hash' is the hash
        ID of the sequence and 'db_name' is the name of the features
        database.

        2. if 'next' is specified and there is an error with the
        sequence or the blat algorithm, redirect to the 'next' URL
        with '/error/' appended. if 'next' is not specified and
        there is an error, just throw the error.

        3. otherwise, redirects to the 'get' view that returns JSON
        array of features.
"""
def post(request):
    assert (request.method == 'POST')
    db_name = request.POST['db']
    sequence = request.POST['sequence']
    try:
        hash = models.Giraffe_Mappable_Model.detect_features(sequence,db_name)
        if 'next' in request.POST:
            u = request.POST['next']
            if u.endswith('/'):
                u = u+hash+'/'+db_name
            else:
                u = u+'/'+hash+'/'+db_name
            return redirect(u)
        return redirect(reverse(get,args=[hash,db_name]))
    except Exception as e:
        # print 'Blat error on sequence: '+str(sequence)
        if 'next' in request.POST:
            print str(e)
            u = request.POST['next']
            if u.endswith('/'):
                u = u+'error/'
            else:
                u = u+'/error/'
            return redirect(u)
        raise e

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

    if db.db_version != sequence.db_version:
        print 'feature list and database out of sync!'
        # feature out of date with database, re gather features
        hash = models.Giraffe_Mappable_Model.detect_features(sequence.sequence,db_name)

    res = []

    # get automated features

    if 'sc' in request.GET:
        features = []
        cutters = {}
        for f in sequence.sequence_feature_set.order_by("start").select_related(
            'feature_db_index',
            'feature_db_index__feature',
            'feature_db_index__feature__type',
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
            'feature_db_index',
            'feature_db_index__feature',
            'feature_db_index__feature__type',
        ):
            res.append(f.to_dict())

    # get annotated features

    for f in sequence.sequence_feature_annotated_set.order_by(
        "start"
    ).select_related('feature_type'):
        res.append(f.to_dict())

    # now sort everything by start

    res.sort(cmp=lambda x,y:cmp(int(x['start']),int(y['start'])))

    res = [len(sequence.sequence),res]

    if 'sequence' in request.GET:
        # also asked for sequence
        res.append(sequence.sequence)

    j = json.JSONEncoder().encode(res)

    if 'jsonp' in request.GET:
        j = request.GET['jsonp']+'('+j+')'
        http_res = HttpResponse(j,mimetype="text/javascript",status=httplib.OK)

    else:
        # technically we should be returning "application/json", but
        # in that case browsers force user to download into a file,
        # and for debugging we want to be able to see the JSON list in
        # browser. looks like most browsers will handle JSON sent back
        # as text/html anyways.
        if request.is_ajax():
            http_res = HttpResponse(j,mimetype="application/json",status=httplib.OK)
        else:
            http_res = HttpResponse(j,status=httplib.OK)

    # we tell browser to cache this; if the sequence change, the hash would
    # change. the only danger is if we re-blat the sequence, in that case the
    # features list cached by browser will be out of date. so client
    # should attach some kind of CGI string to invalidate cache.
    http_res['Cache-Control'] = 'max-age=2592000'
    return http_res



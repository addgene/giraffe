from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings
from django.shortcuts import redirect

import django.contrib.admin
django.contrib.admin.autodiscover()

import views

urlpatterns = patterns('',
    (r'^admin/', include(django.contrib.admin.site.urls)),
    (r'^blat/', include('giraffe.blat.urls')),

    url(r'^test/analyze/(\w+)/(\w+)/?$', views.test_analyze, name='analyzer'),
    (r'^test/draw/(\w+)/(\w+)/?$', views.test_draw),

    (r'^/test/$', direct_to_template, { 'template' : 'test/post.html' }),
    (r'^/?$', direct_to_template, { 'template' : 'test/post.html' }),
)

if settings.DEBUG:
    urlpatterns += patterns('',
        (r'^static/blat/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_BLAT_ROOT}),
        (r'^static/analyze/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_ANALYZE_ROOT}),
        (r'^static/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_ROOT}),
    )


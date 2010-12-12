from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from django.conf import settings
from django.shortcuts import redirect

import django.contrib.admin
django.contrib.admin.autodiscover()

urlpatterns = patterns('',
    (r'^admin/', include(django.contrib.admin.site.urls)),
    (r'^blat/', include('giraffe.blat.urls')),
)

if settings.DEBUG:
    urlpatterns += patterns('',
        (r'^headers/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_HEADERS_ROOT}),
        (r'^images/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_IMAGES_ROOT}),
        (r'^static/(?P<path>.*)$', 'django.views.static.serve',
         {'document_root': settings.STATIC_ROOT}),
    )


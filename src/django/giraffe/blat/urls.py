from django.conf.urls.defaults import *
import views

urlpatterns = patterns('',
    url(r'^/?(\w+)/(\w+)/?$', views.get, name='blat-get'),
    url(r'^/?$', views.post, name='blat-post'),
)


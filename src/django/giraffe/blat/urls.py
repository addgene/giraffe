from django.conf.urls.defaults import *
import views

urlpatterns = patterns('',
    (r'^/?(\w+)/(.+)/?$', views.get),
    (r'^/?(\w+)/?$', views.post),
)


from django.conf.urls.defaults import *
import views

urlpatterns = patterns('',
    (r'^/?(\w+)/(\w+)/?$', views.get),
    (r'^/?$', views.post),
)


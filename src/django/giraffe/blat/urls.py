from django.conf.urls.defaults import *
import views

urlpatterns = patterns('',
    (r'^/?(\w+)/(\w+)/draw/?$', views.draw),
    (r'^/?(\w+)/(\w+)/?$', views.get),
    (r'^/?$', views.post),
)


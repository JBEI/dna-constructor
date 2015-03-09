from django.conf.urls.defaults import *
from django.views.generic import DetailView
from constructor.models import Node

urlpatterns=patterns('',
    url(r'^$', 'constructor.views.index', name='home'),
    url(r'^construct/images/(?P<path>.*)$', 'django.views.static.serve', {'document_root': 'constructor/static/media'}),
    url(r'^construct/ajax/query/(?P<pk>\d+)', 'constructor.views.query'),
    url(r'^construct/ajax/export/(?P<pk>\d+)', 'constructor.views.export'),
    url(r'^construct/ajax/build/', 'constructor.views.build'),
    url(r'^construct/ajax/shrink/', 'constructor.views.shrink'),
    url(r'^construct/ajax/protocols/', 'constructor.views.protocols'),
    url(r'^construct/ajax/report/', 'constructor.views.report'),
    url(r'upload', 'constructor.views.upload'),
    url(r'^construct/(?P<gv_id>\w+).gv', 'constructor.views.gvreturn'),
    url(r'^construct/$', 'constructor.views.construct'),
    url(r'^construct/#', 'constructor.views.construct'),
    url(r'^copyright', 'constructor.views.copyright'),
    url(r'^disclaimer', 'constructor.views.disclaimer')
)

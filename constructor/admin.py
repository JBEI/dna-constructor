from constructor.models import *
from django.contrib import admin

class NodeAdmin(admin.ModelAdmin):
    list_display = ('protocolID', 'name')
    list_display_links = list_display
    
class ProtocolAdmin(admin.ModelAdmin):
    list_display = ('id', 'name')
    list_display_links = list_display
    
class SequenceVariableAdmin(admin.ModelAdmin):
    list_display = ('protocol', 'name')
    
class PrimerAdmin(admin.ModelAdmin):
    list_display = ['sequence']

admin.site.register(Node, NodeAdmin)
admin.site.register(Protocol, ProtocolAdmin)
admin.site.register(SequenceVariable, SequenceVariableAdmin)
admin.site.register(Primer, PrimerAdmin)

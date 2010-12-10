from django.core.exceptions import PermissionDenied
from django.contrib import admin
import models

class FeatureAdmin(admin.ModelAdmin):
    def sequence(obj):
        if len(obj.sequence) < 32:
            return obj.sequence+' ('+str(len(obj.sequence))+')'
        return obj.sequence[0:32]+'... ('+str(len(obj.sequence))+')'

    def db(obj):
        s = []
        for n in obj.feature_database_set.all():
            s.append(n.name)
        return ', '.join(s)
    db.short_description = 'Database'

    list_display = ('name','type',sequence,'cut_after','last_modified',db)
    list_filter = ('type',)
    search_fields = ('name',)

    def delete_view(self, request, object_id, extra_context=None):
        raise PermissionDenied

admin.site.register(models.Feature, FeatureAdmin)



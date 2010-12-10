from django.db import models


class Feature_Type(models.Model):
    type = models.CharField(max_length=64)

    def __unicode__(self):
        return self.type

    class Meta:
        verbose_name = "Feature Type"


class Feature(models.Model):
    type = models.ForeignKey(Feature_Type)
    name = models.CharField(max_length=32,unique=True,db_index=True)
    sequence = models.TextField()
    cut_after = models.PositiveIntegerField(null=True,blank=True)
    last_modified = models.DateTimeField(auto_now=True,db_index=True)


class Feature_Database(models.Model):
    name = models.CharField(max_length=64,unique=True)
    features = models.ManyToManyField(Feature)
    last_built = models.DateTimeField(null=True,blank=True)


class Sequence(models.Model):
    sequence = models.TextField()
    hash = models.CharField(max_length=64,db_index=True)
    modified = models.DateTimeField(auto_now=True)
    feature_database = models.ForeignKey(Feature_Database)


class Sequence_Feature(models.Model):
    sequence = models.ForeignKey(Sequence)
    feature = models.ForeignKey(Feature)
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    clockwise = models.BooleanField()




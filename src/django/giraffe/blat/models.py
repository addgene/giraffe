from django.db import models


class Feature_Type(models.Model):
    type = models.CharField(max_length=64)


class Feature(models.Model):
    type = models.ForeignKey(Feature_Type)
    name = models.CharField(max_length=32,db_index=True)
    sequence = models.TextField()
    cut_site = models.PositiveIntegerField(null=True)
    hash = models.CharField(max_length=64,db_index=True)


class Feature_Database(models.Model):
    name = models.CharField(max_length=64)
    features = models.ManyToManyField(Feature)
    last_built = models.DateTimeField()


class Sequence(models.Model):
    sequence = models.TextField()
    hash = models.CharField(max_length=64,db_index=True)
    modified = models.DateTimeField(auto_now=True)
    feature_database = models.ForeignKey(Feature_Database)


class Sequence_Feature(models.Model):
    sequence = models.FroeignKey(Sequence)
    feature = models.ForeignKey(Feature)
    start = models.PositiveInteger()
    end = models.PositiveInteger()
    clockwise = models.BooleanField()




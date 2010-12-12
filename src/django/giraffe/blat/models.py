from django.db import models
import hashlib
import re


class Sequence(models.Model):
    sequence = models.TextField()
    hash = models.CharField(max_length=64,db_index=True)
    modified = models.DateTimeField(auto_now=True)
    db = models.ForeignKey('Feature_Database')

    class Meta:
        unique_together = (("db","hash"),)

    def save(self):
        self.sequence = Sequence.strip(self.sequence)
        self.hash = hashlib.sha1(self.sequence.lower()).hexdigest()
        return super(Sequence,self).save()

    @staticmethod
    def strip(sequence):
        sequence = re.sub(r'\s', '', sequence)
        return sequence


class Sequence_Feature(models.Model):
    sequence = models.ForeignKey(Sequence)
    feature = models.ForeignKey('Feature')
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    clockwise = models.BooleanField()

    def to_dict(self):
        return {
            "feature_id" : self.feature_id,
            "start" : self.start,
            "end" : self.end,
            "clockwise" : self.clockwise
        }


class Feature_Type(models.Model):
    type = models.CharField(max_length=64)

    def __unicode__(self):
        return self.type

    class Meta:
        verbose_name = "Feature Type"


class Feature(models.Model):
    type = models.ForeignKey(Feature_Type)
    name = models.CharField(max_length=32,db_index=True)
    sequence = models.TextField()
    hash = models.CharField(max_length=64)
    cut_after = models.PositiveIntegerField(null=True,blank=True)
    last_modified = models.DateTimeField(auto_now=True,db_index=True)

    def save(self):
        self.sequence = Sequence.strip(self.sequence)
        self.hash = hashlib.sha1(self.sequence.lower()).hexdigest()
        return super(Feature,self).save()

    def __unicode__(self):
        return self.name

    class Meta:
        unique_together = (("name","hash"),)


class Feature_Database(models.Model):
    name = models.CharField(max_length=64,unique=True)
    features = models.ManyToManyField(Feature)
    last_built = models.DateTimeField(null=True,blank=True)


class Feature_DB_Index(models.Model):
    db = models.ForeignKey(Feature_Database)
    feature_index = models.PositiveIntegerField(db_index=True)
    feature = models.CharField(max_length=32)
    antisense = models.BooleanField()

    class Meta:
        unique_together = (("db","feature_index"),)




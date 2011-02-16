from django.db import models
from django.db import utils
import datetime
import hashlib
import re


class Giraffe_Mappable_Model(models.Model):
    """
    This is an abstract class for other apps to use with their models;
    this class requires the giraffe.blat app being deployed with the
    app using this class.
    """

    class Meta:
        abstract = True

    sequence = models.TextField(null=True,blank=True)
    sequence_giraffe_id = models.CharField(max_length=64,null=True,blank=True)
    sequence_giraffe_time = models.DateTimeField(null=True,blank=True)

    @staticmethod
    def detect_features(sequence,db_name):
        import frags.features
        import orfs
        db = Feature_Database.objects.get(name=db_name)
        # run blat algorithm to automatically detect features
        s = frags.features.blat(db,sequence)
        # detect ORFs
        orfs.detect_orfs(s)
        return s.hash

    def giraffe_ready(self,db_name='default',force=False,save=True):
        if not self.sequence:
            self.sequence_giraffe_id = ''
            self.sequence_giraffe_time = None
            if save:
                self.save()
            return
        run = force
        if not force:
            (s,h) = Sequence.clean_and_hash(self.sequence)
            if h != self.sequence_giraffe_id:
                run = True
        if run:
            s = Giraffe_Mappable_Model.detect_features(self.sequence,db_name)
            self.sequence_giraffe_id = s
            self.sequence_giraffe_time = datetime.datetime.now()
            if save:
                self.save()
    giraffe_ready.alters_data = True

    def save(self,run_giraffe_ready=True):
        if run_giraffe_ready:
            self.giraffe_ready(save=False)
        return super(Giraffe_Mappable_Model,self).save()
    save.alters_data = True


class Sequence_Feature_Base(models.Model):
    class Meta:
        abstract = True

    sequence = models.ForeignKey('Sequence')
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    clockwise = models.BooleanField()

    def to_dict(self):
        d = {
            "start" : self.start,
            "end" : self.end,
            "clockwise" : self.clockwise,
        }
        return d


class Sequence_Feature(Sequence_Feature_Base):
    feature_db_index = models.ForeignKey('Feature_DB_Index')

    @property
    def feature(self):
        return self.feature_db_index.feature

    # Gene variant info
    subset_start = models.PositiveIntegerField(default=0)
    subset_end = models.PositiveIntegerField(default=0)
    is_variant = models.BooleanField(default=False)
    has_gaps = models.BooleanField(default=False)

    class Meta:
        ordering = ['start','end']

    def to_dict(self):
        d = super(Sequence_Feature,self).to_dict()
        d['feature'] = self.feature.name
        d['feature_id'] = self.feature.id
        d['type_id'] = self.feature.type_id
        d['show_feature'] = self.feature_db_index.show_feature

        # Include cut position
        if d["type_id"] == Feature_Type.ENZYME:
            if d["clockwise"]:
                cp = d["start"] + (self.feature.cut_after - 1)
            else:
                cp = d["end"] - (self.feature.cut_after - 1)
            d["cut"] = cp
        # Include gene variant name modifications
        elif d["type_id"] == Feature_Type.GENE: 
            if self.is_variant:
                d["feature"] = d["feature"] + " (variant)";
            elif self.has_gaps:
                d["feature"] = d["feature"] + " (w/ gaps)";
            elif self.subset_end > 0:
                if d["clockwise"]:
                    d["feature"] = "%s (%d - %d)" % (d["feature"], self.subset_start,
                            self.subset_end);
                else:
                    d["feature"] = "%s (%d - %d)" % (d["feature"], self.subset_end,
                            self.subset_start);

        return d


class Sequence(models.Model):
    @staticmethod
    def strip(sequence):
        if sequence:
            sequence = re.sub(r'\s', '', sequence)
        else:
            sequence = ''
        return sequence

    @staticmethod
    def clean_and_hash(sequence):
        """
        Takes a sequence, returns a tuple, the cleaned sequence and
        the sequence hash.
        """
        sequence = Sequence.strip(sequence)
        hash = hashlib.sha1(sequence.lower()).hexdigest()
        return (sequence,hash)

    sequence = models.TextField()
    hash = models.CharField(max_length=64,db_index=True)
    modified = models.DateTimeField(auto_now=True)
    db = models.ForeignKey('Feature_Database')

    class Meta:
        unique_together = (("db","hash"),)

    def save(self):
        (self.sequence,self.hash) = Sequence.clean_and_hash(self.sequence)
        try:
            super(Sequence,self).save()
        except utils.IntegrityError as e:
            s = Sequence.objects.get(hash=self.hash,db=self.db)
            self.id = s.id
    save.alters_data = True

    def clear_features(self,feature_type=None):
        a = { 'sequence' : self }
        if feature_type:
            a['feature__type'] = feature_type
        Sequence_Feature.objects.filter(**a).delete()
    clear_features.alters_data = True 

    def clear_annotated_features(self,feature_type=None):
        a = { 'sequence' : self }
        if feature_type:
            a['feature_type'] = feature_type
        Sequence_Feature_Annotated.objects.filter(**a).delete()
    clear_annotated_features.alters_data = True 

    def clear_orf_features(self):
        Sequence_Feature_Annotated.objects.filter(
            feature_type=Feature_Type.ORF,sequence=self
        ).delete()
        Sequence_Feature_Annotated.objects.filter(
            orf_annotated__isnull=False,sequence=self
        ).delete()
    clear_orf_features.alters_data = True 


class Feature_Type(models.Model):
    type = models.CharField(max_length=64)

    # Feature type ID constants
    (FEATURE,    PROMOTER,   PRIMER,
     ENZYME,     GENE,       ORIGIN,
     REGULATORY, TERMINATOR, EXACT_FEATURE,
     ORF) = range(1, 11)

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
    features = models.ManyToManyField(Feature, through='Feature_In_Database')
    last_built = models.DateTimeField(null=True,blank=True)


class Feature_In_Database(models.Model):
    feature = models.ForeignKey(Feature)
    feature_database = models.ForeignKey(Feature_Database)
    show_feature = models.BooleanField(default=True)


class Feature_DB_Index(models.Model):
    db = models.ForeignKey(Feature_Database)
    feature_index = models.PositiveIntegerField(db_index=True)
    feature = models.ForeignKey(Feature)
    antisense = models.BooleanField()
    show_feature = models.BooleanField(default=True)

    class Meta:
        unique_together = (("db","feature_index"),)


class Sequence_Feature_Annotated(Sequence_Feature_Base):
    feature_name = models.CharField(max_length=64)
    feature_type = models.ForeignKey(Feature_Type)

    orf_frame = models.PositiveIntegerField(null=True)
    orf_annotated = models.ForeignKey('Sequence_Feature_Annotated',null=True)

    def to_dict(self):
        d = super(Sequence_Feature_Annotated,self).to_dict()
        d['feature'] = self.feature_name
        d['type_id'] = self.feature_type_id
        if self.orf_frame:
            d['orf_frame'] = self.orf_frame
        return d


#  import_features.py
#  Author:  John Furr
#  Date:  11/1/2012
#
#  This script is used to read feature files from 
#  blat/fixtures/features/ 
#  and import them into the Features Tables
#
#  Each Feature has a ForeignKey back to the Feature_Database
#
#  The management script has two required options
#  --db=  can be 'default', 'afire', or 'none'
#  --file=  Should be one of the feature files in blat/fixtures/features/ 
#
# class Feature(models.Model):
#    type = models.ForeignKey(Feature_Type)
#    name = models.CharField(max_length=32,db_index=True)
#    db   = models.ForeignKey('Feature_Database')
#    sequence = models.TextField()
#    hash = models.CharField(max_length=64)
#    cut_after = models.PositiveIntegerField(null=True,blank=True)
#    last_modified = models.DateTimeField(auto_now=True,db_index=True)
#    last_built = models.DateTimeField(null=True,blank=True)
#
# class Feature_Database(models.Model):
#    features = models.ManyToManyField(Feature, through='Feature_In_Database')
#    name = models.CharField(max_length=64,unique=True)
#    db_version = models.CharField(max_length=64)
#
#
# class Feature_DB_Index(models.Model):
#    db = models.ForeignKey(Feature_Database)
#    feature = models.ForeignKey(Feature)
#    feature_index = models.PositiveIntegerField(db_index=True)
#    antisense = models.BooleanField()
#    show_feature = models.BooleanField(default=True)

## Python specific imports
import os
import sys
from optparse import make_option
import re
import hashlib
import datetime

## Django specific imports
from django.core.management import BaseCommand

## Addgene/giraffe specific imports
from giraffe.blat.models import Feature, Feature_Type, Feature_Database
from giraffe.settings import FEATURES as feature_dir


class Command(BaseCommand):

    help = u'''\
        Import new Features into blat_features
        '''

    option_list = BaseCommand.option_list + (
        make_option(
            '--db',
            action='store',
            dest='db',
            default='',
            help='The Feature Database (default, afire, none)'
        ),
        make_option(
            '--all',
            action='store_true',
            dest='all',
            default=False,
            help='Run all feature files for the given database. ex. blat/fixtures/features/default/*.feature'
        ),
        make_option(
            '--file',
            action='store',
            dest='file',
            default='',
            help='The Feature File  fixtures/features/*.feature'
        ),
    )


    def get_feature_database(self, db_name):
        '''
            There are currently 3 options that we use for this process
            default, afire, none
        '''
        if db_name == 'none':
            return False

        try:
            fdb = Feature_Database.objects.get(name=db_name)
        except:
            print "Creating new Feature_Database: ", db_name 
            fdb = Feature_Database()
            fdb.name = db_name
            fdb.db_version = hashlib.sha1(str(datetime.datetime.now())).hexdigest()
            fdb.save()

        return fdb


    def create_enzyme_feature(self, line, fdb):
        """ 
            Enzyme have two possible format options
            E:Name,1/5 Sequence
            E*:Name,1/5 Sequence
            
            The 1 above is the cut_after point
            I'd assume the 5 is the cut_before point???
            E specifies that it's an enzyme
            E* is a special enzyme?????  Need to finish documenting this

            Example Explanation from Jason
            E:EcoRI,1/5 GAATTC
            1 and 5 are cut points depending upon whether or
            not the sequence is translated via the 3' end or the 5' end
        """
        m = re.match('^E:(\w+),(\d)\/(\d) (.+)$',line)
        show_feature=False
        if not m:
            m = re.match('^E\*:(\w+),(\d)\/(\d) (.+)$',line)
            show_feature=True

        name = m.group(1)
        cut_after = int(m.group(2))
        sequence = m.group(4)


        ## Create a new feature or get the current feature from the database
        f = Feature()
        if fdb:
            f.db = fdb
        f.type = Feature_Type.objects.get(type='Enzyme')
        f.name = name
        f.sequence = sequence
        f.cut_after = cut_after
        f.show_feature = show_feature

        try:
            f.save()
        except:
            f = Feature.objects.get(name=name,hash=f.hash)

        ## If the feature database was set then we 
        ## set the feature field of this database to teh new feature
        #if fdb:
        #    m = Feature_In_Database()
        #    m.feature = f
        #    m.feature_database = fdb
        #    m.show_feature = show_feature
        #    m.save() 
    
    def get_feature_type(self, n):
        if n == "F":   return Feature_Type.objects.get(type='Feature')
        if n == "G":   return Feature_Type.objects.get(type='Gene')
        if n == "P":   return Feature_Type.objects.get(type='Promoter')
        if n == "O":   return Feature_Type.objects.get(type='Origin')
        if n == "R":   return Feature_Type.objects.get(type='Regulatory')
        if n == "T":   return Feature_Type.objects.get(type='Terminator')
        if n == "f":   return Feature_Type.objects.get(type='ExactFeature')
        if n == "S":   return Feature_Type.objects.get(type='Primer')
        if n == "E":   return Feature_Type.objects.get(type='Enzyme')
        if n == "E*":  return Feature_Type.objects.get(type='Enzyme')

    def create_other_feature(self, line, fdb):
        """ 
            ALl other features come in a standard format
            F:Name  Sequence
            
            F can be, F, G, P, O, R, T, E, E*, f, S
            See get_feature_types() above
        """

        #Don't process blank lines or comments
        if not line.strip():  return
        if line[0] == '#': return

        # Match for F:Name Sequence
        m = re.match('^(\w+):(\S+) (.+)$',line)
        if not m:
            return

        n = m.group(1)
        name = m.group(2)
        sequence = m.group(3)

        f = Feature()
        f.type = self.get_feature_type(n)
        if fdb:
            f.db = fdb
        f.name = name
        f.sequence = sequence
        if fdb: f.show_feature = True
        else: f.show_feature = False
        
        try:
            f.save()
        except Exception as err:
            print "Exception: ", str(err)
            f = Feature.objects.get(name=name,hash=f.hash)

        #if fdb:
        #    m = Feature_In_Database()
        #    m.feature = f
        #    m.feature_database = fdb
        #    m.show_feature = True
        #    m.save()
 
 
    def handle(self, *args, **options):
        ##Get or Create a new feature_database
        ## Current options are default, afire, none
        fdb = self.get_feature_database( options['db'] )    


        ## Now parse the *.feautre files found in blat/fixtures/features/
        ## There are three line format types
        ##
        ## 2 Enzyme Formatting options
        ## E:Name,1/5 sequence
        ## E*:Name,1/5 sequence
        ##
        ## Other feature type lines
        ## F:Name Sequence
        ## F can be one of the following
        ## F: Feature
        ## G: Gene
        ## P: Promoter
        ## O: Origin
        ## R: Regulatory
        ## T: Terminal
        ## E: Enzyme ??
        ## f: exact Feature
        ## S: Primer
        if options['all']:
            feature_files = os.listdir(feature_dir+"/"+options['db'])

            for file in feature_files:
                file_name = feature_dir+"/"+options['db'] + "/" + file
                print "Installing: ", file_name
                with open(file_name, 'r') as f:
                    for line in f.readlines():
                        if line[0] == 'E':
                            self.create_enzyme_feature(line, fdb)
                        else:
                            self.create_other_feature(line, fdb)
        
        elif options['file']:
            with open(options['file'], 'r') as f:
                for line in f.readlines():
                    if line[0] == 'E':
                        self.create_enzyme_feature(line, fdb)
                    else:
                        self.create_other_feature(line, fdb)
        




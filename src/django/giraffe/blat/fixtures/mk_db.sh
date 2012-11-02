#!/bin/sh

## Use the -W ignore::UserWarning to hide
## "The virtualenv distutils package at %s appears to be in the same location as the system distutils?")
PYTHON="python -W ignore::UserWarning"
GIRAFFE="/Users/johnfurr/Git_Tower/giraffe/src/django/giraffe"
FIXTURE="$GIRAFFE/blat/fixtures/"
hname=`hostname`
if [ "$hname" = "carbon" ]
then
	PYTHON='/srv/addgene/bin/python -W ignore::UserWarning'
    GIRAFFE="/home/benjie/git/giraffe/src/django/giraffe"
    FIXTURE="$GIRAFFE/blat/fixtures/"
fi

## This is the older method.  It's left here for refernce sake at this point
## Once the code has been fully converted these lines can be deleted
#$PYTHON import_features.py default < features/generic.features
#$PYTHON import_features.py default < features/generic.primers
#$PYTHON import_features.py default < features/tb.enzymes
#$PYTHON import_features.py default < features/fp.features
#$PYTHON create_frag_db.py default > ../frags/data/default.data
#
#$PYTHON import_features.py afire < features/af.unc
#$PYTHON import_features.py afire < features/af.custom
#$PYTHON import_features.py afire < features/af.features
#$PYTHON import_features.py afire < features/af.enzymes
#$PYTHON create_frag_db.py afire > ../frags/data/afire.data
#
#$PYTHON import_features.py none < features/all.enzymes

## Hold Current Directory
pwd=$PWD

cd $GIRAFFE

#############################################################################################
####                   Default Database Creation                                      
echo "Importing blat/fixtures/features/generic.features to Feature_DataBase(default)"
$PYTHON manage.py import_features --db=default --file=blat/fixtures/features/generic.features

echo "Importing blat/fixtures/features/generic.primers to Feature_DataBase(default)"
$PYTHON manage.py import_features --db=default --file=blat/fixtures/features/generic.primers

echo "Importing blat/fixtures/features/tb.enzymes to Feature_DataBase(default)"
$PYTHON manage.py import_features --db=default --file=blat/fixtures/features/tb.enzymes

echo "Importing blat/fixtures/features/fp.features to Feature_DataBase(default)"
$PYTHON manage.py import_features --db=default --file=blat/fixtures/features/fp.features

cd $FIXTURE
echo "Creating Frag DB default"
$PYTHON create_frag_db.py default > ..//frags/data/default.data
cd $GIRAFFE
##############################################################################################


##############################################################################################
####               aFire Database Creation **The afire db is never used**               ######   
echo "Importing blat/fixtures/features/af.unc to Feature_DataBase(afire)"
$PYTHON manage.py import_features --db=afire --file=blat/fixtures/features/af.unc

echo "Importing blat/fixtures/features/af.custom to Feature_DataBase(afire)"
$PYTHON manage.py import_features --db=afire --file=blat/fixtures/features/af.custom

echo "Importing blat/fixtures/features/af.features to Feature_DataBase(afire)"
$PYTHON manage.py import_features --db=afire --file=blat/fixtures/features/af.features
echo "Importing blat/fixtures/features/af.enzymes to Feature_DataBase(afire)"
$PYTHON manage.py import_features --db=afire --file=blat/fixtures/features/af.enzymes

cd $FIXTURE
echo "Creating Frag DB afire"
$PYTHON create_frag_db.py afire > ../frags/data/afire.data
cd $GIRAFFE
##############################################################################################

##############################################################################################
##                    Import all.enayzmes no database
#echo "Importing blat/fixtures/features/all.enzymes to Feature_DataBase(none)"
#$PYTHON manage.py import_features --db=none --file=blat/fixtures/features/all.enzymes
##############################################################################################

## Now send the user back to were the came from
cd $pwd

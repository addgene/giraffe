#!/bin/sh

## Use the -W ignore::UserWarning to hide
## "The virtualenv distutils package at %s appears to be in the same location as the system distutils?")

PYTHON='/srv/addgene/bin/python -W ignore::UserWarning'
GIRAFFE="/home/benjie/git/giraffe/src/django/giraffe"
FIXTURE="$GIRAFFE/blat/fixtures/"

hname=`hostname`
if [ "$hname" = "John-Furrs-iMac.local" ]
then
    PYTHON="python -W ignore::UserWarning"
    GIRAFFE="/Users/johnfurr/Git_Tower/giraffe/src/django/giraffe"
    FIXTURE="$GIRAFFE/blat/fixtures/"
fi


## Save users current working directory
pwd=$PWD

## Head to the main GIRAFFE directory
cd $GIRAFFE

# Import default features
echo "Importing blat/fixtures/features/default/*.feature"
$PYTHON manage.py import_features --db=default --all

## Import afire features (not ever used)
echo "Importing blat/fixtures/features/afire/*.feature"
$PYTHON manage.py import_features --db=afire --all



##                  Create data files
cd $FIXTURE

echo "Creating Frag DB default"
$PYTHON create_frag_db.py default > ..//frags/data/default.data
echo "Creating Frag DB afire...never used" # never used
$PYTHON create_frag_db.py afire > ../frags/data/afire.data

cd $GIRAFFE


# echo"Import none features (Never used so commented out for new format"
#$PYTHON manage.py import_features --db=none --file=blat/fixtures/features/all.enzymes

## Now send the user back to were the came from
cd $pwd

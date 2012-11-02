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

## Save users current working directory
pwd=$PWD

## Head to the main GIRAFFE directory
cd $GIRAFFE

#############################################################################################
####                   Default Database Creation                                      
echo "Importing blat/fixtures/features/default/*.feature"
$PYTHON manage.py import_features --db=default --all

cd $FIXTURE
echo "Creating Frag DB default"
$PYTHON create_frag_db.py default > ..//frags/data/default.data
cd $GIRAFFE
##############################################################################################


##############################################################################################
####               aFire Database Creation **The afire db is never used**               ######   
echo "Importing blat/fixtures/features/afire/*.feature"
$PYTHON manage.py import_features --db=afire --all

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

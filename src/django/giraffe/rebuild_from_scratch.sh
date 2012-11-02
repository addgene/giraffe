#!/usr/bin/env/bash

echo "You are about to drop all the blat tables adn rebuild them"
echo "On production this will destroy all plasmid map"
echo "It will take about 10 hours to rebuild"

read  -p "Are you sure you want to continue? yes|no > " answer

if [ $answer != "yes" ]
then
    echo "Exiting"
    exit
fi

# This will drop the entire blat table
## Calls python manage.py syncdb when it's done
./drop_blat.sh

## Now create all the feature files
## This also writes the files: 
#   --blat/frags/data/afire.data 
#   --blat/frags/data/default.data
blat/fixtures/mk_db.sh

#Change to the addgene-core lims directory
# Currently Assumes you are at the top of the giraffe tree....above giraffe/blat"
cd ../../../../addgene-core/src/django/lims/

## Then run refresh_blat_sequence on all the sequences
## This management command all takes options:
##      --vdb               # Update vector database only
##      --available-items   # Update available catalog items only
##      --all               # Update all plasmid maps
python manage.py refresh_blat_sequence --all

#!/bin/sh

PYTHON='/srv/addgene/bin/python'

$PYTHON import_features.py default < features/generic.features
$PYTHON import_features.py default < features/generic.primers
$PYTHON import_features.py default < features/tb.enzymes
$PYTHON import_features.py default < features/fp.features
$PYTHON create_frag_db.py default > ../frags/data/default.data

$PYTHON import_features.py afire < features/af.unc
$PYTHON import_features.py afire < features/af.custom
$PYTHON import_features.py afire < features/af.features
$PYTHON import_features.py afire < features/af.enzymes
$PYTHON create_frag_db.py afire > ../frags/data/afire.data

$PYTHON import_features.py none < features/all.enzymes


#!/bin/sh

python import_features.py default < features/generic.features
python import_features.py default < features/generic.primers
python import_features.py default < features/tb.enzymes
python import_features.py default < features/fp.features
python create_frag_db.py default > ../frags/data/default.data

python import_features.py afire < features/af.unc
python import_features.py afire < features/af.custom
python import_features.py afire < features/af.features
python import_features.py afire < features/af.enzymes
python create_frag_db.py afire > ../frags/data/afire.data

python import_features.py none < features/all.enzymes


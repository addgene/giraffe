#!/usr/bin/env bash


echo "drop table blat_feature" | mysql -u root -ppassword addgene
echo "drop table blat_feature_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_db_index" | mysql -u root -ppassword addgene
echo "drop table blat_feature_in_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_type" | mysql -u root -ppassword addgene
echo "drop table blat_sequence" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature_annotated" | mysql -u root -ppassword addgene

echo "drop table blat_feature" | mysql -u root -ppassword addgene
echo "drop table blat_feature_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_db_index" | mysql -u root -ppassword addgene
echo "drop table blat_feature_in_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_type" | mysql -u root -ppassword addgene
echo "drop table blat_sequence" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature_annotated" | mysql -u root -ppassword addgene

echo "drop table blat_feature" | mysql -u root -ppassword addgene
echo "drop table blat_feature_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_db_index" | mysql -u root -ppassword addgene
echo "drop table blat_feature_in_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_type" | mysql -u root -ppassword addgene
echo "drop table blat_sequence" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature_annotated" | mysql -u root -ppassword addgene

echo "drop table blat_feature" | mysql -u root -ppassword addgene
echo "drop table blat_feature_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_db_index" | mysql -u root -ppassword addgene
echo "drop table blat_feature_in_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_type" | mysql -u root -ppassword addgene
echo "drop table blat_sequence" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature_annotated" | mysql -u root -ppassword addgene

echo "drop table blat_feature" | mysql -u root -ppassword addgene
echo "drop table blat_feature_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_db_index" | mysql -u root -ppassword addgene
echo "drop table blat_feature_in_database" | mysql -u root -ppassword addgene
echo "drop table blat_feature_type" | mysql -u root -ppassword addgene
echo "drop table blat_sequence" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature" | mysql -u root -ppassword addgene
echo "drop table blat_sequence_feature_annotated" | mysql -u root -ppassword addgene


python manage.py syncdb

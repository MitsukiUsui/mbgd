target=${1}
./create_database.sh ${target}
./myphylophlan.sh ${target}
./query_lookup.py ${target}

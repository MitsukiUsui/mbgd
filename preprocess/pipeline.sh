target=${1}
#target="bifidobacterium_animalis"

./extract_strain.py ${target}
echo ""

./format_cluster.py ${target}
echo ""

# can be parallelized
./add_family_gene.py ${target}
echo ""

# can not be parallelized so far because of output collision
./distribute_family.py ${target}
echo ""

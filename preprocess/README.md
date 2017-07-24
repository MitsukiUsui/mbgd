## details

### flow
0. wget.sh
 * data acquisition
0. split.sh
 * divides `mbgd_2016-01_gene` by genomes into `/gene/${strain}.gene` files
0. add_family_cluster.py
 * adds unique "family" column to Ortholog Table
0. add_family_gene.py
 * adds "family" column to designated `.gene` files
0. distribute_family.py
 * distributes geneseq & proteinseq by its family


### util
* convert2gff.py
 * converts .gene to .gff file
* sample_cluster.py
 * samples Ortholog Cluster Table to subsets.

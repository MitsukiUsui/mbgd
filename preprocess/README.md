## details

### flow
* wget.sh
  * data acquisition
* split.sh
  * divides `mbgd_2016-01_gene` by genomes into `/gene/${strain}.gene` files
* add_family_cluster.py
  * adds unique "family" column to Ortholog Table
* add_family_gene.py
  * adds "family" column to designated `.gene` files
* distribute_family.py
  * distributes geneseq & proteinseq by its family


### util
* convert2gff.py
  * converts .gene to .gff file
* sample_cluster.py
  * samples Ortholog Cluster Table to subsets.

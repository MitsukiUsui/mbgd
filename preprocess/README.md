## details

### flow
#### only once
0. wget.sh
    * data acquisition
0. split.sh
    * divides `mbgd_2016-01_gene` by genomes into `/gene/${strain}.gene` files
0. add_family_cluster.py
    * adds unique "family" column to Ortholog Table
0. sample_strain.py
    * creates strain.lst 
0. add_family_gene.py
    * adds "family" column to designated `.gene` files
0. distribute_family.py
    * distributes geneseq & proteinseq by its family

#### every time with new strain.lst
0. format_cluster.py
    * add lineage, num_query column
    * convert to csv
    * sample rows if specified

### util
* convert2gff.py
    * converts .gene to .gff file

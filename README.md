# What's this?
Find frameshifted de novo gene.

## Workflow
Follow this flow from top to bottom. For further information, please refer to README on each sub directories.

### preprocess
* Data acquisition and preprocessing

### blastn
* blastn taxonomy restricted gene family to genomes which do not have the family.

### analize
* analize

## Data
MBGD Ortholog Cluster Table with corresponding complete genomes & their annotations.

### cluster.tab
|ClusterID|HomClusterID|Size|...genomes...|Gene|...params...|Description|
|:--:|:--:|:--|:--|:--|:--|:--|
|1|1|#genes|...gene ids...|family name (if exists)|...|description of family (if exists)|

### mbgd_2016-01_gene
annotation file. "from1" & "to1" columns corresponds to [first, last] position.



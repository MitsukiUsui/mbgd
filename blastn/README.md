## details
### flow

0. create_database.sh
0. ~~create_query.py~~
    * preprocess/distribute_family.pyが代用
0. exec.py
    * それぞれのfamilyを持たない株に対して、対応するqueryでblastn検索を行う。
    * UGEを使って並列処理するため、内部で`sge_blastn_args.sh`を呼び出している。
        * 引数として与えられたstrainとfamilyの組み合わせに対して、対応するdbとqueryを用いてblastn検索を行うjobをUGEにqsubする。
        * 内部で`blastn.sh`を呼び出している。
            * 引数として与えられたstrainとfamilyの組み合わせに対して、対応するdbとqueryを用いてblastn検索を行う。

### util
* summarize_blastn.py
    * output 2 csvs
        * count.csv: same data structure as cluster, each cell represent the number of hits
        * all.csv: mere concatenation of blast tabular output with their columns name
* convert2bed.py
    * output .bed of blastn hit

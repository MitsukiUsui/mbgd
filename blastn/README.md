## details

* create_database.sh
  * strain.lstに記載されている株に対してmakeblastdbし、ループでdbを構築する。
* extract_family.sh
  * それぞれのfamilyごとに対応する配列を抜き出し、一つのファイルにqueryとしてまとめる。
  * fastaファイルからの配列の抜き出しは`fatt extract`を呼び出すことで行なった。
* exec.py
  * それぞれのfamilyを持たない株に対して、対応するqueryでblastn検索を行う。
  * UGEを使って並列処理するため、内部で`sge_blastn_args.sh`を呼び出している。
* sge_blastn_args.sh
  * 引数として与えられたstrainとfamilyの組み合わせに対して、対応するdbとqueryを用いてblastn検索を行うjobをUGEにqsubする。
  * 内部で`blastn.sh`を呼び出している。
* blastn.sh
  * 引数として与えられたstrainとfamilyの組み合わせに対して、対応するdbとqueryを用いてblastn検索を行う。

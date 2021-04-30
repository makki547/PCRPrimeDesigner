# PCRPrimeDesigner
プラスミドからのサブクローニングなどで使用するPCRプライマーのデザイン支援スクリプトです。
Python3で動作します、またpandasが必要です。

# 使い方
## スクリプト以外に必要なもの
PCRで用いるDNAテンプレートのFASTAファイルが必要です。

## リポジトリ内のスクリプトの説明
* `make_primer_candidates.py` 入力したシーケンスファイルについて指定した位置から指定長内の塩基配列を取り出し、TmとGC含量、最後がG or Cとなっているかをリストアップします
* `get_sequence_position_from_unique_query.py` 入力したシーケンスファイルについて、ある塩基配列もしくはアミノ酸配列の位置を取得します
* `genomics_tools.py` 上記2つのスクリプトで使用する関数などをまとめたものです、通常、こちらを直接使用することはありません

## DNAテンプレート中のプライマーとして使う配列の位置を取得する
ここでは、Lysozyme.fastaをテンプレートとし、Lysozymeのアミノ酸配列F<sup>11</sup>LPLAA...AWVAWR<sup>130</sup>をサブクローニングするケースを例に説明します。
まず、DNA`get_sequence_position_from_unique_query.py`を使用して、該当する配列の位置を取得します。N末端側については、以下のコマンドで取得します;
```
$ python get_sequence_position_from_unique_query.py lysozyme.fasta FLPLA
```

C末端については以下のコマンドで取得します(-rオプションを指定することで、右端の位置を取得します);
```
$ python get_sequence_position_from_unique_query.py lysozyme.fasta WVAWR -r
```

上記の例では、いずれもアミノ酸配列を用いて検索していますが、塩基配列を用いても検索可能です。アミノ酸配列を指定した場合、翻訳フレームの3パターン全てに対して翻訳を実行・検索します。クエリーに対して翻訳フレームおよびその位置がユニークなときのみ位置を出力します。

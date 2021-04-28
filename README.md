# PCRPrimeDesigner
プラスミドからのサブクローニングなどで使用するPCRプライマーのデザイン支援スクリプトです。
Python3で動作します、またpandasが必要です。

# 使い方
## スクリプト以外に必要なもの
PCRで用いるDNAテンプレートのFASTAファイルが必要です。

## リポジトリ内のスクリプトの説明
* ``
* `get_sequence_position_from_unique_query.py`

## DNAテンプレート中のプライマーとして使う配列の位置を取得する
ここでは、Lysozyme.fastaをテンプレートとし、Lysozymeのアミノ酸配列F<sup>11</sup>LPLAA...AWVAWR<sup>130</sup>をサブクローニングするケースを例に説明します。
まず、DNA`get_sequence_position_from_unique_query.py`を使用します。
```
$ python get_sequence_position_from_unique_query.py lysozyme.fasta 
```

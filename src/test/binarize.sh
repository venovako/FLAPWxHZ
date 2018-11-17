#!/bin/bash
DIR=`dirname $0`
$DIR/binvec.exe $1 > $1.txt
sort -g $1.txt > $1.srt
unset DIR

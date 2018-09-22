#!/bin/bash
if [ -z "$1" ]
then
	echo "$0 FN"
else
	ln -fv $1.YY $1.Y
	ln -fv $1.WW $1.W
	ln -fv $1.JJ $1.J
fi

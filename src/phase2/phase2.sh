#!/bin/bash
if [ -z "$1" ]
then
    echo "$0 FN"
else
    mv -fv $1.YY $1.Y
    mv -fv $1.WW $1.W
    mv -fv $1.JJ $1.J
fi

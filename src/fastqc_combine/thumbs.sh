#!/bin/bash

for i in $1/*
do
	j=`basename $i .png`
	k=`dirname $i`
	/global/mb/amw/soft/ImageMagick-7.0.5-3/install/bin/convert -contrast -thumbnail 110 "$i" $k/thumb.$j.png
done

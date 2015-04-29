#!/bin/sh
if [ -z $1 ]; then
	echo "Need File Name"
else
	NAME=$1;
	if [ $(ls | grep $NAME.lis) ]; then
		rm $NAME.lis
	fi
	ls -d -1 $PWD/*.* | grep "\.root" >> $NAME.lis
fi



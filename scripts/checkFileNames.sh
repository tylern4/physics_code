#!/bin/sh
#
# while each loop to read through file.cut from v3
# 	while each loop to read through file.cut from v2
#		if name in file.cut.v2 == file.cut.v3 then echo file.cut.v3 >> same.txt
#		else file.cut.v3 >> unique_v3.txt
#	end while
# end while
#

V3_FILE_NAME=$1
V2_FILE_NAME=$2


while read LINE1
do
    V3_NAME=$LINE1
    while read LINE2
		do
    		V2_NAME=$LINE2

    		if [ "$V2_NAME" = "$V3_NAME" ]; then
    		  	echo "$V3_NAME:$V2_NAME" >> same.txt
    		fi
	done < $V2_FILE_NAME
done < $V3_FILE_NAME

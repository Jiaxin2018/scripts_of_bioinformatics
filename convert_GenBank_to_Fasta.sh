#!/bin/bash
# Usage: ./convert_GenBank_to_Fasta.sh FastaFile

if [ "$1" == "-h" -o "$1" == "--help" ]
then
        echo "Usage: ./convert_GenBank_to_Fasta.sh FastaFile"
        exit 0
elif [ $# -eq 0 ]
then
        seq_flag=0
        while read line
        do
                # find accession and save its value
                if echo $line |grep "^ACCESSION" >/dev/null
                then
                        acc=`echo "$line" |cut -d' ' -f4-`
			continue
                # find definition and save its value
                elif echo $line |grep "^DEFINITION" >/dev/null
                then
                        def=`echo "$line" |cut -d' ' -f3-`
			continue
                # determine the end of a sequence
                elif echo "$line" |grep "^\/\/" >/dev/null
                then
                        seq_flag=0
			continue
                # output the processed sequence
                elif [ $seq_flag -eq 1 ]
                then
                        echo "$line" |perl -pe "s/[ \d]//g" |tr a-z A-Z
			continue
                # determine the begining of a sequence and output the annotation line
                elif echo "$line" |grep "^ORIGIN" >/dev/null
                then
                        seq_flag=1
                        echo ">${acc} ${def}"
                fi
        done
else
        echo "Wrong arguments!"
        exit 1
fi
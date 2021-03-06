# !/bin/bash
#	DENSITY="300"
#	QUALITY="80"
#while getopts d:q: opts; do
#	case ${opts} in
#		d) DENSITY=${OPTARG} ;;
#		q) QUALITY=${OPTARG} ;;
#	esac
#done
col1="\$1"
col2="\$2"
for file in *.csv 
do 
	out_file="$(echo $file | sed 's/csv/dat/g')"
	sed 's/,/\t/g' $file > $out_file
done
rm *.csv

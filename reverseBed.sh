#!/bin/bash -l

while [ $# -gt 0 ]; do
case "$1" in
  pull)
	;;
  -i*)
	if [[ "$1" != *=* ]]; then shift; fi
	INPUT_FILE="${1#*=}"
	;;
  -c*)
    if [[ "$1" != *=* ]]; then shift; fi
	CENTER="${1#*=}"
	;;
  -h|--help)
	echo -e "NAME\n\t reverseBed.sh - Inverts the strand of a bed (or gtf) annotation. Hugo Guillen, 2020."
	echo -e "usage:\n\t reverseBed.sh [options]"
    echo -e "\t -i <bed/gtf> \n\t\t Input file."
    echo -e "\t -c <center> \n\t\t Center to reverse the bed. Omit to use the midpoint of the input annotation."
	exit 0
	;;
  *)
	>&2 echo "ERROR: Invalid argument. Please run with -h to see valid parameters."
	exit 1
	;;
esac
shift
done
if ! [ -f "${INPUT_FILE}" ]; then
  echo "# Input file not found. Exiting.";
  exit 1
fi

COLS=$(head -n 1 $INPUT_FILE | awk 'BEGIN{FS="\t"}{print NF}')

case "${INPUT_FILE##*.}" in 
  "bed")
    COL_X=2;COL_Y=3;
    if [[ $COLS -ge 6 ]]; then COL_S=6; else COL_S=0; fi
    ;;
  "gtf")
    COL_X=4;COL_Y=5;COL_S=7;
    ;;
esac

if [ -z "$CENTER" ]; then
  CENTER=$(bedtools groupby -i $INPUT_FILE -g 1 -c $COL_X,$COL_Y -o min,max |\
              awk 'BEGIN{FS=OFS="\t"}END{print ($2+$3)/2}')
fi
awk -v X="$COL_X" -v Y="$COL_Y" -v S="$COL_S" -v C="$CENTER" '
  BEGIN{FS=OFS="\t"}
  {x=$X; y=$Y; $X=2*C-y; $Y=2*C-x; if(S>0){($S=="+")?$S="-":$S="+"}; print $0}
' $INPUT_FILE | sort -k1,1 -k"$COL_X","$COL_X"n


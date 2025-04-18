#!/bin/zsh

## $# = the number of arguments
while (( $# )); do
	key="$1"
	case $key in
##		pattern)
##			sentence
##		;;
		-s | --subj)
			subj="$2"
		;;
	esac
	shift ##takes one argument
done
nn=`printf "%02d" $subj`

dir_work="/home/ROBARTS/skim2764/tsclient/skim2764/BIDS"
dir_data="/srv/diedrichsen/data/SeqSpatialSupp_fMRI/BIDS"

parallel -j1 mkdir -p -m 755 "$dir_work/sub-{2}{1}/{3}" ::: "$nn" ::: 'S' 'R' ::: 'anat' 'fmap' 'func'
parallel -j1 cp -ru "$dir_data/sub-{2}{1}/{3}/*" "$dir_work/sub-{2}{1}/{3}/" ::: $nn ::: 'S' 'R' ::: 'anat' 'fmap' 'func'

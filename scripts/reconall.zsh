#!/bin/zsh

## ======================================= ##

## $# = the number of arguments
while (( $# )); do
	key="$1"
	case $key in
##		pattern)
##			sentence
##		;;
		-s | --subject)
			subj="$2"
		;;
	esac
	shift ##takes one argument
done

## ======================================= ##
dir_root="/mnt/f/SeqSpatialSupp_fMRI"
dir_raw="$dir_root/anatomicals/$subj"
dir_FreeSurfer="$dir_root/freesurf"
## ======================================= ##

## I/O path, same as above, following earlier steps
if [[ ! -d $dir_FreeSurfer ]]; then
	mkdir -p -m 755 $dir_FreeSurfer
fi

if [[ ! -d $dir_FreeSurfer/$subj ]]; then
	## FS function
	recon-all \
		-sid $subj \
		-sd $dir_FreeSurfer \
		-i $dir_raw/${subj}_anatomical.nii \
		-all -cw256
fi

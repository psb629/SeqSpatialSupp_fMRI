#!/bin/zsh

##############################################################
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
##############################################################
dir_root='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
dir_anat=$dir_root/anatomicals
##############################################################
dir_FreeSurfer="$dir_root/FreeSurfer"
if [[ ! -d $dir_FreeSurfer ]]; then
	mkdir -p -m 755 $dir_FreeSurfer
fi

## FS function
recon-all \
	-sid	$subj \
	-sd		$dir_FreeSurfer \
	-i		$dir_anat/$subj/${subj}_anatomical.nii \
	-all
	

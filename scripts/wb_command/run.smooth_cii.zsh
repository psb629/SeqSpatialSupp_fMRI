#!/bin/zsh

##############################################################
gs='n'
axis='COLUMN'
## $# = the number of arguments
while (( $# )); do
	key="$1"
	case $key in \
##		pattern)
##			sentence
##		;;
		-i | --input)
			input="$2"
		;;
		-a | --axis)
			axis="$2"
		;;
		-h | --help)
			cat <<EOF
===========================================
Options:
	-i | --input
	-a | --axis (default: column)
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
##############################################################
axis=${axis:u}
fname=${input:t}
dname=${input:h}
##############################################################
case "$(uname)" in
	Darwin)
		dir_root='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
	;;
	Linux)
		dir_root='/mnt/f/SeqSpatialSupp_fMRI'
	;;
esac
dir_fs=$HOME/github/fs_LR_32
##############################################################
wb_command \
	-cifti-smoothing $input \
	2 2 $axis \
	$dname/smooth.$fname \
	-left-surface $dir_fs/fs_LR.32k.L.midthickness.surf.gii \
	-right-surface $dir_fs/fs_LR.32k.R.midthickness.surf.gii


 #wb_command \
 #	-cifti-smoothing bold.dtseries.nii 7 7 COLUMN bold_roi_smooth.dtseries.nii \
 #    -fwhm \
 #    -left-surface subject.L.midthickness.surf.gii \
 #    -right-surface subject.R.midthickness.surf.gii \
 #    -cifti-roi roi_mask.dscalar.nii


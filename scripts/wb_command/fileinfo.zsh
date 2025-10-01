#!/bin/zsh

##############################################################
gs='n'
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
		-h | --help)
			cat <<EOF
===========================================
Options:
	-i | --input
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
##############################################################
fname=${input:t}
dname=${input:h}
prefix=${fname:r}
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
dir_surf=$dir_root/surfaceWB
##############################################################
wb_command -file-information $input > $dir_root/ROI/$prefix.txt

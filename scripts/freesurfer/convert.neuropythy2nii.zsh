#!/bin/zsh

##############################################################
## $# = the number of arguments
while (( $# )); do
	key="$1"
	case $key in \
		-s | --subject)
			subj="$2"
		;;
		-h | --hemisphere)
			hem="$2"
		;;
		--help)
			cat <<EOF
===========================================
Options:
	-s | --subject
	-h | --hemisphere
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
hem=${hem:l}
##############################################################
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
dir_fs=$dir_root/FreeSurfer
dir_work=$dir_fs/$subj
##############################################################
export SUBJECTS_DIR=$dir_fs
##############################################################
mri_surf2vol \
	--o $dir_work/mri/${hem}h.benson14_varea.vol.nii \
	--subject $subj \
	--so $dir_work/surf/${hem}h.white \
	$dir_work/surf/${hem}h.benson14_varea.mgz

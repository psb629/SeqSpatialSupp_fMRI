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
		-s | --subject)
			subj="$2"
		;;
		-h | --hem)
			hem="$2"
		;;
		-gs | --glmsingle)
			gs="$2"
		;;
		-g | --glm)
			glm="$2"
		;;
		-m | --map)
			map="$2"
		;;
		--help)
			cat <<EOF
===========================================
Options:
	-s | --subject
	-h | --hem
	-g | --glm
	-m | --map
	-gs | --glmsingle (default: n)
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
glm1d=`printf "%1d\n" $glm`
hem=${hem:u}
##############################################################
dir_root='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
dir_fs=$dir_root/FreeSurfer
dir_fs=$HOME/github/fs_LR_32
##############################################################
case $gs in \
	'yes' | 'y')
		dir_work=$dir_root/GLMsingle/glm_${glm1d}/surfaceWB/$subj
	;;
	'no' | 'n')
		dir_work=$dir_root/surfaceWB/glm_${glm1d}/$subj
	;;
	*)
		dir_work=$dir_root/glm_${glm1d}
	;;
esac
##############################################################
wb_command \
	-cifti-smoothing $dir_work/$subj.$hem.glm_${glm1d}.$map.dscalar.nii \
	2 2 COLUMN \
	$dir_work/$subj.$hem.glm_${glm1d}.${map}_smooth.dscalar.nii \
    -left-surface $dir_fs/fs_LR.32k.L.midthickness.surf.gii \
    -right-surface $dir_fs/fs_LR.32k.R.midthickness.surf.gii


 #wb_command \
 #	-cifti-smoothing bold.dtseries.nii 7 7 COLUMN bold_roi_smooth.dtseries.nii \
 #    -fwhm \
 #    -left-surface subject.L.midthickness.surf.gii \
 #    -right-surface subject.R.midthickness.surf.gii \
 #    -cifti-roi roi_mask.dscalar.nii


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
		--glmsingle)
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
	-g | -glm
	-m | --map
	--glmsingle
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
glm1d=`printf "%1d\n" $gs`
hem=${hem:l}
##############################################################
dir_root='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
dir_surf=$dir_root/surfaceWB
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
	-cifti-smoothing bold.dtseries.nii 2 2 COLUMN bold_smooth.dtseries.nii \
    -left-surface subject.L.midthickness.surf.gii \
    -right-surface subject.R.midthickness.surf.gii

wb_command \
	-cifti-smoothing bold.dtseries.nii 7 7 COLUMN bold_roi_smooth.dtseries.nii \
    -fwhm \
    -left-surface subject.L.midthickness.surf.gii \
    -right-surface subject.R.midthickness.surf.gii \
    -cifti-roi roi_mask.dscalar.nii

case $map in \
	'beta')
		wb_command \
			-cifti-create-dense-timeseries $dir_work/$subj.$hem.glm_${glm1d}.$map.dtseries.nii \
	    	-left-metric $dir_work/$subj.$hem.glm_${glm1d}.$map.func.gii \
		    -timestep $TR -timestart $onset
	;;	
	*)
		wb_command \
			-cifti-create-dense-scalar $dir_work/$subj.$hem.glm_${glm1d}.$map.dscalar.nii \
	    	-left-metric $dir_work/$subj.$hem.glm_${glm1d}.$map.func.gii
	;;
esac


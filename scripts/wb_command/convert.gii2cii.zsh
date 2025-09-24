#!/bin/zsh

##############################################################
gs='n'
type_='dscalar'
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
		-s | --subject)
			subj="$2"
		;;
		-t | --type)
			type_="$2"
		;;
		-h | --help)
			cat <<EOF
===========================================
Options:
	-i | --input
	-s | --subject
	-t | --type
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
hem=${hem:u}
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
dir_fs=$HOME/github/fs_LR_32
dir_surf=$dir_root/surfaceWB
##############################################################
TR=1
onset=3
##############################################################
case $type_ in \
	'dseries')
		wb_command \
			-cifti-create-dense-timeseries $input \
	    	-left-metric $dir_surf/$subj/$subj.L.glm_${glm1d}.$map.func.gii \
	    	-right-metric $dir_surf/$subj/$subj.R.glm_${glm1d}.$map.func.gii \
		    -timestep $TR -timestart $onset
	;;	
	'dscalar')
		wb_command \
			-cifti-create-dense-scalar $input \
	    	-left-metric $dir_surf/$subj/$subj.L.glm_${glm1d}.$map.func.gii \
	    	-right-metric $dir_surf/$subj/$subj.R.glm_${glm1d}.$map.func.gii \
	;;
esac
echo " Created $fname"

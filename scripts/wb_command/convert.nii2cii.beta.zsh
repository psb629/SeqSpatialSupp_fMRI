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
		-s | --subject)
			subj="$2"
		;;
		-h | --help)
			cat <<EOF
===========================================
Options:
	-i | --input
	-s | --subject
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
L_MID=$dir_fs/fs_LR.32k.L.midthickness.surf.gii
R_MID=$dir_fs/fs_LR.32k.R.midthickness.surf.gii
L_WHITE=$dir_root/surfaceWB/$subj/$subj.L.white.32k.surf.gii
L_PIAL=$dir_root/surfaceWB/$subj/$subj.L.pial.32k.surf.gii
R_WHITE=$dir_root/surfaceWB/$subj/$subj.R.white.32k.surf.gii
R_PIAL=$dir_root/surfaceWB/$subj/$subj.R.pial.32k.surf.gii
##############################################################
## Transform: NIFTI -> surface metric
echo " # $subj: NIFTI -> surface metric (GIFTI)"
dir_tmp=$dname/tmp_metric
if [ ! -d $dir_tmp ]; then
	mkdir -p -m 755 $dir_tmp

	cd $dname
	for nii in beta_*.nii; do
		base=${nii%.nii}
		# 좌반구
	    wb_command -volume-to-surface-mapping \
	        "$nii" \
	        "$L_MID" \
	        "$dir_tmp/$base.L.func.gii" \
	        -ribbon-constrained "$L_WHITE" "$L_PIAL"
	
	    # 우반구
	    wb_command -volume-to-surface-mapping \
	        "$nii" \
	        "$R_MID" \
	        "$dir_tmp/$base.R.func.gii" \
	        -ribbon-constrained "$R_WHITE" "$R_PIAL"	
	done
	echo " # Done!"
fi
##############################################################
## Merge
echo " # Merging the GIFTIs"
cd $dir_tmp
wb_command -metric-merge betas.L.func.gii \
	$(for f in beta_*.L.func.gii; do echo "-metric $f"; done)
wb_command -metric-merge betas.R.func.gii \
	$(for f in beta_*.R.func.gii; do echo "-metric $f"; done)
##############################################################
## Converting: GITFI -> CIFTI
echo " # Creating CIFTI"
wb_command -cifti-create-dense-scalar $dname/betas.L.dscalar.nii \
	-left-metric $dir_tmp/betas.L.func.gii \

wb_command -cifti-create-dense-scalar $dname/betas.R.dscalar.nii \
	-right-metric $dir_tmp/betas.R.func.gii \
##############################################################
rm -r $dir_tmp

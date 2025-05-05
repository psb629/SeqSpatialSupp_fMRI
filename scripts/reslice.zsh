#!/bin/zsh

## ======================================= ##
align_surf=(1 1 1)
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
		-a | --align_surf)
			align_surf="$2"
		;;
	esac
	shift ##takes one argument
done

## ======================================= ##
dir_root="/mnt/f/SeqSpatialSupp_fMRI"
dir_surf="$dir_root/freesurf/$subj"
dir_wb="$dir_root/surfaceWB/$subj"
dir_atlas="/home/sungbeenpark/github/surfAnalysis/standard_mesh"
## ======================================= ##
hemisphere=(1 2)
hem=('lh' 'rh');
Hem=('L' 'R'); 
surf_files=('white' 'pial' 'inflated')
curv_files=('curv' 'sulc' 'area') 
resolution=32
anafile="$dir_surf/mri/brain.mgz"
## ======================================= ##
dir_output="$dir_wb"
if [[ ! -d $dir_output ]]; then
	mkdir -p -m 755 $dir_output
fi

## Transform of voxels in 256x256 image to surface vertices
Mvox2surf=($(mri_info --vox2ras-tkr "$anafile" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'))
Mvox2surf="[${(j:, :)Mvox2surf}]"

## Transform of voxel to subject space
Mvox2space=($(mri_info --vox2ras "$anafile" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'))
Mvox2space="[${(j:, :)Mvox2space}]"

## calculate the transformation matrix
Msurf2space=($(python -c "
import numpy as np
Mvox2surf = np.array($Mvox2surf).reshape(4,4)
 #print(Mvox2surf)
Mvox2space = np.array($Mvox2space).reshape(4,4)
 #print(Mvox2space)
Msurf2space = Mvox2space @ np.linalg.inv(Mvox2surf)
print(' '.join(map(str, Msurf2space.flatten())))
"))
Msurf2space="[${(j:, :)Msurf2space}]"
## ======================================= ##
## process hemispheres
for oo in $hemisphere
{
	H=$Hem[$oo]
	h=$hem[$oo]

	# Convert registration sphere
	reg_sphere="$dir_surf/surf/${h}.sphere.reg.surf.gii"
	if [[ ! -f $reg_sphere ]]; then
		mris_convert $dir_surf/surf/${h}.sphere.reg $reg_sphere
	fi
	
	atlas="$dir_atlas/resample_fsaverage/fs_LR-deformed_to-fsaverage.$H.sphere.${resolution}k_fs_LR.surf.gii"

	# Process surface files
	cnt=1
	for ss in $surf_files
	{
		fname="$dir_surf/surf/$h.$ss.surf.gii"
		output="$dir_output/$subj.$H.$ss.${resolution}k.surf.gii"
 #		if [[ -f $atlas ]]; then
 #			echo "true"
 #		else
 #			echo "false"
 #		fi

		mris_convert $dir_surf/surf/$h.$ss $fname
		wb_command -surface-resample $fname $reg_sphere $atlas BARYCENTRIC $output

		if [[ $align_surf[$cnt] -eq 1 ]]; then
			# Apply affine transform using Python script
			python3 -c "
from os.path import join
import numpy as np
import nibabel as nb

fname = join('/mnt/f/SeqSpatialSupp_fMRI/surfaceWB/$subj/$subj.$H.$ss.${resolution}k.surf.gii')
surf = nb.load(fname)
darray = surf.darrays[0]
M = np.array($Msurf2space).reshape(4,4)
 #print(M)
verts = np.hstack(
    [darray.data, np.ones((darray.data.shape[0],1))]
) @ M.T
darray.data = verts[:,:3]
nb.save(surf,fname)
			"
		fi
		cnt=$((cnt + 1))
	}

	# Process curvature files
	for ii in {1..3}
	{
		cc=$curv_files[$ii]
		ss=$surf_files[$ii]
		ccgii="$dir_surf/surf/$h.$cc.shape.gii"
		ssgii="$dir_surf/surf/$h.$ss.surf.gii"
		output="$dir_output/$subj.$H.$cc.${resolution}k.shape.gii"

		mris_convert -c $dir_surf/surf/$h.$cc $dir_surf/surf/$h.$ss $fname
		wb_command -metric-resample $fname $reg_sphere $atlas BARYCENTRIC $output
	}	
}

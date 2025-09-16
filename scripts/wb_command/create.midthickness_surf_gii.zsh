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
		-h | --help)
			cat <<EOF
===========================================
Options:
	-s | --subject
===========================================
EOF
		exit
		;;
	esac
	shift ##takes one argument
done
hem=${hem:l}
##############################################################
dir_root='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
dir_fs=$dir_root/FreeSurfer
##############################################################
# 경로 설정
dir_work="$dir_fs/$subj"
dir_surf="$dir_work/surf"
dir_tmp="$dir_work/tmp_midthickness"
mkdir -p $dir_tmp

# 좌/우반구 처리 루프
for HEMI in lh rh; do
	# 중간 표면 생성 (좌표 평균)
	mris_avg_surfaces \
		"$dir_tmp/$HEMI.graymid" \
		"$dir_surf/$HEMI.white" \
		"$dir_surf/$HEMI.pial" \
		0.5

	# GIFTI 형식으로 변환
	mris_convert \
		"$dir_tmp/$HEMI.graymid" \
		"$dir_tmp/$HEMI.midthickness.surf.gii"
done

echo "✅ Midthickness surfaces created:"
echo "  $dir_tmp/lh.midthickness.surf.gii"
echo "  $dir_tmp/rh.midthickness.surf.gii"
rm -rf $dir_tmp

#!/bin/zsh

## $# = the number of arguments
while (( $# )); do
	key="$1"
	case $key in
##		pattern)
##			sentence
##		;;
		-s | --subject)
			## string
			subj="$2"
		;;
	esac
	shift ##takes one argument
done
## ======================================= ##
dir_work="/Users/sungbeenpark/Downloads/sub-$subj"
## ======================================= ##
fname="$dir_work/sslb2_$subj.dat"
## ======================================= ##
sed -i '' 's/\.000//g' $fname

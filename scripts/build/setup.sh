#!/bin/bash

# this script sets up $LD_LIBRARY_PATH and .so links for all external libs;
# this allows to run BCL from a checkout outside of the Meilerlab.
#
# heinzes1, 2013-2014

if [[ $0 != ${SHELL} ]]; then echo "Please source script with: source $0"; exit; fi

# determine basic parameters
git_command="git"
svn_command="svn"
command -v git >/dev/null 2>&1 || git_command=""
command -v svn >/dev/null 2>&1 || svn_command=""

git_bcl_root=""
svn_bcl_root=""
if [[ -n ${git_command} ]]; then git_bcl_root=`git rev-parse --show-toplevel 2>/dev/null`/; fi
if [[ -n ${svn_command} ]]; then svn_bcl_root=`svn info | grep "Working Copy Root Path: " | sed s'/Working Copy Root Path: //'`; fi

if [[ -n "${git_bcl_root}" ]]; then
	bcl_root="${git_bcl_root}"
elif [[ -n "${svn_bcl_root}" ]]; then
	bcl_root="${svn_bcl_root}"
else
	echo "Could not determine svn or git workspace root. Exiting."
	exit
fi

arch="x86_64"
lib_prefix="${bcl_root}extern/Linux2/${arch}/"

# create dir list
dir_list=( $(ls -1 ${lib_prefix}) )
version_list=()
list_size=${#dir_list[@]}
index=0
while [ "$index" -lt "$list_size" ]; do
	dir="${lib_prefix}${dir_list[$index]}/"
	this_dir_version_list=( $(ls -v ${dir}) ) # natural version number sort
	last_version=${this_dir_version_list[${#this_dir_version_list[@]} - 1]}
	version_list+=( ${last_version} )
	let "index=$index + 1"
done

# alternatively set lists manually
#dir_list=    (ati  bzip2  freeocl  mysql   mysqlpp  zlib )
#version_list=(2.5  1.0.6  0.3.6    5.6.14  3.2.0    1.2.8)
#list_size=${#version_list[@]}

# assemble ld_library_path value and add if not already added
bcl_ld_library_path=""

index=0
while [ "$index" -lt "$list_size" ]; do
	if [[ -z "${bcl_ld_library_path}" ]]; then
		bcl_ld_library_path="${lib_prefix}${dir_list[$index]}/${version_list[$index]}/lib/"
	else
		bcl_ld_library_path="${lib_prefix}${dir_list[$index]}/${version_list[$index]}/lib/:${bcl_ld_library_path}"
	fi
	let "index=$index + 1"
done

case "${LD_LIBRARY_PATH}" in
	*${bcl_ld_library_path}*) echo "LD_LIBRARY_PATH already set. Skipping."   ;;
	*)                        echo "Setting LD_LIBRARY_PATH."
					if [[ -z "${LD_LIBRARY_PATH}" ]]; then
															export LD_LIBRARY_PATH="${bcl_ld_library_path}"
														else
															export LD_LIBRARY_PATH="${bcl_ld_library_path}:${LD_LIBRARY_PATH}"
														fi;;
esac

# create .so links
index=0
while [ "$index" -lt "$list_size" ]; do
	dir="${lib_prefix}${dir_list[$index]}/${version_list[$index]}/lib/"
	for file in ${dir}*.so; do
		 [ ! -e ${file} ] && continue # if file does not exist, skip

		 objdump_command="objdump -x ${file}"
		 soname=`${objdump_command}|grep SONAME|awk '{print $2}'`

		 this_dir=`dirname ${file}`
		 this_file=`basename ${file}`

		 # if soname == this_file, skip
		 if [ "${this_file}" == "${soname}" ]; then
			 echo "Skipping to create link $soname -> $this_file in $this_dir; names are identical."
			 continue
		 fi

		 echo "Creating link $soname -> $this_file in $this_dir"
		 cd $this_dir
		 ln -sf "$this_file" "$soname"
		 cd - > /dev/null
	done
	let "index=$index + 1"
done

#objdump -x extern/Linux2/x86_64/ati/2.5/lib/libOpenCL.so |grep SONAME|awk '{print $2}'

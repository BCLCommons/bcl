#!/bin/sh

# source the common shell script
source `dirname $0`/common.sh

####################
# build for linux 64
####################
cd ${BCL_ROOT}

# fix line endings for histograms and opencl kernels
# This may not be strictly necessary for proper functioning of the bcl, but if users browse their histogram or opencl
# folder, it is best if the files are formatted properly for their system.  When getline is used read something in but 
# is not in a while loop, it could cause incorrect behavior if files from different systems are read in. 
dos2unix ${BCL_ROOT}/histogram/* ${BCL_ROOT}/opencl_kernels/* ${BCL_ROOT}/model/*/*.model ${BCL_ROOT}/model/*/*.info ${BCL_ROOT}/model/*/*.descriptor ${BCL_ROOT}/model/*/*.result 2> /dev/null
# change the distcc compiler wrappers for Linux to the last version of gcc compiled on CentOS5, to make the binaries
# CentOS5 compatible
BUILD_PATH=build/linux64_release/
rm -rf $BUILD_PATH
mkdir -p $BUILD_PATH
cd ${BUILD_PATH}
cmake -D ADDITIONAL_COMPILER_OPTIONS="-D__USE_XOPEN2K8" \
      -D CMAKE_TOOLCHAIN_FILE=../../cmake/toolchain/x86_64-unknown-linux-gnu.cmake \
      -D BCL_LICENSE=OFF \
      -D CMAKE_BUILD_TYPE=Release \
      -D BCL_EXTERNALS="pthreads;bzip2;zlib;freeocl" \
      ../../
pump make -j120 
find . -name '*.sh' -delete
rm -rf ./_CPack_Packages/

# note on flags:
# CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE just means that we want a single shell scripts for each component
# CPACK_COMPONENTS_ALL            tells cpack what to install, provided CPACK_ARCHIVE_COMPONENT_INSTALL is set to 1
# CPACK_INSTALL_CMAKE_PROJECTS    takes 4 values, the build path, project name, component, root directory
# Note that the setting of CPACK_COMPONENTS_ALL is unused unless the CPACK_INSTALL_CMAKE_PROJECTS component value is set 
# to something other than ALL 
cpack -G "STGZ;TGZ" \
      -D CPACK_COMPONENTS_ALL="${DEFAULT_COMPONENTS};pthreads;freeocl;freeocl_LICENSE;BclOpenCLKernels" \
      -D CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE=ON \
      -D CPACK_INSTALL_CMAKE_PROJECTS="${BCL_ROOT}/${BUILD_PATH};bcl-project;BclReleaseAll;/" --verbose

make package_source

file_content="Load_File('/net/$HOST${BCL_ROOT}/${BUILD_PATH}bcl-${BCL_VERSION}-Linux-x86_64.sh')"
srcfile_content="Load_File('/net/$HOST${BCL_ROOT}/${BUILD_PATH}bcl-${BCL_VERSION}-Source.tar.gz')"
vmajor=`echo $BCL_VERSION | awk -F. '{print $1}'`
vminor=`echo $BCL_VERSION | awk -F. '{print $2}'`
vpatch=`echo $BCL_VERSION | awk -F. '{print $3}'`
echo "insert into bcl_downloads(os_arch,v_major,v_minor,v_patch,compile_date,file,filename) values('Linux x86_64 (64-bit)',$vmajor,$vminor,$vpatch,CURRENT_DATE(),$file_content,'bcl-${BCL_VERSION}-Linux-x86_64.sh') ON DUPLICATE KEY UPDATE compile_date=VALUES(compile_date), file=VALUES(file);" > ./mysql.cmd.txt
echo "insert into bcl_downloads(os_arch,v_major,v_minor,v_patch,compile_date,file,filename) values('Source tar.gz',$vmajor,$vminor,$vpatch,CURRENT_DATE(),$srcfile_content,'bcl-${BCL_VERSION}-Source.tar.gz') ON DUPLICATE KEY UPDATE compile_date=VALUES(compile_date), file=VALUES(file);" >> ./mysql.cmd.txt
chmod a+wx *.sh *.tar.gz mysql.cmd.txt
ssh carbon mysql --defaults-file=$HOME/.my.conf bclweb_sym < /net/$HOST/${BCL_ROOT}/${BUILD_PATH}/mysql.cmd.txt
chmod go-wx *.sh *.tar.gz mysql.cmd.txt

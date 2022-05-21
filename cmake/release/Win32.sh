#!/bin/sh

# source the common shell script
source `dirname $0`/common.sh

##################
# build for win 32
##################

cd ${BCL_ROOT}
BUILD_PATH=build/win32_release/

# need to change line endings on all text files that go into the release
# fix line endings for histograms and opencl kernels and models
# This may not be strictly necessary for proper functioning of the bcl, but if users browse their histogram or opencl
# folder, it is best if the files are formatted properly for their system.  When getline is used read something in but 
# is not in a while loop, it could cause incorrect behavior if files from different systems are read in. 
find ${BCL_ROOT} ${BCL_ROOT}/scripts/build/* \( \( -not -name build -not -name .svn -not -name .git -not -name extern \) -o -prune \) -type f -exec file \{\} \; | egrep '(ASCII|C\+\+| C |text)' | tr -d ':' | awk '{print $1}' | xargs -n10 unix2dos
find ${BCL_ROOT} ${BCL_ROOT}/scripts/build/*  \( \( -not -name build -not -name .svn -not -name .git -not -name extern -not -name example_files \) -o -prune \) -type f -exec file \{\} \; | egrep '( compressed )' | tr -d ':' | awk '{print $1}' > /tmp/compressed_files.txt
grep '.gz$' /tmp/compressed_files.txt | xargs -n1 gunzip
grep '.gz$' /tmp/compressed_files.txt | sed 's/\.gz$//' | xargs -n1 unix2dos
grep '.gz$' /tmp/compressed_files.txt | sed 's/\.gz$//' | xargs -n1 gzip
grep '.bz2$' /tmp/compressed_files.txt | xargs -n1 bunzip2
grep '.bz2$' /tmp/compressed_files.txt | sed 's/\.bz2$//' | xargs -n1 unix2dos
grep '.bz2$' /tmp/compressed_files.txt | sed 's/\.bz2$//' | xargs -n1 bzip2

rm -rf $BUILD_PATH
mkdir -p $BUILD_PATH
cd ${BUILD_PATH}
cmake -D CMAKE_TOOLCHAIN_FILE=../../cmake/toolchain/i686-w64-mingw32.cmake \
      -D BCL_LICENSE=OFF \
      -D CMAKE_BUILD_TYPE=Release \
      -D BCL_EXTERNALS="bzip2;zlib" \
      ../../
pump make -j120 

# see notes about flags in Linux64.sh
cpack -G "NSIS;ZIP" \
      -D CPACK_INSTALL_CMAKE_PROJECTS="${BCL_ROOT}/${BUILD_PATH};bcl-project;BclReleaseAll;/" \
      -D CPACK_COMPONENTS_ALL="${DEFAULT_COMPONENTS}" \
      -D CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE=1 \
      -D CPACK_MONOLITHIC_INSTALL=1 

# revert line endings from dos
find ${BCL_ROOT} ${BCL_ROOT}/scripts/build/*  \( \( -not -name build -not -name .svn -not -name .git -not -name extern \) -o -prune \) -type f -exec file \{\} \; | egrep '(ASCII|C\+\+| C |text)' | tr -d ':' | awk '{print $1}' | xargs -n10 dos2unix
grep '.gz$' /tmp/compressed_files.txt | xargs -n1 gunzip
grep '.gz$' /tmp/compressed_files.txt | sed 's/\.gz$//' | xargs -n1 dos2unix
grep '.gz$' /tmp/compressed_files.txt | sed 's/\.gz$//' | xargs -n1 gzip -n 
grep '.bz2$' /tmp/compressed_files.txt | xargs -n1 bunzip2
grep '.bz2$' /tmp/compressed_files.txt | sed 's/\.bz2$//' | xargs -n1 dos2unix
grep '.bz2$' /tmp/compressed_files.txt | sed 's/\.bz2$//' | xargs -n1 bzip2
svn revert -R model histogram

file_content="Load_File('/net/$HOST${BCL_ROOT}/${BUILD_PATH}bcl-${BCL_VERSION}-Windows-x86.exe')"
vmajor=`echo $BCL_VERSION | awk -F. '{print $1}'`
vminor=`echo $BCL_VERSION | awk -F. '{print $2}'`
vpatch=`echo $BCL_VERSION | awk -F. '{print $3}'`
echo "insert into bcl_downloads(os_arch,v_major,v_minor,v_patch,compile_date,file,filename) values('Windows x86 (32-bit)',$vmajor,$vminor,$vpatch,CURRENT_DATE(),$file_content,'bcl-${BCL_VERSION}-Windows-x86.exe') ON DUPLICATE KEY UPDATE compile_date=VALUES(compile_date), file=VALUES(file);" > ./mysql.cmd.txt
chmod a+wx *.exe mysql.cmd.txt
ssh carbon mysql --defaults-file=$HOME/.my.conf bclweb_sym < /net/$HOST/${BCL_ROOT}/${BUILD_PATH}/mysql.cmd.txt
chmod go-wx *.exe mysql.cmd.txt

cd ${BCL_ROOT}

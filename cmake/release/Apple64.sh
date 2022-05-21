#!/bin/sh

# source the common shell script
source `dirname $0`/common.sh

# fix line endings for histograms and opencl kernels.  
# This may not be strictly necessary for proper functioning of the bcl, but if users browse their histogram or opencl
# folder, it is best if the files are formatted properly for their system.  When getline is used read something in but 
# is not in a while loop, it could cause incorrect behavior if files from different systems are read in.  
dos2unix ${BCL_ROOT}/histogram/* ${BCL_ROOT}/opencl_kernels/* ${BCL_ROOT}/model/*/* 2> /dev/null

#################
# build for apple
#################
cd ${BCL_ROOT}
OPERATING_SYSTEM=Darwin
BUILD_PATH=build/apple64_release/
rm -rf $BUILD_PATH
mkdir -p $BUILD_PATH
cd ${BUILD_PATH}
cmake -D CMAKE_TOOLCHAIN_FILE=../../cmake/toolchain/x86_64-apple-darwin13.cmake \
      -D BCL_LICENSE=OFF \
      -D CMAKE_BUILD_TYPE=Release \
      -D BCL_EXTERNALS="pthreads;bzip2;zlib;opencl" ../../
pump make -j120

find . -name '*.sh' -delete
rm -rf ./_CPack_Packages/ dmg *.jpg

# create a normal .sh script, used only to collect the files in this case
# it is unnecessary to include pthreads since the executable will link against the system library anyway
# see notes about flags in Linux64.sh
cpack \
  -G "STGZ" \
  -D CPACK_INSTALL_CMAKE_PROJECTS="${BCL_ROOT}/${BUILD_PATH};bcl-project;BclReleaseAll;/" \
  -D CPACK_COMPONENTS_ALL="${DEFAULT_COMPONENTS};pthreads;BclOpenCLKernels;" \
  -D CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE=1
  
# call the installer; hit y to both prompts; note, when running this script, it is still necessary to hit space a couple
# times to scroll through the license
/bin/echo -e "y\ny\n" | ./bcl-$BCL_VERSION-Darwin-x86_64.sh

# long term goal: have cmake generate a bundle on mac os systems
# It appears the only way to do this is to run cpack and cmake and make from a mac os machine
# but then we end up having to compile and link from the mac, which is both slow and inconsistent
# with the cross compiler and linker used on linux
# Short term: make a dmg for Mac with the proper image, also change executable names 
mv ./bcl-$BCL_VERSION-Darwin-x86_64/*.exe ./bcl-$BCL_VERSION-Darwin-x86_64/bcl
mkdir dmg
mv bcl-$BCL_VERSION-Darwin-x86_64 dmg/BioChemicalLibrary-$BCL_VERSION

# copy over image and write a script to make the dmg
cp $BCL_ROOT/documentation/image/BCL_m.jpg ./
echo "#!/bin/sh" > ./make_dmg_script.sh
chmod a+x ./make_dmg_script.sh
echo "rm -f bcl-$BCL_VERSION.dmg dmg/Applications" >> ./make_dmg_script.sh
echo "ln -s /Applications dmg/" >> ./make_dmg_script.sh
echo "/dors/meilerlab/apps/Darwin8/noarch/shell_scripts/create-dmg/create-dmg --window-size 500 255 --background BCL_m.jpg --icon-size 96 --volname BCL\\ \\-\\ The\\ BioChemical\\ Library --icon Applications 333 128 --icon BioChemicalLibrary-$BCL_VERSION 166 128 bcl-$BCL_VERSION.dmg dmg" >> ./make_dmg_script.sh

#rsync everything from the local build directory over to blue, since the local drive will not be accessible from the mac
LOCAL_PWD=`echo "$PWD" | sed 's:^/hd0/:/home/:' | sed 's:/ssd./:/home/:'`
mkdir -p $LOCAL_PWD
rsync -a --delete --exclude "**/.svn/" --exclude "**/.git/" --exclude "**/.metadata/" $PWD/ $LOCAL_PWD/
echo "Waiting a minute for the file system to notice the changes on titanium"
sleep 62
ssh titanium "cd $LOCAL_PWD; ./make_dmg_script.sh"

cp $LOCAL_PWD/bcl-$BCL_VERSION.dmg $PWD/

make package_source

file_content="Load_File('/net/$HOST${BCL_ROOT}/${BUILD_PATH}bcl-${BCL_VERSION}.dmg')"
srcfile_content="Load_File('/net/$HOST${BCL_ROOT}/${BUILD_PATH}bcl-${BCL_VERSION}-Source.zip')"
vmajor=`echo $BCL_VERSION | awk -F. '{print $1}'`
vminor=`echo $BCL_VERSION | awk -F. '{print $2}'`
vpatch=`echo $BCL_VERSION | awk -F. '{print $3}'`
echo "insert into bcl_downloads(os_arch,v_major,v_minor,v_patch,compile_date,file,filename) values('Mac OS X 10.4 or greater',$vmajor,$vminor,$vpatch,CURRENT_DATE(),$file_content,'bcl-${BCL_VERSION}.dmg') ON DUPLICATE KEY UPDATE compile_date=VALUES(compile_date), file=VALUES(file);" > ./mysql.cmd.txt
echo "insert into bcl_downloads(os_arch,v_major,v_minor,v_patch,compile_date,file,filename) values('Source zip',$vmajor,$vminor,$vpatch,CURRENT_DATE(),$srcfile_content,'bcl-${BCL_VERSION}-Source.zip') ON DUPLICATE KEY UPDATE compile_date=VALUES(compile_date), file=VALUES(file);" >> ./mysql.cmd.txt
chmod a+wx *.zip *.dmg mysql.cmd.txt
ssh carbon mysql --defaults-file=$HOME/.my.conf bclweb_sym < /net/$HOST/${BCL_ROOT}/${BUILD_PATH}/mysql.cmd.txt
chmod og-wx *.zip *.dmg mysql.cmd.txt

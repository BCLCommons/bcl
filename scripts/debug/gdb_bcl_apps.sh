#!/bin/sh

# determine the path to the bcl directory so that this script can be started from anywhere
EXE_CMD=`echo $0 | sed 's:^\([^./]\):./\1:' | sed 's:/gdb_bcl_apps.sh::'`
BCL_PATH=""
if [ "$EXE_CMD" == "." ]
  then
  # case where the script was started from the script folder
  BCL_PATH=".."
else
  # the script was started outside the script folder
  BCL_PATH=`echo $EXE_CMD | sed 's:/scripts/debug::'`
fi

# set breakpoints for bad pointer casts, asserts, etc.
ASSERT_LINE_NUMBER=`grep -n "Assert::Exit" $BCL_PATH/source/util/bcl_util_assert.cpp | sed 's/^\([^:]*\).*/\1/g'`
BAD_POINTER_CAST_LINE_NUMBER=`grep -n BCL_Message $BCL_PATH/source/util/bcl_util_ptr_interface.cpp | sed 's/^\([^:]*\).*/\1/g'`
BREAKPOINT1="$BCL_PATH/source/util/bcl_util_assert.cpp:$ASSERT_LINE_NUMBER"
BREAKPOINT2="$BCL_PATH/source/util/bcl_util_ptr_interface.cpp:$BAD_POINTER_CAST_LINE_NUMBER"

# make a gdb command file with the breakpoint statements
GDB_CMD_FILE=./.gdb_apps_debug.txt
echo -e "break $BREAKPOINT1\nbreak $BREAKPOINT2" > $GDB_CMD_FILE

# start gdb
gdb -x $GDB_CMD_FILE $BCL_PATH/build/linux64_debug/bin/bcl-apps-static.exe 
rm $GDB_CMD_FILE

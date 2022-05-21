#!/bin/sh
OLD_DIR=$PWD

# determine the path to the bcl directory so that this script can be started from anywhere
EXE_CMD=`echo $0 | sed 's:^\([^./]\):./\1:' | sed 's:/gdb_bcl_examples.sh::'`
BCL_PATH=""
if [ "$EXE_CMD" == "." ]
  then
  # case where the script was started from the script folder
  BCL_PATH=".."
else
  # the script was started outside the scripts folder
  BCL_PATH=`echo $EXE_CMD | sed 's:/scripts::' | sed 's:/debug::'`
fi

# set breakpoints for bad pointer casts, asserts, and failed examples
cd $BCL_PATH
BAD_EXAMPLE_LINE_NUMBER=`grep -n "++m_NumberErrors" $PWD/example/example.cpp | sed 's/^\([^:]*\).*/\1/g'`
ASSERT_LINE_NUMBER=`grep -n "Assert::Exit" $PWD/source/util/bcl_util_assert.cpp | sed 's/^\([^:]*\).*/\1/g'`
BAD_POINTER_CAST_LINE_NUMBER=`grep -n BCL_Message $PWD/source/util/bcl_util_ptr_interface.cpp | sed 's/^\([^:]*\).*/\1/g'`
BREAKPOINT1="$PWD/source/util/bcl_util_assert.cpp:$ASSERT_LINE_NUMBER"
BREAKPOINT2="$PWD/example/example.cpp:$BAD_EXAMPLE_LINE_NUMBER"
BREAKPOINT3="$PWD/source/util/bcl_util_ptr_interface.cpp:$BAD_POINTER_CAST_LINE_NUMBER"

cd -

# make a gdb command file with the breakpoint statements
GDB_CMD_FILE=./.gdb_example_debug.txt
echo -e "break $BREAKPOINT1\nbreak $BREAKPOINT2\nbreak $BREAKPOINT3" > $GDB_CMD_FILE

#start gdb
gdb -x $GDB_CMD_FILE $BCL_PATH/build/linux64_debug/bin/bcl-example-static.exe
rm $GDB_CMD_FILE
cd $OLD_DIR

#!/bin/bash
###############
#
# Author: Jeff Mendenhall
#
#############
help()
{
  echo "Usage: mimic_class OldClassName NewClassName [-ehcp]"
  echo "Without any of ech specified, automatically tries to mimic the h, hpp, cpp, fwd.hh, and example for the class"
  echo "-e mimic the example of OldClassName"
  echo "-c mimic the source of OldClassName"
  echo "-h mimic the header of OldClassName"
  echo "-p mimic the hpp of OldClassName"
  exit 1
}

if [ $# -lt 2 ]
  then
  help
fi

BCL_PATH="./"
SCRIPT_HOME="./scripts/code"

$SCRIPT_HOME/backup_class.sh $2

OLDNAME=$1
OLDCLASS=${OLDNAME##*::}
OLDNAMESPACE=${OLDNAME%::*}
NEWNAME=$2
NEWCLASS=${NEWNAME##*::}
NEWNAMESPACE=${NEWNAME%::*}

OLDFILENAME=`$SCRIPT_HOME/bcl_filename_for_scoped_class $OLDNAMESPACE$OLDCLASS`
NEWFILENAME=`$SCRIPT_HOME/bcl_filename_for_scoped_class $NEWNAMESPACE$NEWCLASS`

echo "Old filename core: $OLDFILENAME"
echo "New filename core: $NEWFILENAME"

DO_EXAMPLE=0
DO_H=0
DO_HPP=0
DO_CPP=0

OLDH="$BCL_PATH/include/$OLDNAMESPACE/bcl_$OLDFILENAME.h"
OLDHPP="$BCL_PATH/include/$OLDNAMESPACE/bcl_$OLDFILENAME.hpp"
OLDCPP="$BCL_PATH/source/$OLDNAMESPACE/bcl_$OLDFILENAME.cpp"
OLDEXAMPLEPATH="$BCL_PATH/example/$OLDNAMESPACE/example_$OLDFILENAME.cpp"
OLDEXAMPLE="$OLDNAMESPACE/example_$OLDFILENAME.cpp"
OLDNAMESPACESC=`echo $OLDNAMESPACE | sed "s/^\([a-z]\)/\U\1/"`
OLDFILEUC="BCL_"`echo $OLDFILENAME | sed "s/\([a-z]\)/\U\1/g"`
NEWH="$BCL_PATH/include/$NEWNAMESPACE/bcl_$NEWFILENAME.h"
NEWHPP="$BCL_PATH/include/$NEWNAMESPACE/bcl_$NEWFILENAME.hpp"
NEWCPP="$BCL_PATH/source/$NEWNAMESPACE/bcl_$NEWFILENAME.cpp"
NEWEXAMPLEPATH="$BCL_PATH/example/$NEWNAMESPACE/example_$NEWFILENAME.cpp"
NEWEXAMPLE="$NEWNAMESPACE/example_$NEWFILENAME.cpp"
NEWNAMESPACESC=`echo $NEWNAMESPACE | sed "s/^\([a-z]\)/\U\1/"`
NEWFILEUC="BCL_"`echo $NEWFILENAME | sed "s/\([a-z]\)/\U\1/g"`

shift 2

while getopts echp opt
do
  case "$opt" in
    e) DO_EXAMPLE=1;;
    h) DO_H=1;;
    c) DO_CPP=1;;
    p) DO_HPP=1;;
    \?) help
  esac
done

if [ $[ $DO_HPP + $DO_EXAMPLE + $DO_H + $DO_CPP ] -eq 0 ]
then
  DO_EXAMPLE=1
  DO_H=1
  DO_HPP=1
  DO_CPP=1
fi 

echo "Old Filename: " $OLDFILENAME
echo "Old Namespace: " $OLDNAMESPACE
echo "Old Class: " $OLDCLASS
echo "New Filename: " $NEWFILENAME
echo "New Namespace: " $NEWNAMESPACE
echo "New Class: " $NEWCLASS
echo "Old filename complete: $OLDH22"


mkdir -p $HOME/tmp/
COMMANDFILE="$HOME/tmp/commands.txt"
DATE_STR=`date "+%b %d, %Y"`
echo -e "s/$OLDFILEUC/$NEWFILEUC/g\ns|$OLDNAMESPACE/$OLDFILENAME|$NEWNAMESPACE/$NEWFILENAME|g\ns/$OLDFILENAME/$NEWFILENAME/g\n" > $COMMANDFILE
if [ "$OLDCLASS" != "" ]; then 
  echo -e "s/$OLDCLASS/$NEWCLASS/g\n" >> $COMMANDFILE
fi
echo -e "s|$OLDEXAMPLE|$NEWEXAMPLE|g\ns/Example$OLDNAMESPACESC/Example$NEWNAMESPACESC/gi\ns/$OLDNAMESPACE::$NEWCLASS/$NEWNAMESPACE::$NEWCLASS/g\ns://! @author .*://! @author $USER:\ns://! @date .*://! @date $DATE_STR:" >> $COMMANDFILE
echo -e "s|namespace $OLDNAMESPACE|namespace $NEWNAMESPACE|gi" >> $COMMANDFILE
echo "Sed commands written to " $COMMANDFILE

echo "mimicing..." 

if [ $DO_EXAMPLE -eq 1 ]; then
  if [ -e $OLDEXAMPLEPATH ]; then
    mkdir -p ${NEWEXAMPLEPATH%/*}
    echo "new example: " $NEWEXAMPLEPATH
    sed -f $COMMANDFILE $OLDEXAMPLEPATH > $NEWEXAMPLEPATH
  else
    echo "no example to mimic (looked at $OLDEXAMPLEPATH)"
  fi
fi
if [ $DO_H -eq 1 ]; then
  if [ -e $OLDH ]; then
    mkdir -p ${NEWH%/*} 
    echo "new header: " $NEWH
    sed -f $COMMANDFILE $OLDH > $NEWH
  else
    echo "no header to mimic (looked at $OLDH)"
  fi
fi
if [ $DO_HPP -eq 1 ]; then
  if [ -e $OLDHPP ]; then
    mkdir -p ${NEWHPP%/*} 
    echo "new hpp: " $NEWHPP
    sed -f $COMMANDFILE $OLDHPP > $NEWHPP
  else
    echo "no hpp to mimic (looked at $OLDHPP)"
  fi
fi
if [ $DO_CPP -eq 1 ]; then
  if [ -e $OLDCPP ]; then
    mkdir -p ${NEWCPP%/*}
    echo "new cpp: " $NEWCPP 
    sed -f $COMMANDFILE $OLDCPP > $NEWCPP
  else
    echo "no cpp to mimic (looked at $OLDCPP)"
  fi
fi

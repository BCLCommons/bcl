#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: backup_class ClassName"
  exit 1
fi

BCLDIR="./"
NAME=$1
CLASS=${NAME##*::}
NAMESPACE=${NAME%::*}
NEWDIR="../backup_$NAMESPACE"

FILENAME=`bcl_filename_for_scoped_class $NAMESPACE$CLASS`

mkdir -p $NEWDIR
if [ -f "$BCLDIR/include/$NAMESPACE/bcl_$FILENAME.h" ]
  then
  cp --backup=numbered $BCLDIR/include/$NAMESPACE/bcl_$FILENAME.h $NEWDIR/
fi
if [ -f "$BCLDIR/include/$NAMESPACE/bcl_$FILENAME.fwd.hh" ]
  then
  cp --backup=numbered $BCLDIR/include/$NAMESPACE/bcl_$FILENAME.fwd.hh $NEWDIR/
fi
if [ -f "$BCLDIR/source/$NAMESPACE/bcl_$FILENAME.cpp" ]
  then
  cp --backup=numbered $BCLDIR/source/$NAMESPACE/bcl_$FILENAME.cpp $NEWDIR/
fi
if [ -f "$BCLDIR/example/$NAMESPACE/example_$FILENAME.cpp" ]
  then
  cp --backup=numbered $BCLDIR/example/$NAMESPACE/example_$FILENAME.cpp $NEWDIR/
fi


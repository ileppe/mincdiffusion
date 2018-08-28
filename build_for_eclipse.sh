#! /bin/sh

VTK=$1
INSTALL=$2

cmake  -G "Eclipse CDT4 - Unix Makefiles" \
       -D "CMAKE_BUILD_TYPE:String=DEBUG" \
       -D "CMAKE_INSTALL_PREFIX:PATH=$INSTALL" \
       -D "VTK_DIR:PATH=$VTK" \
       $(dirname $0)

echo "*******************************"
echo
echo run \"make\" to compile
echo run \"make install\" to install into target directory
echo
echo "******************************"

#! /bin/sh

VTK=$1
LIBMINC=$2
INSTALL=$3


if [ "x$VTK" == "x" -o "x$LIBMINC" == "x" ];then
echo Usage $0 VTK_DIR LIBMINC_BASE_DIR [INSTALL_DIR]
exit 1
fi

cmake  -D "CMAKE_BUILD_TYPE:STRING=Release" \
       -D "CMAKE_INSTALL_PREFIX:PATH=$INSTALL" \
       -D "VTK_DIR:PATH=$VTK" \
       -D "LIBMINC_DIR:PATH=${LIBMINC}/lib" \
       $(dirname $0)

echo "*******************************"
echo
echo run \"make\" to compile
echo run \"make install\" to install into target directory
echo
echo "******************************"

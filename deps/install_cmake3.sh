#!/bin/bash
#!/bin/bash

ver="cmake-3.0.2"
wget http://cmake.org/files/v3.0/"$ver".tar.gz
mkdir -v ./cmake_install


tar xfz ./"$ver".tar.gz
install_prefix=`pwd`/cmake_install
echo "installing to $install_prefix"

cd ./"$ver"
./bootstrap --prefix="$install_prefix"
make -j 4
make install

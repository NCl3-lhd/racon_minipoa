racon-minipoa
==
To replace SPOA with minipoa in Racon:
```
# download Racon
wget https://github.com/lbcb-sci/racon/releases/download/1.4.13/racon-v1.4.13.tar.gz
tar -zxvf racon-v1.4.13.tar.gz && cd racon-v1.4.13

# download minipoa
git clone https://github.com/NCl3-lhd/minipoa.git

# replace SPOA with minipoa:
cp ../minipoa -r vendor/
cp ../CMakeLists.txt CMakeLists.txt
cp ../polisher.cpp src/polisher.cpp
cp ../polisher.hpp src/polisher.hpp
cp ../window.cpp src/window.cpp
cp ../window.hpp src/window.hpp

# build racon-minipoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
rm -rf eloss *.mod ./build
mkdir build
cd build
cmake ..
make
cp eloss ../
cd ../
./eloss

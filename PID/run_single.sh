cd build 
rm -rf *
cmake ../
make
cd ..
root -l run_single.cpp

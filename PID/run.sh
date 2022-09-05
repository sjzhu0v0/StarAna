cd build 
rm -rf *
cmake ../
make 
cd ../
root -l run.cpp

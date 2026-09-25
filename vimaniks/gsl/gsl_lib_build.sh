
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
tar -xf gsl-latest.tar.gz

cd ./gsl-2.8

./configure --prefix=$(dirname ${PWD})/local  CFLAGS=-fPIC CXXFLAGS=-fPIC --enable-shared

make && make install
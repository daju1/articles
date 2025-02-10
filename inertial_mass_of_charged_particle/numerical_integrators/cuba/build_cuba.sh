wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
tar zxvf Cuba-4.2.tar.gz
cd Cuba-4.2
./configure --prefix=$PWD/../local
make
make install
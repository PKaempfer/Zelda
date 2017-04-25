#make clean
autoheader
aclocal
automake --add-missing
autoconf
./configure --prefix $(pwd)
make
make install

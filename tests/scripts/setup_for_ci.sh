pip install rosetta_finder -U
sudo apt-get install perl -y

wget https://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2 -O parallel.tar.bz2
tar -xf parallel.tar.bz2; rm -f parallel.tar.bz2;cd parallel-20*
./configure && make && sudo make install



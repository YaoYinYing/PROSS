pip install rosetta_finder -U
pip install biopython

apt-get update -y; apt-get install perl wget build-essential -y

wget https://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2 -O parallel.tar.bz2

tar -xf parallel.tar.bz2; rm -f parallel.tar.bz2;
cd $(ls |grep '^parallel-20' |head -1); ./configure && make && make install
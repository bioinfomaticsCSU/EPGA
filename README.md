EPGA
====

De Novo Assembler

EPGA is more efficient for paired-read libraries of read length shorter than 40bp and coverage larger than 100.

1)Installing.

EPGA is written C++ and therefore will require a machine with GNU C++ pre-installed.

Create a main directory (eg:EPGA). Copy all source code to this directory.

Run command line: g++ main.cpp -o epga -lpthread

2)Running.

./epga libraryName insertsize sd libraryName1 insertsize1 sd1 kmerLength threadNumber

libraryName is paired-end reads file name in Fasta format (eg: read1.fa). Paired-end reads are in single file. Left mate read and right mate read stored in file one by one.(For mate-paired reads, please transform them to paired-end reads.)

insertsize is the sequence fragment length about paired-end reads (eg: 500).

sd is the standard deviation of insertsize.

kmerLength is one integer shorter than read length which is used for building De Bruijn graph.

threadNumber is thread number of program.


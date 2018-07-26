CXX=g++
SAMTOOLS_ROOT=vendor/samtools-0.1.19

FLAGS=-O2 -fopenmp
LFLAGS=-lm -L${SAMTOOLS_ROOT} -lbam -lz -lpthread
IFLAGS=-I${SAMTOOLS_ROOT}
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	        $(CXX) $(CFLAGS) $(CXXFLAGS) $(FLAGS) $(IFLAGS) -c $< -o $@

all: samtools msisensor

samtools:
	        $(MAKE) -C ${SAMTOOLS_ROOT}

msisensor: $(OBJS)
	        $(CXX) $^ $(CFLAGS) $(CXXFLAGS) $(FLAGS) $(LFLAGS) -o $@ 

clean:
	        rm -f *.o msisensor
			        $(MAKE) -C ${SAMTOOLS_ROOT} clean

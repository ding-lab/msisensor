CXX=g++
SAMTOOLS_ROOT=vendor/samtools-0.1.19

CXXFLAGS=-O2 -fopenmp
CXXFLAGS+=-I${SAMTOOLS_ROOT}
CXXLDFLAGS=-lm -L${SAMTOOLS_ROOT} -lbam -lz -lpthread
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	    $(CXX) $(CXXFLAGS) -c $< -o $@

all: samtools msisensor

samtools:
	    $(MAKE) -C ${SAMTOOLS_ROOT}

msisensor: $(OBJS)
	    $(CXX) $^ $(CXXFLAGS) $(CXXLDFLAGS) -o $@ 

clean:
	    rm -f *.o msisensor

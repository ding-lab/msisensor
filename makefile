CXX=g++
SAMTOOLS_ROOT=vendor/samtools-0.1.19


CXXFLAGS+=    -O2 -fopenmp
LDFLAGS+=    -L${SAMTOOLS_ROOT}
LIBS+=    -lm -lbam -lz -lpthread
INCLUDES+=    -I${SAMTOOLS_ROOT}
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	        $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: samtools msisensor

samtools:
	        $(MAKE) -C ${SAMTOOLS_ROOT}

msisensor: $(OBJS)
	        $(CXX) $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ 

clean:
	        rm -f *.o msisensor
			        $(MAKE) -C ${SAMTOOLS_ROOT} clean


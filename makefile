CC=g++

FLAGS=-O2 -fopenmp
CFLAGS=-O2 -fopenmp

#FLAGS=-g -fopenmp
#CFLAGS=-g -fopenmp

#SAMTOOLS_ROOT=/home/scbniu/software/samtools-1.8/
#HTSLIB_ROOT=/home/scbniu/software/samtools-1.8/htslib-1.8/
FLAGS+=-I${SAMTOOLS_ROOT}
FLAGS+=-I${HTSLIB_ROOT}
LFLAGS=-lm -L${SAMTOOLS_ROOT} -L${HTSLIB_ROOT} -lbam -lhts -lz -lpthread
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic
OBJS= $(patsubst %,%.o,$(SOURCE))

all: check-samtools check-htslib msisensor

%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@

check-samtools:
    ifndef SAMTOOLS_ROOT
        $(error SAMTOOLS_ROOT is undefined)
    endif

check-htslib:
    ifndef HTSLIB_ROOT
        $(error HTSLIB_ROOT is undefined)
    endif

msisensor: $(OBJS)
	$(CC) $^ $(CFLAGS) $(LFLAGS) -o $@ 
clean:
	rm -f *.o msisensor



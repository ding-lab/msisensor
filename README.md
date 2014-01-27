MSIsensor
===========
MSIsensor is a C++ program for automatically detecting somatic and germline variants at microsatellite regions. It computes length distributions of microsatellites per site in paired tumor and normal sequence data, subsequently using these to statistically compare observed distributions in both samples. Comprehensive testing indicates MSIsensor is an efficient and effective tool for deriving MSI status from standard tumor-normal paired sequence data.

Usage
-----

        Version 0.2
        Usage:  msisensor <command> [options]

Key commands:

        scan            scan homopolymers and miscrosatelites
        msi             msi scoring

msisensor scan [options]:
       
       -d   <string>   reference genome sequences file, *.fasta format
       -o   <string>   output homopolymer and microsatelites file

       -l   <int>      minimal homopolymer size, default=5
       -c   <int>      context length, default=5
       -m   <int>      maximal homopolymer size, default=50
       -s   <int>      maximal length of microsatellite, default=5
       -r   <int>      minimal repeat times of microsatellite, default=3
       -p   <int>      output homopolymer only, 0: no; 1: yes, default=0
       
       -h   help
 
msisensor msi [options]:

       -d   <string>   homopolymer and microsatellites file
       -n   <string>   normal bam file ( bam index file is needed )
       -t   <string>   tumor  bam file ( bam index file is needed )
       -o   <string>   output distribution file

       -e   <string>   bed file, to select a few resions
       -f   <double>   FDR threshold for somatic sites detection, default=0.05 
       -r   <string>   choose one region, format: 1:10000000-20000000
       -l   <int>      mininal homopolymer size, default=5
       -p   <int>      mininal homopolymer size for distribution analysis, default=10
       -m   <int>      maximal homopolymer size for distribution analysis, default=50
       -q   <int>      mininal microsatellites size, default=3
       -s   <int>      mininal number of repeats in microsatellites for distribution analysis, default=5
       -w   <int>      maximal microsatellites size for distribution analysis, default=40
       -u   <int>      span size around window for extracting reads, default=500
       -b   <int>      threads number for parallel computing, default=1
       -x   <int>      output homopolymer only, 0: no; 1: yes, default=0
       -y   <int>      output microsatellite only, 0: no; 1: yes, default=0
       
       -h   help

Install
-------
The Makefile assumes that you have the samtools source code in an environment variable `$SAMTOOLS_ROOT`. 

you don't know what that means, then simply follow these steps from any directory that you have permissions to write into:
Install some prerequisite packages if you are using Debian or Ubuntu:

    sudo apt-get install git libbam-dev zlib1g-dev

If you are using Fedora, CentOS or RHEL, you'll need these packages instead:

    sudo yum install git samtools-devel zlib-devel

Download the samtools-0.1.19 from SOURCEFORGE (http://sourceforge.net/projects/samtools/files/samtools/0.1.19):

    tar jxf samtools-0.1.19.tar.bz2
    cd samtools-0.1.19
    make
    export SAMTOOLS_ROOT=$PWD

Clone the msisensor repos, and build the `msisensor` binary:

    git clone https://github.com/ding-lab/msisensor.git
    cd msisensor
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv msisensor /usr/local/bin/

We have also provided pre-build binary distributions for Linux x86_64 and Mac OS X in directory: ./binary

    msisensor_Linux_x86_64: for Linux x86_64
    msisensor_Mac_OS_X    : for Mac OS X


Example
-------
1. Scan microsatellites from reference genome:
  
        msisensor scan -d referen.fa -o microsatellites.list

2. Msi scorring: 

        msisensor msi -d microsatellites.list -n normal.bam -t tumor.bam -e bed.file -o output.prefix -l 1 -q 1 -b 2

   Note: normal and tumor bam index files are needed in the same directory as bam files 

Output
-------
There will be one microsatellite list output in "scan" step. Msi scorring step will give 4 output files based on given output prefix:
        
        output.prefix
        output.prefix_dis
        output.prefix_germline
        output.prefix_somatic

1. microsatellites.list: microsatellite list output ( columns with *_binary means: binary conversion of DNA bases based on A=00, C=01, G=10, and T=11 )

        chromosome      location        repeat_unit_length     repeat_unit_binary    repeat_times    left_flank_binary     right_flank_binary      repeat_unit_bases      left_flank_bases       right_flank_bases
        1       10485   4       149     3       150     685     GCCC    AGCCG   GGGTC
        1       10629   2       9       3       258     409     GC      CAAAG   CGCGC
        1       10652   2       2       3       665     614     AG      GGCGC   GCGCG
        1       10658   2       9       3       546     409     GC      GAGAG   CGCGC
        1       10681   2       2       3       665     614     AG      GGCGC   GCGCG

2. output.prefix: msi score output

        Total_Number_of_Sites   Number_of_Somatic_Sites %
        640     75      11.72

3. output.prefix_dis: read count distribution (N: normal; T: tumor)

        1 10529896 CTTTC 15[T] GAGAC
        N: 0 0 0 0 0 0 0 1 0 0 8 9 1 7 17 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        T: 0 0 0 0 0 0 0 0 0 1 19 14 17 9 32 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 

4. output.prefix_somatic: somatic sites detected ( FDR: false discovery rate ) 

        chromosome   location        left_flank     repeat_times    repeat_unit_bases    right_flank      difference      P_value    FDR     rank
        1       16200729        TAAGA   10      T       CTTGT   0.55652 2.8973e-15      1.8542e-12      1
        1       75614380        TTTAC   14      T       AAGGT   0.82764 5.1515e-15      1.6485e-12      2
        1       70654981        CCAGG   21      A       GATGA   0.80556 1e-14   2.1333e-12      3
        1       65138787        GTTTG   13      A       CAGCT   0.8653  1e-14   1.6e-12 4
        1       35885046        TTCTC   11      T       CCCCT   0.84682 1e-14   1.28e-12        5
        1       75172756        GTGGT   14      A       GAAAA   0.57471 1e-14   1.0667e-12      6
        1       76257074        TGGAA   14      T       GAGTC   0.66023 1e-14   9.1429e-13      7
        1       33087567        TAGAG   16      A       GGAAA   0.53141 1e-14   8e-13   8
        1       41456808        CTAAC   14      T       CTTTT   0.76286 1e-14   7.1111e-13      9

5. output.prefix_germline: germline sites detected
    
        chromosome   location        left_flank     repeat_times    repeat_unit_bases    right_flank      genotype
        1       1192105 AATAC   11      A       TTAGC   5|5
        1       1330899 CTGCC   5       AG      CACAG   5|5
        1       1598690 AATAC   12      A       TTAGC   5|5
        1       1605407 AAAAG   14      A       GAAAA   1|1
        1       2118724 TTTTC   11      T       CTTTT   1|1


Test sample
-------
We provided one small sized sample data (tumor and matched normal bam files) for user to try msi scoring step.
It is very simple to run this test using sample data:

        cd ./test
        bash run.sh

Contact
-------
Please contact Beifang Niu by bniu@genome.wustl.edu and Kai Ye by kye@genome.wustl.edu if you have any questions.


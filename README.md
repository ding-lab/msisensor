MSIsensor
===========
MSIsensor is a c++ program for automatically detecting somatic microsatellite changes. It computes length distributions of microsatellites per site in paired tumor and normal sequence data, subsequently using these to statistically compare observed distributions in both samples. Comprehensive testing indicates MSIsensor is an efficient and effective tool for deriving MSI status from standard tumor-normal paired sequence data.

Usage
-----

        Version 0.1
        Usage:  msisensor <command> [options]

Key commands:

        scan            scan homopolymers and miscrosatelites
        msi             msi scoring

msisensor scan [options]:
       
       -d   <string>   reference genome sequences file, *.fasta format
       -o   <string>   output homopolymer and microsatelittes file

       -l   <int>      minimal homopolymer size, default=5
       -c   <int>      context length, default=5
       -m   <int>      maximal homopolymer size, default=50
       -s   <int>      maximal length of microsate, default=5
       -r   <int>      minimal repeat times of microsate, default=3
       -p   <int>      output homopolymer only, 0: no; 1: yes, default=0
       
       -h   help
 
msisensor msi [options]:

       -d   <string>   homopolymer and microsates file
       -n   <string>   normal bam file
       -t   <string>   tumor  bam file
       -o   <string>   output distribution file

       -e   <string>   bed file
       -f   <double>   FDR threshold for somatic sites detection, default=0.05 
       -r   <string>   choose one region, format: 1:10000000-20000000
       -l   <int>      mininal homopolymer size, default=5
       -p   <int>      mininal homopolymer size for distribution analysis, default=10
       -m   <int>      maximal homopolymer size for distribution analysis, default=50
       -q   <int>      mininal microsates size, default=3
       -s   <int>      mininal microsates size for distribution analysis, default=5
       -w   <int>      maximal microstaes size for distribution analysis, default=40
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

Clone the samtools and msisensor repos, and build the `msisensor` binary:

    git clone https://github.com/samtools/samtools.git
    export SAMTOOLS_ROOT=$PWD/samtools
    git clone https://github.com/ding-lab/msisensor.git
    cd msisensor
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv msisensor /usr/local/bin/

Example
-------
1. Scan microsatellites from reference genome:
  
        msisensor scan -d referen.fa -o microsatellites.list

2. Msi scorring: 

        msisensor msi -d microsatellites.list -n normal.bam -t tumor.bam -e bed.file -o output.prefix -l 1 -q 1 -b 2


Output
-------
There will be one microsatellite list output in "scan" step. 
Msi scorring step will give 4 output files based on given output prefix:
        output.prefix
        output.prefix_dis
        output.prefix_germline
        output.prefix_somatic

1. microsatellites.list: microsatellite list output 

        chromosome      location        site_length     site_content    repeat_times    front_flank     tail_flank      site_bases      front_flank_bases       tail_flank_bases
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

4. output.prefix_somatic: somatic sites detected
  
        chromosome   location        front_flank     repeat_times    site_content    tail_flank      difference      P_value
        1       10357206        TTGAA   17      T       ACTTC   0.45670 0.00014
        1       11140610        TCTGG   11      A       CACAC   0.80855 0.00000
        1       11156045        ACATC   15      T       GAGAC   0.75281 0.00001
        1       12368705        GAGTG   15      T       GAGAT   0.51139 0.00000
        1       16200729        TAAGA   10      T       CTTGT   0.55652 0.00000
        1       16245610        AAGGG   10      T       GCATA   0.75928 0.00000

5. output.prefix_germline: germline sites detected
    
        chromosome   location        front_flank     repeat_times    site_content    tail_flank      xxx|xxxx
        1       1192105 AATAC   11      A       TTAGC   5|5
        1       1330899 CTGCC   5       AG      CACAG   5|5
        1       1598690 AATAC   12      A       TTAGC   5|5
        1       1605407 AAAAG   14      A       GAAAA   1|1
        1       2118724 TTTTC   11      T       CTTTT   1|1


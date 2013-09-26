MSIsensor
===========
Homopolymer and microsatellite analysis using bam files

Usage
-----

        Version 0.1
        Usage:  msisensor <command> [options]

Key commands:

        scan            scan homopolymers and miscrosatelites
        msi             msi scoring

This tool was originally designed to do homopolymer and microsatellites analysis. 


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
    git clone https://github.com/Beifang/msisensor.git
    cd msisensor
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv msisensor /usr/local/bin/

xxx


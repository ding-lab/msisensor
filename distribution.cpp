/*
 * distribution.cpp for MSIsensor
 * Copyright (c) 2013 Beifang Niu && Kai Ye WUGSC All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// System header files
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

// Static function declaration
#include "param.h"
#include "polyscan.h"
#include "distribution.h"
#include "utilities.h"
#include "sample.h"

// branch 
#include "cmds.h"

Param paramd;
PolyScan polyscan;
Sample sample;

std::string homoFile;
std::string bamListFile;
std::string normalBam;
std::string tumorBam;
std::string bedFile;
std::string disFile;

std::ifstream finH;
std::ifstream finM;
std::ifstream finB;
std::ofstream foutD;

std::string one_region;

void DisUsage(void) {
    std::cerr<<"\nUsage:  msisensor msi [options] \n\n"
        <<"       -d   <string>   homopolymer and microsates file\n"
        <<"       -n   <string>   normal bam file\n"
        <<"       -t   <string>   tumor  bam file\n"
        <<"       -o   <string>   output distribution file\n\n"

        <<"       -e   <string>   bed file, optional\n"
        <<"       -f   <double>   FDR threshold for somatic sites detection, default="<<paramd.fdrThreshold<<"\n"
        <<"       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default="<<paramd.covCutoff<<"\n"
        <<"       -r   <string>   choose one region, format: 1:10000000-20000000\n"    
        <<"       -l   <int>      mininal homopolymer size, default="<<paramd.MininalHomoSize<<"\n"
        <<"       -p   <int>      mininal homopolymer size for distribution analysis, default="<<paramd.MininalHomoForDis<<"\n"
        <<"       -m   <int>      maximal homopolymer size for distribution analysis, default="<<paramd.MaxHomoSize<<"\n"

        <<"       -q   <int>      mininal microsates size, default="<<paramd.MinMicrosate<<"\n"
        <<"       -s   <int>      mininal microsates size for distribution analysis, default="<<paramd.MinMicrosateForDis<<"\n"
        <<"       -w   <int>      maximal microstaes size for distribution analysis, default="<<paramd.MaxMicrosateForDis<<"\n"

        <<"       -u   <int>      span size around window for extracting reads, default="<<paramd.DisSpan<<"\n"
        <<"       -b   <int>      threads number for parallel computing, default="<<paramd.numberThreads<<"\n"
        <<"       -x   <int>      output homopolymer only, 0: no; 1: yes, default="<<paramd.HomoOnly<<"\n"
        <<"       -y   <int>      output microsatellite only, 0: no; 1: yes, default="<<paramd.MicrosateOnly<<"\n"
        <<"       \n"
        <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}

int dGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': homoFile = rgv[++i]; break;
            case 'n': normalBam  = rgv[++i]; break;
            case 't': tumorBam  = rgv[++i]; break;
            case 'o': disFile  = rgv[++i]; break;
            case 'e': bedFile  = rgv[++i]; break;
            case 'r': one_region = rgv[++i]; break;
            case 'f': paramd.fdrThreshold  = atof(rgv[++i]); break;
            case 'c': paramd.covCutoff = atoi(rgv[++i]); break;
            case 'l': paramd.MininalHomoSize = atoi(rgv[++i]); break;
            case 'p': paramd.MininalHomoForDis = atoi(rgv[++i]); break;
            case 'u': paramd.DisSpan = atoi(rgv[++i]); break;
            case 'm': paramd.MaxHomoSize = atoi(rgv[++i]); break;
            case 'q': paramd.MinMicrosate = atoi(rgv[++i]); break;
            case 's': paramd.MinMicrosateForDis = atoi(rgv[++i]); break;
            case 'w': paramd.MaxMicrosateForDis = atoi(rgv[++i]); break;
            case 'b': paramd.numberThreads = atoi(rgv[++i]); break;
            case 'x': paramd.HomoOnly= atoi(rgv[++i]); break;
            case 'y': paramd.MicrosateOnly = atoi(rgv[++i]); break;
            break;
            case 'h':DisUsage();
            case '?':DisUsage();    
        }
    }
    return i;
}

int HomoAndMicrosateDisMsi(int argc, char *argv[]) {
    if (argc == 1) DisUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout <<"Start at:  "<<Curr_Time() << std::endl;

    int noptions = dGetOptions(argc, argv);
    // process user defined region
    if (!one_region.empty()) {  
        if (!polyscan.ParseOneRegion(one_region)) {
            std::cerr<<"fatal error: Please give correct defined region format (-r) \n";
            exit(1);
        }
        polyscan.ifUserDefinedRegion = true;
    } else {
        polyscan.ifUserDefinedRegion = false;
    }
    // reading bed file if is exist
    finB.open(bedFile.c_str());
    if (finB) {
        std::cout << "loading bed regions ..." << std::endl;
        polyscan.LoadBeds(finB);
        polyscan.BedFilterorNot();
    }

    // load bam files
    polyscan.LoadBams( normalBam, tumorBam );

    // check homo/microsate file
    finH.open(homoFile.c_str());
    if (!finH) {
        std::cerr<<"fatal error: failed to open homopolymer and microsatellites file\n";
        exit(1);
    }
    std::cout << "loading homopolymer and microsatellite sites ..." << std::endl;
    polyscan.LoadHomosAndMicrosates(finH);
    finH.close();
    //polyscan.TestHomos();
    polyscan.SplitWindows();
    //polyscan.TestWindows();
    std::cout << "\nTotal loading windows:  " << polyscan.totalWindowsNum << " \n\n";
    std::cout << "\nTotal loading homopolymer and microsatellites:  " << polyscan.totalHomosites << " \n\n";

    // change code to one sample
    polyscan.GetHomoDistribution(sample, disFile);

    std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";

    return 0;
}



/*
 * scan.cpp for MSIsensor
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
#include <omp.h>

#include "param.h"
#include "refseq.h"
#include "scan.h"
#include "utilities.h"
#include "cmds.h"

Param param;
RefSeq ref;

std::string ref_file;
std::string out_file;
std::ifstream fin_d;
std::ofstream fout;

void ScanUsage(void) {
    std::cerr<<"\nUsage:  msisensor scan [options] \n\n"
        <<"       -d   <string>   reference genome sequences file, *.fasta format\n"
        <<"       -o   <string>   output homopolymer and microsatelittes file\n\n"
        <<"       -l   <int>      minimal homopolymer size, default="<<param.MininalHomoSize<<"\n"
        <<"       -c   <int>      context length, default="<<param.ContextLength<<"\n"
        <<"       -m   <int>      maximal homopolymer size, default="<<param.MaxHomoSize<<"\n"
        <<"       -s   <int>      maximal length of microsate, default="<<param.MaxMicrosate<<"\n"
        <<"       -r   <int>      minimal repeat times of microsate, default="<<param.Repeats<<"\n"
        <<"       -p   <int>      output homopolymer only, 0: no; 1: yes, default="<<param.HomoOnly<<"\n"
        <<"       \n"
        <<"       -h   help\n\n"
        << std::endl;
    exit(1);
}

int mGetOptions(int rgc, char *rgv[]) {
    int i;
    for (i=1; i<rgc; i++) {
        if (rgv[i][0] != '-') return i;
        switch(rgv[i][1]) {
            case 'd': ref_file = rgv[++i]; break;
            case 'o': out_file = rgv[++i]; break;
            case 'l': param.MininalHomoSize = atoi(rgv[++i]); break;
            case 'c': param.ContextLength = atoi(rgv[++i]); break;
            case 'm': param.MaxHomoSize = atoi(rgv[++i]); break;
            case 's': param.MaxMicrosate = atoi(rgv[++i]); break;
            case 'r': param.Repeats = atoi(rgv[++i]); break;
            case 'p': param.HomoOnly = atoi(rgv[++i]); break;
            break;
            case 'h':ScanUsage();
            case '?':ScanUsage();    
        }
    }
    return i;
}

int HomoAndMicrosateScan(int argc, char *argv[]) {
    if (argc == 1) ScanUsage();
    for (int i=0; i<argc; i++) {
        std::cout <<argv[i]<<' ';
    }
    Initial_Time();
    std::cout <<"Start at:  "<<Curr_Time()<< std::endl;
    int noptions = mGetOptions(argc, argv);
    // check refseq file
    fin_d.open(ref_file.c_str());
    if (!fin_d) {
        std::cerr<<"fatal error: failed to open ref file\n";
        exit(1);
    }
    // output calling results
    fout.open(out_file.c_str());
    if (!fout) {
        std::cerr <<"failed to open file: "<<out_file<< std::endl;
        exit(1);
    }
    ref.PouroutHeader(fout);
    // reading refseq and count homo sites
    ref.ScanHomoAndMicrosate(fin_d);
    ref.PouroutBuffer(fout);
    std::cout << "\nTotal time consumed:  " << Cal_AllTime() << " secs\n\n";
    fout.close();
    fin_d.close();

    return 0;
}


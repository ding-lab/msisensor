/*
 * window.cpp for MSIsensor 
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

#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>

#include "window.h"
#include "param.h"
#include "polyscan.h"

extern Param paramd;
extern PolyScan polyscan;
extern char homo_code[];
extern char uhomo_code[];
extern bit8_t alphabet[];

Window::Window() 
    : _start(0)
    , _end(0)
    , _chr("")
    , _siteCount(0)
    , _startSite(NULL)
    , _endSite(NULL)
{
    //xxxx
    //
};

Window::~Window() {
    //xxx
};

void Window::InitialDisW() {
    HomoSite *p = NULL;
    for (unsigned short i=0; i<_siteCount; i++) {
        p = _startSite + i;
        p->InitialDis();
    }
};

void Window::ClearDis() {
    HomoSite *p = NULL;
    for (unsigned short i=0; i<_siteCount; i++) {
        p = _startSite + i;
        p->ReleaseMemory();
    }
};

void Window::OutputDisW() {
    HomoSite *p = NULL;
    for (unsigned short i=0; i<_siteCount; i++) {
        p = _startSite + i;
        p->OutputDis();
    }
};

void Window::PouroutDisW(Sample &oneSample) {
    HomoSite *p = NULL;
    for (unsigned short i=0; i<_siteCount; i++) {
        p = _startSite + i;
        p->PouroutDis(oneSample);
    }
};

void Window::DisGenotypingW(Sample &oneSample) {
    HomoSite *p = NULL;
    for (unsigned short i=0; i<_siteCount; i++) {
        p = _startSite + i;
        p->DisGenotyping(oneSample);
    }
};

// change start
void Window::ChangeStart() {
    if ((_start - MAX_SPAN_SIZE) < 0) {
        _start = 0;
    } else { _start -= MAX_SPAN_SIZE; }
}

void Window::GetDistribution(std::vector <SPLIT_READ> &readsInWindow) {
    for (unsigned short j=0; j<polyscan.totalBamPairsNum; j++) {
        // normal
        if (!polyscan.totalBamPairs[j].normal_bam.empty()) {
            // extract reads
            LoadReads(readsInWindow, polyscan.totalBamPairs[j].normal_bam.c_str());
            ScanReads(readsInWindow, j, false);
            readsInWindow.clear();
        }
        // tumor
        if (!polyscan.totalBamPairs[j].tumor_bam.empty()) {
            // extract reads
            LoadReads(readsInWindow, polyscan.totalBamPairs[j].tumor_bam.c_str());
            ScanReads(readsInWindow, j, true);
            readsInWindow.clear();
        }
    }
}

void Window::LoadReads( std::vector <SPLIT_READ> &readsInWindow, 
                        const std::string bam ) {
    std::string tag = "";
    if (!bam.empty()) {
        // extract reads
        ReadInBamReads(bam.c_str(), _chr, _start, _end, readsInWindow, tag);
    }
}

void Window::ScanReads( const std::vector <SPLIT_READ> &readsInWindow, 
                        unsigned short bamIndex, 
                        bool isTumor) {

    // openmp parallel 
    omp_set_num_threads( paramd.numberThreads );
#pragma omp parallel for
    for (unsigned short i=0; i<_siteCount; i++) {
        HomoSite *p = _startSite + i;
        unsigned long tsize = readsInWindow.size();
        for (unsigned long j=0; j<tsize; j++) {
            if ( readsInWindow[j].Mapped ) {
                if ( (readsInWindow[j].MatchedRelPos < p->lowcut) || (readsInWindow[j].MatchedRelPos > p->highcut) ) continue; 
            }
            unsigned short tCount = DoOneRead(readsInWindow[j].ReadSeq, p);
            if ( (tCount > 0) && (tCount < paramd.s_dispots) ) {
                if (isTumor) {
                    p->tumorDis[bamIndex][tCount-1]++;
                } else {
                    p->normalDis[bamIndex][tCount-1]++;
                }
            } else {
                // don't scan reverse if mapped
                if ( readsInWindow[j].Mapped ) continue; 
                // reverse
                std::string tStr = readsInWindow[j].ReadSeq;
                ReverseComplement(tStr);
                tCount = DoOneRead(tStr, p);
                if ( (tCount > 0) && (tCount < paramd.s_dispots) ) {
                    if (isTumor) {
                        p->tumorDis[bamIndex][tCount-1]++;
                    } else {
                        p->normalDis[bamIndex][tCount-1]++;
                    }
                }
            }
        }
    }
}

unsigned short Window::DoOneRead( const std::string &oneRead, const HomoSite *p ) {
    std::string::size_type startPos = 0;
    unsigned short count = 0;
    while (std::string::npos != (startPos = oneRead.find(p->fbases, startPos))) {
        //std::cout<<startPos<<std::endl;
        count = 0;
        std::string::size_type tstart0 = startPos+p->fbases.length();
        std::string::size_type tstart  = tstart0;
        while ( tstart0 == (tstart = oneRead.find(p->bases, tstart)) ) {
            count++;
            tstart += p->bases.length();
            tstart0 = tstart;
        }
        // if get one
        if ( (p->typeLen == 1) && (count >= paramd.MininalHomoSize)
             ||
             (p->typeLen >  1) && (count >= paramd.MinMicrosate)
           ) {
            tstart = tstart0;
            if (tstart == (tstart0 = oneRead.find(p->ebases, tstart0))) {
                return count;
            }
        }
        startPos++;
    }
    return 0;
}

// reverse complement string
void Window::ReverseComplement(std::string &theWord) {       
    char tempChar;
    unsigned int t_uint0;
    unsigned int t_uint1;    
    for (int i=0; i<theWord.length()/2; i++) {
        tempChar = theWord[i];
        t_uint0 = alphabet[tempChar];
        t_uint1 = alphabet[theWord[theWord.length()-i-1]];
        theWord[i] = uhomo_code[t_uint1];
        theWord[theWord.length()-i-1] = uhomo_code[t_uint0];
    }
    if (theWord.length()%2) {
        tempChar = theWord[theWord.length()/2];
        t_uint0 = alphabet[tempChar];
        theWord[theWord.length()/2] = uhomo_code[t_uint0];
    }
}  


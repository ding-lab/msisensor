/*
 * param.h for MSIsensor
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

#ifndef _PARAM_H_
#define _PARAM_H_

#include <string>

// maximal length of 
// microsatellite
#define MAX_MICROSATE_LEN 8
#define MAX_FLANK_REGION 8
#define MAX_TRANSFER_LINE_LENGTH 100
#define MAX_WINDOW 1000000
#define MAX_SPAN_SIZE 1000
#define MAX_READ_LENGTH 110

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;

typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;

class Param {
public:
    Param();
    ~Param();

    int max_dbseq_size; 
    int append_dbseq_size; 
    int bufSize;
    // Homo sites
    int MininalHomoSize;
    int ContextLength;
    // Microsate
    unsigned int MaxMicrosate;
    unsigned int Repeats;

    unsigned int MinMicrosate;
    unsigned int MinMicrosateForDis;
    unsigned int MaxMicrosateForDis;

    // filtering
    int HomoOnly;
    int MicrosateOnly;
    int ncpu;
    int chains;
    // Homo sites
    int MininalHomoForDis;
    int MaxHomoSize;
    int SpanSize;
    int DisSpan;
    // Thread number
    unsigned int numberThreads;
    // statistic var
    unsigned int s_dispots; 
    unsigned int PercentPairs;
    unsigned int PercentPairsNumber;
    unsigned int HomoCoverage;
    // window size
    unsigned int windowSize;

    // genotyping 
    unsigned int covCutoff;
    double fdrThreshold;

};

#endif //_PARAM_H_


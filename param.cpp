
/*
 * param.cpp for MSIsensor
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

#include "param.h"
#include <iostream>
#include <bitset>

bit8_t homo_alphabet[256];
void inital_homo_phabet() { 
    // convert char 'A' to 1
    for (int i=0; i<256; i++){
        homo_alphabet[i] = 0;
    }
    homo_alphabet['a'] = homo_alphabet['A'] = 1;
}
// Aa=0 
// Cc=1 
// Gg=2 
// Tt=3 
// N=4
//
bit8_t alphabet[256];
void initalphabet() { 
    // convert unknown char as 4
    for (int i=0; i<256; i++) { 
        alphabet[i] = 4;
    }
    alphabet['a'] = alphabet['A'] = 0;
    alphabet['c'] = alphabet['C'] = 1;
    alphabet['g'] = alphabet['G'] = 2;
    alphabet['t'] = alphabet['T'] = 3;
}
bit8_t rev_alphabet[256];
void initrevalphabet() { 
    // complementary chain 
    for (int i=0; i<256; i++) {
        rev_alphabet[i] = 4;
    }
    rev_alphabet['a'] = rev_alphabet['A'] = 3;
    rev_alphabet['c'] = rev_alphabet['C'] = 2;
    rev_alphabet['g'] = rev_alphabet['G'] = 1;
    rev_alphabet['t'] = rev_alphabet['T'] = 0;
}

// for homopolymer
char homo_code[5] = {'A', 'C', 'G', 'T', 'N'};
char uhomo_code[5] = {'T', 'G', 'C', 'A', 'N'};
char chain_flag[2] = {'+', '-'};
char nt_code[4] = {'A', 'C', 'G', 'T'};
char revnt_code[4] = {'T', 'G', 'C', 'A'};

Param::Param()
    : bufSize(500000)
    , ncpu(1)
    , chains(0)
    , max_dbseq_size(300000000) //300Mb
    , append_dbseq_size(0x1000000) //16Mb
    , MininalHomoSize(5)
    , ContextLength(5)
    , MaxHomoSize(50)
    , SpanSize(500)
    , MininalHomoForDis(10)
    // microsate
    , MinMicrosate(3)
    , MinMicrosateForDis(5)
    , MaxMicrosateForDis(40)

    , DisSpan(500)
    , HomoOnly(0)
    , MicrosateOnly(0)
    , s_dispots(100)
    , MaxMicrosate(5)
    , Repeats(3)
    , numberThreads(1)
    , PercentPairs(20)
    , PercentPairsNumber(2)
    , HomoCoverage(40)
    , windowSize(500000) // window size (default 50k)
    , covCutoff( 20 )
    , Normalization(0)
    , fdrThreshold( 0.05 ) 
    , comentropyThreshold( 1 ) // for tumor only
{
    inital_homo_phabet();
    initalphabet();
    initrevalphabet();
};

Param::~Param() {
    //xxx
};


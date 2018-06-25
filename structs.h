
/*
 * structures.h for MSIsensor
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


#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

// user defined region
struct UserDefinedRegion {
    UserDefinedRegion()
        : chr("")
        , start(-1)
        , end(-1)
    {
        //xxxx
    }
    std::string chr;
    int start;
    int end;
};

// bed region 
struct BedRegion {
    BedRegion()
        : start(0)
        , end(0)
    {
        //xxxx
    }
    int start; 
    int end;
};

// bed regions located on one chromosome
struct BedChr {
    BedChr()
        : chr("")
    {
        //xxx
    }
    std::string chr;
    std::vector< BedRegion > regions_list;
};

// genotype by Kai
struct Genotype {
    Genotype() {
        //xxx
        GT[0] = -1;
        GT[1] = -1;
        WithGenotype = false;
    }
    short GT[2];
    bool WithGenotype;
};

// Bam file pairs
struct BamPairs {
    BamPairs()
        : sName("")
        , normal_bam("")
        , tumor_bam("")
    {
        //xxx
    }
    std::string sName;
    // bam files
    std::string normal_bam;
    std::string tumor_bam;

};

// Bam file tumors
struct BamTumors {
    BamTumors()
        : sName("")
        , tumor_bam("")
    {
        //xxx
    }
    std::string sName;
    // bam files
    std::string tumor_bam;

};

#endif //_STRUCTS_H_


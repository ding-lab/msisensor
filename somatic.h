
/*
 * somatic.h for MSIsensor
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

#ifndef _SOMATIC_H_
#define _SOMATIC_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// homopolymer site
class SomaticSite {
public:
    SomaticSite();
    ~SomaticSite();

    // homo or microsate 
    // A/C/G/T/AC/AGC
    //
    std::string chr;

    // location
    int location;

    // repeat times 
    unsigned short length;

    // homo or microsate
    // content 
    std::string bases;
    // front kmer
    std::string fbases;
    // tail kmer
    std::string ebases;
    // difference between
    // normal and tumor 
    double diff;
    double pValue;
    double FDR;
    unsigned short rank;
    bool somatic;
    
    // output content
    void PourOut(); 
    // sorting based on p-value
    // 
    bool operator < (const SomaticSite& rhs) const { return pValue < rhs.pValue; }

    protected:
        // xxx
};

#endif //_SOMATIC_H_


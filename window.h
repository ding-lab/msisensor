
/*
 * window.h for MSIsensor
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

#ifndef _WINDOW_H_
#define _WINDOW_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"
#include "homo.h"
#include "bamreader.h"
#include "sample.h"

// homopolymer site
class Window {
public:
    Window();
    ~Window();

    int _start;
    int _end;
    unsigned short _siteCount;
    std::string _chr;
    HomoSite *_startSite;
    HomoSite *_endSite;

    void InitialDisW(); 
    void OutputDisW();
    void ClearDis();
    void ChangeStart(); 
    void GetDistribution(std::vector <SPLIT_READ> &readsInWindow);
    void LoadReads(std::vector <SPLIT_READ> &readsInWindow, const std::string bam);
    void ScanReads(const std::vector <SPLIT_READ> &readsInWindow, unsigned short bamIndex, bool isTumor);
    void ReverseComplement(std::string &theWord);
    unsigned short DoOneRead(const std::string &oneRead, const HomoSite *p);
    void PouroutDisW(Sample &oneSample);
    void DisGenotypingW(Sample &oneSample);

protected:
    //xxxxxx
};

#endif //_WINDOW_H_


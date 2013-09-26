
/*
 * polyscan.h for MSIsensor
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


#ifndef _POLYSCAN_H_
#define _POLYSCAN_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"
#include "structs.h"

#include "homo.h"
#include "window.h"
#include "sample.h"

class PolyScan
{
    public:

         PolyScan();
        ~PolyScan();

        // user defined region
        UserDefinedRegion region_one;
        bool ifUserDefinedRegion;
        void eliminate(const char ch, std::string & str);
        bool ParseOneRegion(const std::string &region);

        // read bed regions
        bool ifUserDefinedBed;
        std::map<std::string, bit16_t> chrMaptoIndex;
        std::vector< BedChr > beds;
        void LoadBeds(std::ifstream &fin);
        void BedFilterorNot();

        // load bam list file 
        std::vector< BamPairs > totalBamPairs;
        //void LoadBams(std::ifstream &fin);
        void LoadBams(const std::string &nBam, const std::string &tBam);
        unsigned int totalBamPairsNum;

        // load homos and microsatellites
        unsigned long totalHomosites; 
        //std::vector< HomoSite * > totalSites;
        std::vector< HomoSite  > totalSites;
        void LoadHomosAndMicrosates(std::ifstream &fin);
        void TestHomos();

        std::vector< HomoSite > homosBuffer;

        // windows
        std::vector< Window > totalWindows;
        void SplitWindows();
        void TestWindows();
        unsigned long totalWindowsNum;

        // distribution 
        void InithializeDistributions();
        void outputDistributions();
        void releaseDistributions(); 
        void GetHomoDistribution( Sample &oneSample, const std::string &prefix );
        
protected:

        // xxxxxx
        // xxxx

};

#endif //_POLYSCAN_H_


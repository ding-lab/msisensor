/*
 * polyscan.cpp for MSIsensor 
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
#include <map>
#include <omp.h>

#include "utilities.h"
#include "polyscan.h"
#include "bamreader.h"
#include "param.h"

extern Param paramd;
extern bit8_t alphabet[];
extern bit8_t rev_alphabet[];
extern char uhomo_code[];
extern char homo_code[];

PolyScan::PolyScan() { 
    homosBuffer.reserve(paramd.bufSize);
    totalSites.reserve(paramd.bufSize);
}

PolyScan::~PolyScan() { 
    totalSites.clear();
}

// eliminates a character from the input string
void PolyScan::eliminate(const char ch, std::string & str){
    size_t eliminateCharPos = str.find(ch);
    while (eliminateCharPos != std::string::npos) {
        str.erase(eliminateCharPos, 1);
        eliminateCharPos = str.find(ch);
    }
}

// Parse one region 
bool PolyScan::ParseOneRegion(const std::string & regionString) {
    size_t separatorPos = regionString.find(":");
    bool correctParse   = false;
    bool m_endDefined   = false;
    bool m_startDefined = false;
    int m_start = -1;
    int m_end   = -1;
    std::string m_targetChromosomeName;
    // found a separator
    if (separatorPos != std::string::npos) {
        m_targetChromosomeName = regionString.substr(0, separatorPos);
        std::string coordinates = regionString.substr(separatorPos + 1);
        // removes the ',' in 1,000 or 1,000,000 that users may add 
        // for readability but wreak havoc with atoi
        eliminate(',', coordinates); 
        size_t startEndSeparatorPos = coordinates.find("-");
        // there are two coordinates
        if (startEndSeparatorPos != std::string::npos) {
            std::string secondPositionStr = coordinates.substr(startEndSeparatorPos + 1);
            m_end = atoi(secondPositionStr.c_str());
            m_endDefined = true;
        }

        m_start = atoi(coordinates.c_str());
        m_startDefined = true;
        if (m_start < 0 || (m_endDefined && (m_end < m_start))) {
            correctParse = false;
        } else { correctParse = true; }
    }
    // no separator found
    else {
        m_targetChromosomeName = regionString;
        correctParse = true;
    }
    // assign values
    region_one.chr   = m_targetChromosomeName;
    region_one.start = m_start;
    region_one.end   = m_end;

    return correctParse;
}

// Loading bed regions
void PolyScan::LoadBeds(std::ifstream &fin) {
    std::string chr;
    std::string line;
    std::string tempChr = "";
    int start;
    int stop;
    int i = -1;
    while (getline(fin, line)){
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> start;
        linestream >> stop;
        /*
        std::cout<< chr <<"\t"
                 << start <<"\t"
                 << stop << "\n";
        */
        if (chr == tempChr) {
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
        } else {
            ++i;
            BedChr tempBedChr;
            tempBedChr.chr = chr;
            beds.push_back(tempBedChr);
            // load mapping
            chrMaptoIndex.insert(std::pair<std::string, bit16_t>(chr, i));
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
            tempChr = chr;
        } 
        linestream.clear();
        linestream.str("");
    } 
}

// loading bam list
// only load one bam file
void PolyScan::LoadBams(const std::string &bam1, const std::string &bam2) {
    BamPairs t_bampair;
    t_bampair.sName = "sample_name";

    if (bam1.find(".bam") != std::string::npos) { 
        t_bampair.normal_bam = bam1;
    } else { std::cerr << "please provide valid format normal bam file ! \n"; exit(0); }
    if (bam2.find(".bam") != std::string::npos) {
        t_bampair.tumor_bam = bam2; 
    } else { std::cerr << "please provide valid format tumor bam file ! \n"; exit(0); }

    // loading
    totalBamPairs.push_back(t_bampair);
    totalBamPairsNum++;
}

// read and load sites
void PolyScan::LoadHomosAndMicrosates(std::ifstream &fin) {
    std::string chr;
    std::string bases;
    std::string fbases;
    std::string ebases;
    std::string line;
    std::string tChr = "";
    // count total loading sites
    //
    totalHomosites = 0;
    int loc;
    bit8_t  siteLength;
    bit16_t tsiteLength;
    bit16_t siteBinary;
    bit16_t siteRepeats;
    bit16_t frontF;
    bit16_t tailF;

    int j = 0;
    BedChr tbedChr;
    BedRegion tbedRegion;
    bit16_t tIndex;
    
    // skip title
    getline(fin, line);
    while (getline(fin, line)) {
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> loc;
        linestream >> tsiteLength;
        linestream >> siteBinary;
        linestream >> siteRepeats;
        linestream >> frontF;
        linestream >> tailF;
        //xxx
        linestream >> bases;
        linestream >> fbases;
        linestream >> ebases;

        // filtering
        if (tsiteLength > 1 && paramd.HomoOnly == 1) continue;
        if (tsiteLength == 1 && paramd.MicrosateOnly == 1) continue;
        if (tsiteLength == 1 && ((siteRepeats < paramd.MininalHomoForDis) || (siteRepeats > paramd.MaxHomoSize)) ) continue;
        if (tsiteLength > 1 && ((siteRepeats < paramd.MinMicrosateForDis) || (siteRepeats > paramd.MaxMicrosateForDis)) ) continue;

        siteLength = tsiteLength & 255;
        // defined one region
        if (ifUserDefinedRegion) {
            if (chr != region_one.chr) {
                continue;
            } else {
                if ( 
                     (loc < region_one.start) 
                     || 
                     ((loc + siteLength * siteRepeats) > region_one.end)
                   ) { continue; }
            }
        }
        // bed filtering
        if (ifUserDefinedBed) {
            // new chr 
            if (tChr != chr) {
                j = 0;
                if (chrMaptoIndex.count(chr) > 0) {
                    tbedChr = beds[chrMaptoIndex[chr]];
                    tChr = tbedChr.chr;
                    tbedRegion = tbedChr.regions_list[j];
                } else { continue; }
            }
            // filtering 
            if (loc < tbedRegion.start) continue;
            if (loc > tbedRegion.end) {
                for (j; j<tbedChr.regions_list.size(); j++) {
                    tbedRegion = tbedChr.regions_list[j];
                    if (loc < tbedRegion.end) { break; }
                }
                if (j >= tbedChr.regions_list.size()) continue;
            }
            if ((loc + siteLength * siteRepeats) > tbedRegion.end) continue;
        }

        // load sites 
        //HomoSite *toneSite = new HomoSite;
        HomoSite toneSite;
        toneSite.chr = chr;
        toneSite.location = loc;
        toneSite.typeLen = siteLength;
        toneSite.homoType = siteBinary;
        toneSite.length = siteRepeats;
        toneSite.frontKmer = frontF;
        toneSite.endKmer = tailF;
        toneSite.bases = bases;
        toneSite.fbases = fbases;
        toneSite.ebases = ebases;

        toneSite.lowcut = ((loc - MAX_READ_LENGTH) > 0) ? (loc - MAX_READ_LENGTH) : 0;
        toneSite.highcut = loc + MAX_READ_LENGTH;

        totalSites.push_back(toneSite);
        totalHomosites++;

        linestream.clear();
        linestream.str("");

    } // end while

}

// bed regions ?
void PolyScan::BedFilterorNot() {
    if (beds.size() > 0) ifUserDefinedBed = true;
}

// test sites loading
void PolyScan::TestHomos() {
    for (unsigned long i=0; i<totalHomosites; i++) {
        HomoSite *toneSite = &totalSites[i];
        std::cout << toneSite->chr<<"\t"
                  << toneSite->location<<"\t"
                  << int(toneSite->typeLen)<<"\t"
                  << toneSite->homoType<<"\t"
                  << toneSite->length<<"\t"
                  << toneSite->frontKmer<<"\t"
                  << toneSite->endKmer<<"\t"
                  << sizeof(*toneSite) <<"\n";
    }
}

// split windows
void PolyScan::SplitWindows() {
    Window oneW;
    HomoSite *second;
    HomoSite *first = &totalSites[0];

    oneW._start = first->location;
    oneW._end = oneW._start;
    oneW._chr = first->chr;
    oneW._startSite = oneW._endSite = &totalSites[0];
    for (int i=1; i< totalHomosites; i++) {
        first = &totalSites[i];
        if ( (first->chr == oneW._chr) 
             && 
             (first->location - oneW._start) < paramd.windowSize) {
            continue;
        }
        oneW._end = totalSites[i-1].location + MAX_SPAN_SIZE;
        oneW._endSite = &totalSites[i-1];
        oneW._siteCount = oneW._endSite - oneW._startSite + 1;
        // record one window
        oneW.ChangeStart();
        totalWindows.push_back(oneW);
        totalWindowsNum++;
        oneW._start = first->location;
        oneW._end = oneW._start;
        oneW._chr = first->chr;
        oneW._startSite = oneW._endSite = &totalSites[i];
    }
    oneW._end = first->location + MAX_SPAN_SIZE;
    oneW._endSite = first;
    oneW._siteCount = oneW._endSite - oneW._startSite + 1;
    // record this window
    oneW.ChangeStart();
    totalWindows.push_back(oneW);
    totalWindowsNum++;
}

// test windows
void PolyScan::TestWindows() {
    Window *oneW;
    for (int i=0; i< totalWindowsNum; i++) {
        oneW = &totalWindows[i];
        std::cout << oneW->_chr <<"\t"
                  << oneW->_start <<"\t"
                  << oneW->_siteCount <<"\t"
                  << oneW->_startSite->chr<<"\t"
                  << oneW->_startSite->location<<"\t"
                  << oneW->_endSite->chr<<"\t"
                  << oneW->_endSite->location<<"\n";
    }
}

// initial distribution
void PolyScan::InithializeDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
    }
}

// release distribution
void PolyScan::releaseDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].ClearDis();
    }
}

// output distribution
void PolyScan::outputDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].OutputDisW();
    }
}

// get distribution 
//void PolyScan::GetHomoDistribution( std::ofstream &fout ) {
void PolyScan::GetHomoDistribution( Sample &oneSample, const std::string &prefix ) {
    oneSample.iniOutput(prefix);
    std::vector< SPLIT_READ > readsInWindow;
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
        totalWindows[i].GetDistribution(readsInWindow);
        totalWindows[i].PouroutDisW(oneSample);
        totalWindows[i].DisGenotypingW(oneSample);
        totalWindows[i].ClearDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" <<  totalWindows[i]._chr << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end << std::endl;
    }
    // FDR
    oneSample.calculateFDR();
    oneSample.pourOutSomaticFDR();
    // MSI score
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}


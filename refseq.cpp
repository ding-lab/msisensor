/*
 * refseq.cpp for MSIsensor 
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

#include "utilities.h"
#include "refseq.h"

extern Param param;

extern bit8_t alphabet[];
extern bit8_t rev_alphabet[];

extern char uhomo_code[];
extern char homo_code[];

extern std::ofstream fout;

RefSeq::RefSeq()
    : totalHomosites(0)
    , totalMicrosates(0)
    , totalSites(0)
{
    // buffer resize
    homosBuffer.reserve(param.bufSize);
};

RefSeq::~RefSeq() { 
   //*** 
};

// output and clear buffer
bool RefSeq::PouroutHeader(std::ofstream &fout) {
    fout << "chromosome"         <<"\t"
         << "location"           <<"\t"
         << "repeat_unit_length" <<"\t"
         << "repeat_unit_binary" <<"\t"
         << "repeat_times"       <<"\t"
         << "left_flank_binary"  <<"\t"
         << "right_flank_binary" <<"\t"
         << "repeat_unit_bases"         <<"\t"
         << "left_flank_bases"  <<"\t"
         << "right_flank_bases"   <<"\n";
}

// output and clear buffer
bool RefSeq::PouroutBuffer(std::ofstream &fout) {
    for (int i=0; i<homosBuffer.size(); i++) {
        fout << homosBuffer[i].chr       <<"\t"
             << homosBuffer[i].location  <<"\t"
             << int(homosBuffer[i].typeLen)   <<"\t"
             << homosBuffer[i].homoType  <<"\t"
             << homosBuffer[i].length    <<"\t"
             << homosBuffer[i].frontKmer <<"\t"
             << homosBuffer[i].endKmer   <<"\t";
        // append readable output
        homosBuffer[i].TransferString();
        fout << homosBuffer[i].transfer << "\n";
    }
    homosBuffer.clear();
    totalSites = 0;
}

// Load next chromosome
ref_loc_t RefSeq::LoadNextSeq(std::ifstream &fin) {
    char c;
    char ch[1000];
    std::string s;
    fin>>c;
    if (fin.eof()) return 0;
    _length = 0;
    // get name
    fin>>_name;
    fin.getline(ch, 1000);
    // get seq
    while (!fin.eof()) {
        fin>>c;  
        if (fin.eof()) break;
        fin.unget();
        if (c == '>') break;
        fin>>s;
        if (_length + s.size() >= param.max_dbseq_size) {
            if (s.size() > param.append_dbseq_size) {
                param.max_dbseq_size += (s.size() + 10);
            } else { 
                param.max_dbseq_size += param.append_dbseq_size; 
            }
            _seq.resize(param.max_dbseq_size);
        }
        copy(s.begin(), s.end(), _seq.begin() + _length);
        _length += s.size();
    }
    return _length;
}

// load one homo/microsate
bit8_t RefSeq::LoadOneSite( const std::string & chr, 
                            const std::string & seq, 
                            int loc, 
                            bit8_t hLen, 
                            bit16_t type, 
                            bit16_t len, 
                            HomoSite & oneSite ) {
    std::string::iterator h  = _seq.begin();
    std::string::iterator p0 = _seq.begin();
    std::string::iterator p1 = _seq.begin();
    // flank region
    bit16_t flankH = 0;
    bit16_t flankT = 0;
    int contextStart = loc - param.ContextLength;
    int contextStop  = loc + len*hLen + param.ContextLength;
    if ((contextStart < 0) || (contextStop >= _length)) return 0;
    for (unsigned char s=0; s<param.ContextLength; s++) {
        p0 = h+contextStart+s;
        p1 = h+loc+len*hLen+s;
        if ((alphabet[*p0] == 4) || (alphabet[*p1] == 4)) { return 0; }
        flankH<<=2;
        flankH|=alphabet[*p0];
        flankT<<=2;
        flankT|=alphabet[*p1];
    }
    oneSite.chr = chr;
    oneSite.location = loc;
    oneSite.typeLen = hLen;
    oneSite.homoType = type;
    oneSite.length = len;
    oneSite.frontKmer = flankH;
    oneSite.endKmer = flankT;

    return 1;
}

// test sites binary
void RefSeq::TestSitesBinary() {
    std::cout << "chromosome"   <<"\t"
              << "location"     <<"\t"
              << "sitelength"   <<"\t"
              << "sitecontent"  <<"\t"
              << "repeattimes"  <<"\t"
              << "front_flank"  <<"\t"
              << "tail_flank"   << std::endl;
    for (int i=0; i<homosBuffer.size(); i++) {
        std::cout << homosBuffer[i].chr       <<"\t"
                  << homosBuffer[i].location  <<"\t"
                  << int(homosBuffer[i].typeLen)   <<"\t"
                  << homosBuffer[i].homoType  <<"\t"
                  << homosBuffer[i].length    <<"\t"
                  << homosBuffer[i].frontKmer <<"\t"
                  << homosBuffer[i].endKmer   << std::endl; 
    }
}

// redable test
void RefSeq::TestSites() {
    for (int i=0; i<homosBuffer.size(); i++) {
        homosBuffer[i].TransferString();
        std::cout << homosBuffer[i].transfer << std::endl;
    }
}

// Scanning ref seqs
void RefSeq::DoScan( int length, 
                     const std::string &name,
                     unsigned int index ) {
    std::string::iterator h   = _seq.begin();
    std::string::iterator p   = _seq.begin();
    std::string::iterator p0  = _seq.begin();
    int location;
    int endLocation;
    int frontLocation;
    unsigned int endLength;
    unsigned int frontLength;
    unsigned int tempChar = 0;
    unsigned int hitLength = 0;
    int i  = 0;
    int i0 = 0;
    while (i < _length) {
        tempChar = alphabet[*p];
        // filter 'N'
        if ( tempChar == 4 ) { 
            ++i;
            ++p;
            continue; 
        }
        // record loc
        p0 = p;
        i0 = i;

        hitLength = 0;
        // homopolymer
        while (alphabet[*p] == tempChar) { 
            ++hitLength;
            ++i; 
            ++p; 
        }

        if (hitLength >= param.MininalHomoSize) {
            // record homopolymer
            HomoSite toneSite;
            if (LoadOneSite(name, _seq, i0, 1, tempChar, hitLength, toneSite)) {
                homosBuffer.push_back(toneSite);
                totalHomosites++;
                totalSites++;
                // pour out from buffer
                if (totalSites == param.bufSize){
                  PouroutBuffer(fout);       
                }
            }
        } else {
            bool ifHit = false;
            bool ifN   = false;
            unsigned short  k = 2;
            for (k=2; k<=param.MaxMicrosate; k++) {
                // traceback
                p = p0;
                i = i0;
                bit16_t s0 = 0;
                bit16_t s1 = 0;
                int m = 0;
                while ((m<k) && (i < _length)) {
                    if (alphabet[*p] == 4) {
                        ifN = true;
                        break;
                    }
                    s0<<=2;
                    s0|=alphabet[*p];
                    m++;
                    i++;
                    p++;
                }
                if ((ifN) || (i >= _length)) break;
                s1 = s0;
                hitLength = 0;
                while (s0 == s1) {
                    hitLength++;
                    if ((ifN) || (i >= _length)) break;
                    m = 0;
                    s0 = 0;
                    while ((m<k) && (i < _length)) {
                        if (alphabet[*p] == 4) {
                            ifN = true;
                            break;
                        }
                        s0<<=2;
                        s0|=alphabet[*p];
                        m++;
                        i++;
                        p++;
                    }
                }
                // record one 
                if (hitLength >= param.Repeats) {
                   // homo only ? no 
                   if (param.HomoOnly == 0) {
                       // load one microsate
                       HomoSite toneSite;
                       if (LoadOneSite(name, _seq, i0, k, s1, hitLength, toneSite)) {
                           homosBuffer.push_back(toneSite);
                           totalMicrosates++;
                           totalSites++;
                           // pour out from buffer
                           if (totalSites == param.bufSize) {
                               PouroutBuffer(fout);      
                           }
                       }
                   }
                   // load one microsate
                   ifHit = true;
                   break;
                }
            } // micro searching end

            // return to head
            if (ifHit) {
                p = p0 + k*hitLength;
                i = i0 + k*hitLength;
            } else {
                p = p0 + 1;
                i = i0 + 1;
            }
        } // end if
    }
}

// Scan homosites and windows
void RefSeq::ScanHomoAndMicrosate(std::ifstream &fin) {
    _seq.resize(param.max_dbseq_size);
    total_num = sum_length = 0;
    unsigned int index = 0;
    _count = 0; 
    while (LoadNextSeq(fin)) {
        // filtering little reference sequences
        if (_length < 20 ) continue;
        RefTitle r;
        r._name =_name;
        r._size = _length;
        title.push_back(r);
        // scan window and homopolymer site
        DoScan(_length, _name, index);
        std::cout << "scanning chomosome "
                  << _name << " done. "
                  << Cal_AllTime() << " secs passed" << std::endl;
        index++;
        total_num++;
        sum_length += _length;
    }
    _seq.clear();
}


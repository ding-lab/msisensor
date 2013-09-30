/*
 * sample.cpp for MSIsensor 
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

#include "sample.h"
#include "param.h"

extern Param paramd;

Sample::Sample()
    : outputPrefix( "test" )
    , output( NULL )
   //, outputPSomatic( NULL )
    , outputSomatic( NULL )
    , outputGermline( NULL )
    , outputDistribution( NULL )
    , numberOfSites( 0 )
    , precisionNumS( 2 )
    , precisionNumL( 5 ) 
    , numberOfDataPoints( 0 )
    , numberOfMsiDataPoints( 0 )
    , numberOftotalSites( 0 )
{
    //xxxx
    output.precision( precisionNumS );
   // outputPSomatic.precision( precisionNumL );
    outputSomatic.precision( precisionNumL );
    outputGermline.precision( precisionNumL );
    outputDistribution.precision( precisionNumL );

};

Sample::~Sample() {
    // xxxxx
};

void Sample::iniOutput( const std::string &gavePrefix ) { 
    if ( !gavePrefix.empty()  ) { outputPrefix = gavePrefix; }
    // init pour out result files
    output.open( outputPrefix.c_str() );
   // outputPSomatic.open( (outputPrefix + "_p_somatic").c_str() );
    outputSomatic.open( (outputPrefix + "_somatic").c_str() );
    outputGermline.open( (outputPrefix + "_germline").c_str() );
    outputDistribution.open( (outputPrefix + "_dis").c_str() );

    //if ( !output || !outputPSomatic || !outputSomatic || !outputGermline || !outputDistribution ) {
    if ( !output || !outputSomatic || !outputGermline || !outputDistribution ) {
        std::cerr <<"failed to open output files to write !"<< std::endl; 
        exit(1);
    }
}

void Sample::pourOutMsiScore() {
    output << "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%" << std::endl;
    output << numberOfDataPoints << "\t" 
           << numberOfMsiDataPoints << "\t" 
           << std::fixed 
           << (numberOfMsiDataPoints / (double)numberOfDataPoints) * 100.0 << std::endl;
}

void Sample::closeOutStream() {
    output.close();
    //outputPSomatic.close();
    outputSomatic.close();
    outputGermline.close();
    outputDistribution.close();
}

// FDR determination
void Sample::calculateFDR() {
    // sorting by p_value
    sort( totalSomaticSites.begin(), totalSomaticSites.end() );
    unsigned short rank = 1;
    // FDR calculation
    for ( std::vector< SomaticSite >::iterator _it = totalSomaticSites.begin(); _it != totalSomaticSites.end(); ++_it ) {
        //_it->PourOut();
        _it->FDR = _it->pValue * numberOfDataPoints / rank ;
        if ( _it->FDR > paramd.fdrThreshold ) { 
            rank++; 
            continue; 
        } else {
            _it->rank = rank;
            _it->somatic = true;
            numberOfMsiDataPoints++;
            rank++;
        }
    }
}

// report somatics && FDR
void Sample::pourOutSomaticFDR() {
    for ( std::vector< SomaticSite >::iterator _it = totalSomaticSites.begin(); _it != totalSomaticSites.end(); ++_it ) {
        if ( ! _it->somatic ) continue;
        outputSomatic << _it->chr << "\t"
                      << _it->location << "\t"
                      << _it->fbases << "\t"
                      << _it->length << "\t"
                      << _it->bases << "\t"
                      << _it->ebases << "\t"
                      << _it->diff << "\t"
                      << _it->pValue << "\t"
                      << _it->FDR << "\t"
                      << _it->rank << "\n";
    }
}

// verbose 
void Sample::VerboseInfo() {
    std::cerr << "\n*** Summary information ***\n\n"
              << "Number of total sites: " 
              << numberOftotalSites << "\n"
              << "Number of sites with enough coverage: "
              << numberOfDataPoints << "\n"
              << "Number of MSI sites: " 
              << numberOfMsiDataPoints << "\n"; 
}


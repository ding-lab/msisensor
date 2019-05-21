/*
 * cmds.cpp for MSIsensor
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

// System header files
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "cmds.h"

#ifndef VERSION
#define VERSION "v0.6"
#endif

int usage(void) {
    std::cerr<<"\n\n"
    << "Program: msisensor (homopolymer and miscrosatelite analysis using bam files)\n"
    << "Version: "<<VERSION<<"\n"
    << "Author: Beifang Niu && Kai Ye\n\n"
    << "Usage:   msisensor <command> [options]\n\n"
    << "Key commands:\n\n"
    << " scan            scan homopolymers and miscrosatelites\n"
    << " msi             msi scoring\n"
    << "\n\n";
    return 1; 
}

int main(int argc, char **argv) {
    try {
        if (argc < 2) {
            return usage();
        }
        if (strcmp(argv[1], "scan") == 0) {
            // scan homopolymer and microsate
            HomoAndMicrosateScan(argc-1, argv+1);
            return 0;
        } else if (strcmp(argv[1], "msi") == 0) {
            // distribution && msi scoring analysis 
            HomoAndMicrosateDisMsi(argc-1, argv+1);
            return 0;
        } else {
            std::cerr<<"ERROR: unrecognized command "<<argv[1]<<"\n";
            return 2;
        }
    } catch (const char *e) {
        std::cerr << e << std::endl;
        return 3;
    }

    return 0;
}


/*
 * bamreader.cpp for MSIsensor
 * Copyright (c) 2013 Beifang Niu && Kai Ye WUGSC All Rights Reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// System header files
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>

// Samtools header files
#include "bam.h"
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "ksort.h"

// Bam header files
#include "bamreader.h"
#include "distribution.h"

// Static function declaration
static int fetch_func_SR (const bam1_t * b1, void *data);

unsigned int g_InWinPlus = 0; 
unsigned int g_InWinMinus = 0; 
unsigned int g_NumReadScanned = 0;
unsigned int g_NumReadInWindow = 0; 
unsigned int g_CloseMappedPlus = 0; 
unsigned int g_CloseMappedMinus = 0; 

//init hash/maps for read pairing on the fly
KSORT_INIT_GENERIC (uint32_t) KHASH_MAP_INIT_STR (read_name, bam1_t *)
struct fetch_func_data_SR {
    fetch_func_data_SR () {
        LeftReads = NULL;
        read_to_map_qual = NULL;
        header = NULL;
        b1_flags = NULL;
        b2_flags = NULL;
        Tag = "";
        InsertSize = 0;
    }
    std::vector < SPLIT_READ > *LeftReads;
    khash_t (read_name) * read_to_map_qual;
    bam_header_t *header;
    flags_hit *b1_flags;
    flags_hit *b2_flags;
    std::string Tag;
    int InsertSize;
};

void parse_flags_and_tags (const bam1_t * b, flags_hit * flags) {
    const bam1_core_t *c = &b->core;
    char xt_code = '\0';
    int mf_code = 0, nm_code = 0, best_hits = 0;
    flags->unique = 0;
    flags->mapped = !(c->flag & BAM_FUNMAP);
    uint8_t *s = bam_aux_get (b, "XT");
    if (s != 0) {
        xt_code = bam_aux2A (s);
        if (xt_code == 'U') {
            flags->unique = 1;
        } else { 
            flags->unique = 0;
        }
    } else { flags->unique = 1; }
    s = NULL;
    s = bam_aux_get (b, "X0");
    if (s != 0) { 
        best_hits = bam_aux2i (s); 
    }
    s = NULL;
    s = bam_aux_get (b, "X1");
    if (s != 0) {
        int sub_hits = bam_aux2i (s);
        if (best_hits + sub_hits == 1) {
            flags->suboptimal = 0;
        } else { flags->suboptimal = 1; }
    } else { flags->suboptimal = 0; }
    if (xt_code == 'M' || mf_code == 130) {
        //short term fix to unset unique if the maq read was s-w mapped. 
        //bwa can't set U and M at once.
        flags->sw = 1;
    } else { flags->sw = 0; }
    s = NULL;
    s = bam_aux_get (b, "NM");
    if (s != 0) {
        nm_code = bam_aux2i (s);
        flags->edits = nm_code;
    } else {
        nm_code = 0;
        flags->edits = nm_code;
    }

    return;
}

/*

static int fetch_func_ALL (const bam1_t * b1, void *data) {
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    khash_t (read_name) * read_to_map_qual = (khash_t (read_name) *) data_for_bam->read_to_map_qual;

    flags_hit *b1_flags = data_for_bam->b1_flags;
    flags_hit *b2_flags = data_for_bam->b2_flags;

    //SPLIT_READ Temp_One_Read;
    const bam1_core_t *b1_core;
    bam1_t *b2;
    bam1_core_t *b2_core;
    b1_core = &b1->core;
    std::string read_name = bam1_qname (b1);
    khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));
    
    if (key == kh_end (read_to_map_qual)) {
        int ret = 0;
        key = kh_put (read_name, read_to_map_qual, strdup (bam1_qname (b1)), &ret);
        kh_value (read_to_map_qual, key) = bam_dup1 (b1);
        return 0;
    } else {
        b2 = bam_dup1 (kh_value (read_to_map_qual, key));
        bam_destroy1 (kh_value (read_to_map_qual, key));
        b2_core = &b2->core;
        //this seems stupid, but in order to manage the read names, necessary
        free ((char *) kh_key (read_to_map_qual, key));
        kh_del (read_name, read_to_map_qual, key);
        std::string c_sequence;
    }
    parse_flags_and_tags (b1, b1_flags);
    parse_flags_and_tags (b2, b2_flags);
    build_record(b1, b2, data, b1_flags);
    build_record(b2, b1, data, b2_flags);
    bam_destroy1 (b2);
    return 0;
}

*/

// new version fetch function for single end_reads
static int fetch_func_ALL (const bam1_t * b1, void *data) {
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    khash_t (read_name) * read_to_map_qual = (khash_t (read_name) *) data_for_bam->read_to_map_qual;

    flags_hit *b1_flags = data_for_bam->b1_flags;

    //SPLIT_READ Temp_One_Read;
    const bam1_core_t *b1_core;
    b1_core = &b1->core;
    parse_flags_and_tags (b1, b1_flags);
    build_record( b1, data, b1_flags );

    return 0;
}

/** 'ReadInRead' reads in reads from Pindel input file. */
bool ReadInBamReads( const char *bam_path, 
                     const std::string & FragName, 
                     unsigned start, 
                     unsigned end, 
                     std::vector < SPLIT_READ > & AllReads, 
                     std::string Tag) { 
    bamFile fp;
    fp = bam_open (bam_path, "r");
    assert (fp);
    bam_index_t *idx;
    idx = bam_index_load (bam_path); // load BAM index
    assert (idx);
    bam_header_t *header = bam_header_read (fp);
    bam_init_header_hash (header);
    assert (header);

    int tid;

    tid = bam_get_tid (header, FragName.c_str ());
    fetch_func_data_SR data;
    data.header = header;
    data.LeftReads = &AllReads;
    data.read_to_map_qual = NULL;
    data.read_to_map_qual = kh_init (read_name);

    flags_hit b1_flags, b2_flags;
    data.b1_flags = &b1_flags;
    data.b2_flags = &b2_flags;
    data.Tag = Tag;
    // std:: cout << " before bam_fetch " << std::endl;
    // give warning and abort if using dif refs
    if (tid == -1) {
        std::cout << "Program aborted: " << std::endl;
        std::cout << "Same reference genome file should be used in both 'msisensor scan' and 'msisensor msi' steps!!!"<<std::endl;
        exit(1);
    }
    bam_fetch (fp, idx, tid, start, end, &data, fetch_func_ALL);
    // std:: cout << " after bam_fetch " << std::endl;
    khint_t key;
    if (kh_size (data.read_to_map_qual) > 0) {
        for (key = kh_begin (data.read_to_map_qual); key != kh_end (data.read_to_map_qual); ++key) {
            if (kh_exist (data.read_to_map_qual, key)) {
                bam_destroy1 (kh_value (data.read_to_map_qual, key));
                free ((char *) kh_key (data.read_to_map_qual, key));
            }
        }
    }
    kh_clear (read_name, data.read_to_map_qual);
    kh_destroy (read_name, data.read_to_map_qual);
    bam_header_destroy (header);
    bam_index_destroy (idx);
    bam_close (fp);

    return true;
}

// new version build_record for single end reads 
void build_record( const bam1_t * current_read, 
                   void *data, 
                   const flags_hit *flag_current_read ) {   

    SPLIT_READ Temp_One_Read;
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    std::string Tag = (std::string) data_for_bam->Tag;

    const bam1_core_t *current_core;
    current_core = &current_read->core;
    // Determine sample name for read.
    // std::string c_sequence;
    uint8_t *s = bam1_seq (current_read);
    for (int i = 0; i <current_core->l_qseq; ++i) {
        Temp_One_Read.ReadSeq.append (1, bam_nt16_rev_table[bam1_seqi (s, i)]);
    }
    //std::cout<<Temp_One_Read.ReadSeq<<"\n";
    if ( !flag_current_read->mapped ) {
       Temp_One_Read.Mapped = false;
    } else {
       Temp_One_Read.MatchedRelPos = current_core->pos;
       Temp_One_Read.Mapped = true;
    }
    //std::cout<< Temp_One_Read.MatchedRelPos <<"\n";
    //
    // load one read
    data_for_bam->LeftReads->push_back(Temp_One_Read);

    return;
}

/*
void build_record( const bam1_t * current_read, 
                   const bam1_t * mate_read, 
                   void *data, 
                   const flags_hit *flag_current_read ) {   

    SPLIT_READ Temp_One_Read;
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    std::string Tag = (std::string) data_for_bam->Tag;

    const bam1_core_t *current_core;
    const bam1_core_t *mate_core;

    current_core = &current_read->core;
    mate_core = &mate_read->core;
    // Determine sample name for read.
    // std::string c_sequence;
    uint8_t *s = bam1_seq (current_read);
    for (int i = 0; i <current_core->l_qseq; ++i) {
        Temp_One_Read.ReadSeq.append (1, bam_nt16_rev_table[bam1_seqi (s, i)]);
    }
    //std::cout<<Temp_One_Read.ReadSeq<<"\n";
    if ( !flag_current_read->mapped ) {
       Temp_One_Read.Mapped = false;
    } else {
       Temp_One_Read.MatchedRelPos = current_core->pos;
       Temp_One_Read.Mapped = true;
    }
    //std::cout<< Temp_One_Read.MatchedRelPos <<"\n";
    //
    // load one read
    data_for_bam->LeftReads->push_back(Temp_One_Read);

    return;
}

*/

int32_t bam_cigar2len (const bam1_core_t * c, const uint32_t * cigar) {
    uint32_t k;
    int32_t l = 0;
    for (k = 0; k < c->n_cigar; ++k) {
        int op = cigar[k] & BAM_CIGAR_MASK;
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            l += cigar[k] >> BAM_CIGAR_SHIFT;
        }
        if (op == BAM_CDEL) { 
            l -= cigar[k] >> BAM_CIGAR_SHIFT; 
        }
   }
   return l;
}

int32_t bam_cigar2mismatch( const bam1_core_t *readCore, const uint32_t *cigar) {
    uint32_t cigarIndex;
    int32_t numberOfMismatches = 0;
    for (cigarIndex = 0; cigarIndex < readCore->n_cigar; ++cigarIndex) {
        int elementType = cigar[cigarIndex] & BAM_CIGAR_MASK;
        if (elementType != BAM_CMATCH) { 
            numberOfMismatches += cigar[cigarIndex] >> BAM_CIGAR_SHIFT; 
        }
    }
    return numberOfMismatches;
}


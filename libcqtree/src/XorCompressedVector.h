/*
 * XorCompressedVector.h
 *
 *  Created on: Jan 15, 2015
 *      Author: diegocaro
 */

#ifndef XORCOMPRESSEDVECTOR_H_
#define XORCOMPRESSEDVECTOR_H_

#include "XorCode.h"

#include <cstdio>
#include <cstdlib>

#include <libcdsBasics.h>
#include <Coder.h>

#include <utility>
#include <vector>

#include <coding_policy.h>

namespace cqtree_static {

class XorCompressedVector: public XorCode {
public:
    XorCompressedVector(MyCompressionCoder *cc, const vector<unsigned> &vi, unsigned sample_rate = 128) : sample_rate_(sample_rate) {
        assert(cc != NULL);

        N_ = vi.size();

        if (N_ == 0) {
            samples_ = 0;
            words_ = 0;
            coded_ = new unsigned[1];
            samples_and_pointers_ = new size_t[samples_*2+2];
            samples_and_pointers_[0] = 0;
            samples_and_pointers_[1] = 0;
            //printf("empty vector\n");
            return;
        }

        unsigned v1,v2;

        // (1) Calculate maximal value of samples and of deltas, and the size for coded_
        //unsigned max_sample_value = 0;
        words_ = 0;

        samples_ = 0;

        size_t jj = 0;

        vector <unsigned> inv;

        auto it = vi.begin();
        v1 = *it;
        for(; it != vi.end(); ++it ) {
            unsigned stream[1280];

            if (++jj%10000 == 0) fprintf(stderr,"First pass: %.2f\r", 100.0*jj/vi.size());
            v2 = *it;

                samples_++;

                inv.clear();
                v1 = v2;
                ++it;
                while(inv.size() < sample_rate_-1) {
                    if (it == vi.end()) { break; }
                    v2 = *it;

                    inv.push_back(v1 ^ v2);
                    v1 = v2;

                    ++it;

                    if (++jj%10000 == 0) fprintf(stderr,"First pass: %.2f\r", 100.0*jj/vi.size());

                }
                --it;

                //printf("encoding %u %u -> %u\n",v1,v2,v1^v2);
                if(inv.size()> 0)
                words_ += cc->Compress(inv.data(),stream,inv.size());

        }
        //fprintf(stderr, "First pass OK\n");
        //printf("bits used -> %lu (%lu words)\n",bits_, (bits_/32 + 1));

        // (2) Write sample values and xor with huffman
        samples_and_pointers_ = new size_t [samples_*2+2];
        coded_ = new unsigned [ words_ ];


        size_t *sp_it = samples_and_pointers_;
        size_t ptr = 0;
        jj = 0;

        it = vi.begin();
        v1 = *it;
        for(; it != vi.end(); ++it) {
            if (++jj%10000 == 0) fprintf(stderr,"Second pass: %.2f\r", 100.0*jj/vi.size() );
            v2 = *it;

                *sp_it = v2; ++sp_it;
                *sp_it = ptr; ++sp_it;

                inv.clear();
                ++it;
                v1 = v2;
                while(inv.size() < sample_rate_  -1) {
                    if (it == vi.end()) {  break; }
                    v2 = *it;
                    inv.push_back(v1 ^ v2);
                    v1 = v2;
                    it++;

                    if (++jj%10000 == 0) fprintf(stderr,"Second pass: %.2f\r", 100.0*jj/vi.size());

                }



                it--;
                //printf("encoding %u %u -> %u\n",v1,v2,v1^v2);
                if(inv.size()> 0)
                ptr += cc->Compress(inv.data(),&coded_[ptr],inv.size());




        }

        *sp_it = 0; ++sp_it; // initialize
        *sp_it = ptr; ++sp_it; // last entry
        //fprintf(stderr, "Second pass OK\n");

    }

    ~XorCompressedVector() {

    }

    virtual size_t getLength() {
        return N_;
    }

    size_t getSize() {
        return sizeof(XorCompressedVector) + (samples_+1)*2*sizeof(size_t) + (words_)*4;
    }

    unsigned at(MyCoder *cc, size_t i) {
        assert(i < N_);

        size_t idx = i/sample_rate_;

        unsigned mod = i - sample_rate_*idx; //idx % sample_rate

        unsigned v1 = samples_and_pointers_[idx*2];
        unsigned v2;
        size_t ptr = samples_and_pointers_[idx*2+1];

        unsigned sizeblock = sample_rate_-1;
        if (idx == samples_-1) { // last block
            sizeblock =  (N_-1) - idx*sample_rate_;
        }

        if(sizeblock > 0) {
            unsigned stream[1280];

            cc->Decompress(&coded_[ptr], stream, sizeblock);

            for(unsigned j = 0; j < mod; j++) {
                v2 = stream[j];
                v1 = v2 ^ v1;
            }
        }

        return v1;
    }

    //decode the block (things between sample rate), including the sample!
    // return the number of items in output
    size_t getBlock(MyCoder *cc, size_t idx, unsigned *output) {
        unsigned stream[1280];

        unsigned *o = output;

        unsigned v1 = samples_and_pointers_[idx*2];
        unsigned v2;

        size_t ptr = samples_and_pointers_[idx*2+1];

        *o = v1; o++;

        unsigned sizeblock = sample_rate_-1;

        if (idx == samples_-1) { // last block
            sizeblock =  (N_-1) - idx*sample_rate_;
        }

        if(sizeblock > 0) {
            cc->Decompress(&coded_[ptr], stream, sizeblock);

            for(unsigned j = 0; j < sizeblock; j++) {
                v2 = stream[j];
                v1 = v2 ^ v1;
                *o = v1; o++;
            }
        }
        return o - output;

    }

    // get values between [lo, hi)
    size_t getRange(MyCoder *cc, size_t lo, size_t hi, unsigned *output) {
        assert(lo < hi);
        assert(hi <= N_);

        size_t idxlo = lo/sample_rate_;
        size_t idxhi = hi/sample_rate_;

        unsigned modlo = lo - sample_rate_*idxlo;

        unsigned *o = output;

        // first block
        o += getBlock(cc,idxlo, o);

        // erase data between [0, modlo)
        std::move(output+modlo, o, output);
        o = o - modlo;

        // inner to last blocks
        for(size_t j = idxlo+1; j <= idxhi; j++) {
            o += getBlock(cc,j,o);
        }

        return hi - lo;
    }


    XorCompressedVector(ifstream &f) {
        uint type = loadValue<uint>(f);
        // TODO:throw an exception!
        if(type!=XORCOMPRESS_SAV) {
          abort();
        }

        loadValue(f,samples_);
        loadValue(f,words_);
        loadValue(f,sample_rate_);
        loadValue(f,N_);

        samples_and_pointers_ = cds_utils::loadValue<size_t>(f,samples_*2+2);
        coded_ = cds_utils::loadValue<unsigned>(f, words_);
    }

    virtual void save(ofstream &f) {
        uint wr = XORCOMPRESS_SAV;
        cds_utils::saveValue(f,wr);

        cds_utils::saveValue(f,samples_);
        cds_utils::saveValue(f,words_);
        cds_utils::saveValue(f,sample_rate_);
        cds_utils::saveValue(f,N_);

        cds_utils::saveValue(f,samples_and_pointers_,samples_*2+2);
        cds_utils::saveValue(f,coded_, words_);
    }

private:
    size_t samples_;
    size_t words_;

    unsigned sample_rate_;
    size_t *samples_and_pointers_;

    unsigned *coded_;
    unsigned N_; //number of elements in the input vector
};

class XorCompressedVectorBuilder: public XorCodeBuilder {
 public:
    XorCompressedVectorBuilder(char *codername, unsigned sample_rate):sample_rate_(sample_rate) {
        sprintf(codername_,"%s",codername);
        cc_ = new MyCompressionCoder(codername);
    }

    virtual XorCodeBuilder* copy() {
        return new XorCompressedVectorBuilder(codername_, sample_rate_);
    }

    virtual ~XorCompressedVectorBuilder() {

    }
    virtual XorCode* build(const vector<unsigned> &vi) {
        return new XorCompressedVector(cc_, vi, sample_rate_);
    }

    // not used
    virtual void updateDict(const std::vector<unsigned> &vi){

    }

    virtual MyCoder* getCoder() {
        return cc_;
    }

 private:
    unsigned sample_rate_;
    MyCompressionCoder *cc_;
    char codername_[255];
};


} /* namespace cqtree_static */

#endif /* XORCOMPRESSEDVECTOR_H_ */

/*
 * XorHuffmanVector.h
 *
 *  Created on: Jan 14, 2015
 *      Author: diegocaro
 */

#ifndef XORHUFFMANVECTOR_H_
#define XORHUFFMANVECTOR_H_

#include "XorCode.h"

#include <cstdio>
#include <cstdlib>

#include <libcdsBasics.h>
#include <Coder.h>

#include <utility>
#include <vector>



using namespace std;
using namespace cds_utils;
using namespace cds_static;
using namespace cqtree_utils;

namespace cqtree_static {

class XorHuffmanVector: public XorCode {
public:
    XorHuffmanVector(MyHuffmanCoder *hc, const vector<unsigned> &vi, unsigned sample_rate = 128) : sample_rate_(sample_rate) {
        assert(hc != NULL);

        N_ = vi.size();

        if (N_ == 0) {
            samples_ = 0;
            bits_huffmancoded_ = 0;
            huffmancoded_ = new unsigned[1];
            samples_and_pointers_ = new size_t[samples_*2+2];
            samples_and_pointers_[0] = 0;
            samples_and_pointers_[1] = 0;
            //printf("empty vector\n");
            return;
        }

        unsigned v1,v2;

        unsigned sample_now = 0;
        auto it = vi.begin();

        // (1) Calculate maximal value of samples and of deltas, and the size for huffmancoded_
        //unsigned max_sample_value = 0;
        bits_huffmancoded_ = 0;
        sample_now = 0;
        samples_ = 0;
        uint stream[10];
        size_t jj = 0;
        for(it = vi.begin(), v1 = *it; it != vi.end(); ++it ) {
            if (++jj%10000 == 0) fprintf(stderr,"First pass: %.2f\r", 100.0*jj/vi.size() );

            v2 = *it;
            if (sample_now == 0) {
                sample_now = sample_rate_;
                //if (max_sample_value < v2) max_sample_value = v2;
                samples_++;
            }
            else {
                //printf("encoding %u %u -> %u\n",v1,v2,v1^v2);
                bits_huffmancoded_ += hc->encode(v1 ^ v2,stream,0); //encoding xor
            }
            v1 = v2;
            --sample_now;
        }
        //fprintf(stderr, "First pass OK\n");
        //printf("bits used -> %lu (%lu words)\n",bits_huffmancoded_, (bits_huffmancoded_/32 + 1));

        // (2) Write sample values and xor with huffman
        samples_and_pointers_ = new size_t [samples_*2+2];
        huffmancoded_ = new unsigned [ (bits_huffmancoded_/32 + 1) ];

        sample_now = 0;

        size_t *sp_it = samples_and_pointers_;
        size_t ptr = 0;
        jj = 0;
        for(it = vi.begin(), v1 = *it; it != vi.end(); ++it) {
            if (++jj%10000 == 0) fprintf(stderr,"Second pass: %.2f\r", 100.0*jj/vi.size() );
            v2 = *it;
            if (sample_now == 0) {
                sample_now = sample_rate_;
                //if (max_sample_value < v2) max_sample_value = v2;

                *sp_it = v2; ++sp_it;
                *sp_it = ptr; ++sp_it;
            }
            else {
                ptr = hc->encode(v1 ^ v2, huffmancoded_,ptr); //encoding xor
            }
            v1 = v2;
             --sample_now;
        }

        *sp_it = 0; ++sp_it; // initialize
        *sp_it = ptr; ++sp_it; // last entry
        //fprintf(stderr, "Second pass OK\n");

    }

    ~XorHuffmanVector() {

    }

    virtual size_t getLength() {
        return N_;
    }

    virtual size_t getSize() {
        return sizeof(XorHuffmanVector) + (samples_+1)*2*sizeof(size_t) + (bits_huffmancoded_/32+1)*4;
    }

    virtual unsigned at(MyCoder *hc,size_t i) {
        assert(i < N_);

        size_t idx = i/sample_rate_;

        unsigned mod = i - sample_rate_*idx; //idx % sample_rate

        unsigned v1 = samples_and_pointers_[idx*2];
        unsigned v2;
        size_t ptr = samples_and_pointers_[idx*2+1];

        for(unsigned j = 0; j < mod; j++) {
             ptr = hc->decode(&v2,huffmancoded_,ptr);
             v1 = v2 ^ v1;
        }

        return v1;
    }

    //decode the block (things between sample rate), including the sample!
    // return the number of items in output
    virtual size_t getBlock(MyCoder *hc,size_t idx, unsigned *output) {
        unsigned *o = output;

        unsigned v1 = samples_and_pointers_[idx*2];
        unsigned v2;

        size_t ptr = samples_and_pointers_[idx*2+1];

        *o = v1; o++;

        unsigned mod = sample_rate_-1;

        if (idx == samples_-1) { // last block
            mod =  (N_-1) - idx*sample_rate_;
        }

        for(unsigned j = 0; j < mod; j++) {
             ptr = hc->decode(&v2,huffmancoded_,ptr);
             v1 = v2 ^ v1;
             *o = v1; o++;
        }
        return o - output;
    }

    // get values between [lo, hi)
    virtual size_t getRange(MyCoder *hc, size_t lo, size_t hi, unsigned *output) {
        assert(lo < hi);
        assert(hi <= N_);

        size_t idxlo = lo/sample_rate_;
        size_t idxhi = hi/sample_rate_;

        unsigned modlo = lo - sample_rate_*idxlo;

        unsigned *o = output;

        // first block
        o += getBlock(hc,idxlo, o);

        // erase data between [0, modlo)
        std::move(output+modlo, o, output);
        o = o - modlo;

        // inner to last blocks
        for(size_t j = idxlo+1; j <= idxhi; j++) {
            o += getBlock(hc,j,o);
        }

        return hi - lo;
    }

    XorHuffmanVector(ifstream &f) {
        uint type = loadValue<uint>(f);
        // TODO:throw an exception!
        if(type!=XORHUFFMAN_SAV) {
          abort();
        }

        loadValue(f,samples_);
        loadValue(f,bits_huffmancoded_);
        loadValue(f,sample_rate_);
        loadValue(f,N_);

        samples_and_pointers_ = cds_utils::loadValue<size_t>(f,samples_*2+2);
        huffmancoded_ = cds_utils::loadValue<unsigned>(f, bits_huffmancoded_/32 + 1);
    }

    virtual void save(ofstream &f) {
        uint wr = XORHUFFMAN_SAV;
        cds_utils::saveValue(f,wr);

        cds_utils::saveValue(f,samples_);
        cds_utils::saveValue(f,bits_huffmancoded_);
        cds_utils::saveValue(f,sample_rate_);
        cds_utils::saveValue(f,N_);

        cds_utils::saveValue(f,samples_and_pointers_,samples_*2+2);
        cds_utils::saveValue(f,huffmancoded_, bits_huffmancoded_/32 + 1);
    }


private:
    size_t samples_;
    size_t bits_huffmancoded_;

    unsigned sample_rate_;
    size_t *samples_and_pointers_;

    unsigned *huffmancoded_;
    unsigned N_; //number of elements in the input vector

};


class XorHuffmanVectorBuilder: public XorCodeBuilder {
 public:
    XorHuffmanVectorBuilder(unsigned sample_rate):sample_rate_(sample_rate) {
        fixedDict = false;
        hc_ = NULL;
    }
    virtual ~XorHuffmanVectorBuilder() {

    }

    virtual XorCodeBuilder* copy() {
        return new XorHuffmanVectorBuilder(sample_rate_);
    }

    virtual XorCode* build(const vector<unsigned> &vi) {
        getCoder();
        return new XorHuffmanVector(hc_, vi, sample_rate_);
    }

    // update the dictionary, with only non-sample values
    virtual void updateDict(const vector<unsigned> &vi) {
        assert(fixedDict == false);

        unsigned sample_now = 0;

        if (vi.size() == 0) return;
        unsigned v1,v2;
        size_t jj = 0;
        auto it = vi.begin();
        for(v1 = *it; it != vi.end(); ++it ) {
            if (++jj%10000 == 0) fprintf(stderr,"Updating dict: %.2f\r", 100.0*jj/vi.size() );

            v2 = *it;
            if (sample_now == 0) {
                sample_now = sample_rate_;
            }
            else {
                dict[v1 ^ v2] += 1;
            }
            v1 = v2;
            --sample_now;
        }
    }

    virtual MyCoder* getCoder() {
        fixedDict = true;

        if (hc_ == NULL) {
            if(dict.size() == 0) dict[0]=1;
            hc_ = new MyHuffmanCoder(dict);
        }

        return hc_;
    }


 private:
    unsigned sample_rate_;
    unordered_map<unsigned,unsigned> dict;
    bool fixedDict;
    MyHuffmanCoder *hc_;
};

}
#endif /* XORHUFFMANVECTOR_H_ */

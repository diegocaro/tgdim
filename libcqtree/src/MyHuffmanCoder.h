/*
 * MyHuffmanCoder.h
 *
 *  Created on: Jan 15, 2015
 *      Author: diegocaro
 */

#ifndef MYHUFFMANCODER_H_
#define MYHUFFMANCODER_H_

#include <Coder.h> //libcds
#include "MyCoder.h"

using namespace cds_static;

namespace cqtree_static {

class MyHuffmanCoder: public MyCoder {
 public:
    MyHuffmanCoder(unordered_map<unsigned,unsigned> &dict) {
        hc = new HuffmanCoder(dict);
    }

    virtual size_t encode(uint symb, uint * stream, size_t pos) const { return hc->encode(symb,stream,pos); };
    virtual size_t decode(uint * symb, uint *stream, size_t pos) const { return hc->decode(symb,stream,pos); };
    virtual size_t maxLength() const { return hc->maxLength(); };
    virtual size_t getSize() const {return hc->getSize(); }

    virtual void save(ofstream & fp) const {
        uint wr = MYHUFF_HDR;
        cds_utils::saveValue(fp,wr);
        hc->save(fp);};

    MyHuffmanCoder (ifstream & fp) {
        uint type = loadValue<uint>(fp);
        if(type != MYHUFF_HDR) {   //throw exception
            abort();
        }
        hc = HuffmanCoder::load(fp);
    }

    virtual size_t Compress(unsigned *input, unsigned *ouput, int n) {
    	UNUSED(input); UNUSED(ouput); UNUSED(n); return 0;
    };
    virtual size_t Decompress(unsigned *input, unsigned *ouput, int n) {
    	UNUSED(input); UNUSED(ouput); UNUSED(n); return 0;
    };

 private:
    HuffmanCoder *hc;
};


} /* namespace cqtree_static */

#endif /* MYHUFFMANCODER_H_ */

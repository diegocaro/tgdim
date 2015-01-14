/*
 * XorCode.h
 *
 *  Created on: Jan 14, 2015
 *      Author: diegocaro
 */

#ifndef XORCODE_H_
#define XORCODE_H_

#include "utils.h"
#include <vector>
// Libcds
#include <libcdsBasics.h>

//CODES to save structures
#define XORHUFFMAN_SAV 2
#define XORCOMPRESS_SAV 3


using namespace cqtree_utils;

namespace cqtree_static {

class XorCode {
 public:
    XorCode() {}
    virtual ~XorCode() {};

    virtual size_t getSize()=0;

    virtual unsigned at(size_t i)=0;

    //decode the block (things between sample rate), including the sample!
    // return the number of items in output
    virtual size_t getBlock(size_t idx, unsigned *output) = 0;

    // get values between [lo, hi)
    // return the number of items in output
    virtual size_t getRange(size_t lo, size_t hi, unsigned *output) = 0;

    virtual size_t getLength() = 0;

    virtual void save(ofstream &f) = 0;

    static XorCode* load(ifstream &f);
};


class XorCodeBuilder {
 public:
    XorCodeBuilder() {}
    virtual ~XorCodeBuilder() {}
    virtual XorCode* build(const std::vector<unsigned> &vi) = 0;
};

} /* namespace cqtree_static */

#include "XorHuffmanVector.h"

#endif /* XORCODE_H_ */

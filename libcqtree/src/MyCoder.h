/*
 * MyCoder.h
 *
 *  Created on: Jan 15, 2015
 *      Author: diegocaro
 */

#ifndef MYCODER_H_
#define MYCODER_H_

#include <Coder.h> //libcds
#include "utils.h"

#define MYHUFF_HDR 1
#define MYCOMPR_HDR 2


namespace cqtree_static {

class MyCoder {

 public:
        virtual ~MyCoder() {}

             virtual size_t encode(uint symb, uint * stream, size_t pos) const = 0;
             virtual size_t decode(uint * symb, uint *stream, size_t pos) const = 0;
             virtual size_t maxLength() const = 0;
             virtual size_t getSize() const = 0;

             virtual void save(ofstream & fp) const = 0;
             static  MyCoder * load(ifstream & fp);

             virtual size_t Compress(unsigned *input, unsigned *ouput, int n) = 0;
             virtual size_t Decompress(unsigned *input, unsigned *ouput, int n) = 0;
};

} /* namespace cqtree_static */
#include "MyHuffmanCoder.h"
#include "MyCompressionCoder.h"

#endif /* MYCODER_H_ */

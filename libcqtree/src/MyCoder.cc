/*
 * MyCoder.cc
 *
 *  Created on: Jan 15, 2015
 *      Author: diegocaro
 */

#include "MyCoder.h"

namespace cqtree_static {

MyCoder * MyCoder::load(ifstream & fp) {
    uint r = cds_utils::loadValue<unsigned>(fp);
    size_t pos = fp.tellg();
    fp.seekg(pos-sizeof(uint));
    switch(r) {
      case MYHUFF_HDR: return new MyHuffmanCoder(fp);
      case MYCOMPR_HDR: return new MyCompressionCoder(fp);
    }
    return NULL;
  }

} /* namespace cqtree_static */

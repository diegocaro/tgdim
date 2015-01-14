/*
 * XorCode.cc
 *
 *  Created on: Jan 14, 2015
 *      Author: diegocaro
 */

#include "XorCode.h"

namespace cqtree_static {

XorCode * XorCode::load(ifstream & fp) {
    uint r = cds_utils::loadValue<unsigned>(fp);
    size_t pos = fp.tellg();
    fp.seekg(pos-sizeof(uint));
    switch(r) {
      case XORHUFFMAN_SAV: return new XorHuffmanVector(fp);
      case XORCOMPRESS_SAV: return NULL;//new MXCompactQtreeFixed(fp);
    }
    return NULL;
  }

} /* namespace cqtree_static */

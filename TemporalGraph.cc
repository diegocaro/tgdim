/*
 * TemporalGraph.cc
 *
 *  Created on: Sep 2, 2014
 *      Author: diegocaro
 */

#include "TemporalGraph.h"



TemporalGraph* TemporalGraph::load(ifstream &fp) {
    uint r = loadValue<uint>(fp);
    size_t pos = fp.tellg();
    fp.seekg(pos-sizeof(uint));
    switch(r) {
      case TG_INTERV: return new IntervalContactGraph(fp);
      case TG_GROWTH: return new GrowingContactGraph(fp);
      case TG_POINT: return new PointContactGraph(fp);
    }
    return NULL;
  }


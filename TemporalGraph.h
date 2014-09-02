/*
 * TemporalGraph.h
 *
 *  Created on: Apr 22, 2014
 *      Author: diegocaro
 */

#ifndef TEMPORALGRAPH_H_
#define TEMPORALGRAPH_H_

#include <sys/types.h>

#include <CompactQtree.h>

enum bitseq {
  eRG,
  eRRR,
  eSD
};

enum TypeGraph {
  kInterval,
  kGrowth,
  kPoint
};

enum typeds {
  ePRBlack,
  ePRWhite,
  eMXDepth,
  eMXFixed,
};

struct opts {
  enum bitseq bs;  //bits for T bitmaps (tree)
  enum bitseq bb;  //bits for B bitmaps (leaves)
  enum typeds ds;  //type of compact qtree
  const char *outfile;

  /* Name of input file.  May be "-".  */
  const char *infile;

  enum TypeGraph typegraph;

  int k1;
  int k2;
  int lk1; //levels for k1
  int lki; //levels for ki
  int lf;  //levels for fixed mx
};


#define TG_INTERV 3
#define TG_GROWTH 4
#define TG_POINT 5


using namespace cqtree_static;



//TODO Add Interval operations
class TemporalGraph {
public:
    
    
	TemporalGraph() {
	    qt_ = NULL;
	    nodes_ = 0;
        edges_ = 0;
        lifetime_ = 0;
        contacts_ = 0;
	}

	static TemporalGraph* load(ifstream &f);
    
	~TemporalGraph() {
	    if (qt_!=NULL) {
	        delete qt_;
	    }
	    qt_ = NULL;
	};

    void setDs(CompactQtree *q){qt_ = q;};

    void setInfo(uint nodes, uint edges, uint lifetime, uint contacts) {
        nodes_ = nodes;
        edges_ = edges;
        lifetime_ = lifetime;
        contacts_ = contacts;
    }

    void setOpts(struct opts &opts){ opts_=opts;};


    virtual void save(ofstream &f)=0;

	//interface to query a temporal graph
	virtual void direct_point(uint u, uint t, uint *res)=0;
	virtual void direct_weak(uint u, uint tstart, uint tend, uint *res)=0;
	virtual void direct_strong(uint u, uint tstart, uint tend, uint *res)=0;

	virtual void reverse_point(uint v, uint t, uint *res)=0;
	virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res)=0;
	virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res)=0;

	virtual int edge_point(uint u, uint v, uint t)=0;
	virtual int edge_weak(uint u, uint v, uint tstart, uint tend)=0;
	virtual int edge_strong(uint u, uint v, uint tstart, uint tend)=0;
	virtual int edge_next(uint u, uint v, uint t)=0;

	virtual unsigned long snapshot(uint t)=0;


protected:
	uint nodes_;
	uint edges_;
	uint lifetime_;
	uint contacts_;
    
    CompactQtree *qt_;
    struct opts opts_;
    vector<Point<uint> > vp;
};



class IntervalContactGraph: public TemporalGraph {
public:
    class dirnei:public _Compare {
     public:
        dirnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[0] <= u_ && lo[0]+nk > u_ && lo[2] <= t_ && lo[3]+nk > t_+1) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[0] <= u_ && hi[0] > u_ && lo[2] <= t_ && hi[3] > t_+1) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    class revnei:public _Compare {
     public:
        revnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[1] <= u_ && lo[1]+nk > u_ && lo[2] <= t_ && lo[3]+nk > t_+1) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[1] <= u_ && hi[1] > u_ && lo[2] <= t_ && hi[3] > t_+1) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    IntervalContactGraph() {};

    ~IntervalContactGraph () {
        if (qt_!=NULL) {
                    delete qt_;
                }
                qt_ = NULL;
    };

    // Load & Save
    IntervalContactGraph(ifstream &f) {
        uint type = loadValue<uint>(f);
     // TODO:throw an exception!
     if(type!=TG_INTERV) {
       abort();
     }

     loadValue(f,nodes_);
     loadValue(f,edges_);
     loadValue(f,lifetime_);
     loadValue(f,contacts_);
     loadValue(f,opts_);
     qt_ = CompactQtree::load(f);
    }


   void save(ofstream &f) {
       uint wr = TG_INTERV;
           saveValue(f,wr);
           saveValue(f,nodes_);
           saveValue(f,edges_);
           saveValue(f,lifetime_);
           saveValue(f,contacts_);
           saveValue(f,opts_);

           qt_->save(f);
       }





    /// Interface

	virtual void direct_point(uint u, uint t, uint *res) {
        vp.clear();
        dirnei a(u, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][1];
                }
        #endif
	}

	virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
		//*res = 0;
		//cqtree->direct_weak(u,tstart,tend,res);
	}
	virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
		//*res = 0;
		//cqtree->direct_strong(u,tstart,tend,res);
	}

	virtual void reverse_point(uint v, uint t, uint *res) {
        vp.clear();
        dirnei a(v, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][0];
                }
        #endif
	}

	virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
		//*res = 0;
		//cqtree->reverse_weak(v,tstart,tend,res);
	}
	virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
		//*res = 0;
		//cqtree->reverse_weak(v,tstart,tend,res);
	}

	virtual int edge_point(uint u, uint v, uint t) {
		//return cqtree->edge_point(u, v, t);
	}

	virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
		//return cqtree->edge_weak(u,v,tstart,tend);
	}
	virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
		//return cqtree->edge_strong(u,v,tstart,tend);
	}
	virtual int edge_next(uint u, uint v, uint t) {
		//return cqtree->edge_next(u,v,t);
	}

	virtual unsigned long snapshot(uint t) {
		//return cqtree->snapshot(t);
	}

};


class GrowingContactGraph: public TemporalGraph {
public:
    class dirnei:public _Compare {
     public:
        dirnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[0] <= u_ && lo[0]+nk > u_ && lo[2] <= t_) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[0] <= u_ && hi[0] > u_ && lo[2] <= t_) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    class revnei:public _Compare {
     public:
        revnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[1] <= u_ && lo[1]+nk > u_ && lo[2] <= t_) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[1] <= u_ && hi[1] > u_ && lo[2] <= t_) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    GrowingContactGraph() {};

    ~GrowingContactGraph () {
        if (qt_!=NULL) {
                    delete qt_;
                }
                qt_ = NULL;
    };

    // Load & Save
    GrowingContactGraph(ifstream &f) {
        uint type = loadValue<uint>(f);
     // TODO:throw an exception!
     if(type!=TG_INTERV) {
       abort();
     }

     loadValue(f,nodes_);
     loadValue(f,edges_);
     loadValue(f,lifetime_);
     loadValue(f,contacts_);
     loadValue(f,opts_);
     qt_ = CompactQtree::load(f);
    }


   void save(ofstream &f) {
       uint wr = TG_INTERV;
           saveValue(f,wr);
           saveValue(f,nodes_);
           saveValue(f,edges_);
           saveValue(f,lifetime_);
           saveValue(f,contacts_);
           saveValue(f,opts_);

           qt_->save(f);
       }





    /// Interface

    virtual void direct_point(uint u, uint t, uint *res) {
        vp.clear();
        dirnei a(u, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][1];
                }
        #endif
    }

    virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->direct_weak(u,tstart,tend,res);
    }
    virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->direct_strong(u,tstart,tend,res);
    }

    virtual void reverse_point(uint v, uint t, uint *res) {
        vp.clear();
        dirnei a(v, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][0];
                }
        #endif
    }

    virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->reverse_weak(v,tstart,tend,res);
    }
    virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->reverse_weak(v,tstart,tend,res);
    }

    virtual int edge_point(uint u, uint v, uint t) {
        //return cqtree->edge_point(u, v, t);
    }

    virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
        //return cqtree->edge_weak(u,v,tstart,tend);
    }
    virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
        //return cqtree->edge_strong(u,v,tstart,tend);
    }
    virtual int edge_next(uint u, uint v, uint t) {
        //return cqtree->edge_next(u,v,t);
    }

    virtual unsigned long snapshot(uint t) {
        //return cqtree->snapshot(t);
    }

};

class PointContactGraph: public TemporalGraph {
public:
    class dirnei:public _Compare {
     public:
        dirnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[0] <= u_ && lo[0]+nk > u_ && lo[2] <= t_) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[0] <= u_ && hi[0] > u_ && lo[2] <= t_) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    class revnei:public _Compare {
     public:
        revnei(uint u, uint t): u_(u), t_(t) {};

    virtual bool operator()(const Point<uint> &lo, const uint nk) const {
               if (lo[1] <= u_ && lo[1]+nk > u_ && lo[2] <= t_) {
                   return true;
               }

               return false;
       }

       virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {

                  if (lo[1] <= u_ && hi[1] > u_ && lo[2] <= t_) {
                      return true;
                  }

                  return false;
          }
        uint u_;
        uint t_;
    };


    PointContactGraph() {};

    ~PointContactGraph () {
        if (qt_!=NULL) {
                    delete qt_;
                }
                qt_ = NULL;
    };

    // Load & Save
    PointContactGraph(ifstream &f) {
        uint type = loadValue<uint>(f);
     // TODO:throw an exception!
     if(type!=TG_INTERV) {
       abort();
     }

     loadValue(f,nodes_);
     loadValue(f,edges_);
     loadValue(f,lifetime_);
     loadValue(f,contacts_);
     loadValue(f,opts_);
     qt_ = CompactQtree::load(f);
    }


   void save(ofstream &f) {
       uint wr = TG_INTERV;
           saveValue(f,wr);
           saveValue(f,nodes_);
           saveValue(f,edges_);
           saveValue(f,lifetime_);
           saveValue(f,contacts_);
           saveValue(f,opts_);

           qt_->save(f);
       }





    /// Interface

    virtual void direct_point(uint u, uint t, uint *res) {
        vp.clear();
        dirnei a(u, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][1];
                }
        #endif
    }

    virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->direct_weak(u,tstart,tend,res);
    }
    virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->direct_strong(u,tstart,tend,res);
    }

    virtual void reverse_point(uint v, uint t, uint *res) {
        vp.clear();
        dirnei a(v, t);
        qt_->range(vp, a);
        res[0] = vp.size();
        #ifndef EXPERIMENTS
                for(size_t i=0; i < vp.size(); i++) {
                    res[i+1] = vp[i][0];
                }
        #endif
    }

    virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->reverse_weak(v,tstart,tend,res);
    }
    virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
        //*res = 0;
        //cqtree->reverse_weak(v,tstart,tend,res);
    }

    virtual int edge_point(uint u, uint v, uint t) {
        //return cqtree->edge_point(u, v, t);
    }

    virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
        //return cqtree->edge_weak(u,v,tstart,tend);
    }
    virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
        //return cqtree->edge_strong(u,v,tstart,tend);
    }
    virtual int edge_next(uint u, uint v, uint t) {
        //return cqtree->edge_next(u,v,t);
    }

    virtual unsigned long snapshot(uint t) {
        //return cqtree->snapshot(t);
    }

};

//template <class C>
//class PointContactGraph: public TemporalGraph {
//private:
//	C *ds;
//
//public:
//	PointContactGraph() {ds=NULL;};
//
//	PointContactGraph( std::vector<Point3d> &points, BitSequenceBuilder *bs) {
//		ds = new C(points,bs);
//	}
//
//	void direct_point(uint u, uint t, uint *res) {
//		*res = 0;
//		ds->template direct_interval<is_weak_point>(u,  t, t+1, res);
//	}
//
//	void direct_weak(uint u, uint tstart, uint tend, uint *res) {
//		*res = 0;
//		ds->template direct_interval<is_weak_point>(u, tstart, tend, res);
//		//FIXME Remove duplicates
//		//qsort(&res[1], *res, sizeof(unsigned int), compare);
//		//remove_duplicates(res);
//	};
//
//
//	void reverse_point(uint v, uint t, uint *res) {
//		*res = 0;
//		ds->template reverse_interval<is_weak_point>(v,  t, t+1, res);
//	}
//
//	void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
//		*res = 0;
//
//		ds->template reverse_interval<is_weak_point>(v, tstart, tend, res);
//		//FIXME Remove duplicates
//		//qsort(&res[1], *res, sizeof(unsigned int), compare);
//		//remove_duplicates(res);
//	};
//
//	int edge_point(uint u, uint v, uint t) {
//		return ds->template edge_interval<is_weak_point>(u, v,  t, t+1);
//	}
//
//	int edge_weak(uint u, uint v, uint tstart, uint tend){
//		return ds->template edge_interval<is_weak_point>(u, v, tstart, tend);
//	};
//
//
//	ulong snapshot(uint t) {
//		Point3d from; // all zeros
//		from[2] = t;
//
//		Point3d to;
//		to[0] = nodes_;
//		to[1] = nodes_;
//		to[2] = t+1;
//
//		ulong items=0;
//
//		ds->range(from, to, items);
//		//printf("items: %lu\n", items);
//		// return the number of active edges_
//		return items;
//
//	};
//
//	//TODO Falta edge_next
////
////	    int edge_next(uint u, uint v, uint t) {
////	        return edge_next(maxvalue_, u, v, t, 0, -1, depth_);
////	    }
////
////	    int edge_next(uint n, uint u, uint v, uint tq, uint a, ulong z, uint l);
////
//};
//
//template <class C>
//class GrowingContactGraph: public TemporalGraph {
//private:
//	C *ds;
//
//public:
//
//	GrowingContactGraph() {ds=NULL;};
//
//	GrowingContactGraph( std::vector<Point3d> &points, BitSequenceBuilder *bs) {
//		ds = new C(points,bs);
//	}
//
//	void direct_point(uint u, uint t, uint *res) {
//		*res = 0;
//		ds->template direct_interval<is_weak_grow>(u,  t, t+1, res);
//	}
//
//	void reverse_point(uint v, uint t, uint *res) {
//		*res = 0;
//		ds->template reverse_interval<is_weak_grow>(v,  t, t+1, res);
//	}
//
//	int edge_point(uint u, uint v, uint t) {
//		return ds->template edge_interval<is_weak_point>(u, v,  t, t+1);
//	}
//
//	ulong snapshot(uint t) {
//		Point3d from; // all zeros
//
//		Point3d to;
//		to[0] = nodes_;
//		to[1] = nodes_;
//		to[2] = t+1;
//
//		ulong items=0;
//
//		ds->range(from, to, items);
//		//printf("items: %lu\n", items);
//		// return the number of active edges_
//		return items;
//
//	};
//
//	//TODO Falta edge_next
//
//	//	    int edge_next(uint u, uint v, uint t) {
//	//	        return edge_next(maxvalue_, u, v, t, 0, -1, depth_);
//	//	    }
//	//
//	//	    int edge_next(uint n, uint u, uint v, uint tq, uint a, ulong z, uint l);
//
//};

#endif /* TEMPORALGRAPH_H_ */

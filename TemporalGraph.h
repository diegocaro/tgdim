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

enum TypeGraph {
  kInterval,
  kGrowth,
  kPoint,
  kIntervalGrowth
};

enum typeds {
  ePRBlack,
  ePRWhite,
  eMXDepth,
  eMXFixed,
  ePRB2Black
};

struct opts {
  BitSequenceBuilder *bs;  //bits for T bitmaps (tree)
  BitSequenceBuilder *bb;  //bits for B bitmaps (leaves)
  BitSequenceBuilder *bc;  //bits for C bitmaps (count leaves)

  enum typeds ds;  //type of compact qtree
  enum TypeGraph typegraph;

  const char *outfile;

  /* Name of input file.  May be "-".  */
  const char *infile;

  /* string parameters*/
  char params_char[255];
  map<string,string> params;

  int k1;
  int k2;
  int F;
  int lk1;  //levels for k1
  int lki;  //levels for ki
  int lf;  //levels for fixed mx
};

#define TG_INTERV 3
#define TG_GROWTH 4
#define TG_POINT 5
#define TG_INTERVPRO 6

using namespace cqtree_static;

void readflags(struct opts *opts, const char *flags);
BitSequenceBuilder* getBSBuilder(string e);
uint *getBitmap(BitSequence *bs);


class UpdateMXCompactQtree: public MXCompactQtree {
    public:
  void updateBitmaps(BitSequenceBuilder *bt) {
      uint *btemp;
      size_t len;

      for(int i = 0; i < depth_; i++ ) {
          fprintf(stderr, "New bitmaps at level %d\n",i);

          btemp = getBitmap(T_[i]);
          len = T_[i]->getLength();


          delete T_[i];
          T_[i] = bt->build(btemp, len);
          delete btemp;
      }
  }
};

class UpdatePRBCompactQtree: public PRBCompactQtree {
    public:
  void updateBitmaps(BitSequenceBuilder *bt,BitSequenceBuilder *bb) {
      uint *btemp;
      size_t len;

      for(int i = 0; i < depth_; i++ ) {
          fprintf(stderr, "New bitmaps at level %d\n",i);

          btemp = getBitmap(T_[i]);
          len = T_[i]->getLength();


          delete T_[i];
          T_[i] = bt->build(btemp, len);
          delete btemp;



          btemp = getBitmap(B_[i]);
          len = B_[i]->getLength();

          delete B_[i];
          B_[i] = bb->build(btemp, len);
          delete btemp;
      }
  }
};

  class UpdatePRB2CompactQtree: public PRB2CompactQtree {
  public:
    void updateBitmaps(BitSequenceBuilder *bt,BitSequenceBuilder *bb, BitSequenceBuilder *bc) {
        uint *btemp;
        size_t len;

        for(int i = 0; i < depth_; i++ ) {
            fprintf(stderr, "New bitmaps at level %d\n",i);

            btemp = getBitmap(T_[i]);
            len = T_[i]->getLength();


            delete T_[i];
            T_[i] = bt->build(btemp, len);
            delete btemp;



            btemp = getBitmap(B_[i]);
            len = B_[i]->getLength();

            delete B_[i];
            B_[i] = bb->build(btemp, len);
            delete btemp;

            btemp = getBitmap(C_[i]);
            len = C_[i]->getLength();

            delete C_[i];
            C_[i] = bb->build(btemp, len);
            delete btemp;
        }
    }
};


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
    if (qt_ != NULL) {
      delete qt_;
    }
    qt_ = NULL;
  }
  ;

  void setDs(CompactQtree *q) {
    qt_ = q;
  }
  ;

  void setInfo(uint nodes, uint edges, uint lifetime, uint contacts) {
    nodes_ = nodes;
    edges_ = edges;
    lifetime_ = lifetime;
    contacts_ = contacts;
  }

//  void setOpts(struct opts &opts) {
//    opts_ = opts;
//  }
//  ;

  void stats() {
    qt_->stats_space();
  }

  void print_leaves() {


    if (dynamic_cast<PRBCompactQtree*>(qt_)) {
        ((PRBCompactQtree* )qt_ )->print_leaves();
    }
    else if (dynamic_cast<PRB2CompactQtree*>(qt_)) {
            ((PRB2CompactQtree* )qt_ )->print_leaves();
        }
    else {
        fprintf(stderr, "Input is not a PRB data structure\n");
    }
  }

  void updateBitmaps(BitSequenceBuilder *bt, BitSequenceBuilder *bb, BitSequenceBuilder *bc) {
      if (dynamic_cast<PRBCompactQtree*>(qt_)) {
          ((UpdatePRBCompactQtree* )qt_ )->updateBitmaps(bt,bb);
      }
      else if (dynamic_cast<PRB2CompactQtree*>(qt_)) {
              ((UpdatePRB2CompactQtree* )qt_ )->updateBitmaps(bt,bb,bc);
          }
      else if (dynamic_cast<MXCompactQtree*>(qt_)) {
         ((UpdateMXCompactQtree* )qt_ )->updateBitmaps(bt);
     }
      else {
          fprintf(stderr, "CompactQtree class cannot be updated to new bitmaps.\n");
          return;
      }
  }

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
  virtual unsigned long contacts()=0;

 protected:
  uint nodes_;
  uint edges_;
  uint lifetime_;
  uint contacts_;

  CompactQtree *qt_;
  //struct opts opts_;
  vector<Point<uint> > vp;
};

class IntervalContactGraph : public TemporalGraph {
 public:
  IntervalContactGraph() {
  }
  ;

  ~IntervalContactGraph() {
    if (qt_ != NULL) {
      delete qt_;
    }
    qt_ = NULL;
  }
  ;

  // Load & Save
  IntervalContactGraph(ifstream &f) {
    uint type = loadValue<uint>(f);
    // TODO:throw an exception!
    if (type != TG_INTERV) {
      abort();
    }

    loadValue(f, nodes_);
    loadValue(f, edges_);
    loadValue(f, lifetime_);
    loadValue(f, contacts_);
    //loadValue(f, opts_);
    qt_ = CompactQtree::load(f);
  }

  void save(ofstream &f) {
    uint wr = TG_INTERV;
    saveValue(f, wr);
    saveValue(f, nodes_);
    saveValue(f, edges_);
    saveValue(f, lifetime_);
    saveValue(f, contacts_);
    //saveValue(f, opts_);

    qt_->save(f);
  }

  /// Interface

  virtual void direct_point(uint u, uint t, uint *res) {
    direct_weak(u, t, t+1, res);
  }

  virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = u;
      from[1] = 0;
      from[2] = 0;
      from[3] = tstart+1;

      to[0] = u+1;
      to[1] = nodes_;
      to[2] = tend;
      to[3] = lifetime_;



#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
        //*res = vp.size();
      for (size_t i = 0; i < vp.size(); i++) {
        res[i + 1 + *res] = vp[i][1];
      }
      *res += vp.size();
#endif
  }
  virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = u;
      from[1] = 0;
      from[2] = tstart;
      from[3] = tstart+1;

      to[0] = u+1;
      to[1] = nodes_;
      to[2] = tend-1;
      to[3] = tend;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
        //*res = vp.size();
      for (size_t i = 0; i < vp.size(); i++) {
        res[i + 1 + *res] = vp[i][1];
      }
      *res += vp.size();
#endif
  }

  virtual void reverse_point(uint v, uint t, uint *res) {
    reverse_weak(v, t, t+1, res);
  }

  virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = 0;
      from[1] = v;
      from[2] = 0;
      from[3] = tstart+1;

      to[0] = nodes_;
      to[1] = v+1;
      to[2] = tend;
      to[3] = lifetime_;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
        //*res = vp.size();
      for (size_t i = 0; i < vp.size(); i++) {
        res[i + 1 + *res] = vp[i][0];
      }
      *res += vp.size();
#endif
  }
  virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = 0;
      from[1] = v;
      from[2] = tstart;
      from[3] = tstart+1;

      to[0] = nodes_;
      to[1] = v+1;
      to[2] = tend-1;
      to[3] = tend;

#ifdef EXPERIMENTS
      *res = qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
        //*res = vp.size();
      for (size_t i = 0; i < vp.size(); i++) {
        res[i + 1 + *res] = vp[i][0];
      }
      *res += vp.size();
#endif
  }

  virtual int edge_point(uint u, uint v, uint t) {
    return edge_weak(u,v,t,t+1);
  }

  virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = u;
      from[1] = v;
      from[2] = 0;
      from[3] = tstart+1;

      to[0] = u+1;
      to[1] = v+1;
      to[2] = tend;
      to[3] = lifetime_;

      qt_->range(from,to,vp,true);
      if (vp.size() > 0) return 1;
      return 0;
  }
  virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
    vp.clear();

    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = u;
      from[1] = v;
      from[2] = 0;
      from[3] = tstart+1;

      to[0] = u+1;
      to[1] = v+1;
      to[2] = tend;
      to[3] = lifetime_;

      qt_->range(from,to,vp,true);
      if (vp.size() > 0) return 1;
      return 0;
  }
  virtual int edge_next(uint u, uint v, uint t) {
    //return cqtree->edge_next(u,v,t);
    return 0;
  }

  virtual unsigned long snapshot(uint t) {
      vp.clear();
    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = 0;
      from[1] = 0;
      from[2] = 0;
      from[3] = t+1;

      to[0] = nodes_;
      to[1] = nodes_;
      to[2] = t+1;
      to[3] = lifetime_;

      return qt_->range(from,to,vp,false);
  }

  virtual unsigned long contacts() {
      vp.clear();
    Point<uint> from(4);
      Point<uint> to(4);

      from[0] = 0;
      from[1] = 0;
      from[2] = 0;
      from[3] = 0;

      to[0] = nodes_;
      to[1] = nodes_;
      to[2] = lifetime_;
      to[3] = lifetime_;

      return qt_->range(from,to,vp,false);
  }


};

class GrowingContactGraph : public TemporalGraph {
 public:
  GrowingContactGraph() {
  }
  ;

  ~GrowingContactGraph() {
    if (qt_ != NULL) {
      delete qt_;
    }
    qt_ = NULL;
  }
  ;

  // Load & Save
  GrowingContactGraph(ifstream &f) {
    uint type = loadValue<uint>(f);
    // TODO:throw an exception!
    if (type != TG_GROWTH) {
      abort();
    }

    loadValue(f, nodes_);
    loadValue(f, edges_);
    loadValue(f, lifetime_);
    loadValue(f, contacts_);
    //loadValue(f, opts_);
    qt_ = CompactQtree::load(f);
  }

  void save(ofstream &f) {
    uint wr = TG_GROWTH;
    saveValue(f, wr);
    saveValue(f, nodes_);
    saveValue(f, edges_);
    saveValue(f, lifetime_);
    saveValue(f, contacts_);
    //saveValue(f, opts_);

    qt_->save(f);
  }

  /// Interface

  virtual void direct_point(uint u, uint t, uint *res) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = u;
      from[1] = 0;
      from[2] = 0;


      to[0] = u+1;
      to[1] = nodes_;
      to[2] = t+1;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
      //*res = vp.size();
          for (size_t i = 0; i < vp.size(); i++) {
            res[i + 1 + *res] = vp[i][1];
          }
          *res += vp.size();
      #endif
  }
  virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
    direct_point(u,tstart,res);
  }

  virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
    direct_point(u,tend,res);
  }

  virtual void reverse_point(uint v, uint t, uint *res) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = 0;
      from[1] = v;
      from[2] = 0;


      to[0] = nodes_;
      to[1] = v+1;
      to[2] = t+1;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
      //*res = vp.size();
          for (size_t i = 0; i < vp.size(); i++) {
            res[i + 1 + *res] = vp[i][0];
          }
          *res += vp.size();
      #endif
  }

  virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
    reverse_point(v,tend,res);
  }
  virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
    reverse_point(v,tstart,res);
  }

  virtual int edge_point(uint u, uint v, uint t) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = u;
      from[1] = v;
      from[2] = 0;


      to[0] = nodes_;
      to[1] = nodes_;
      to[2] = t+1;

      qt_->range(from,to,vp,true);
      if (vp.size() > 0 ) return 1;
      return 0;
  }

  virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
    return edge_point(u,v,tend);
  }
  virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
    return edge_point(u,v,tstart);
  }
  virtual int edge_next(uint u, uint v, uint t) {
    //return cqtree->edge_next(u,v,t);
    return 0;
  }

  virtual unsigned long snapshot(uint t) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = 0;
      from[1] = 0;
      from[2] = 0;


      to[0] = nodes_;
      to[1] = nodes_;
      to[2] = t+1;

      return qt_->range(from,to,vp,0);
  }

  virtual unsigned long contacts() {
      vp.clear();

      Point<uint> from(3);
        Point<uint> to(3);

        from[0] = 0;
        from[1] = 0;
        from[2] = 0;


        to[0] = nodes_;
        to[1] = nodes_;
        to[2] = lifetime_;

      return qt_->range(from,to,vp,false);
  }

};

class PointContactGraph : public TemporalGraph {
 public:

  PointContactGraph() {
  }
  ;

  ~PointContactGraph() {
    if (qt_ != NULL) {
      delete qt_;
    }
    qt_ = NULL;
  }
  ;

  // Load & Save
  PointContactGraph(ifstream &f) {
    uint type = loadValue<uint>(f);
    // TODO:throw an exception!
    if (type != TG_POINT) {
      abort();
    }

    loadValue(f, nodes_);
    loadValue(f, edges_);
    loadValue(f, lifetime_);
    loadValue(f, contacts_);
    //loadValue(f, opts_);
    qt_ = CompactQtree::load(f);
  }

  void save(ofstream &f) {
    uint wr = TG_POINT;
    saveValue(f, wr);
    saveValue(f, nodes_);
    saveValue(f, edges_);
    saveValue(f, lifetime_);
    saveValue(f, contacts_);
    //saveValue(f, opts_);

    qt_->save(f);
  }

  /// Interface

  virtual void direct_point(uint u, uint t, uint *res) {
    direct_weak(u,t,t+1,res);
  }

  virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = u;
      from[1] = 0;
      from[2] = tstart;


      to[0] = u+1;
      to[1] = nodes_;
      to[2] = tend;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
      //*res = vp.size();
          for (size_t i = 0; i < vp.size(); i++) {
            res[i + 1 + *res] = vp[i][0];
          }
          *res += vp.size();
      #endif
  }
  virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
    if (tstart+1 == tend) {
      direct_point(u,tstart,res);
    }
    else {
      *res = 0;
    }

  }

  virtual void reverse_point(uint v, uint t, uint *res) {
    reverse_weak(v,t,t+1,res);
  }

  virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = 0;
      from[1] = v;
      from[2] = tstart;


      to[0] = nodes_;
      to[1] = v+1;
      to[2] = tend;

#ifdef EXPERIMENTS
      *res += qt_->range(from,to,vp,false);
#else
      qt_->range(from,to,vp,true);
      //*res = vp.size();
          for (size_t i = 0; i < vp.size(); i++) {
            res[i + 1 + *res] = vp[i][1];
          }
          *res += vp.size();
      #endif
  }
  virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
    if (tstart+1 == tend) {
      reverse_point(v,tstart,res);
    }
    else {
      *res = 0;
    }
  }

  virtual int edge_point(uint u, uint v, uint t) {
   return edge_weak(u,v,t,t+1);
  }

  virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = u;
      from[1] = v;
      from[2] = tstart;


      to[0] = u+1;
      to[1] = v+1;
      to[2] = tend;

      qt_->range(from,to,vp,true);
      if (vp.size() > 0) return 1;
      else return 0;
  }
  virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
    if (tstart+1 == tend) {
      return edge_point(u,v,tstart);
    }
    return 0;
  }
  virtual int edge_next(uint u, uint v, uint t) {
    //return cqtree->edge_next(u,v,t);
    return 0;
  }

  virtual unsigned long snapshot(uint t) {
    vp.clear();

    Point<uint> from(3);
      Point<uint> to(3);

      from[0] = 0;
      from[1] = 0;
      from[2] = t;


      to[0] = nodes_;
      to[1] = nodes_;
      to[2] = t+1;

     return qt_->range(from,to,vp,false);

  }

  virtual unsigned long contacts() {
      vp.clear();

      Point<uint> from(3);
        Point<uint> to(3);

        from[0] = 0;
        from[1] = 0;
        from[2] = 0;


        to[0] = nodes_;
        to[1] = nodes_;
        to[2] = lifetime_;

       return  qt_->range(from,to,vp,false);

  }

};

// @TODO this can be improved..., is ugly the way I mixed up the classes
// to get an itnerval and growing graph
class IntervalContactGraphImproved : public TemporalGraph {
 public:
    IntervalContactGraphImproved() {
  }
  ;

  ~IntervalContactGraphImproved() {
      if (past_!=NULL) delete past_;
      if (curr_!=NULL) delete curr_;
      past_=NULL;
      curr_=NULL;
  }
  ;

  // Load & Save
  IntervalContactGraphImproved(ifstream &f) {
    uint type = loadValue<uint>(f);
    // TODO:throw an exception!
    if (type != TG_INTERVPRO) {
      abort();
    }

    loadValue(f, nodes_);
    loadValue(f, edges_);
    loadValue(f, lifetime_);
    loadValue(f, contacts_);
    //loadValue(f, opts_);

    past_ = new IntervalContactGraph(f);
    curr_ = new GrowingContactGraph(f);

    qt_ = NULL;
  }

  void save(ofstream &f) {
    uint wr = TG_INTERVPRO;
    saveValue(f, wr);
    saveValue(f, nodes_);
    saveValue(f, edges_);
    saveValue(f, lifetime_);
    saveValue(f, contacts_);
    //saveValue(f, opts_);

    past_->save(f);
    curr_->save(f);
  }

 void setGraphs(IntervalContactGraph *past, GrowingContactGraph *curr) {
    past_ = past;
    curr_ = curr;

}

  /// Interface

  virtual void direct_point(uint u, uint t, uint *res) {
    past_->direct_point(u,t,res);
    curr_->direct_point(u,t,res);
  }

  virtual void direct_weak(uint u, uint tstart, uint tend, uint *res) {
      past_->direct_weak(u,tstart,tend,res);
      curr_->direct_weak(u,tstart,tend,res);
  }
  virtual void direct_strong(uint u, uint tstart, uint tend, uint *res) {
      past_->direct_strong(u,tstart,tend,res);
      curr_->direct_strong(u,tstart,tend,res);
  }

  virtual void reverse_point(uint v, uint t, uint *res) {
      past_->reverse_point(v,t,res);
      curr_->reverse_point(v,t,res);
  }

  virtual void reverse_weak(uint v, uint tstart, uint tend, uint *res) {
    past_->reverse_weak(v,tstart,tend,res);
    curr_->reverse_weak(v,tstart,tend,res);
  }
  virtual void reverse_strong(uint v, uint tstart, uint tend, uint *res) {
    past_->reverse_strong(v,tstart,tend,res);
    curr_->reverse_strong(v,tstart,tend,res);
  }

  virtual int edge_point(uint u, uint v, uint t) {
    return (past_->edge_point(u,v,t) || curr_->edge_point(u,v,t));
  }

  virtual int edge_weak(uint u, uint v, uint tstart, uint tend) {
      return (past_->edge_weak(u,v,tstart,tend) || curr_->edge_weak(u,v,tstart,tend));
  }
  virtual int edge_strong(uint u, uint v, uint tstart, uint tend) {
      return (past_->edge_strong(u,v,tstart,tend) || curr_->edge_strong(u,v,tstart,tend));
  }
  virtual int edge_next(uint u, uint v, uint t) {
    //return cqtree->edge_next(u,v,t);
    return 0;
  }

  virtual unsigned long snapshot(uint t) {
      return past_->snapshot(t) + curr_->snapshot(t);
  }

  virtual unsigned long contacts() {
      return past_->contacts() + curr_->contacts();
  }

  IntervalContactGraph *past_;
  GrowingContactGraph *curr_;

};


#endif /* TEMPORALGRAPH_H_ */

/*
 * MXCompactQtreeFixed.h
 *
 *  Created on: Apr 15, 2014
 *      Author: diegocaro
 */

#ifndef MXCOMPACTQTREEFIXED_H_
#define MXCOMPACTQTREEFIXED_H_

#include "CompactQtree.h"
#include "utils.h"

// Libcds
#include <libcdsBasics.h>
#include <BitSequenceBuilder.h>

#include <vector>

using namespace cqtree_utils;
using namespace cds_static;

namespace cqtree_static {

//const uint _K = 2;  //fixed K in K^n trees :)



//typename int or long
//template<typename T>
class MXCompactQtreeFixed:public CompactQtree {
  struct Node {
      size_t lo, hi;
      int level;
  };

    struct Less {
            Less(const MXCompactQtreeFixed& cl) : c(cl) {
//                printf("depth: %d\n",c.depth_);
//                for(int i = 0; i < c.depth_; i++) {
//                    printf("k_[%d] = %d\n",i,c.k_[i]);
//                    printf("nk_[%d] = %d\n",i,c.nk_[i]);
//                }


            }
            bool operator () (const Point<uint> &x, const Point<uint> &y) {
                for(int i=0; i < c.depth_; i++) {
                    for(int j=c.rangedim_by_level_[i].first; j < c.rangedim_by_level_[i].second; j++) {
                        int a=(x[j]/c.nk_[i])%c.k_[i];
                        int b=(y[j]/c.nk_[i])%c.k_[i];

                        if (a!=b) return (a < b);
                    }
                }
                return false;
            }

            const MXCompactQtreeFixed& c;
        };

 public:
    MXCompactQtreeFixed() {
        __setdefaultvalues();
    }

    MXCompactQtreeFixed(std::vector<Point<uint> > &points, BitSequenceBuilder *bs, BitSequenceBuilder *bb, int k1=2, int k2=2, int max_level_k1=0, int max_levels_ki=0, int max_levels_fixed=1);

    virtual ~MXCompactQtreeFixed() {
        for(size_t i=0; i < T_.size(); i++) {
                    delete T_[i];
                }
                T_.clear();

                delete B_;

                for(int i=0; i < num_dims_; i++) {
                    delete leaves_[i];
                }

                delete [] leaves_;

                __setdefaultvalues();

    }

    virtual void stats_space() const;

    void range(Point<uint> &p, size_t z, int level, Point<uint> &from, Point<uint> &to, vector<Point<uint> > &vpall,size_t &items, bool pushval);
    virtual size_t range(Point<uint> &from, Point<uint> &to,vector<Point<uint> > &vpall, bool pushval=true){
      Point<uint> p(num_dims_);
      size_t items=0;
      range(p, -1, -1, from,to,vpall,items,pushval);
      return items;
    }


    void all(Point<uint> p, size_t z, int level, vector<Point<uint> > &vpall);
    virtual void all(vector<Point<uint> > &vpall) {
        Point<uint> p(num_dims_);
        all(p, -1, -1, vpall);
    }

    // fast calculation of morton code, just when k1=k2=2
    inline int code_k_eq_2(const Point<uint> &p, int level) const {
       int r=0; for(int i=0; i < num_dims_; i++) {
         r |= (( (p[i] >> (depth_+levels_fx_ -level-1)) & 1 ) << (num_dims_-i-1) ); // ((p[i]/2^(depth_-level-1)%2)^(num_dims_-i-1)
      }
      return r;
      }

    // Return the morton code, considering different arity on each level of the tree
    int code(const Point<uint> &p, int level) const {
      assert(level < depth_);


        int r = 0;
        int c;

        if (false == is_interleaved_) {
            if ((levels_k1_ == 0 && k2_ == 2 ) || (k1_==2 && k2_==2)) {
                return code_k_eq_2(p,level);
            }

          for (int i = 0; i < num_dims_; i++) {
            c = (p[i] / nk_[level]) % k_[level];
            r += c * kpower_per_level_dim_[level*num_dims_+i];
          }
        }
        else {
            //caso k1 y k2, con k2 interleaved

          for (int i = rangedim_by_level_[level].first;
                  i < rangedim_by_level_[level].second; i++) {
              c = (p[i] / nk_[level]) % k_[level];
              r += c * kpower_per_level_dim_[level*num_dims_+i];
          }
        }

        return r;
    }


    MXCompactQtreeFixed(ifstream & f) {
      uint type = loadValue<uint>(f);
                  // TODO:throw an exception!
                  if(type!=MXQFIX_SAV) {
                    abort();
                  }
        loadValue(f,levels_k1_);
        loadValue(f,levels_k2_);
        loadValue(f,levels_ki_);
        loadValue(f,levels_fx_);

        loadValue(f,k1_);
        loadValue(f,k2_);
        loadValue(f,ki_);

        loadValue(f,depth_);
        loadValue(f,maxvalue_);

        loadValue(f,num_dims_);

        loadValue(f,max_children_);
        loadValue(f,is_interleaved_);

       loadVector(f,k_);
       loadVector(f,children_);
       loadVector(f,nk_);
       loadVector(f,rangedim_by_level_);
       loadVector(f, dims_);
       loadVector(f, nodes_by_level_);
       loadVector(f,kpower_per_level_dim_);

        for(int i = 0; i < depth_; i++) {
            T_.push_back(BitSequence::load(f));
            //printf("s T_[%d].size() = %lu %p\n", i,T_[i]->getLength(),T_[i]);
        }

         B_ = BitSequence::load(f);

         leaves_ = new Array*[num_dims_];
         for(int i = 0; i < num_dims_; i++) {
             leaves_[i] = new Array(f);
         }


        }

    virtual void save(ofstream &f) const {
      uint wr = MXQFIX_SAV;
            cds_utils::saveValue(f,wr);
        cds_utils::saveValue(f,levels_k1_);
        cds_utils::saveValue(f,levels_k2_);
        cds_utils::saveValue(f,levels_ki_);
        cds_utils::saveValue(f,levels_fx_);

        cds_utils::saveValue(f,k1_);
        cds_utils::saveValue(f,k2_);
        cds_utils::saveValue(f,ki_);

        cds_utils::saveValue(f,depth_);
        cds_utils::saveValue(f,maxvalue_);

        cds_utils::saveValue(f,num_dims_);

        cds_utils::saveValue(f,max_children_);
        cds_utils::saveValue(f,is_interleaved_);

        saveVector(f,k_);
        saveVector(f,children_);
        saveVector(f,nk_);
        saveVector(f,rangedim_by_level_);
        saveVector(f, dims_);
        saveVector(f, nodes_by_level_);
        saveVector(f, kpower_per_level_dim_);


        for(int i = 0; i < depth_; i++) {
            T_[i]->save(f);
        }

        B_->save(f);

          for(int i = 0; i < num_dims_; i++) {
                      leaves_[i]->save(f);
                  }
    }


 protected:
    void __setdefaultvalues() {
        depth_ = 0;
        max_children_ = 0;
        maxvalue_ = 0;
        levels_ki_ = 0;
        is_interleaved_ = false;


        leaves_ = NULL;
        B_ = NULL;

        T_.clear();
    }

    size_t rank(const std::vector<Point<uint> > &vp, int key, uint level,
                long lo, long hi) const {
        if (hi < lo)
            return lo;  //this is because we are using size_t instead of int as indices
        long mid = lo + (hi - lo) / 2;

        //printf("%u %u\n", vp[mid].getMorton(level), code(vp[mid], depth_-level-1));

        //printf("lo: %u - hi: %u\n", lo, hi);
        // printf("vp_[%u].get(%u) = %u\n", mid, level, vp_[mid].get(level));
        if (code(vp[mid],level) == key
                && (mid == 0 || (mid > 0 && code(vp[mid-1],level) != key))) {
            //printf("found in %u\n", mid);
            return mid;
        } else if (key <= code(vp[mid],level)) {
            return rank(vp, key, level, lo, mid - 1);
        } else {  // if (key > keys[mid])
            return rank(vp, key, level, mid + 1, hi);
        }
    }


    void get_stats(const std::vector<Point<uint> > &vp);

    void create(const std::vector<Point<uint> > &vp,
                BitSequenceBuilder *bs, BitSequenceBuilder *bb);
    void build(const std::vector<Point<uint> > &vp,
                    BitSequenceBuilder *bs, BitSequenceBuilder *bb);

    int levels_k1_;
    int levels_k2_; //int virtual_depth_k2_;
    int levels_ki_; //int virtual_depth_ki_;// levels interleaved
    int levels_fx_; //level fixed

    int k1_;
    int k2_;
    int ki_; //k used for interleaved, for now is fixed to 2

    int depth_;
    uint maxvalue_;

    int num_dims_; //number of dimensions

    vector<int> k_;
    vector<int> children_; //children per node, per level
    int max_children_;

    vector<uint> nk_; //k^d per level (o cuanto aporta cada nivel a cada dimension)

    bool is_interleaved_;

    vector<pair<int,int> > rangedim_by_level_; // range of dimensions used by level
    vector<int> dims_; //number of dimensions processed by level

    vector<size_t> nodes_by_level_; // nodes by level, or, bitmap size per level :D
    vector<BitSequence*> T_;  // tree bitmaps, one per level

    BitSequence *B_;  // bitmap to count points per gray node
    Array **leaves_; //encoded path for each leaf

    // to replace mypow and get code values
       vector<uint> kpower_per_level_dim_;
};

} /* namespace cqtree_static */

#endif /* MXCOMPACTQTREE_H_ */

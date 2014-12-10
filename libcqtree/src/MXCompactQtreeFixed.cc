/*
 * MXCompactQtreeFixed.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: diegocaro
 */

#include "MXCompactQtreeFixed.h"
#include "utils.h"
#include "debug.h"
#include "bitvector.h"
#include <algorithm> //sort
#include <vector>
#include <queue>
#include <climits>

namespace cqtree_static {

MXCompactQtreeFixed::MXCompactQtreeFixed(vector<Point<uint> > &vp,
                                         BitSequenceBuilder *bs, BitSequenceBuilder *bb, int k1, int k2, int max_level_k1, int max_levels_ki, int max_levels_fixed)
: levels_fx_(max_levels_fixed), k1_(k1), k2_(k2) {

    __setdefaultvalues();  //set default values for this class

    //preconditions
    assert(max_levels_fixed > 0);

    if (max_levels_ki > 0) {
        assert(max_levels_ki %2 == 0); //number of levels interleaved must be even

        is_interleaved_ = true;
        ki_ = 2; //by default k for interleaved section is 2

        levels_ki_ = max_levels_ki;
    }


    Point<uint> pmax_dim(max(vp));
    num_dims_ = pmax_dim.num_dims();

    uint maxval;
    maxval = max(pmax_dim);

    assert((uint)max_levels_fixed < mylog(2,maxval));

    //depth_ = ceil(log((double) maxval) / log((double) _K));

    uint max_k1=0;
    //max_k1 = (uint) pow((double) k1, (double) max_level_k1);
    max_k1 = mypow(k1, max_level_k1);
    LOG("max_k1: %u\n",max_k1);
    levels_k1_ = max_level_k1;
    levels_k2_ = 0;
    LOG("resto1: %u\n",maxval/max_k1+1);

    // k1 must have available levels
    assert(max_k1 < maxval);

    uint max_fixed = mypow(2,max_levels_fixed); //max value stored on deeper levels
    LOG("max_fixed: %u\n",max_fixed);
    if (false == is_interleaved_) {

        // the number of levels availables must be at least the number of levels for fixed depth
        int lvl_availables_ = mylog(2, maxval/max_k1+1);
        assert(lvl_availables_ >= max_levels_fixed);

        //levels using k2_ are automatically calculated (taking account the fixed depth)
        levels_k2_ = mylog(k2, maxval/max_k1/max_fixed+1);

        maxvalue_ = mypow(k1_,levels_k1_) * mypow(k2_,levels_k2_) * mypow(2,max_levels_fixed);
        //virtual_depth_k2_ = mylog(k2, maxvalue_);
        //printf("virtual_depth_k2_: %d\n",virtual_depth_k2_);

        depth_ = levels_k1_ + levels_k2_;

        dims_.insert(dims_.end(), depth_, num_dims_);
        rangedim_by_level_.insert(rangedim_by_level_.end(), depth_, make_pair(0, num_dims_));


        uint m=maxvalue_;
        for(int i=0; i < levels_k1_; i++) {
            m = m / k1_;
            nk_.push_back(m);
            //virtual_depth_.push_back(mylog(k1_,m+1));
        }

        for(int i = 0; i < levels_k2_; i++) {
            m = m / k2_;
            nk_.push_back(m);
            //virtual_depth_.push_back(mylog(k2_,m+1));
        }
        assert(m == max_fixed);


        for(int i=0; i < levels_k1_; i++) {
            for(int j=0; j < num_dims_; j++) {
              kpower_per_level_dim_.push_back(mypow(k1_, num_dims_ - j - 1));
            }
        }

        for(int i = 0; i < levels_k2_; i++) {
            for(int j=0; j < num_dims_; j++) {
              kpower_per_level_dim_.push_back(mypow(k2_, num_dims_ - j - 1));
            }
        }
    }
    else {


        // the number of levels availables must be at least the number of levels for ki_
        int lvl_availables_ki = mylog(ki_, maxval/max_k1/max_fixed+1);
        assert(lvl_availables_ki >= levels_ki_/2);

        int max_ki = mypow(ki_, levels_ki_/2);

//        printf("ki_: %d\n",ki_);
//        printf("l/2: %d\n", levels_ki_);
//        printf("max_ki: %d",max_ki);
//        printf("resto2: %u\n", (maxval/max_k1)/max_ki/max_fixed);

        //levels using k2_ are automatically calculated (taking account the fixed depth and interleaved levels)
        levels_k2_ = mylog(k2, (maxval/max_k1)/max_ki/max_fixed+1);
        LOG("levels for k2: %d\n", levels_k2_);

        maxvalue_ = mypow(k1_,levels_k1_) * mypow(k2_,levels_k2_) * mypow(ki_, levels_ki_/2) * mypow(2,max_levels_fixed);;

        depth_ = levels_k1_ + levels_k2_ + levels_ki_;

        dims_.insert(dims_.end(), levels_k1_, num_dims_);
        dims_.insert(dims_.end(), levels_k2_, num_dims_);

        rangedim_by_level_.insert(rangedim_by_level_.end(), levels_k1_, make_pair(0, num_dims_));
        rangedim_by_level_.insert(rangedim_by_level_.end(), levels_k2_, make_pair(0, num_dims_));

        int cdim = ceil((float)num_dims_/ki_);
        //int fdim = floor((float)num_dims_/2);

        LOG("cdim: %d\n",cdim);
        for(int i = 0; i < levels_ki_; i++) {
            if (i%2 == 0) {
                dims_.push_back(cdim);
                rangedim_by_level_.push_back(make_pair(0, cdim));
            }
            else {
                dims_.push_back(num_dims_-cdim);
                rangedim_by_level_.push_back(make_pair(cdim, num_dims_));
            }
        }


        uint m=maxvalue_;
        for(int i=0; i < levels_k1_; i++) {
            m = m / k1_;
            nk_.push_back(m);
            //virtual_depth_.push_back(mylog(k1_,m));
        }

        for(int i=0; i < levels_k2_; i++) {
            m = m / k2_;
            nk_.push_back(m);
            //   virtual_depth_.push_back(mylog(k2_,m));
        }

        for(int i = 0; i < levels_ki_; i++) {
            if (i%2 == 0) {
                m = m / ki_;
            }
            nk_.push_back(m);
            //virtual_depth_.push_back(mylog(ki_,m));
        }

        assert(m == max_fixed);


        for(int i=0; i < levels_k1_; i++) {
          for(int j=0; j < num_dims_; j++) {
            kpower_per_level_dim_.push_back(mypow(k1_, num_dims_ - j - 1));
          }
        }

        for(int i=0; i < levels_k2_; i++) {
          for(int j=0; j < num_dims_; j++) {
            kpower_per_level_dim_.push_back(mypow(k2_, num_dims_ - j - 1));
          }
        }

        for(int i = 0; i < levels_ki_; i++) {
          if (i%2 == 0) {

            for(int j=0; j < cdim; j++) {
              kpower_per_level_dim_.push_back(mypow(ki_, cdim - j - 1));
            }

            for(int j=cdim; j < num_dims_; j++) {
              kpower_per_level_dim_.push_back(0);
            }
          }
          else {
            for(int j=0; j < cdim; j++) {
              kpower_per_level_dim_.push_back(0);
            }

            for(int j=0; j < num_dims_-cdim; j++) {
              kpower_per_level_dim_.push_back(mypow(ki_, num_dims_-cdim - j - 1));
            }
          }
        }

    }

    LOG("depth: %d\n", depth_);
    LOG("levelsk2: %d\n", levels_k2_);
    k_.insert(k_.begin(), levels_k1_, k1);
    k_.insert(k_.end(), levels_k2_, k2);
    k_.insert(k_.end(), levels_ki_, ki_);

    max_children_ = 0;
    for(int i =0; i < depth_; i++) {
        //children_.push_back((uint) pow((double) k_[i], (double) pmax_dim.num_dims()));
        children_.push_back(mypow(k_[i], dims_[i]));

        max_children_ = max(max_children_,children_[i]);


        LOG("children[%d]: %d\n",i, children_[i]);
        LOG("dims[%d]: %d\n",i,dims_[i]);
        LOG("rangedim_by_level[%d]: [%d,%d)\n",i,rangedim_by_level_[i].first,rangedim_by_level_[i].second);
        //printf("virtual_depth_[%d]: %d\n",i, virtual_depth_[i]);
        LOG("nk_[%d]: %d\n",i, nk_[i]);
    }

    LOG("depth: %d\n",depth_);

    //childs_ = (uint) pow((double) _K, (double) pmax_dim.num_dims());
    //maxvalue_ = (uint) pow((double) _K, (double) depth_);
    //maxvalue_ = mypow(k1_,levels_k1_) * mypow(k2_,levels_k2_);
    LOG("max value: %u\n",maxvalue_);

    printf("Sorting...\n");
    // we need to sort, this is why we cannot set "vector<T> &vp" to be const.
    if ((levels_k1_ == 0 && k2_ == 2 ) || (k1_==2 && k2_==2)) {
        bool is_sorted = true;
        for (size_t i = 1; i < vp.size(); i++) {
            if (Point<uint>::cmpmorton(vp[i - 1], vp[i]) == false) {
                is_sorted = false;
                break;
            }
        }

        if (!is_sorted) {
            std::sort(vp.begin(), vp.end(), Point<uint>::cmpmorton);  //sorting using morton only for k=2
        }
        else {
            printf("Sequence is sorted!\n");
        }
    }
    else {
    // but as we are using different values of k and more modifications...
      std::sort(vp.begin(), vp.end(),Less(*this));
    }
    printf("Done!\n");

    // remove duplicated elements
    vector<Point<uint> >::iterator last = unique(vp.begin(), vp.end());
    vp.erase(last, vp.end());


//    get_stats(vp);
//
//    size_t treebits=0;
//    for(int i = 0; i < depth_; i++) {
//        printf("T[%d]: %lu\n", i, nodes_by_level_[i]);
//        treebits += nodes_by_level_[i];
//    }
//    printf("Total size T: %lu\n", treebits);
//
//    create(vp, bs, bb);

    items_ = vp.size();

    build(vp,bs,bb);
    size_t treebits=0;
    for(int i = 0; i < depth_; i++) {
        printf("T[%d]: %lu\n", i, nodes_by_level_[i]);
        treebits += nodes_by_level_[i];
    }
    printf("Total size T: %lu\n", treebits);
}

void MXCompactQtreeFixed::get_stats(const std::vector<Point<uint> > &vp) {
    struct Node z;
    queue<struct Node> q;
    size_t r[max_children_ + 1];

    z.level = 0;
    z.lo = 0;
    z.hi = vp.size();
    q.push(z);

    vector<uint> entries;
    entries.insert(entries.begin(), depth_, 0);

    nodes_by_level_.insert(nodes_by_level_.begin(), depth_, 0);

    //int ktreenodes=0;
    while (!q.empty()) {
        z = q.front();
        q.pop();

        r[0] = z.lo;
        r[children_[z.level]] = z.hi;
        //ktreenodes+=1;
        nodes_by_level_[z.level] += children_[z.level];

        for (int j = 1; j < children_[z.level]; j++)
            //r[j] = rank(vp, j, depth_- z.level - 1, z.lo, z.hi - 1);
            r[j] = rank(vp, j, z.level, z.lo, z.hi - 1);

        //printf("r: %d %d %d %d\n", r[0], r[1], r[2], r[3]);

        for (int j = 1; j < children_[z.level] + 1; j++) {
            if (r[j] - r[j - 1] > 0) {
                struct Node t;
                t.lo = r[j - 1];

                t.hi = r[j];

                t.level = z.level + 1;

                if (z.level < depth_ - 1)
                    q.push(t);
                //printf("set bit %u\n", (ktreenodes-1)*childs_ + j-1);

            }
        }
    }
}

void MXCompactQtreeFixed::create(const std::vector<Point<uint> > &vp,
                                 BitSequenceBuilder *bs, BitSequenceBuilder *bb) {
    struct Node z;
    queue<struct Node> q;
    size_t r[max_children_ + 1];

    //int leaves[depth_];
    //for(int i=0; i < depth_; i++) leaves[i]=0;
    //int childs[depth_];
    //for(int i=0; i < depth_; i++) childs[i]=0;

    z.level = 0;
    z.lo = 0;
    z.hi = vp.size();

    q.push(z);

    //ulong ktreenodes = 0;

    vector<uint*> treebits;
    for(int i = 0; i < depth_; i++) {
        treebits.push_back(new uint[nodes_by_level_[i] / W + 1]);
        fill_n(treebits[i], nodes_by_level_[i] / W + 1, 0); //fill with zeroes
    }

    uint *gapbits = new uint[vp.size() / W + 1];
    fill_n(gapbits, vp.size() / W + 1, 0); //fill with zeroes

    vector<size_t> curr_nodes;
    curr_nodes.insert(curr_nodes.begin(), depth_, 0);

    while (!q.empty()) {
        z = q.front();
        q.pop();

        //ktreenodes += 1;

        r[0] = z.lo;
        r[children_[z.level]] = z.hi;

        for (int j = 1; j < children_[z.level]; j++)
            //r[j] = rank(vp, j, depth_- z.level - 1, z.lo, z.hi - 1);
            r[j] = rank(vp, j, z.level, z.lo, z.hi - 1);

        //printf("r: %d %d %d %d %d\n", r[0], r[1], r[2], r[3], r[4]);

        for (int j = 1; j < children_[z.level] + 1; j++) {
            if (r[j] - r[j - 1] > 0) {
                struct Node t;
                t.lo = r[j - 1];

                t.hi = r[j];

                t.level = z.level + 1;

                if (z.level < depth_ - 1)
                    q.push(t);

                //at leaf level
                if(z.level == depth_ -1) {
                    //printf("set bit B_ %lu\n", t.hi-1);
                    cds_utils::bitset(gapbits, t.hi-1);
                }

                //printf("set bit T %lu\n", (ktreenodes-1)*childs_ + j-1);
                //printf("set bit T_[%d] %lu\n", z.level, curr_nodes[z.level] + j - 1);
                cds_utils::bitset(treebits[z.level], curr_nodes[z.level] + j - 1);
            }
        }

        curr_nodes[z.level] += children_[z.level];

    }

    //building bitmaps for T and cleaning
    for(int i = 0; i < depth_; i++) {
        T_.push_back( bs->build(treebits[i], nodes_by_level_[i]));
        delete [] treebits[i];
    }


    B_ = bb->build(gapbits, vp.size());
    delete [] gapbits;

    uint maxval = nk_[depth_-1];
    leaves_ = new Array*[num_dims_];

    for(int i = 0; i < num_dims_; i++) {
        leaves_[i] = new Array(vp.size(), maxval-1);
    }

    for(int j = 0; j < num_dims_; j++) {
        for(size_t i=0; i < vp.size(); ++i) {
            leaves_[j]->setField(i, vp[i][j] % maxval);
        }
    }

}




void MXCompactQtreeFixed::build(const std::vector<Point<uint> > &vp,
                                 BitSequenceBuilder *bs, BitSequenceBuilder *bb) {
    struct Node z;
    queue<struct Node> q;
    size_t r[max_children_ + 1];

    //int leaves[depth_];
    //for(int i=0; i < depth_; i++) leaves[i]=0;
    //int childs[depth_];
    //for(int i=0; i < depth_; i++) childs[i]=0;

    z.level = 0;
    z.lo = 0;
    z.hi = vp.size();

    q.push(z);

    //ulong ktreenodes = 0;

    bitvector *bt = new bitvector(children_[0]);

    uint *gapbits = new uint[vp.size() / W + 1];
    fill_n(gapbits, vp.size() / W + 1, 0); //fill with zeroes

    vector<size_t> curr_nodes;
    curr_nodes.insert(curr_nodes.begin(), depth_, 0);

    int curr_level = 0;
    while (!q.empty()) {
        z = q.front();
        q.pop();

        //ktreenodes += 1;

        r[0] = z.lo;
        r[children_[z.level]] = z.hi;

        if (curr_level != z.level) {
            curr_level = z.level;
            T_.push_back( bs->build((uint*)bt->data(), bt->length()));

            delete bt;

            bt = new bitvector(T_[curr_level-1]->countOnes()*children_[curr_level]);
        }

        for (int j = 1; j < children_[z.level]; j++)
            //r[j] = rank(vp, j, depth_- z.level - 1, z.lo, z.hi - 1);
            r[j] = rank(vp, j, z.level, z.lo, z.hi - 1);

        //printf("r: %d %d %d %d %d\n", r[0], r[1], r[2], r[3], r[4]);

        for (int j = 1; j < children_[z.level] + 1; j++) {
            if (r[j] - r[j - 1] > 0) {
                struct Node t;
                t.lo = r[j - 1];

                t.hi = r[j];

                t.level = z.level + 1;

                if (z.level < depth_ - 1)
                    q.push(t);

                //at leaf level
                if(z.level == depth_ -1) {
                    //printf("set bit B_ %lu\n", t.hi-1);
                    cds_utils::bitset(gapbits, t.hi-1);
                }

                //printf("set bit T %lu\n", (ktreenodes-1)*childs_ + j-1);
                //printf("set bit T_[%d] %lu\n", z.level, curr_nodes[z.level] + j - 1);
                bt->bitset(curr_nodes[z.level] + j - 1);
            }
        }

        curr_nodes[z.level] += children_[z.level];

    }


    T_.push_back( bs->build((uint*)bt->data(), bt->length()));
    delete bt;

    B_ = bb->build(gapbits, vp.size());
    delete [] gapbits;

    uint maxval = nk_[depth_-1];
    leaves_ = new Array*[num_dims_];

    for(int i = 0; i < num_dims_; i++) {
        leaves_[i] = new Array(vp.size(), maxval-1);
    }

    for(int j = 0; j < num_dims_; j++) {
        for(size_t i=0; i < vp.size(); ++i) {
            leaves_[j]->setField(i, vp[i][j] % maxval);
        }
    }

    //updates nodes_by_level_
    for(int i=0; i < depth_; i++) {
        nodes_by_level_.push_back(curr_nodes[i]);
    }
}


//
//void MXCompactQtreeFixed::range(uint n, Point<uint> from, Point<uint> to, Point<uint> p, size_t z, int level, size_t &items) {
//    if (level != -1 && (p+n <= from || p > to) ) {
//        return;
//    }
//
//    size_t y = 0;
//
//    //printf("checking T[%d] %lu\n",level,z);
//
//
//    if (level == depth_-1) {
//        if (T_[level]->access(z) == 0) return;
//
//        size_t lo,hi;
//
//        y = T_[level]->rank1(z-1);
//        printf("y: %lu\n",y);
//
//        if (y == 0) {
//            lo = 0;
//            hi = B_->select1(y+1);
//        }
//        else {
//            lo = B_->select1(y)+1;
//            hi = B_->select1(y+1);
//        }
//
//        printf("lo: %lu\n",lo);
//        printf("hi: %lu\n",hi);
//
//        Point<uint> c(num_dims_);
//
//        for(size_t j=lo; j <= hi; j++) {
//            printf("p[%lu]:", items);
//            for(int i=0; i < num_dims_; i++) {
//                c[i] = p[i] + leaves_[i]->getField(j);
//                printf(" %u", c[i]);
//            }
//            printf("\n");
//
//            items++;
//        }
//        return;
//
//    }
//
//    else if (level == -1 || T_[level]->access(z) == 1) {
//        if( level == -1 ) y = 0;
//        else y = T_[level]->rank1(z-1) * children_[level+1];
//        // else y = T_[level]->rank1(z-1)*mypow(k_[level], p.num_dims());
//
//
//
//        /*// simple code, for the simple mx quadtree
//         * for (uint u =0 ; u < _K; u++) { c[0] = p[0] + u*nk;
//            for (uint v=0; v < _K; v++) { c[1] = p[1] + v*nk;
//                for( uint a=0; a < _K; a++) { c[2] = p[2] + a*nk;
//                    range(nk, from ,to, c, y + u*(1u << 2) + v*(1u << 1) + a, level+1, items);
//                }
//            }
//        }
//         */
//
//        /* or this simple code */
//        /*
//        Point<uint> c(p);
//        for(int i=0; i < children_[level+1]; i++) {
//            for(int j=0; j < num_dims_; j++) {
//                //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);
//
//                c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
//            }
//
//            range(nk, from ,to, c, y + i, level+1, items);
//        }
//         */
//
//        uint nk=nk_[level+1];
//        printf("nk[%d]: %u\n", level+1, nk_[level+1]);
//
//        if (false == is_interleaved_) {
//            //uint nk=n/k_[level+1];
//
//            Point<uint> c(p);
//            for(int i=0; i < children_[level+1]; i++) {
//                for(int j=0; j < num_dims_; j++) {
//                    //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);
//
//                    //c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
//                  c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
//                }
//
//                range(nk, from ,to, c, y + i, level+1, items);
//            }
//        }
//        else {
//            for(int i=0; i < children_[level+1]; i++) {
//                Point<uint> c(p);
//                for(int j=rangedim_by_level_[level+1].first; j < rangedim_by_level_[level+1].second; j++) {
//                    //printf("updating dim: %d\n", rangedim_by_level_[level+1].first + j);
//                    //c[rangedim_by_level_[level+1].first+j] += nk * ((i/mypow(k_[level+1], dims_[level+1]-j-1))%k_[level+1]);
//                    c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
//                }
//
//                printf("c: %u %u %u\n",c[0],c[1],c[2]);
//
//                range(nk, from ,to, c, y + i, level+1, items);
//            }
//
//
//
//        }
//
//    }
//}




void MXCompactQtreeFixed::all(Point<uint> p, size_t z, int level, vector<Point<uint> > &vpall) {
    size_t y = 0;

    //printf("checking T[%d] %lu\n",level,z);


    if (level == depth_-1) {
        if (T_[level]->access(z) == 0) return;

        size_t lo,hi;

        y = T_[level]->rank1(z-1);

        if (y == 0) {
            lo = 0;
            hi = B_->select1(y+1);
        }
        else {
            lo = B_->select1(y)+1;
            hi = B_->select1(y+1);
        }

        Point<uint> c(num_dims_);

        for(size_t j=lo; j <= hi; j++) {

            for(int i=0; i < num_dims_; i++) {
                c[i] = p[i] + leaves_[i]->getField(j);

            }
            vpall.push_back(c);
            if (vpall.size()%100000 == 0) {
                fprintf(stderr, "Progress: %.2f%% \r", (float)vpall.size()/items_*100);
            }

        }
        return;

    }

    else if (level == -1 || T_[level]->access(z) == 1) {
        if( level == -1 ) y = 0;
        else y = T_[level]->rank1(z-1) * children_[level+1];

        uint nk=nk_[level+1];


        if (false == is_interleaved_) {
            //uint nk=n/k_[level+1];

            Point<uint> c(p);
            for(int i=0; i < children_[level+1]; i++) {
                for(int j=0; j < num_dims_; j++) {
                    //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);

                    //c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
                  c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }

                all(c, y + i, level+1, vpall);
            }
        }
        else {
            for(int i=0; i < children_[level+1]; i++) {
                Point<uint> c(p);
                for(int j=rangedim_by_level_[level+1].first; j < rangedim_by_level_[level+1].second; j++) {
                    //printf("updating dim: %d\n", rangedim_by_level_[level+1].first + j);
                    //c[rangedim_by_level_[level+1].first+j] += nk * ((i/mypow(k_[level+1], dims_[level+1]-j-1))%k_[level+1]);
                    c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }
                all(c, y + i, level+1, vpall);
            }



        }

    }
}





void MXCompactQtreeFixed::range(Point<uint> &p, size_t z, int level, Point<uint> &from, Point<uint> &to, vector<Point<uint> > &vpall, size_t &items, bool pushval) {
  if (level != -1 ) {
    Point<uint> hi(p);
    if (false == is_interleaved_) {
      hi += nk_[level];
    }
    else {

      if (rangedim_by_level_[level].first == 0) {

        for (int j = rangedim_by_level_[level].first;
            j < rangedim_by_level_[level].second; j++) {
          hi[j] += nk_[level];
        }

        for (int j = rangedim_by_level_[level].second; j < num_dims_; j++) {
          hi[j] += (level==0)?maxvalue_:nk_[level - 1];
        }
      } else {

        //all are the same
        for (int i = 0; i < num_dims_; i++) {
          hi[i] += nk_[level];
        }

      }

    }
    if (!(p < to && from < hi)) {
      return;
    }
  }

    size_t y = 0;

    //printf("checking T[%d] %lu\n",level,z);


    if (level == depth_-1) {
        if (T_[level]->access(z) == 0) return;

        size_t lo,hi;

        y = T_[level]->rank1(z-1);

        if (y == 0) {
            lo = 0;
            hi = B_->select1(y+1);
        }
        else {
            lo = B_->select1(y)+1;
            hi = B_->select1(y+1);
        }

        Point<uint> c(num_dims_);

        for(size_t j=lo; j <= hi; j++) {

            for(int i=0; i < num_dims_; i++) {
                c[i] = p[i] + leaves_[i]->getField(j);

            }

            if (from <= c && c < to) {
                if(pushval) vpall.push_back(c);
                items++;
            }

        }
        return;

    }

    else if (level == -1 || T_[level]->access(z) == 1) {
        if( level == -1 ) y = 0;
        else y = T_[level]->rank1(z-1) * children_[level+1];

        uint nk=nk_[level+1];


        if (false == is_interleaved_) {
            //uint nk=n/k_[level+1];

            Point<uint> c(p);
            for(int i=0; i < children_[level+1]; i++) {
                for(int j=0; j < num_dims_; j++) {
                    //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);

                    //c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
                  c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }

                range(c, y + i, level+1,from,to, vpall,items,pushval);
            }
        }
        else {
            for(int i=0; i < children_[level+1]; i++) {
                Point<uint> c(p);
                for(int j=rangedim_by_level_[level+1].first; j < rangedim_by_level_[level+1].second; j++) {
                    //printf("updating dim: %d\n", rangedim_by_level_[level+1].first + j);
                    //c[rangedim_by_level_[level+1].first+j] += nk * ((i/mypow(k_[level+1], dims_[level+1]-j-1))%k_[level+1]);
                    c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }
                range(c, y + i, level+1, from,to, vpall,items,pushval);
            }



        }

    }
}

void MXCompactQtreeFixed::stats_space() const {  
  size_t treebits=0;
  printf(" l  k  d  k^d\t");
  printf("%10s\t%10s\t%10s\t%10s\t%s\n", "|T|","T.1s","T.1s/|T|", "T.bytes","T.cratio");
  for(int i = 0; i < depth_; i++) {
      printf("%2d %2d %2d %3d\t", i, k_[i],dims_[i],children_[i]);
      printf("%10lu\t%10lu\t%.2f\t%10lu\t%.2f\n", T_[i]->getLength(), T_[i]->countOnes(),1.0*T_[i]->countOnes()/T_[i]->getLength() ,T_[i]->getSize(), 1.0*T_[i]->getSize()*8/T_[i]->getLength());
      treebits += T_[i]->getSize();
  }
  printf("T space (MBytes): %0.1f\n", 1.0*treebits/1024/1024);
  size_t A=0;
  for(int j=0; j < num_dims_; j++) {
    A += leaves_[j]->getSize();
  }
  printf("A space (MBytes): %0.1f\n", 1.0*A/1024/1024);
  
  printf("Total space (MBytes): %0.1f\n", 1.0*treebits/1024/1024+1.0*A/1024/1024);
  
}

}




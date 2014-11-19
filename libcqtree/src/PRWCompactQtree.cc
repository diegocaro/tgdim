/*
 * PRWCompactQtreeHeight.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: diegocaro
 */

#include "PRWCompactQtree.h"
#include "utils.h"
#include "debug.h"

#include <algorithm> //sort
#include <vector>
#include <queue>
#include <climits>

namespace cqtree_static {

PRWCompactQtree::PRWCompactQtree(vector<Point<uint> > &vp,
                                 BitSequenceBuilder *bs, BitSequenceBuilder *bb, int k1, int k2, int max_level_k1, int max_levels_ki)
: k1_(k1), k2_(k2) {

    __setdefaultvalues();  //set default values for this class

    //preconditions
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

    //depth_ = ceil(log((double) maxval) / log((double) _K));

    uint max_k1=0;
    //max_k1 = (uint) pow((double) k1, (double) max_level_k1);
    max_k1 = mypow(k1, max_level_k1);

    levels_k1_ = max_level_k1;
    levels_k2_ = 0;
    LOG("resto1: %u\n",maxval/max_k1+1);



    if (false == is_interleaved_) {

        //levels using k2_ are automatically calculed
        if (max_k1 < maxval) {
            //depth_ += ceil(log((double) maxval/max_k1) / log((double) k2));
            levels_k2_ = mylog(k2, maxval/max_k1+1);
        }

        maxvalue_ = mypow(k1_,levels_k1_) * mypow(k2_,levels_k2_);
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
        // k1 must have available levels
        assert(max_k1 < maxval);

        // the number of levels availables must be at least the number of levels for ki_
        int lvl_availables_ki = mylog(ki_, maxval/max_k1+1);
        assert(lvl_availables_ki >= levels_ki_/2);

        int max_ki = mypow(ki_, levels_ki_/2);
        LOG("ki_: %d\n",ki_);
        LOG("l/2: %d\n", levels_ki_);
        LOG("max_ki: %d",max_ki);
        LOG("resto2: %u\n", (maxval/max_k1)/max_ki);

        levels_k2_ = mylog(k2, (maxval/max_k1)/max_ki+1);
        LOG("levels for k2: %d\n", levels_k2_);

        maxvalue_ = mypow(k1_,levels_k1_) * mypow(k2_,levels_k2_) * mypow(ki_, levels_ki_/2);

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


    get_stats(vp);

    size_t treebits=0;
    for(int i = 0; i < depth_; i++) {
        printf("T[%d]: %lu\n", i, nodes_by_level_[i]);
        treebits += nodes_by_level_[i];
    }
    printf("Total size T: %lu\n", treebits);

    size_t leafbits=0;
    for(int i = 0; i < depth_; i++) {
        printf("B[%d]: %lu\n", i, graywhite_[i]);
        leafbits += graywhite_[i];
    }
    printf("Total size B: %lu\n", leafbits);

    for(int i = 0; i < depth_; i++) {
        printf("leaves at level %d: %lu\n",i,leaves_per_level_[i]);
    }

    create(vp, bs, bb);
}

void PRWCompactQtree::get_stats(const std::vector<Point<uint> > &vp) {
    struct Node z;
    queue<struct Node> q;
    size_t r[max_children_ + 1];

    z.level = 0;
    z.lo = 0;
    z.hi = vp.size();
    q.push(z);

    graywhite_.insert(graywhite_.begin(), depth_, 0);

    nodes_by_level_.insert(nodes_by_level_.begin(), depth_, 0);

    leaves_per_level_.insert(leaves_per_level_.begin(), depth_, 0);

    vector<size_t> setbits;
    setbits.insert(setbits.begin(), depth_, 0);


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

                if (z.level < depth_ - 1 && r[j]-r[j-1] == 1) {
                    //black node found! (a leaf!)
                    leaves_per_level_[z.level] += 1;
                }
                else {
                    if (z.level < depth_ - 1) q.push(t);
                    setbits[z.level] += 1;
                }



            }
        }
    }

    for(int i=0; i < depth_-1; i++) {
        // number of nodes with one or zero leaves
        graywhite_[i] = nodes_by_level_[i]-setbits[i];
    }
}

void PRWCompactQtree::create(const std::vector<Point<uint> > &vp,
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

    vector<uint*> leafbits;
    for(int i = 0; i < depth_; i++) {
        leafbits.push_back(new uint[graywhite_[i] / W + 1]);
        fill_n(leafbits[i], graywhite_[i] / W + 1, 0); //fill with zeroes
    }

    vector<size_t> setbits;
    setbits.insert(setbits.begin(), depth_, 0);

    vector<size_t> curr_nodes;
    curr_nodes.insert(curr_nodes.begin(), depth_, 0);

    vector<size_t> curr_pos;
    curr_pos.insert(curr_pos.begin(), depth_, 0);



    leaves_ = new Array**[depth_];
    uint maxval, maxvalb;
    for(int i = 0; i < depth_-1; i++) {
        maxval = nk_[i];


        leaves_[i] = new Array* [num_dims_];

        if (false == is_interleaved_) {

            for(int j = 0; j < num_dims_; j++) {

                leaves_[i][j] = new Array(leaves_per_level_[i], maxval-1);
            }
        }
        else {

            if (i == 0) {
                maxvalb = maxvalue_;
            }
            else {
                maxvalb = nk_[i-1];
            }

            //k1 and k2 levels
            //all dimensions have the same remaining info
            if (i < levels_k1_ + levels_k2_) {
                for(int j=0; j < num_dims_; j++) {
                    leaves_[i][j] = new Array(leaves_per_level_[i], maxval-1);
                }
            }
            else {
                if (rangedim_by_level_[i].first == 0) {
                    //less info on first levels
                    if (i != depth_-2) {
                        for(int j=rangedim_by_level_[i].first; j < rangedim_by_level_[i].second; j++) {
                            leaves_[i][j] = new Array(leaves_per_level_[i], maxval-1);
                        }
                    } else {
                        for(int j=rangedim_by_level_[i].first; j < rangedim_by_level_[i].second; j++) {
                            leaves_[i][j] = new Array(0UL,0);
                        }
                    }

                    for(int j=rangedim_by_level_[i].second; j < num_dims_; j++) {
                        leaves_[i][j] = new Array(leaves_per_level_[i], maxvalb-1);
                    }
                }
                else {
                    //all dimensions have the same remaining info
                    for(int j = 0; j < num_dims_; j++) {
                        leaves_[i][j] = new Array(leaves_per_level_[i], maxval-1);
                    }
                }
            }
        }
    }

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

                if (z.level < depth_ - 1 && r[j]-r[j-1] == 1) {
                    //black node found!

                    // printf("level %u (bits) cut here %u %u\n", t.l, t.lo, t.hi);
                    //printf("B set bit %u\n", setbits);


                    //                    printf("set bit B_[%d] %lu\n", z.level, curr_nodes[z.level] - setbits[z.level] + j - 1);
                    cds_utils::bitset(leafbits[z.level], curr_nodes[z.level] - setbits[z.level] + j - 1);

                    Point<uint> p(vp[r[j-1]]);
                    uint maxval = nk_[z.level];
                    //                    printf("maxval: %u at level %d\n",maxval,z.level);

                    if (false == is_interleaved_) {
                        for(int i=0; i < num_dims_; i++) {
                            //                        printf("saving value %u at dim %i at pos %lu -> %u\n",p[i] % maxval,i,curr_pos[z.level],p[i]);
                            leaves_[z.level][i]->setField(curr_pos[z.level], p[i] % maxval);
                        }
                    }
                    else {
                        uint maxvalb;
                        if (z.level == 0) {
                            maxvalb = maxvalue_;
                        }
                        else {
                            maxvalb = nk_[z.level-1];
                        }

                        //k1 and k2 levels
                        if (z.level < levels_k1_ + levels_k2_) {
                            for(int i=0; i < num_dims_; i++) {
                                //                                               printf("level %d storing %u (dima %d) %d\n",z.level, p[i] % maxval,i,p[i]);
                                leaves_[z.level][i]->setField(curr_pos[z.level], p[i] % maxval);
                            }
                        }
                        else {

                            if (0 == rangedim_by_level_[z.level].first) {
                                //                                                 printf("range: [%d,%d)\n",rangedim_by_level_[z.level].first,rangedim_by_level_[z.level].second);

                                //less info on first levels
                                if (z.level != depth_-2) {
                                    for(int i=0; i < rangedim_by_level_[z.level].second; i++) {
                                        //                                                     printf("level %d storing %u (dimb %d) %d\n",z.level, p[i] % maxval,i,p[i]);
                                        leaves_[z.level][i]->setField(curr_pos[z.level], p[i] % maxval  );
                                    }
                                }

                                for(int i=rangedim_by_level_[z.level].second; i < num_dims_; i++) {
                                    //                                                     printf("level %d storing %u (dimc %d) %d\n",z.level, p[i] % maxvalb,i,p[i]);
                                    leaves_[z.level][i]->setField(curr_pos[z.level], p[i] % maxvalb);
                                }
                            }
                            else {
                                for(int i=0; i < num_dims_; i++) {
                                    //                                                     printf("level %d storing %u (dimd %d) %d\n",z.level, p[i] % maxval,i,p[i]);
                                    leaves_[z.level][i]->setField(curr_pos[z.level], p[i] % maxval);
                                }

                            }
                        }



                    }

                    curr_pos[z.level] += 1;
                }
                else {

                    if (z.level < depth_ - 1) q.push(t);

                    //printf("set bit T %lu\n", (ktreenodes-1)*childs_ + j-1);
                    //                                      printf("set bit T_[%d] %lu\n", z.level, curr_nodes[z.level] + j - 1);
                    cds_utils::bitset(treebits[z.level], curr_nodes[z.level] + j - 1);
                    setbits[z.level] += 1;
                }


            }
        }

        curr_nodes[z.level] += children_[z.level];

    }

    //building bitmaps for T and cleaning
    for(int i = 0; i < depth_; i++) {
        T_.push_back( bs->build(treebits[i], nodes_by_level_[i]));
        delete [] treebits[i];
    }

    //building bitmaps for B and cleaning
    for(int i = 0; i < depth_; i++) {
        B_.push_back( bb->build(leafbits[i], graywhite_[i]));
        delete [] leafbits[i];
    }
}
//
//void PRWCompactQtree::range(uint n, Point<uint> from, Point<uint> to, Point<uint> p, size_t z, int level, size_t &items) {
//    if (level != -1 && (p+n <= from || p > to) ) {
//        return;
//    }
//
//    size_t y = 0;
//
//    printf("checking T[%d] %lu\n",level,z);
//
//
//    if (level == depth_-1) {
//
//        if (T_[level]->access(z) == 0) return;
//        //printf("T_[%d] bit %d is set\n",level, z);
//        printf("p[%lu]:", items);
//        for(int i=0; i < p.num_dims(); i++) {
//            printf(" %u", p[i]);
//        }
//        printf("\n");
//        items++;
//
//        return;
//    }
//    else if (level > -1 &&  T_[level]->access(z) == 0) {
//      y = T_[level]->rank0(z)-1;
//
//      if (B_[level]->access(y) == 1) {
//                // node is a leaf
//        size_t b = B_[level]->rank1(y)-1;
//        Point<uint> c(p);
//
//                    if (false == is_interleaved_) {
//                        for(int i=0; i < num_dims_; i++) {
//                                      c[i] = p[i] + leaves_[level][i]->getField(b);
//                                   }
//                    }
//                    else {
//                        //k1 and k2 levels
//                      if (level < levels_k1_ + levels_k2_) {
//                          for(int i=0; i < num_dims_; i++) {
//                              c[i] = p[i] + leaves_[level][i]->getField(b);
//                        }
//                      }
//                      else {
//
//                          if (0 == rangedim_by_level_[level].first) {
//
//
//                              //less info on first levels
//                              if (level != depth_-2) {
//                                  for(int i=0; i < rangedim_by_level_[level].second; i++) {
//                                      c[i] = p[i] + leaves_[level][i]->getField(b);
//                                  }
//                              }
//
//                              for(int i=rangedim_by_level_[level].second; i < num_dims_; i++) {
//                                  c[i] = p[i] + leaves_[level][i]->getField(b);
//                              }
//                          }
//                          else {
//                              for(int i=0; i < num_dims_; i++) {
//                                  c[i] = p[i] + leaves_[level][i]->getField(b);
//                              }
//
//                          }
//                      }
//                    }
//
//                  printf("p[%lu]:", items);
//                  for(int i=0; i < num_dims_; i++) {
//                      printf(" %u", c[i]);
//                   }
//                   printf("\n");
//                   items++;
//
//                return;
//              }
//
//    }
//    else if (level == -1 || T_[level]->access(z) == 1) {
//      if( level == -1 ) y = 0;
//      else y = T_[level]->rank1(z-1) * children_[level+1];
//
//        /*// simple code, for the simple mx quadtree
//         * for (uint u =0 ; u < _K; u++) { c[0] = p[0] + u*nk;
//            for (uint v=0; v < _K; v++) { c[1] = p[1] + v*nk;
//                for( uint a=0; a < _K; a++) { c[2] = p[2] + a*nk;
//                    range(nk, from ,to, c, y + u*(1u << 2) + v*(1u << 1) + a, level+1, items);
//                }
//            }
//        }
//        */
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
//        */
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
//
//                  c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
//                }
//                printf("c: %u %u %u\n",c[0],c[1],c[2]);
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
//



void PRWCompactQtree::all(Point<uint> p, size_t z, int level, vector<Point<uint> > &vpall) {
    size_t y = 0;

    //printf("checking T[%d] %lu\n",level,z);


    if (level == depth_-1) {

        if (T_[level]->access(z) == 0) return;
        vpall.push_back(p);

        return;
    }
    else if (level > -1 &&  T_[level]->access(z) == 0) {
        y = T_[level]->rank0(z)-1;

        if (B_[level]->access(y) == 1) {
            //          printf("z: %lu\n",z);
            //          printf("y[%d]: %lu\n",level,y);
            //          printf("Checking T[%lu]\n",z);
            // node is a leaf
            size_t b = B_[level]->rank1(y)-1;
            //        printf("%u %u %u %u\n",p[0],p[1],p[2],p[3]);
            Point<uint> c(p);
            //        printf("b: %lu\n",b);
            if (false == is_interleaved_) {
                for(int i=0; i < num_dims_; i++) {
                    c[i] = p[i] + leaves_[level][i]->getField(b);
                }
                //                                    printf("c: %u %u %u %u\n",c[0],c[1],c[2],c[3]);
            }
            else {
                //k1 and k2 levels
                if (level < levels_k1_ + levels_k2_) {
                    for(int i=0; i < num_dims_; i++) {
                        c[i] = p[i] + leaves_[level][i]->getField(b);
                    }
                }
                else {

                    if (0 == rangedim_by_level_[level].first) {


                        //less info on first levels
                        if (level != depth_-2) {
                            for(int i=0; i < rangedim_by_level_[level].second; i++) {
                                c[i] = p[i] + leaves_[level][i]->getField(b);
                            }
                        }

                        for(int i=rangedim_by_level_[level].second; i < num_dims_; i++) {
                            c[i] = p[i] + leaves_[level][i]->getField(b);
                        }
                    }
                    else {
                        for(int i=0; i < num_dims_; i++) {
                            c[i] = p[i] + leaves_[level][i]->getField(b);
                        }

                    }
                }
            }

            vpall.push_back(c);

            return;
        }

    }
    else if (level == -1 || T_[level]->access(z) == 1) {
        if( level == -1 ) y = 0;
        else y = T_[level]->rank1(z-1) * children_[level+1];


        uint nk=nk_[level+1];
        // printf("nk[%d]: %u\n", level+1, nk_[level+1]);

        if (false == is_interleaved_) {
            //uint nk=n/k_[level+1];

            Point<uint> c(p);
            for(int i=0; i < children_[level+1]; i++) {
                for(int j=0; j < num_dims_; j++) {
                    //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);

                    //c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
                    c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }
                // printf("c: %u %u %u\n",c[0],c[1],c[2]);
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

                //printf("c: %u %u %u\n",c[0],c[1],c[2]);

                all(c, y + i, level+1, vpall);
            }



        }

    }
}




void PRWCompactQtree::range(Point<uint> &p, size_t z, int level, Point<uint> &from, Point<uint> &to, vector<Point<uint> > &vpall, size_t &items, bool pushval) {
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


  //    if (level != -1 && !comp(p, nk_[level])) {
//        return;
//    }
    size_t y = 0;

    //printf("checking T[%d] %lu\n",level,z);


    if (level == depth_-1) {

        if (T_[level]->access(z) == 0) return;
            if (pushval) vpall.push_back(p);
            items++;

        return;
    }
    else if (level > -1 &&  T_[level]->access(z) == 0) {
        y = T_[level]->rank0(z)-1;

        if (B_[level]->access(y) == 1) {
            //          printf("z: %lu\n",z);
            //          printf("y[%d]: %lu\n",level,y);
            //          printf("Checking T[%lu]\n",z);
            // node is a leaf
            size_t b = B_[level]->rank1(y)-1;
            //        printf("%u %u %u %u\n",p[0],p[1],p[2],p[3]);
            Point<uint> c(p);
            //        printf("b: %lu\n",b);
            if (false == is_interleaved_) {
                for(int i=0; i < num_dims_; i++) {
                    c[i] = p[i] + leaves_[level][i]->getField(b);
                }
                //                                    printf("c: %u %u %u %u\n",c[0],c[1],c[2],c[3]);
            }
            else {
                //k1 and k2 levels
                if (level < levels_k1_ + levels_k2_) {
                    for(int i=0; i < num_dims_; i++) {
                        c[i] = p[i] + leaves_[level][i]->getField(b);
                    }
                }
                else {

                    if (0 == rangedim_by_level_[level].first) {


                        //less info on first levels
                        if (level != depth_-2) {
                            for(int i=0; i < rangedim_by_level_[level].second; i++) {
                                c[i] = p[i] + leaves_[level][i]->getField(b);
                            }
                        }

                        for(int i=rangedim_by_level_[level].second; i < num_dims_; i++) {
                            c[i] = p[i] + leaves_[level][i]->getField(b);
                        }
                    }
                    else {
                        for(int i=0; i < num_dims_; i++) {
                            c[i] = p[i] + leaves_[level][i]->getField(b);
                        }

                    }
                }
            }

            if (from <= c && c < to) {
                if(pushval) vpall.push_back(c);
                items++;
            }
            return;
        }

    }
    else if (level == -1 || T_[level]->access(z) == 1) {
        if( level == -1 ) y = 0;
        else y = T_[level]->rank1(z-1) * children_[level+1];


        uint nk=nk_[level+1];
        // printf("nk[%d]: %u\n", level+1, nk_[level+1]);

        if (false == is_interleaved_) {
            //uint nk=n/k_[level+1];

            Point<uint> c(p);
            for(int i=0; i < children_[level+1]; i++) {
                for(int j=0; j < num_dims_; j++) {
                    //c[j] = p[j] + nk * ((i/(1 << (c.num_dims()-j-1)) )%k_[level+1]);

                    //c[j] = p[j] + nk * ((i/mypow(k_[level+1], num_dims_-j-1))%k_[level+1]);
                    c[j] = p[j] + nk * ((i/ kpower_per_level_dim_[num_dims_*(level+1)+j]  )%k_[level+1]);
                }
                // printf("c: %u %u %u\n",c[0],c[1],c[2]);
                range(c, y + i, level+1, from, to, vpall,items, pushval);
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

                //printf("c: %u %u %u\n",c[0],c[1],c[2]);

                range(c, y + i, level+1, from, to, vpall,items, pushval);
            }



        }

    }
}
/*
 template <typename T>
 void PRWCompactQtree<T>::save(ofstream & f) const
 {
 uint wr = MXQTRE_SAV;
 saveValue(f,wr);

 saveValue< PRWCompactQtree<T> >(f, this, 1);
 T_->save(f);
 }

 template <typename T>
 PRWCompactQtree<T>* PRWCompactQtree<T>::load(ifstream & f) {
 assert(f.good());
 uint type = loadValue<uint>(f);
 if(type!=MXQTRE_SAV) {	 // throw exception
 abort();
 }
 PRWCompactQtree<T> * ret = loadValue< PRWCompactQtree<T> >(f, 1);
 ret->T_ = BitSequence::load(f);

 return ret;
 }
 */

} /* namespace cqtree_static */




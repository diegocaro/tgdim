/*
 * bitvector.h
 *
 *  Created on: Aug 28, 2014
 *      Author: diegocaro
 */

#ifndef BITVECTOR_H_
#define BITVECTOR_H_

#include <vector>

namespace cqtree_utils {
const unsigned int W_ = 32;
class bitvector {
 public:
    bitvector() {
        length_ = 0;
    }

    // build vector of a given size
    bitvector(size_t size) {
        reserve(size);
        length_ = size;
    }

    ~bitvector() {

    }

    // reserve spaces and increase capacity
    void reserve(size_t s) {
//        printf("BitVector: increasing size\n");
        data_.reserve(s/W_+1);
        data_.insert(data_.end(), s/W_+1 - data_.size(), 0);
    }

    /** sets bit p in data */
    void bitset(size_t p) {
        if (p >= length_) {
            //if current data_ size need to be updated

            if (p/W_+1 < data_.size()) {
                printf("reserving more data...\n");
                reserve(2*p);
            }
            length_ = p+1;
        }

        data_[p/W_] |= (1<<(p%W_));
    }

    const unsigned int *data() const {
        return data_.data();
    }

    size_t length() const {
        return length_;
    }


 private:
    size_t length_;
    std::vector<unsigned int> data_;

};

} /* namespace cqtree_utils */

#endif /* BITVECTOR_H_ */

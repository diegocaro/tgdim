/*
 * MyCompressionCoder.h
 *
 *  Created on: Jan 15, 2015
 *      Author: diegocaro
 */

#ifndef MYCOMPRESSIONCODER_H_
#define MYCOMPRESSIONCODER_H_

#include <coding_policy.h>
#include "MyCoder.h"

namespace cqtree_static {

#define POLICYSIZE 255

class MyCompressionCoder : public MyCoder {
 public:
    MyCompressionCoder(char *policy) {
        policy_ = new char[POLICYSIZE];

        sprintf(policy_,"%s",policy);
        hc = new CodingPolicy(CodingPolicy::kPosition);
        hc->LoadPolicy(policy);
    }

    virtual size_t encode(uint symb, uint * stream, size_t pos) const {
    	UNUSED(symb); UNUSED(stream); UNUSED(pos); return 0;
    };
    virtual size_t decode(uint * symb, uint *stream, size_t pos) const {
    	UNUSED(symb); UNUSED(stream); UNUSED(pos); return 0;
    };
    virtual size_t maxLength() const { return 0; };
    virtual size_t getSize() const {return sizeof(MyCompressionCoder); }

    virtual void save(ofstream & fp) const {
        uint wr = MYCOMPR_HDR;
        cds_utils::saveValue(fp,wr);
        cds_utils::saveValue(fp,policy_,POLICYSIZE);
    };

    MyCompressionCoder (ifstream & fp) {
        uint type = loadValue<uint>(fp);
        if(type != MYCOMPR_HDR) {   //throw exception
            abort();
        }
        policy_ = cds_utils::loadValue<char>(fp,POLICYSIZE);
        hc = new CodingPolicy(CodingPolicy::kPosition);
        hc->LoadPolicy(policy_);
    }

    virtual size_t Compress(unsigned *input, unsigned *output, int n) {
        return hc->Compress(input,output,n);
    };
    virtual size_t Decompress(unsigned *input, unsigned *output, int n) {
        return hc->Decompress(input,output,n);
    }

 private:
    CodingPolicy *hc;
    char *policy_;
};

} /* namespace cqtree_static */

#endif /* MYCOMPRESSIONCODER_H_ */

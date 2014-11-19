#include "utils.h"
#include "MXCompactQtreeFixed.h"
#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace cqtree_utils;
using namespace cqtree_static;
using namespace cds_static;

namespace {

	class MXFTest : public ::testing::Test {
		protected:
			virtual void SetUp() {
				char inputfile[]="test/test6.input";
				
				FILE *f = fopen(inputfile, "r");
				if (f == NULL) {
					fprintf(stderr, "Test file '%s' not found.\n", inputfile);
					exit(1);
				}
				
				size_t nodes, edges, maxtime, contacts;
				    fscanf(f, "%lu %lu %lu %lu", &nodes, &edges, &maxtime, &contacts);

				 
				    Point<uint> c(4);
				    while(EOF != fscanf(f,"%u %u %u %u", &c[0], &c[1], &c[2], &c[3] )) {
				        vp_.push_back(c);
				    }
				
		  }
		
		  void AllPointsTester(int k1=2, int k2=2, int max_level_k1=0, int max_level_ki=0, int max_levels_fixed=1) {
		      printf("checking k1=%d k2=%d lk1=%d lki=%d fd=%d\n",k1,k2,max_level_k1,max_level_ki,max_levels_fixed);
		  	MXCompactQtreeFixed a(vp_, new BitSequenceBuilderRG(20),new BitSequenceBuilderRG(20),k1,k2,max_level_k1,max_level_ki,max_levels_fixed);
			vector<Point<uint> > vpall;
			    a.all(vpall);

			    ASSERT_EQ(vp_.size(), vpall.size());
			    for(size_t i=0; i < vp_.size(); i++) {
			        EXPECT_EQ(vp_[i], vpall[i]);
			    }
		  }

		  void LoadSave(int k1=2, int k2=2, int max_level_k1=0, int max_level_ki=0, int max_levels_fixed=1) {
		    MXCompactQtreeFixed a(vp_, new BitSequenceBuilderRG(20),new BitSequenceBuilderRG(20),k1,k2,max_level_k1,max_level_ki,max_levels_fixed);

		      ofstream f;
		      f.open("/tmp/aaa",ios::binary);
		      a.save(f);
		      f.close();

		      ifstream fp;
		      fp.open("/tmp/aaa",ios::binary);
		      MXCompactQtreeFixed b(fp);
		      fp.close();


		      vector<Point<uint> > vpall;
		      b.all(vpall);

		      ASSERT_EQ(vp_.size(), vpall.size());
		      for(size_t i=0; i < vp_.size(); i++) {
		          EXPECT_EQ(vp_[i], vpall[i]);
		      }
		  }

		  
		vector<Point<uint> > vp_;
		
	};
	
	TEST_F(MXFTest, AllPointsSimple) {
		AllPointsTester(2,2,0,0,1);
	}
	
	TEST_F(MXFTest, K1Equals4) {
	    int k1=4;
	    int k2=2;
	    Point<uint> pmax_dim(max(vp_));
	    uint maxval;
	    maxval = max(pmax_dim);

	    uint l4 = mylog(k1, maxval) / 2;

	    for(uint i=0; i < l4; i++) {
	        AllPointsTester(k1,k2,i,0,3);
	    }
	}

    TEST_F(MXFTest, Ki) {
        int k1=2;
        int k2=2;
        Point<uint> pmax_dim(max(vp_));
        uint maxval;
        maxval = max(pmax_dim);

        uint li = mylog(2, maxval)-8;

        for(uint i=0; i < li; i++) {
            AllPointsTester(k1,k2,1,i*2,2);
        }
    }

    TEST_F(MXFTest, AllMods) {
        int k1=4;
        int k2=2;
        Point<uint> pmax_dim(max(vp_));
        uint maxval;
        maxval = max(pmax_dim);

        uint l4 = mylog(k1, maxval);
        uint li = mylog(2, maxval);

        for(uint i=0; i < l4; i++) {
            for(uint j=0; j < 2*(l4-i); j+=2) {
                for(uint k=1; k < li-j-2*i; k++) {
                    AllPointsTester(k1,k2,i,j,k);
                }
            }
        }
    }

    TEST_F(MXFTest, LoadSaveBasic) {
        LoadSave();
    }

    TEST_F(MXFTest, LoadSaveAllMods) {
      int k1=4;
      int k2=2;
      Point<uint> pmax_dim(max(vp_));
      uint maxval;
      maxval = max(pmax_dim);

      uint l4 = mylog(k1, maxval);
      uint li = mylog(2, maxval);

      for(uint i=0; i < l4; i++) {
          for(uint j=0; j < 2*(l4-i); j+=2) {
              for(uint k=1; k < li-j-2*i; k++) {
                  LoadSave(k1,k2,i,j,k);
              }
          }
      }
    }

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

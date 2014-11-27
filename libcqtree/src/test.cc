/*
 * test.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: diegocaro
 */


#include "utils.h"
//#include "MXCompactQtree.h"
//#include "MXCompactQtreeFixed.h"
//#include "PRBCompactQtree.h"
//#include "PRWCompactQtree.h"
#include "PRB2CompactQtree.h"

using namespace cqtree_static;
using namespace cqtree_utils;

using namespace cds_static;

using namespace std;

class a:public _Compare {
virtual bool operator()(const Point<uint> &lo, const Point<uint> &hi) const {
           if (lo[0] <= 16 && hi[0] >= 16) {
               return true;
           }
           return false;
   }
} Less;

int main() {

    size_t nodes, edges, maxtime, contacts;
    scanf("%lu %lu %lu %lu", &nodes, &edges, &maxtime, &contacts);

    vector<Point<uint> > vp;
    Point<uint> c(4);
    while(EOF != scanf("%u %u %u %u", &c[0], &c[1], &c[2], &c[3] )) {
        vp.push_back(c);
    }

    assert(contacts == vp.size());
    BitSequenceBuilderRG rg(20);
    //PRWCompactQtree a(vp, &rg, &rg, 4,2,1,0);
    //MXCompactQtree a(vp, &rg, 4,2,0,0);
    //MXCompactQtreeFixed a(vp, &rg, &rg, 4,2,0,0);
    //PRBCompactQtree a(vp, &rg, &rg, 4,2,1,0);
    PRB2CompactQtree a(vp, &rg, &rg, 2,2,0,0);
    size_t items=0;

    vector<Point<uint> > vpall;
    a.all(vpall);

    printf("vp size: %lu\n",vp.size());
    printf("ans size: %lu\n",vpall.size());
    assert(vp.size() == vpall.size());
    for(size_t i=0; i < vp.size(); i++) {
        
        printf("a: %u %u %u %u\n",vpall[i][0],vpall[i][1],vpall[i][2],vpall[i][3]);
        printf("b: %u %u %u %u\n",vp[i][0],vp[i][1],vp[i][2],vp[i][3]);
        assert(vp[i] == vpall[i]);
    }


    for(size_t i=0; i < vp.size(); i++) {
            //assert(vp[i] == vpall[i]);
            printf("b: %u %u %u %u\n",vp[i][0],vp[i][1],vp[i][2],vp[i][3]);
        }




        vpall.clear();
printf("vpallsize: %lu\n",vpall.size());
a.range(vpall, Less);
for(size_t i=0; i < vpall.size(); i++) {
    printf("answer: ");
    vpall[i].print();

}

/*
    printf("saving a\n");
    ofstream f;
f.open("/tmp/aaa",ios::binary);
    a.save(f);
    f.close();

printf("loading b\n");
    ifstream fp;
    fp.open("/tmp/aaa",ios::binary);
    MXCompactQtree *b= new MXCompactQtree(fp);
    fp.close();

printf("a pointer: %p\n",&a);
    printf("b pointer: %p\n",b);


printf("removing b:\n");
delete b;

printf("removing a:\n");


*/

//    Point<uint> from(4);
//    from[0]=0;from[1]=0;from[2]=0;from[3]=0;
//    Point<uint> to(4);
//    to[0]=nodes;to[1]=nodes;to[2]=maxtime;to[3]=maxtime;
//    a.range(from, to, items);
//

//	vector<Point<uint> > vp;
//
//	Point<uint> p(3);
//	p[0] = 1;
//	p[1] = 2;
//	p[2] = 3;
//
//	vp.push_back(p);
//	p[0] = 2;
//	p[1] = 4;
//	p[2] = 6;
//
//	vp.push_back(p);
//
//	p[0] = 3;
//	p[1] = 5;
//	p[2] = 7;
//
//	vp.push_back(p);



//	MXCompactQtree a(vp, new BitSequenceBuilderRG(20),4,2,1,2);
//
//	Point<uint> from(3);
//	from[0]=0;from[1]=0;from[2]=0;
//	Point<uint> to(3);
//	to[0]=8;to[1]=8;to[2]=8;
//
//	size_t items=0;
//	a.range(from, to, items);

//    PRWCompactQtree a(vp, new BitSequenceBuilderRG(20),4,2,1,0);
//
//    Point<uint> from(3);
//    from[0]=0;from[1]=0;from[2]=0;
//    Point<uint> to(3);
//    to[0]=8;to[1]=8;to[2]=8;
//
//    size_t items=0;
//    a.range(from, to, items);




//    vector<Point<uint> > vp2;
//
//    Point<uint> q(2);
//    q[0] = 1;
//    q[1] = 2;
//
//    vp2.push_back(q);
//    q[0] = 2;
//    q[1] = 4;
//
//    vp2.push_back(q);
//
//    q[0] = 3;
//    q[1] = 5;
//
//    vp2.push_back(q);
//
//    MXCompactQtree b(vp2, new BitSequenceBuilderRG(20));
//    from = q;
//    to = q;
//    from[0]=0;from[1]=0;
//    to[0]=8;to[1]=8;
//    items=0;
//    b.range(from, to, items);


	/*
	PRCompactQtreeBlack<Point3d> a(vp, new BitSequenceBuilderRG(20));

	vector<Point4d> vp2;
	Point4d q;
	q[0] = 1;
	q[1] = 2;
	q[2] = 3;
	q[3] = 4;
	vp2.push_back(q);
	q[0] = 2;
	q[1] = 4;
	q[2] = 6;
	q[3] = 8;
	vp2.push_back(q);

	q[0] = 3;
	q[1] = 5;
	q[2] = 7;
	q[3] = 9;
	vp2.push_back(q);
	vp2.push_back(p);

	PRCompactQtreeBlack<Point4d> b(vp2, new BitSequenceBuilderRG(20));

	ofstream fileout;
	fileout.open("/tmp/test4d", ios::binary);
	b.save(fileout);
	fileout.close();

	ifstream filein;
	filein.open("/tmp/test4d", ios::binary);
	CompactQtree<Point4d> *c = PRCompactQtreeBlack<Point4d>::load(filein);
	filein.close();
	*/

	return 0;
}

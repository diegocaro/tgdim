#ifndef PARAMS_H_
#define PARAMS_H_ 

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

#endif
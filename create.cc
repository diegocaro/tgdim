//
// Maximo numero de vertexes: 2^32-1
// Maximo numero de timepoints 2^32-1
//

// TODO: mejorar localidad en arreglo de hojas, de tal manera que todas las componentes esten cerca..
// e_[4*i+0] e_[4*i+1] e_[4*i+2] e_[4*i+3] -> parece que eso no mejora :P

#include <cstdio>
#include <cmath>
#include <climits>
#include <cstring>

#include <vector>
#include <queue>

#include <algorithm>

// Libcds
#include <BitSequence.h>
#include <BitSequenceBuilder.h>

// mandatory
//#include "TemporalGraph.h"
#include <Point.h>
#include <CompactQtree.h>
//#include <MXCompactQtreeFixed.h>
//#include <PRBCompactQtree.h>
//#include <PRWCompactQtree.h>
// extra
#include "debug.h"

using namespace cds_static;
using namespace cqtree_static;

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
  char *outfile;

  enum TypeGraph typegraph;

  int k1;
  int k2;
  int lk1; //levels for k1
  int lki; //levels for ki
  int lf;  //levels for fixed mx
};

int readopts(int argc, char **argv, struct opts *opts) {
  int o;

  int dsflag = 0;
  // Default options
  opts->bs = eRG;
  opts->bb = eRG;
  opts->typegraph = kInterval;

  opts->k1 = 4;
  opts->k2 = 2;

  opts->lk1 = 0;
  opts->lki = 0;

  opts->lf = 0;

  while ((o = getopt(argc, argv, "t:b:s:g:z:k:x:c:f:")) != -1) {
    switch (o) {
      case 't':
        if (strcmp(optarg, "RG") == 0) {
          INFO("Using RG for T bitmaps");
          opts->bs = eRG;
        } else if (strcmp(optarg, "RRR") == 0) {
          INFO("Using RRR for T bitmaps");
          opts->bs = eRRR;
        } else if (strcmp(optarg, "SD") == 0) {
          INFO("Using SDarray for T bitmaps");
          opts->bs = eSD;
        }
        break;
      case 'b':
        if (strcmp(optarg, "RG") == 0) {
          INFO("Using RG for B bitmaps");
          opts->bb = eRG;
        } else if (strcmp(optarg, "RRR") == 0) {
          INFO("Using RRR for B bitmaps");
          opts->bb = eRRR;
        } else if (strcmp(optarg, "SD") == 0) {
          INFO("Using SDarray for B bitmaps");
          opts->bb = eSD;
        }
        break;
      case 's':
        dsflag = 1;
        if (strcmp(optarg, "PRW") == 0) {
          INFO("Using PRWhite");
          opts->ds = ePRWhite;
        } else if (strcmp(optarg, "PRB") == 0) {
          INFO("Using PRBlack");
          opts->ds = ePRBlack;
        } else if (strcmp(optarg, "MXD") == 0) {
          INFO("Using MXDepth");
          opts->ds = eMXDepth;
        } else if (strcmp(optarg, "MXF") == 0) {
          INFO("Using MXFixed");
          opts->ds = eMXFixed;
        } else {
          dsflag = 0;
        }
        break;
      case 'g':
        if (strcmp(optarg, "I") == 0) {
          INFO("Interval-contact Temporal Graph");
          opts->typegraph = kInterval;
        } else if (strcmp(optarg, "P") == 0) {
          INFO("Point-contact Temporal Graph");
          opts->typegraph = kPoint;
        } else if (strcmp(optarg, "G") == 0) {
          INFO("Growing Temporal Graph");
          opts->typegraph = kGrowth;
        }
        break;
      case 'z':
        opts->k1 = atoi(optarg);
        break;
      case 'k':
        opts->k2 = atoi(optarg);
        break;
      case 'x':
        opts->lk1 = atoi(optarg);
        break;
      case 'c':
        opts->lki = atoi(optarg);
        break;
      case 'f':
        opts->lf = atoi(optarg);
        break;
      default: /* '?' */
        break;
    }
  }

  if (optind >= argc || (argc - optind) < 1 || dsflag == 0
      || (opts->lf == 0 && opts->ds == eMXFixed)) {
    fprintf(stderr,
        "%s -s {MXD,MXF,PRB,PRW} [-g I,P,G] [-t RG,RRR,SD] [-b RG,RRR,SD] <outputfile> \n",
        argv[0]);
    fprintf(stderr, "Expected data structure (-s):\n");
    fprintf(stderr, "\tMXD for MatriX Quadtree (automatic depth)\n");
    fprintf(stderr, "\tMXF for MatriX Quadtree Fixed Depth\n");
    fprintf(stderr, "\tPRB for Point Region Quadtree Leaves as Black Nodes\n");
    fprintf(stderr, "\tPRW for Point Region Quadtree Leaves as White Nodes\n");

    fprintf(stderr, "\nExpected data structure flags:\n");
    fprintf(stderr,  "\t -z set k1\n");
    fprintf(stderr,  "\t -k set k2\n");
    fprintf(stderr,  "\t -x set lk1\n");
    fprintf(stderr,  "\t -c set lki\n");
    fprintf(stderr,  "\t -f set the number of fixed levels (MXF only)\n");


    fprintf(stderr, "\nExpected type of graph (-g):\n");
    fprintf(stderr, "\tI for Interval-contact Temporal Graph\n");
    fprintf(stderr, "\tP for Point-contact Temporal Graph\n");
    fprintf(stderr, "\tG for Growing Temporal Graph\n");



    fprintf(stderr, "\nExpected argument after options\n");
    fprintf(stderr, "\t<outputfile> destination file\n");
    exit(EXIT_FAILURE);
  }

  opts->outfile = argv[optind];

  return optind;

}

void printsettings(struct opts *opts) {
  switch (opts->typegraph) {
    case kGrowth:
      printf("G Growing Temporal Graph (3d)\n");
      break;
    case kInterval:
      printf("I Interval Temporal Graph (4d)\n");
      break;
    case kPoint:
      printf("P Point Temporal Graph (3d)\n");
      break;

  }

  switch (opts->ds) {
    case eMXDepth:
      printf("MXD compact data structure\n");
      break;
    case eMXFixed:
      printf("MXF compact data structure\n");
      break;
    case ePRBlack:
      printf("PRB compact data structure\n");
      break;
    case ePRWhite:
      printf("PRW compact data structure\n");
      break;
  }

  printf("k1: %d\n",opts->k1);
  printf("k2: %d\n",opts->k2);
  printf("lk1: %d\n",opts->lk1);
  printf("lki: %d\n",opts->lki);
  printf("lf: %d\n",opts->lf);

  switch (opts->bb) {
    case eRG:
      printf("RG for T bits\n");
      break;
    case eRRR:
      printf("RRR for T bits\n");
      break;
    case eSD:
      printf("SDarray for T bits\n");
      break;
  }

  switch (opts->bs) {
    case eRG:
      printf("RG for B bits\n");
      break;
    case eRRR:
      printf("RRR for B bits\n");
      break;
    case eSD:
      printf("SDarray for B bits\n");
      break;
  }

}

int main(int argc, char *argv[]) {
  uint nodes, edges, lifetime, contacts;

  struct opts opts;
  readopts(argc, argv, &opts);

  printsettings(&opts);

  // Reading input
  INFO("Reading input...");
  scanf("%u %u %u %u", &nodes, &edges, &lifetime, &contacts);
  LOG("nodes: %u", nodes);
  LOG("edges: %u", edges);
  LOG("maxtime: %u", maxtime);
  LOG("contacts: %u", contacts);




  size_t readcontacts = 0;
  vector<Point<uint> > vp;

  if (opts.typegraph == kInterval) {
    //4dim data
    Point<uint> c(4);
    while(EOF != scanf("%u %u %u %u", &c[0], &c[1], &c[2], &c[3] )) {
      readcontacts++;
      if (readcontacts%10000==0)fprintf(stderr, "Reading data: %.2f%% \r", (float)readcontacts/contacts*100);

      assert(c[0] < nodes);
      assert(c[1] < nodes);
      assert(c[2] < lifetime);
      assert(c[3] < lifetime);

      vp.push_back(c);
    }
  }
  else {
    //3dim data
    Point<uint> c(3);
    uint dummy;
    while(EOF != scanf("%u %u %u %u", &c[0], &c[1], &c[2], &dummy )) {
      readcontacts++;
      if (readcontacts%10000==0)fprintf(stderr, "Reading data: %.2f%% \r", (float)readcontacts/contacts*100);

      assert(c[0] < nodes);
      assert(c[1] < nodes);
      assert(c[2] < lifetime);

      vp.push_back(c);
    }
  }

  assert(readcontacts == contacts);

  BitSequenceBuilder *bs=NULL;
  switch(opts.bs) {
  case eRG:
	  bs = new BitSequenceBuilderRG(20); // by default, 5% of extra space for bitmaps
	  break;
  case eRRR:
	  bs = new BitSequenceBuilderRRR(32); // DEFAULT_SAMPLING for RRR is 32
	  break;
  case eSD:
      bs = new BitSequenceBuilderSDArray(); // ?
      break;
  }

  BitSequenceBuilder *bb=NULL;
  switch(opts.bb) {
  case eRG:
    bb = new BitSequenceBuilderRG(20); // by default, 5% of extra space for bitmaps
    break;
  case eRRR:
    bb = new BitSequenceBuilderRRR(32); // DEFAULT_SAMPLING for RRR is 32
    break;
  case eSD:
      bb = new BitSequenceBuilderSDArray(); // ?
      break;
  }


  CompactQtree *cq;

  ofstream file;
  LOG("Saving graph file in '%s'", opts.outfile);
  file.open(opts.outfile, ios::binary);
  switch(opts.ds) {
    case ePRBlack:
      cq = new PRBCompactQtree(vp,bs,bb,opts.k1,opts.k2,opts.lk1,opts.lki);
      break;
    case ePRWhite:
      cq = new PRWCompactQtree(vp,bs,bb,opts.k1,opts.k2,opts.lk1,opts.lki);
      break;
    case eMXDepth:
      cq = new MXCompactQtree(vp,bs,opts.k1,opts.k2,opts.lk1,opts.lki);
      break;
    case eMXFixed:
      cq = new MXCompactQtreeFixed(vp,bs,bb,opts.k1,opts.k2,opts.lk1,opts.lki,opts.lf);
      break;
  }

  //cq->setInfo(nodes,edges,lifetime,contacts);
  cq->save(file);

  file.close();

  delete bb;
  delete bs;

  return 0;

}

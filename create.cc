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
#include <point.h>
#include <CompactQtree.h>
//#include <MXCompactQtreeFixed.h>
//#include <PRBCompactQtree.h>
//#include <PRWCompactQtree.h>
// extra
#include "debug.h"

#include "TemporalGraph.h"

using namespace cds_static;
using namespace cqtree_static;

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

  printf("Reading input file '%s'\n",opts->infile);

}
int readflags(struct opts *opts, char *flags) {
  // "k1,k2,lk1,lkf,lf"
  vector<string> f;
  tokenize(flags,f,',');

  if (f.size() < 4) {
    return -1; //error
  }

  opts->k1 = atoi(f[0].c_str());
  opts->k2 = atoi(f[1].c_str());
  opts->lk1 = atoi(f[2].c_str());
  opts->lki = atoi(f[3].c_str());

  if (f.size() == 5) {
      opts->lf = atoi(f[4].c_str());
    }

  return f.size();
}

int readopts(int argc, char **argv, struct opts *opts) {
  int o;

  int fflags = 0;
  int dsflag = 0;
  // Default options
  opts->infile = "-";

  opts->bs = eRG;
  opts->bb = eRG;
  opts->typegraph = kInterval;

  opts->k1 = 4;
  opts->k2 = 2;

  opts->lk1 = 0;
  opts->lki = 0;

  opts->lf = 1;

  while ((o = getopt(argc, argv, "t:b:s:g:f:i:")) != -1) {
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
      case 'f':
        fflags = readflags(opts,optarg);
        break;
      case 'i':
        opts->infile = optarg;
        break;
      default: /* '?' */
        break;
    }
  }

  if (optind >= argc || (argc - optind) < 1 || dsflag == 0 || fflags == -1
      || (opts->lf == 0 && opts->ds == eMXFixed)) {
    fprintf(stderr,
        "%s -s {MXD,MXF,PRB,PRW} [-f k1:k2:lk1:lki:lf] [-g I,P,G] [-t RG,RRR,SD] [-b RG,RRR,SD] [-i <inputfile>] <outputfile> \n",
        argv[0]);
    fprintf(stderr, "Expected data structure (-s):\n");
    fprintf(stderr, "\tMXD for MatriX Quadtree (automatic depth)\n");
    fprintf(stderr, "\tMXF for MatriX Quadtree Fixed Depth\n");
    fprintf(stderr, "\tPRB for Point Region Quadtree Leaves as Black Nodes\n");
    fprintf(stderr, "\tPRW for Point Region Quadtree Leaves as White Nodes\n");

    fprintf(stderr, "\nExpected data structure flags (-f k1,k2,lk1,lki,lf):\n");
    fprintf(stderr,  "\t lk1 set the number of levels using k1\n");
    fprintf(stderr,  "\t lk1 set the number of levels using half dimensions\n");
    fprintf(stderr,  "\t lf set the number of fixed levels (MXF only)\n");

    fprintf(stderr, "\nExpected type of graph (-g):\n");
    fprintf(stderr, "\tI for Interval-contact Temporal Graph\n");
    fprintf(stderr, "\tP for Point-contact Temporal Graph\n");
    fprintf(stderr, "\tG for Growing Temporal Graph\n");

    fprintf(stderr, "\nExpected input file -i input file can be set to '-' to read stdin\n");

    fprintf(stderr, "\nExpected argument after options\n");
    fprintf(stderr, "\t<outputfile> destination file\n");
    exit(EXIT_FAILURE);
  }

  opts->outfile = argv[optind];

  return optind;

}



int main(int argc, char *argv[]) {
  uint nodes, edges, lifetime, contacts;

  struct opts opts;
  readopts(argc, argv, &opts);

  printsettings(&opts);



  FILE *infile;
  if ( strcmp(opts.infile,"-") == 0 ) {
    infile = stdin;
  }
  else {
    infile = fopen(opts.infile, "r");
  }

  // Reading input
  INFO("Reading input...");
  fscanf(infile,"%u %u %u %u", &nodes, &edges, &lifetime, &contacts);
  LOG("nodes: %u", nodes);
  LOG("edges: %u", edges);
  LOG("maxtime: %u", maxtime);
  LOG("contacts: %u", contacts);


  size_t readcontacts = 0;
  vector<Point<uint> > vp;

  if (opts.typegraph == kInterval) {
    //4dim data
    Point<uint> c(4);
    while(EOF != fscanf(infile,"%u %u %u %u", &c[0], &c[1], &c[2], &c[3] )) {
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
    while(EOF != fscanf(infile,"%u %u %u %u", &c[0], &c[1], &c[2], &dummy )) {
      readcontacts++;
      if (readcontacts%10000==0)fprintf(stderr, "Reading data: %.2f%% \r", (float)readcontacts/contacts*100);

      assert(c[0] < nodes);
      assert(c[1] < nodes);
      assert(c[2] < lifetime);

      vp.push_back(c);
    }
  }

  fclose(infile);
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


  //First check points...
  vector<Point<uint> > vpall;
  cq->all(vpall);

  if (vp.size() != vpall.size()) {
    fprintf(stderr, "Error: data from data structure doesnt match the input\n");
    abort();
  }
  for(size_t i=0; i < vp.size(); i++) {
    if (i%10000==0)fprintf(stderr, "Checking data: %.2f%% \r", (float)i/contacts*100);
    assert(vp[i] == vpall[i]);
  }


  TemporalGraph *tg;

  switch(opts.typegraph) {
      case kInterval:
          tg = new IntervalContactGraph();
      break;
  }



  ofstream file;
  LOG("Saving graph file in '%s'", opts.outfile);
  
  fprintf(stderr, "Saving data structure\n");
  file.open(opts.outfile, ios::binary);
  
  tg->setDs(cq);
  tg->setInfo(nodes,edges,lifetime,contacts);
  
  tg->save(file);

  file.close();

  delete bb;
  delete bs;

  delete tg;

  return 0;

}
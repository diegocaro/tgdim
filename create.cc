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
    case kIntervalPro:
          printf("R Interval Pro 4d+3d\n");
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
    case ePRB2Black:
       printf("PRB2 compact data structure\n");
       break;
    case ePRWhite:
      printf("PRW compact data structure\n");
      break;
  }

//  printf("params_char = '%s'\n",opts->params_char);
//  for(map<string,string>::iterator it=opts->params.begin(); it != opts->params.end(); ++it) {
//      cout << "params[" << it->first << "] = " << it->second << endl;
//  }
  printf("T: %s\n",opts->params["T"].c_str());
  printf("B: %s\n",opts->params["B"].c_str());
  printf("C: %s\n",opts->params["C"].c_str());

  printf("k1: %d\n",opts->k1);
  printf("k2: %d\n",opts->k2);

  printf("lk1: %d\n",opts->lk1);
  printf("lki: %d\n",opts->lki);
  printf("lf: %d\n",opts->lf);

  printf("F: %d\n", opts->F);


  printf("Reading input file '%s'\n",opts->infile);

}

int readopts(int argc, char **argv, struct opts *opts) {
    int o;

    int dsflag = 0;

    // Default options
    opts->infile = "-";

    const char *def_params = "T:RG5,B:RG5,C:RG5,k1:4,k2:2,lk1:0,lki:0,F:2,lf:1";
    readflags(opts, def_params);

    opts->typegraph = kInterval;

    while ((o = getopt(argc, argv, "s:g:f:i:")) != -1) {
        switch (o) {
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
                } else if (strcmp(optarg, "PRB2") == 0) {
                    INFO("Using PRB2Black");
                    opts->ds = ePRB2Black;
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
                } else if (strcmp(optarg, "R") == 0) {
                    INFO("Interval Pro");
                    opts->typegraph = kIntervalPro;
                }
                break;
            case 'f':
                readflags(opts, optarg);
                break;
            case 'i':
                opts->infile = optarg;
                break;
            default: /* '?' */
                break;
        }
    }

    opts->bs = getBSBuilder(opts->params["T"]);
    opts->bb = getBSBuilder(opts->params["B"]);
    opts->bc = getBSBuilder(opts->params["C"]);

    if (optind >= argc || (argc - optind) < 1 || dsflag == 0 || opts->bs == NULL || opts->bb == NULL || opts->bc == NULL
            || (opts->lf == 0 && opts->ds == eMXFixed)) {
        fprintf(stderr,
                "%s -s {MXD,MXF,PRB,PRW,PRB2} [-f k1,k2,lk1,lki,lf,F] [-g I,P,G] [-i <inputfile>] <outputfile> \n",
                argv[0]);
        fprintf(stderr, "Expected data structure (-s):\n");
        fprintf(stderr, "\tMXD for MatriX Quadtree (automatic depth)\n");
        fprintf(stderr, "\tMXF for MatriX Quadtree Fixed Depth\n");
        fprintf(stderr,
                "\tPRB for Point Region Quadtree Leaves as Black Nodes\n");
        fprintf(stderr,
                "\tPRB2 for Point Region Quadtree Leaves as Black Nodes (variable)\n");
        fprintf(stderr,
                "\tPRW for Point Region Quadtree Leaves as White Nodes\n");

        fprintf(stderr,
                "\nExpected data structure flags (-f k1:2,k2:2,lk1:0,lki:1,lf:2,F:2,T:RG,B:RG,C:RG):\n");
        fprintf(stderr, "\t lk1 set the number of levels using k1\n");
        fprintf(stderr,
                "\t lk1 set the number of levels using half dimensions\n");
        fprintf(stderr, "\t lf set the number of fixed levels (MXF only)\n");
        fprintf(stderr, "\t F set the maximum number of leaves (PRB2 only)\n");
        fprintf(stderr, "\t T,B and C can be of type  {RG,RRR,SD,SRRR15}\n");

        fprintf(stderr, "\nExpected type of graph (-g):\n");
        fprintf(stderr, "\tI for Interval-contact Temporal Graph\n");
        fprintf(stderr, "\tP for Point-contact Temporal Graph\n");
        fprintf(stderr, "\tG for Growing Temporal Graph\n");
        fprintf(stderr, "\tR for Interval Pro Graph\n");

        fprintf(stderr,
                "\nExpected input file -i input file can be set to '-' to read stdin\n");

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

  vector<Point<uint> > vpcurr;
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
  else if (opts.typegraph == kIntervalPro) {
      //4dim data
          Point<uint> c(4);
          Point<uint> c3(3);
          while(EOF != fscanf(infile,"%u %u %u %u", &c[0], &c[1], &c[2], &c[3] )) {
            readcontacts++;
            if (readcontacts%10000==0)fprintf(stderr, "Reading data: %.2f%% \r", (float)readcontacts/contacts*100);

            assert(c[0] < nodes);
            assert(c[1] < nodes);
            assert(c[2] < lifetime);
            assert(c[3] < lifetime);

            if (c[3] == lifetime-1) {
                c3[0]=c[0];
                c3[1]=c[1];
                c3[2]=c[2];

                vpcurr.push_back(c3);
            }
            else {
                vp.push_back(c);
            }
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

  CompactQtree *cq;

  switch(opts.ds) {
    case ePRBlack:
      cq = new PRBCompactQtree(vp,opts.bs,opts.bb,opts.k1,opts.k2,opts.lk1,opts.lki);
      break;
    case ePRB2Black:
          cq = new PRB2CompactQtree(vp,opts.bs,opts.bb,opts.bc,opts.k1,opts.k2,opts.F,opts.lk1,opts.lki);
          break;
//    case ePRWhite:
//      cq = new PRWCompactQtree(vp,bs,bb,opts.k1,opts.k2,opts.lk1,opts.lki);
//      break;
    case eMXDepth:
      cq = new MXCompactQtree(vp,opts.bs,opts.k1,opts.k2,opts.lk1,opts.lki);
      break;
    case eMXFixed:
      cq = new MXCompactQtreeFixed(vp,opts.bs,opts.bb,opts.k1,opts.k2,opts.lk1,opts.lki,opts.lf);
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
    if (i%1000==0)fprintf(stderr, "Checking data: %.2f%% \r", (float)i/contacts*100);
    //assert(vp[i] == vpall[i]);
    //printf("%d %d %d %d %d %d %d %d\n",vp[i][0],vpall[i][0],vp[i][1],vpall[i][1],vp[i][2],vpall[i][2],vp[i][3],vpall[i][3]);
    if (vp[i] != vpall[i]) {
        fprintf(stderr,"Construction failed\n");
        abort();
    }
  }

  TemporalGraph *tg;

  switch(opts.typegraph) {
      case kInterval:
        tg = new IntervalContactGraph();
      break;
      case kGrowth:
        tg = new GrowingContactGraph();
        break;
      case kPoint:
        tg = new PointContactGraph();
        break;
      case kIntervalPro:
          tg = new IntervalContactGraphImproved();
        break;
  }

    if (opts.typegraph == kIntervalPro) {
        CompactQtree *cqcurr;

        switch(opts.ds) {
          case ePRBlack:
              cqcurr = new PRBCompactQtree(vpcurr,opts.bs,opts.bb,opts.k1,opts.k2,opts.lk1,opts.lki);
            break;
          case ePRB2Black:
              cqcurr = new PRB2CompactQtree(vpcurr,opts.bs,opts.bb,opts.bc,opts.k1,opts.k2,opts.F,opts.lk1,opts.lki);
                break;
      //    case ePRWhite:
      //      cq = new PRWCompactQtree(vp,bs,bb,opts.k1,opts.k2,opts.lk1,opts.lki);
      //      break;
          case eMXDepth:
              cqcurr = new MXCompactQtree(vpcurr,opts.bs,opts.k1,opts.k2,opts.lk1,opts.lki);
            break;
          case eMXFixed:
              cqcurr = new MXCompactQtreeFixed(vpcurr,opts.bs,opts.bb,opts.k1,opts.k2,opts.lk1,opts.lki,opts.lf);
            break;
        }




        //First check points...
        vector<Point<uint> > vpallcurr;
        cqcurr->all(vpallcurr);

        if (vpcurr.size() != vpallcurr.size()) {
          fprintf(stderr, "Error: data from data structure doesnt match the input (current graph) \n");
          abort();
        }
        for(size_t i=0; i < vpcurr.size(); i++) {
          if (i%1000==0)fprintf(stderr, "Checking data: %.2f%% \r", (float)i/contacts*100);
          //assert(vp[i] == vpall[i]);
          //printf("%d %d %d %d %d %d %d %d\n",vp[i][0],vpall[i][0],vp[i][1],vpall[i][1],vp[i][2],vpall[i][2],vp[i][3],vpall[i][3]);
          if (vpcurr[i] != vpallcurr[i]) {
              fprintf(stderr,"Construction failed (current graph)\n");
              abort();
          }
        }


        IntervalContactGraph *tgpast = new IntervalContactGraph();
        GrowingContactGraph *tgcurr = new GrowingContactGraph();

        tgpast->setInfo(nodes,edges,lifetime,vp.size());
        tgcurr->setInfo(nodes,edges,lifetime,vpcurr.size());

        tgpast->setDs(cq);
        tgcurr->setDs(cqcurr);

        ((IntervalContactGraphImproved *)tg)->setGraphs(tgpast,tgcurr);


    }
    else {
      tg->setDs(cq);
    }

    ofstream file;
         LOG("Saving graph file in '%s'", opts.outfile);

         fprintf(stderr, "Saving data structure\n");
         file.open(opts.outfile, ios::binary);


         tg->setInfo(nodes,edges,lifetime,contacts);

         tg->save(file);

         file.close();

         delete tg;



  return 0;

}

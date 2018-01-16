#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "util/exception.h"

struct options {
	//-> required parameter
	char input[65532];
	char peer[65532];
	char output[65532];
	//-> key parameter
	int radius;
	int level;
	float scale0;
	int neib;
	//-> vice parameter
	int verbose;
	int test;
	int mode;
	int kmer;
};

inline int GetOpts(int argc, char **argv, options* opts_) {

    static struct option longopts[] = {
	{ "help",            no_argument,            NULL,              'h' },
	{ "input",     required_argument,      NULL,    			  'i' },
	{ "peer",    	     required_argument,      NULL,      	'p' },
	{ "output",      	required_argument,      NULL,              'o' },
	{ "radius",      required_argument,      NULL,              'r' },
	{ "level",      required_argument,      NULL,              'l' },
	{ "scale",      required_argument,      NULL,              's' },
	{ "verbose",         no_argument,            NULL,              'v' },
	{ "test",            no_argument,            NULL,              't' },
	{ "mode",            no_argument,            NULL,              'm' },
	{ "kmer",            no_argument,            NULL,              'k' },
	{ NULL,              0,                      NULL,               0 }
    };

    if((argc !=5 && argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17 && argc != 19 && argc != 21 ) && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1) {
	EX_TRACE("----------- cwDTW_genome ---------- \n");
	EX_TRACE("version v0.03 (NOV 4 2017) \n");
	EX_TRACE("-------------------------------------------------------------\n");
	EX_TRACE("required:\n");
        EX_TRACE("[-i INPUT GENOME 1][-p INPUT GENOME 2][-o OUTPUT] \n");
	EX_TRACE("optional:\n");
	EX_TRACE("([-r RADIUS])([-l LEVEL])([-s SCALE])([-R NEIBOR])\n");
//	EX_TRACE("([-v verbose])([-t test])([-m mode]) \n");
	EX_TRACE("-------------------------------------------------------------\n");
	EX_TRACE("[note]: by default, r=15, l=4, s=sqrt(2), R=25 \n");
//	EX_TRACE("[note]: by default, r=50, l=3, s=sqrt(2), v=0, t=0, m=0 \n");
	EX_TRACE("        for more detailed description, type '-h'  \n");
        return -1;
    }

    int ch;
    while((ch = getopt_long(argc, argv, "hi:p:o:r:l:s:R:v:t:m:k:", longopts, NULL))!= -1) {
        switch(ch) {

        case '?':
	{
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;
	}

        case ':':
	{
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;
	}

        case 'h':
	{
            EX_TRACE("----------- cwDTW_genome ---------- \n"
                     "version v0.03 (NOV 4 2017) \n"
                     "-------------------------------------------------------------\n"
                     "required:\n"
                     "[-i INPUT GENOME 1][-p INPUT GENOME 2][-o OUTPUT] \n"
                     "optional:\n"
                     "([-r RADIUS])([-l LEVEL])([-s SCALE])([-R NEIBOR])\n"
//                     "([-v verbose])([-t test])([-m mode]) \n"
                     "-------------------------------------------------------------\n"
                     "**** required: ******\n"
                     "INPUT GENOME 1: the first input genome, such as ATCG...;\n"
                     "INPUT GENOME 2: the second input genome, such as TCGA...;\n"
                     "OUTPUT:   genome alignment; if not specified, then no output be generated;  \n"
                     "**** key parameters: ******\n"
                     "RADIUS:   warp search radius (default 15);\n"
                     "LEVEL:    sampling level in continous wavelet (default 4);\n"
                     "SCALE:    base scale in continous wavelet (default sqrt(2));\n"
                     "NEIBOR:   final boundary to align two DNA sequences (default 25);\n");
//                     "**** vice parameters: ******\n"
//                     "verbose:  0 for NO screenout message, 1 for screenout (default 0);\n"
//                     "test:     test mode. 0 not_use; 1 equal_ave; 2 peak_ave; 3 FastDTW (default 0) \n"
//                     "mode:     bound mode. 0 block_mode; 1 diagonol_mode (default 0) \n"
//                     "kmer:     kmer pore model. 0 for 5mer; 1 for 6mer (default 0) \n"
            return -1;
	}

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'p':
        {
            std::istringstream iss(optarg);
            iss >> opts_->peer;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
		
	case 'r':
        {
            std::istringstream iss(optarg);
            iss >> opts_->radius;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
		
	case 'l':
        {
            std::istringstream iss(optarg);
            iss >> opts_->level;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
		
	case 's':
        {
            std::istringstream iss(optarg);
            iss >> opts_->scale0;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'R':
        {
            std::istringstream iss(optarg);
            iss >> opts_->neib;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

	case 'v':
	{
             std::istringstream iss(optarg);
             iss >> opts_->verbose;
             if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
             }
	}
	break;

	case 't':
	{
             std::istringstream iss(optarg);
             iss >> opts_->test;
             if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
             }
	}
	break;

	case 'm':
	{
             std::istringstream iss(optarg);
             iss >> opts_->mode;
             if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
             }
	}

        case 'k':
        {
             std::istringstream iss(optarg);
             iss >> opts_->kmer;
             if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
             }
        }


        case 0:
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif


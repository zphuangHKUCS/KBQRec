////
//  main.cpp
//   EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//


#include <iostream>
#include "EQFG.h"
#include "yagoReader.h"
#include "HIN_Graph.h"
#include <time.h>
#include <cstring>
#include <fstream>

using namespace std;

//#define MAC
#define LINUX


#define QUERYFOLLOWPATH "/Users/huangzhipeng/research/term-query-graph-master/data/AOL/01_QueryFollow.txt"
#define ENTITYFOLLOWPATH "/Users/huangzhipeng/research/term-query-graph-master/data/AOL/01_EntityFollow.txt"
#define QUERYENTITYPATH "/Users/huangzhipeng/research/term-query-graph-master/data/AOL/01_queryEntity.txt"
#define GRAPHDATAPATH "/Users/huangzhipeng/research/term-query-graph-master/data/AOL/graph-data/"
#define QUERYCOUNTPATH "/Users/huangzhipeng/research/term-query-graph-master/data/AOL/01_queryCount.txt"

#define DATADIR "/Users/huangzhipeng/XcodeWorkPlace/EQFG/data/AOL01/"


bool buildQueryLog(string graphDataDir, string queryEntityPath, string queryCountPath, string outPath)
{
    EQFG eqfg(graphDataDir, queryEntityPath, queryCountPath);
    cerr << "query log building complete!" << endl;
    eqfg.saveToFiles(outPath);
    return true;
}
bool trainMetaPath(string queryLogDir, string outPath)
{
    EQFG eqfg(queryLogDir);
    
    eqfg.trainMetaPaths();
    eqfg.saveMetaPaths(outPath);
    return true;
}
bool trainAllMetaPath(string queryLogDir, string outDir)
{
    EQFG eqfg(queryLogDir);
    eqfg.trainMetaPaths_all(outDir);
    return true;
}
void rec_QFG(string queryLogDir, string queryPath)
{
	EQFG eqfg(queryLogDir);
	eqfg.rec_QFG_fromFile(queryPath);
}
void rec_EQFG_t(string queryLogDir, string queryPath, string cachePath)
{
        EQFG eqfg(queryLogDir);
        eqfg.rec_EQFG_t_fromFile(queryPath, cachePath);
}


void rec_EQFG(string queryLogDir, string queryPath, string cachePath)
{
	EQFG eqfg(queryLogDir);
	eqfg.rec_EQFG_fromFile(queryPath, cachePath);
}

void rec_P(string queryLogDir,string metapathPath, string cachePath, string queryPath)
{
	EQFG eqfg(queryLogDir);
	eqfg.loadYagoByDefault();
	eqfg.rec_P_fromFile(queryPath, cachePath, metapathPath);
}
void rec_QFGP(string queryLogDir,string metapathPath, string cachePath, string queryPath, string ratio)
{
	EQFG eqfg(queryLogDir);
	eqfg.loadYagoByDefault();
	eqfg.rec_QFGP_fromFile(queryPath, cachePath, metapathPath, atof(ratio.c_str()));
}
void rec_EQFGP(string queryLogDir,string metapathPath, string cachePath, string queryPath, string ratio)
{
        EQFG eqfg(queryLogDir);
        eqfg.loadYagoByDefault();
        eqfg.rec_EQFGP_fromFile(queryPath, cachePath, metapathPath, atof(ratio.c_str()));
}
void rec_QFGPEQFG(string queryLogDir, string cachePath, string queryPath, string ratio)
{
        EQFG eqfg(queryLogDir);
        //eqfg.loadYagoByDefault();
        eqfg.rec_QFGPEQFG_fromFile(queryPath, cachePath, atof(ratio.c_str()));
}

void rec_QTFGP(string queryLogDir,string metapathPath, string cachePath, string queryPath, string QTFGpath, string ratio)
{
        EQFG eqfg(queryLogDir);
        eqfg.loadYagoByDefault();
        eqfg.rec_QTFGP_fromFile(queryPath, cachePath, metapathPath, QTFGpath, atof(ratio.c_str()));
}

void rec_QTFGDQGP(string queryLogDir,string metapathPath, string cachePath, string queryPath, string QTFGpath, string DQGPath, string ratio)
{
        EQFG eqfg(queryLogDir);
        eqfg.loadYagoByDefault();
        eqfg.rec_QTFG_DQGP_fromFile(queryPath, cachePath, metapathPath, QTFGpath, DQGPath, atof(ratio.c_str()));
}
void rec_DQG(string queryLogDir, string cachePath, string queryPath)
{
	EQFG eqfg(queryLogDir);
	//eqfg.loadYagoByDefault();
	//eqfg.rec_QTFGP_fromFile(queryPath, cachePath, metapathPath, QTFGpath, atof(ratio.c_str()));
	eqfg.loadDQG(queryLogDir);
	
	eqfg.rec_DQG_fromFile(queryPath, cachePath);
}

void entityExpand_exp(string querylogPath, string queryEntityPath, string metapathPath)
{
	EQFG eqfg(querylogPath);
	eqfg.loadYagoByDefault();
	eqfg.expand_exp(queryEntityPath, metapathPath);
}



void printUsage(const char * argv[])
{
	cout << "Usage:" << endl;
	cout << argv[0] << " -b graphDataDir queryEntityPath queryCountPath outPath\tfor building the query log." << endl;
	cout << argv[0] << " -t queryLogDir outPath\tfor training the meta path." << endl;
	cout << argv[0] << " -r QFG queryLogDir queryPath\tfor recommanding through QFG." << endl;
	cout << argv[0] << " -r EQFG quertLogDir queryPath cachePath\tfor recommanding through EQFG only. (the queries should be after entity linking)" << endl;
	cout << argv[0] << " -r EQFG_t queryLogDir queryPath cachePath\tfor recommanding through EQFG. (use center piece method)" << endl;
	cout << argv[0] << " -r P queryLogDir metapathPath cachePath queryPath\tfor recommanding through YAGO only." << endl;
	cout << argv[0] << " -r EQFGP_t queryLogDir metapathPath cachePath queryPath\tfor recommanding through EQFGP.(use center piece method)" << endl;
	cout << argv[0] << " -r QFGP queryLogDir metapathPath cachePath queryPath ratio\tfor recommanding through QFGP." << endl;
	cout << argv[0] << " -r EQFGP queryLogDir metapathPath cachePath queryPath ratio\tfor recommanding through EQFGP." << endl;
	cout << argv[0] << " -r QFGPEQFG queryLogDir cachePath queryPath ratio\tfor recommanding through QFG+EQFG." << endl;
	cout << argv[0] << " -r QTFGP queryLogDir metapathPath cachePath queryPath QTFGresultPath ratio\tfor recommanding through EQFGP." << endl;
	cout << argv[0] << " -r DQG queryLogDir cachePath queryPath \tfor recommanding through DQG." << endl;
	cout << argv[0] << " -r QGFGDQGP queryLogDir metapathPath cachePath queryPath QTFGresultPath DQGresultPath ratio\tfor recommanding through QTFGDQGP." << endl;
	cout << argv[0] << " -e EE queryLogDir queryEntityPath metapathPath\t for experiment on entity expansion." << endl;
	cout << argv[0] << " -e PPR queryLogDir metapathPath cachePath queryPath QTFGresultPath ratio \t for experiments PPR" << endl;
	cout << argv[0] << " -T queryLogDir outDir\tfor training all the heuristics." << endl;
	cout << argv[0] << " -c testpath result1 result2 alpha \t for combining two results with alpha." << endl;
	cout << argv[0] << " -C testpath result1 result2 result3 alpha beta \t for combining three results." << endl;
    cout << argv[0] << " -j testpath result1 result2 \t for joining two results with piority." << endl;
    cout << argv[0] << " -J testpath result1 result2 result3 \t for combining three results with piority." << endl;
}

#ifdef LINUX
int main(int argc, const char * argv[]) {
	if(argc < 2){
		printUsage(argv);
		return -1;
	}else if(argc == 6 && argv[1][1] == 'b'){
		buildQueryLog(argv[2], argv[3], argv[4], argv[5]);
		return 0;
	}else if(argc == 4 && argv[1][1] == 't'){
		trainMetaPath(argv[2], argv[3]);
		return 0;
	}else if(argc == 4 && argv[1][1] == 'T'){
		trainAllMetaPath(argv[2], argv[3]);
		return 0;
	}else if(argv[1][1] == 'r'){
		if(strcmp(argv[2], "QFG") == 0){
			rec_QFG(argv[3], argv[4]);
			return 0;
		}
		else if(strcmp(argv[2], "EQFG") == 0){
			rec_EQFG(argv[3], argv[4], argv[5]);
		}else if(strcmp(argv[2], "EQFG_t") == 0){
            rec_EQFG_t(argv[3], argv[4], argv[5]);
        }else if(strcmp(argv[2], "P") == 0){
			rec_P(argv[3],argv[4], argv[5], argv[6]);
        }else if(strcmp(argv[2], "QFGP") == 0){
			rec_QFGP(argv[3],argv[4], argv[5], argv[6], argv[7]);
		}else if(strcmp(argv[2], "EQFGP") == 0){
			rec_EQFGP(argv[3],argv[4], argv[5], argv[6], argv[7]);
		}else if(strcmp(argv[2], "QFGPEQFG") == 0){
			rec_QFGPEQFG(argv[3],argv[4], argv[5], argv[6]);
		}else if(strcmp(argv[2], "QTFGP") == 0){
			rec_QTFGP(argv[3],argv[4], argv[5], argv[6], argv[7], argv[8]);
		}else if(strcmp(argv[2], "DQG") == 0){
			rec_DQG(argv[3], argv[4], argv[5]);
		}else if(strcmp(argv[2], "QTFGDQGP") == 0){
			rec_QTFGDQGP(argv[3],argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
		}
	}else if(argv[1][1] == 'e'){
		// for experiment
		if(strcmp(argv[2], "EE") == 0){
			entityExpand_exp(argv[3], argv[4], argv[5]);
		}	
		else if(strcmp(argv[2], "PPR") == 0){
			rec_QTFGP(argv[3],argv[4], argv[5], argv[6], argv[7], argv[8]);
		}
	}else if(argv[1][1] == 'c'){
		EQFG::combineTwoResult(argv[2], argv[3], argv[4], atof(argv[5]));
    }else if(argv[1][1] == 'j'){
        EQFG::joinTwoResult(argv[2], argv[3], argv[4]);
    }else if(argv[1][1] == 'C'){
		EQFG::combineThreeResult(argv[2], argv[3], argv[4], argv[5], atof(argv[6]), atof(argv[7]));
	}
	else{
		printUsage(argv);
		return -1;
	}
    return 0;
}
#endif
#ifdef MAC

int main()
{
    
    return 0;
}

#endif

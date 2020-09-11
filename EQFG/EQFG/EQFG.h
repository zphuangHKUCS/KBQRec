//
//  EQFG.h
//  EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//

#ifndef __EQFG__EQFG__
#define __EQFG__EQFG__

#include <stdio.h>
#include <cstring>
#include <vector>
#include <map>
#include <iostream>

#include "HIN_Graph.h"
#include "Meta-Structure.h"
#include "pruned_landmark_labeling.h"


#define EQFG_NUMOF_E 10
#define PPR_E_ALPHA 0.15
#define PPR_Q_ALPHA 0.15

#define MAXDISTANCE_T 2
#define PPR_EPS_Q 0.001
#define PPR_EPS_E 0.001


using namespace std;

class RecNode
{
public:
    int qid_;
    double w_;
    RecNode(int qid, double w);
    
    bool operator < (const RecNode & r);
};

class EQFG_Edge
{
public:
    int sid_, eid_;
    double w_;
    EQFG_Edge(int id1, int id2, double w);
    EQFG_Edge(const EQFG_Edge & e);
};

class EQFG_QNode
{
public:
    int id_;
    int count_;
    vector<EQFG_Edge> toQueryEdges_; // outedges
    //vector<EQFG_Edge> inQueryEdges_; // inedges
    vector<EQFG_Edge> toEntityEdges_;
    EQFG_QNode(int id);
};


class EQFG_ENode
{
public:
    int id_;
    
    vector<EQFG_Edge> toQueryEdges_;
    vector<EQFG_Edge> toEntityEdges_;
    //vector<EQFG_Edge> inEntityEdges_;
    EQFG_ENode(int id);
};

class EQFG
{
private:
	double ratio_;
    bool entityShortestPath(int id1, int id2, int dist, vector<int> & pathNodes);
    bool entityAllShortestPaths(int id1, int id2, int dist, vector< vector< vector<int> > > &pathNodes);
    void getMetaPath(const vector<int> & pathNodes, vector<int> & nodeTypes, vector<int> & edgeTypes);
    void insertMetaPath(int eid, const vector<int> & nodeTypes, const vector<int> & edgeTypes, double w);
    void insertMetaPathToOneList(int eid, const vector<int> & nodeTypes, const vector<int> & edgeTypes, double w, map<int, map<int, double> > & PathList);
    void selectTwoTypeEntity(int t1, int t2);
    int getYAGOEID(int eid);
    
    string cachePath_EQFG_;
    
    vector<int> expandEntity(const vector<int> & eids, int k = EQFG_NUMOF_E); // expand the entity set as EQFG
	
    void PPR_E(int eid, map<int, double> & ret, double alpha = PPR_E_ALPHA);
    void PPR_Q(int eid, map<int, double> & ret, double alpha = PPR_Q_ALPHA);
    string metapath2string(const vector<int> & metapath);
    vector<int> string2metapath(string s);

	
    double timeForPPR, timeForExp;
    
	
public:
    HIN_Graph hin_;
    PrunedLandmarkLabeling<> pll_;
    
    map<string, int> query2id_;
    map<string, int> entity2id_;
    map<string, int> doc2id_;
	
	map<int, map<int, double>> entity2docPro_;
	map<int, map<int, double>> doc2entityPro_;
	map<int, map<int, double>> doc2queryPro_;
	
    vector<EQFG_QNode> QNodes_;
    vector<EQFG_ENode> ENodes_;
    
    vector<string> queries_;
    vector<string> entities_;
	vector<string> documents_;
    
    map<int, map<int, double> > MetaPathsLists_; // for each class, we save a set of meta paths
    map<string, int> metapath_2id_;
    vector<vector<int> > metapaths_;
    
    vector<int> YAGOEID;
    
    EQFG(string queryFollowPath, string queryEntityPath, string queryCountPath);
    EQFG(string indexPath);
    void saveToFiles(string dirPath);
    
    vector<RecNode> rec_QFG(string query); // the most original one: query flow graph
    vector<RecNode> rec_EQFG(int qid, vector<int> eids); // EQFG
    vector<RecNode> rec_EQFG_t(int qid, vector<int> eids); // EQFG but centerpiece
    vector<RecNode> rec_P(int qid, vector<string> entities); // Yago only
	
	vector<RecNode> rec_QTFGP(int qid, vector<string> entities, map<string, double> & rec); // QFG Plus yago
	vector<RecNode> rec_QFGP(int qid, vector<string> entities); // QFG Plus yago
	vector<RecNode> rec_EQFGP(int qid, vector<string> entities); // yago
	vector<RecNode> rec_DQG(int qid, vector<string> entities); // DQG
    vector<RecNode> rec_QFGPEQFG(int qid, vector<string> entities);

    void rec_EQFGP_t_fromFile(string queryEntityPath);
    // should not use this name
	void rec_P_fromFile(string queryEntityPath, string cachePath, string metapathPath);
    void rec_EQFG_fromFile(string queryEntityPath, string cachePath);
    void rec_EQFG_t_fromFile(string queryEntityPath, string cachePath);
    void rec_QFG_fromFile(string queryPath);
	
	void rec_QFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, double ratio);
	void rec_EQFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, double ratio);
	void rec_QFGPEQFG_fromFile(string queryEntityPath, string cachePath, double ratio);
	void rec_QTFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, string QTFGPath, double ratio);
	
	void rec_QTFG_DQGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, string QTFGPath, string DQGPath, double ratio);
	
	void rec_DQG_fromFile(string queryEntityPath, string cachePath);
    
    void showRecNodes(string query, const vector<RecNode> & v, int k = 5);
    
	map<int, double> expandEntity_sp(const vector<string> & entities, int k = EQFG_NUMOF_E); // expand the entity set with shortest path
	
    void loadYagoByDefault();
	void loadDQG(string qlDir);
    void trainMetaPaths();
	void trainMetaPaths_all(string ourDir);
    void saveMetaPaths(string path);
    void saveMetaPathsList(string path, map<int, map<int, double> > & PathList);
    void loadMetaPaths(string path);
	
	void expand_exp(string queryEntityPath, string metapathPath);
    
	static void combineTwoResult(string testpath, string resultpath1, string resultpath2, double alpha);
	static void combineThreeResult(string testpath, string resultpath1, string resultpath2, string resultpath3, double alpha, double beta);

	static void joinTwoResult(string testpath, string resultpath1, string resultpath2);
};



#endif /* defined(__EQFG__EQFG__) */

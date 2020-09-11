//
//  EQFG.cpp
//  EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//

#include "EQFG.h"
#include "BoundHeap.h"
#include "SimCalculator.h"
#include "Meta-Structure.h"
#include <cstring>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/heap/fibonacci_heap.hpp>

#include "Tools.h"

using namespace std;

bool pairIntDoubleLargerCMP(const pair<int, double>& a1, const pair<int, double>& a2)
{
    return a1.second > a2.second;
}

bool pairStringDoubleLargerCMP(const pair<string, double>& a1, const pair<string, double>& a2)
{
    return a1.second > a2.second;
}

EQFG_QNode::EQFG_QNode(int id): id_(id)
{

}
EQFG_ENode::EQFG_ENode(int id): id_(id)
{

}
EQFG_Edge::EQFG_Edge(int id1, int id2, double w): sid_(id1), eid_(id2), w_(w)
{

}
EQFG_Edge::EQFG_Edge(const EQFG_Edge & e)
{
    sid_ = e.sid_;
    eid_ = e.eid_;
    w_ = e.w_;
}

RecNode::RecNode(int id, double w):qid_(id), w_(w)
{

}


bool RecNode::operator<(const RecNode &r)
{
    return w_ < r.w_;
}
bool lessRecNode(const RecNode & m1, const RecNode & m2) {
    return m1.w_ > m2.w_;
}

void EQFG::showRecNodes(string query, const vector<RecNode> &v, int k)
{
    cout << query;
    for(int i = 0; i < v.size() && i < k; ++i){
        cout << '\t' << queries_[v[i].qid_] << '\t' << v[i].w_;
    }
    cout << endl;
}

void EQFG::loadDQG(string qlDir)
{
	// entity2docPro_
	map<int, map<int, double>> entity2docCount;
	cerr << "start loading the DQG." << endl;
    string temps = qlDir + "doc2que_count.txt";
    ifstream doc2queIn(temps.c_str(), ios::in);
	string line;
    while(getline(doc2queIn, line)){
		//cerr << line << endl;
        vector<string> strs = split(line, "\t");
        if(query2id_.find(strs[1]) == query2id_.end()){
			//cerr << strs[1] << "\t not found!" << endl;
			continue;
		}	
		int qid = query2id_[strs[1]];
		int count = atoi(strs[2].c_str());
		if (doc2id_.find(strs[0]) == doc2id_.end()){
			doc2id_[strs[0]] = documents_.size();
			documents_.push_back(strs[0]);
		}
		int docid = doc2id_[strs[0]];
		//if (doc2queryPro_.find(docid) == doc2queryPro_.end()){
		//	doc2queryPro_[docid] = map<int, double>();
		//}
		//doc2queryPro_[docid][qid] = count;
		for(int i = 0; i < QNodes_[qid].toEntityEdges_.size(); ++i){
			int eid = QNodes_[qid].toEntityEdges_[i].sid_;
			double w = QNodes_[qid].toEntityEdges_[i].w_ * count;
			if(entity2docCount.find(eid) == entity2docCount.end()){
				entity2docCount[eid] = map<int, double>();
			}
			entity2docCount[eid][docid] = w;
		}
		
		// dealing with doc2entityPro_
		if (doc2entityPro_.find(docid) == doc2entityPro_.end()){
			doc2entityPro_[docid] = map<int,double>();
		}
		for(int i = 0; i < QNodes_[qid].toEntityEdges_.size(); ++i){
			int eid = QNodes_[qid].toEntityEdges_[i].sid_;
			double w = (0.0 + count) / QNodes_[qid].toEntityEdges_.size();
			if(doc2entityPro_[docid].find(eid) == doc2entityPro_[docid].end()){
				doc2entityPro_[docid][eid] = 0.0;
			}
			doc2entityPro_[docid][eid] += w;
		}
    }
    doc2queIn.close();
	
	//for(map<int, map<int, double>>::iterator i = doc2queryPro_.begin(); i != doc2queryPro_.end(); ++i){
		//int docid = i->first;
		//double total = 0.0;
		//for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
		//	total += j->second;
		//}
		//for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
		//	j->second /= total;
		//}
	//}
	
	for(map<int, map<int, double>>::iterator i = doc2entityPro_.begin(); i != doc2entityPro_.end(); ++i){
		int docid = i->first;
		double total = 0.0;
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			total += j->second;
		}
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			j->second /= total;
		}
	}
	
	/*
	for(map<int, map<int, double>>::iterator i = doc2queryPro_.begin(); i != doc2queryPro_.end(); ++i){
		cerr << documents_[i->first];
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			cerr << '\t' << queries_[j->first];
		}
		cerr << endl;
	}
	
	*/
	//cerr << "size of doc2query: " << doc2queryPro_.size() << endl;
	
	for(map<int, map<int, double>>::iterator i = entity2docCount.begin(); i != entity2docCount.end(); ++i){
		int eid = i->first;
		double total = 0.0;
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			total += j->second;
		}
		if (entity2docPro_.find(eid) == entity2docPro_.end()){
			entity2docPro_[eid] = map<int, double>();
		}
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			int docid = j->first;
			entity2docPro_[eid][docid] = j->second / total;
		}
	}
	
	cerr << "size of entity2doc: " << entity2docPro_.size() << endl;
	
	
	for(map<int, map<int, double>>::iterator i = doc2entityPro_.begin(); i != doc2entityPro_.end(); ++i){
		BoundHeap entity2docHeap(5);
		map<int, double> tempMap;
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			entity2docHeap.push(make_pair(j->first, j->second));
		}
		while(entity2docHeap.size() > 0){
			pair<int, double> p = entity2docHeap.pop();
			tempMap[p.first] = p.second;
		}
		i->second = tempMap;
	}
	for(map<int, map<int, double>>::iterator i = entity2docPro_.begin(); i != entity2docPro_.end(); ++i){
		BoundHeap entity2docHeap(5);
		map<int, double> tempMap;
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			entity2docHeap.push(make_pair(j->first, j->second));
		}
		while(entity2docHeap.size() > 0){
			pair<int, double> p = entity2docHeap.pop();
			tempMap[p.first] = p.second;
		}
		i->second = tempMap;
	}
	
	
	
	/*
	for(map<int, map<int, double>>::iterator i = entity2docPro_.begin(); i != entity2docPro_.end(); ++i){
		cerr << entities_[i->first];
		for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
			cerr << '\t' << documents_[j->first];
		}
		cerr << endl;
	}
	*/
	//cin.get();
}

EQFG::EQFG(string graph_dataPath, string queryEntityPath, string queryCountPath)
{	
    if(graph_dataPath[graph_dataPath.size() -1] != '/')
	graph_dataPath += "/";
    string line;
    // read the query2id file
    string query2idPath = graph_dataPath + "query-id-map.tsv";
    cerr << "start reading " << query2idPath << endl;
    ifstream query2idIn(query2idPath.c_str(), ios::in);
    int maxID = -1;
    while(getline(query2idIn, line)){
	//cerr << line << endl;
        vector<string> strs = split(line, "\t");
        int id = atoi(strs[1].c_str());
        if(id > maxID)
            maxID = id;
        string query = strs[0];
        query2id_[query] = id;
    }
    query2idIn.close();
    queries_ = vector<string>(maxID + 1);
    for(map<string, int>::iterator i = query2id_.begin(); i != query2id_.end(); ++i){
        queries_[i->second] = i->first;
    }
    // add QNodes
    for(int i = 0; i < queries_.size(); ++i){
        QNodes_.push_back(EQFG_QNode(i));
    }
    // add Q->Q edges and weights
    
    string QQEdgePath = graph_dataPath + "query-rewrite-matrix.tsv";
    //cerr << "start reading " << QQEdgePath << endl;
    ifstream edgeIn(QQEdgePath.c_str(), ios::in);
    while(getline(edgeIn, line)){
        vector<string> strs = split(line, "\t");
        int sid = atoi(strs[0].c_str());
        for(int i = 1; i < strs.size(); i += 2){
            int eid = atoi(strs[i].c_str());
            double w = atof(strs[i+1].c_str());
            EQFG_Edge tempEdge(sid, eid, w);
            QNodes_[sid].toQueryEdges_.push_back(tempEdge);
            //QNodes_[eid].inQueryEdges_.push_back(tempEdge);
        }
    }
    edgeIn.close();
    cerr << "size of some_Query:\t" << QNodes_.size() << endl;
    // read the queryEntity file
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    cerr << "start reading " << queryEntityPath << endl;
    while(getline(queryEntityIn, line)){
        vector<string> strs = split(line, "\t");
        string query = strs[0];
        if(query2id_.find(query) == query2id_.end()){
            int id = QNodes_.size();
            query2id_[query] = id;
            QNodes_.push_back(EQFG_QNode(id));
            queries_.push_back(query);
        }
        string entity = strs[1];
        if(entity2id_.find(entity) == entity2id_.end()){
            int id = ENodes_.size();
            entity2id_[entity] = id;
            ENodes_.push_back(EQFG_ENode(id));
            entities_.push_back(entity);
        }
        int qid = query2id_[query];
        int eid = entity2id_[entity];
        EQFG_Edge tempEdge(eid, qid, 1.0);
        ENodes_[eid].toQueryEdges_.push_back(tempEdge);
        QNodes_[qid].toEntityEdges_.push_back(tempEdge);
    }
    queryEntityIn.close();
    cerr << "size of all_Query:\t" << QNodes_.size() << endl;
    cerr << "start reading " << queryCountPath << endl;
    // set the weight of query-entity edges
    ifstream queryCountIn(queryCountPath.c_str(), ios::in);
    while(getline(queryCountIn, line)){
        vector<string> strs = split(line, "\t");
        int qid = query2id_[strs[0]];
        int count = atoi(strs[1].c_str());
        QNodes_[qid].count_ = count;
    }
    queryCountIn.close();
    
    for(int i = 0; i < ENodes_.size(); ++i){
        int totalCount = 0;
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            int qid = ENodes_[i].toQueryEdges_[j].eid_;
            totalCount += QNodes_[qid].count_;
        }
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            int qid = ENodes_[i].toQueryEdges_[j].eid_;
            ENodes_[i].toQueryEdges_[j].w_ = (double)(QNodes_[qid].count_) / totalCount;
        }
    }
    cerr << "start setting the weight of e2e edges." << endl;
    // set the weight of entity-entity edges
    map<pair<int, int>, double> tempMap;
    for(int i = 0; i < QNodes_.size(); ++i){
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            int qid1 = i;
            int qid2 = QNodes_[i].toQueryEdges_[j].eid_;
            if(QNodes_[qid1].toEntityEdges_.size() == 0 || QNodes_[qid2].toEntityEdges_.size() == 0)
                continue;
            for(int k1 = 0; k1 < QNodes_[qid1].toEntityEdges_.size(); ++k1){
                for(int k2 = 0; k2 < QNodes_[qid2].toEntityEdges_.size(); ++k2){
                    int eid1 = QNodes_[qid1].toEntityEdges_[k1].sid_;
                    int eid2 = QNodes_[qid2].toEntityEdges_[k2].sid_;
                    EQFG_Edge tempEdge(eid1, eid2, 0.0);
                    bool flag = true;
                    for(int kk = 0; kk < ENodes_[eid1].toEntityEdges_.size(); ++ kk){
                        if(ENodes_[eid1].toEntityEdges_[kk].eid_ == tempEdge.eid_){
                            flag = false;
                            break;
                        }
                    }
                    if(flag){
                        ENodes_[eid1].toEntityEdges_.push_back(tempEdge);
                    }
                    
                    tempMap[make_pair(eid1, eid2)] = 1.0;
                }
            }
            
        }
    }
    for(int i = 0; i < QNodes_.size(); ++i){
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            int qid1 = i;
            int qid2 = QNodes_[i].toQueryEdges_[j].eid_;
            if(QNodes_[qid1].toEntityEdges_.size() == 0 || QNodes_[qid2].toEntityEdges_.size() == 0)
                continue;
            double tobetimes = 1.0 - QNodes_[i].toQueryEdges_[j].w_ / (QNodes_[qid1].toEntityEdges_.size() * QNodes_[qid2].toEntityEdges_.size());
            for(int k1 = 0; k1 < QNodes_[qid1].toEntityEdges_.size(); ++k1){
                for(int k2 = 0; k2 < QNodes_[qid2].toEntityEdges_.size(); ++k2){
                    int eid1 = QNodes_[qid1].toEntityEdges_[k1].sid_;
                    int eid2 = QNodes_[qid2].toEntityEdges_[k2].sid_;
                    tempMap[make_pair(eid1, eid2)] *= tobetimes;
                }
            }
        }
    }
    for(map<pair<int, int>, double>::iterator i = tempMap.begin(); i != tempMap.end(); ++i){
        i->second = 1.0 - i->second;
    }
    for(int i = 0; i < ENodes_.size(); ++i){
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            int id2 = ENodes_[i].toEntityEdges_[j].eid_;
            ENodes_[i].toEntityEdges_[j].w_ = tempMap[make_pair(i, id2)];
        }
    }
    cerr << "end building the graph." << endl;
}

void EQFG::saveToFiles(string dirPath)
{
	cerr << "start saveing the query log to " << dirPath << endl;
    ofstream query2idOut(dirPath + "query2id.txt", ios::out);
    for(int i = 0; i < queries_.size(); ++i){
        query2idOut << queries_[i] << '\t' << i << endl;
    }
    query2idOut.close();
    
    ofstream entity2idOut(dirPath + "entity2id.txt", ios::out);
    for(int i = 0; i < entities_.size(); ++i){
        entity2idOut << entities_[i] << '\t' << i << endl;
    }
    entity2idOut.close();
    
    ofstream query2query_w_out(dirPath + "query2query_w.txt", ios::out);
    for(int i = 0; i < QNodes_.size(); ++i){
        query2query_w_out << i;
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            query2query_w_out << '\t' << QNodes_[i].toQueryEdges_[j].eid_ << '\t' << QNodes_[i].toQueryEdges_[j].w_;
        }
        query2query_w_out << endl;
    }
    query2query_w_out.close();
    
    ofstream entity2query_w_out(dirPath + "entity2query_w.txt", ios::out);
    for(int i = 0; i < ENodes_.size(); ++i){
        entity2query_w_out << i;
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            entity2query_w_out << '\t' << ENodes_[i].toQueryEdges_[j].eid_ << '\t' << ENodes_[i].toQueryEdges_[j].w_;
        }
        entity2query_w_out << endl;
    }
    entity2query_w_out.close();
    
    ofstream entity2entity_w_out(dirPath + "entity2entity_w.txt", ios::out);
    for(int i = 0; i < ENodes_.size(); ++i){
        entity2entity_w_out << i;
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            entity2entity_w_out << '\t' << ENodes_[i].toEntityEdges_[j].eid_ << '\t' << ENodes_[i].toEntityEdges_[j].w_;
        }
        entity2entity_w_out << endl;
    }
    entity2entity_w_out.close();
}

EQFG::EQFG(string indexPAth)
{
    string line;
    cerr << "start loading the query nodes." << endl;
    string temps = indexPAth + "query2id.txt";
    ifstream query2idIn(temps.c_str(), ios::in);
    while(getline(query2idIn, line)){
        vector<string> strs = split(line, "\t");
        QNodes_.push_back(EQFG_QNode(queries_.size()));
        query2id_[strs[0]] = queries_.size();
        queries_.push_back(strs[0]);
    }
    query2idIn.close();
    cerr << "start loading the entity nodes." << endl;
    temps = indexPAth + "entity2id.txt";
    ifstream entity2idIn(temps.c_str(), ios::in);
    while(getline(entity2idIn, line)){
        vector<string> strs = split(line, "\t");
        ENodes_.push_back(EQFG_ENode(entities_.size()));
        entity2id_[strs[0]] = entities_.size();
        entities_.push_back(strs[0]);
    }
    entity2idIn.close();
    cerr << "start loading the query2query edges." << endl;
    string tempPath = indexPAth + "query2query_w.txt";
    ifstream query2query_w_in(tempPath.c_str(), ios::in);
    while(getline(query2query_w_in, line)){
        vector<string> strs = split(line, "\t");
        int sid = atoi(strs[0].c_str());
        for(int i = 1; i < strs.size(); i += 2){
            EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i+1].c_str()));
            QNodes_[sid].toQueryEdges_.push_back(tempEdge);
        }
    }
    query2idIn.close();
    cerr << "start loading the entity2query edges." << endl;
    tempPath = indexPAth + "entity2query_w.txt";
    ifstream entity2queryIn(tempPath.c_str(), ios::in);
    while(getline(entity2queryIn, line)){
        vector<string> strs = split(line, "\t");
        int sid = atoi(strs[0].c_str());
        for(int i = 1; i < strs.size(); i += 2){
            EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i+1].c_str()));
            ENodes_[sid].toQueryEdges_.push_back(tempEdge);
            QNodes_[tempEdge.eid_].toEntityEdges_.push_back(tempEdge);
        }
    }
    entity2queryIn.close();
    cerr << "start loading the entity2entity edges" << endl;
    tempPath = indexPAth + "entity2entity_w.txt";
    ifstream entity2entityIn(tempPath.c_str(), ios::in);
    while(getline(entity2entityIn, line)){
        vector<string> strs = split(line, "\t");
        int sid = atoi(strs[0].c_str());
		double sum = 0.0;
        for(int i = 1; i < strs.size(); i += 2){
            EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i+1].c_str()));
			sum += atof(strs[i+1].c_str());
            ENodes_[sid].toEntityEdges_.push_back(tempEdge);
        }
		for(int i = 0; i < ENodes_[sid].toEntityEdges_.size(); i ++){
            ENodes_[sid].toEntityEdges_[i].w_ /= sum;
        }
    }
    entity2entityIn.close();
	cerr << "end of building the graph." << endl;
}

vector<RecNode> EQFG::rec_QFG(string query)
{
    vector<RecNode> ret;
    if(query2id_.find(query) == query2id_.end())
        return ret;
    int qid = query2id_[query];
    for(int i = 0; i < QNodes_[qid].toQueryEdges_.size(); ++i){
        ret.push_back(RecNode(QNodes_[qid].toQueryEdges_[i].eid_, QNodes_[qid].toQueryEdges_[i].w_));
    }
    sort(ret.begin(), ret.end(), lessRecNode);
    return ret;
}

void EQFG::PPR_Q(int eid, map<int, double> &ret, double alpha)
{
	
	stringstream ss;
	ss << cachePath_EQFG_ << eid << ".txt";
	string cacheP;
	ss >> cacheP;
    //string cacheP = cachePath_EQFG_;
	//cacheP += eid;
	//cacheP += ".txt";
    int hit = access(cacheP.c_str(), 0);
//	if(false){    
	if(hit == 0){
        ifstream in(cacheP.c_str(), ios::in);
        int id;
        double d;
        while(in >> id >> d){
            ret[id] = d;
        }
        in.close();
        return;
    }
    
    BoundHeap activeInk(999999);
    
    for(int i = 0; i < ENodes_[eid].toQueryEdges_.size(); ++i){
        if(queries_[ENodes_[eid].toQueryEdges_[i].eid_] == "-")
            continue;
        activeInk.push(make_pair(ENodes_[eid].toQueryEdges_[i].eid_, ENodes_[eid].toQueryEdges_[i].w_));
    }
    while(activeInk.size() > 0){
        pair<int, double> p = activeInk.pop();
        if(queries_[p.first] == "-")
            continue;
        if(ret.find(p.first) == ret.end()){
            ret[p.first] = 0.0;
        }
        ret[p.first] += p.second * alpha;
        int degree = QNodes_[p.first].toQueryEdges_.size();
        for(int j = 0; j < degree; ++j){
            int newQid = QNodes_[p.first].toQueryEdges_[j].eid_;
            double newW = (1.0 - alpha) * p.second * QNodes_[p.first].toQueryEdges_[j].w_;
		if(newW < PPR_EPS_Q)
			continue;
            activeInk.push(make_pair(newQid, newW));
        }
    }
    
    ofstream out(cacheP.c_str(), ios::out);
    for(map<int, double>::iterator i = ret.begin(); i != ret.end(); ++i){
        out << i->first << '\t' << i->second << endl;
    }
    out.close();
	
}
struct compare_node
{
    bool operator()(const pair<int, double> & n1, const pair<int, double> & n2) const
    {
        return n1.second > n2.second;
    }
};
void EQFG::PPR_E(int eid, map<int, double> &ret, double alpha)
{
	cerr << "in PPR_E" << endl;
    BoundHeap activeInk(999999);

    activeInk.push(make_pair(eid, 1.0));
    
    while(activeInk.size() > 0){
        pair<int, double> p = activeInk.pop();
        if(ret.find(p.first) == ret.end()){
            ret[p.first] = 0.0;
        }
        ret[p.first] += p.second * alpha;
        if(p.second < PPR_EPS_E){
            continue;
        }
        int degree = ENodes_[p.first].toEntityEdges_.size();
        for(int j = 0; j < degree; ++j){
            int newEid = ENodes_[p.first].toEntityEdges_[j].eid_;
        	activeInk.push(make_pair(newEid, (1.0 - alpha) * p.second * ENodes_[p.first].toEntityEdges_[j].w_));    
	//activeInk.push(make_pair(newEid, (1.0 - alpha) * p.second / degree));
        }
    }
	cerr << "out PPRE" << endl;
}

vector<int> EQFG::expandEntity(const vector<int> &eids, int k)
{
    vector<int> ret;
    map<int, double> totalWeight;
    
    for(int i = 0; i < eids.size(); ++i){
        map<int, double> tempWeight;
        PPR_E(eids[i], tempWeight);
        for(map<int, double>::iterator iter = tempWeight.begin(); iter != tempWeight.end(); ++iter){
            if(totalWeight.find(iter->first) == totalWeight.end()){
                totalWeight[iter->first] = 0.0;
            }
            totalWeight[iter->first] += iter->second;
        }
    }
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < k && i < results.size(); ++i){
        ret.push_back(results[i].first);
    }
    return ret;
}

map<int, double> EQFG::expandEntity_sp(const vector<string> &entities, int k)
{
	map<int, double> totalWeight;
	vector<int> Yeids;
	for(int i = 0; i < entities.size(); ++i){
		if(hin_.key2id_.find(entities[i]) != hin_.key2id_.end()){
			Yeids.push_back(hin_.key2id_[entities[i]]);
		}
	}
	for(int i = 0; i < Yeids.size(); ++i){
		//static void PCRW(int src, int dst, Meta_Paths & m_path, int p_id, map<int, double> & id_sim, vector< vector<pair<int, double> > > & pathInstances);
		map<int, double> metapathWeight;
		for(int j = 0; j < hin_.nodes_[Yeids[i]].types_id_.size(); ++j){
			int tid = hin_.nodes_[Yeids[i]].types_id_[j];
			if(MetaPathsLists_.find(tid) == MetaPathsLists_.end())
				continue;
			for(map<int, double>::iterator iter = MetaPathsLists_[tid].begin(); iter != MetaPathsLists_[tid].end(); ++iter){
				if(metapathWeight.find(iter->first) == metapathWeight.end())
					metapathWeight[iter->first] = 0.0;
				metapathWeight[iter->first] += iter->second;
			}
		}
		Meta_Paths tempmetapath(hin_);
		for(map<int, double>::iterator iter = metapathWeight.begin(); iter != metapathWeight.end(); ++iter){
			tempmetapath.weights_.push_back(iter->second);
			tempmetapath.linkTypes_.push_back(metapaths_[iter->first]);
		}
		map<int, double> id_sim;
		SimCalculator::calSim_PCRW(Yeids[i], -1, tempmetapath, id_sim);
		for(map<int, double>::iterator iter = id_sim.begin(); iter != id_sim.end(); ++iter){
			string tempS = "http://en.wikipedia.org/wiki/";
			tempS += hin_.nodes_[iter->first].key_.substr(1);
			tempS = tempS.substr(0, tempS.size() - 1);
			if(entity2id_.find(tempS) == entity2id_.end()) continue;
			int eid = entity2id_[tempS];
			if(totalWeight.find(eid) == totalWeight.end())
				totalWeight[eid] = 0.0;
			totalWeight[eid] += iter->second;
		}
	}
	vector<pair<int, double>> results;
	for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
		results.push_back(*iter);
	}
	sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
	map<int, double> ret;
	
	// we only keep the top k (10) most important entities
	for(int i = 0; i < k && i < results.size(); ++i){
		ret[results[i].first] = results[i].second;
	}
	return ret;
}

// recommend for a labeled query with EQFG
vector<RecNode> EQFG::rec_EQFG(int qid, vector<int> eids)
{
    vector<RecNode> ret;
	cerr << "start to expand"<< endl;
    vector<int> expanded_eids = expandEntity(eids);
	cerr << "end of expand" << endl;
    map<int, double> totalWeight;
    for(int i = 0; i < expanded_eids.size(); ++i){
        map<int, double> tempWeight;
        PPR_Q(expanded_eids[i], tempWeight);
        for(map<int, double>::iterator iter = tempWeight.begin(); iter != tempWeight.end(); ++iter){
            if(totalWeight.find(iter->first) == totalWeight.end()){
                totalWeight[iter->first] = 0.0;
            }
            totalWeight[iter->first] += iter->second;
        }
    }
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

// recommend for a labeled query with EQFG but centerpice
vector<RecNode> EQFG::rec_EQFG_t(int qid, vector<int> eids)
{
    vector<RecNode> ret;
    vector<int> expanded_eids = expandEntity(eids);
    map<int, double> totalWeight;
    for(int i = 0; i < expanded_eids.size(); ++i){
        map<int, double> tempWeight;
        PPR_Q(expanded_eids[i], tempWeight);
        for(map<int, double>::iterator iter = tempWeight.begin(); iter != tempWeight.end(); ++iter){
            if(totalWeight.find(iter->first) == totalWeight.end()){
                totalWeight[iter->first] = 1.0;
            }
            totalWeight[iter->first] *= iter->second;
        }
    }
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

vector<RecNode> EQFG::rec_QTFGP(int qid, vector<string> entities, map<string, double> & rec)
{
    vector<RecNode> ret;
    clock_t t1 = clock();
    map<int, double> expanded_eids = expandEntity_sp(entities);
    clock_t t2 = clock();
    timeForExp += (t2 - t1);
	double tempsum = 0;
	for(map<int, double>::iterator i = expanded_eids.begin(); i != expanded_eids.end(); ++i){
		tempsum += i->second;
	}
	for(map<int, double>::iterator i = expanded_eids.begin(); i != expanded_eids.end(); ++i){
		i->second /= tempsum;
	}
    map<int, double> totalWeight;
	double sumOfQTFG = 0;
	// initialize by the QFG results
	for(map<string, double>::iterator i = rec.begin(); i != rec.end(); ++i){
		if(query2id_.find(i->first) != query2id_.end()){
			int id = query2id_[i->first];
			if(totalWeight.find(id) == totalWeight.end())
				totalWeight[id] = 0.0;
			sumOfQTFG += i->second;
			totalWeight[id] += i->second;
		}
	}
	for(map<int, double>::iterator i = totalWeight.begin(); i != totalWeight.end(); ++i){
		i->second /= sumOfQTFG;
	}
	
	// revised by the YAGO result
	double totalTime = 0;
	int times = 0;
	for(map<int, double>::iterator iter = expanded_eids.begin(); iter != expanded_eids.end(); ++iter){
        map<int, double> tempWeight;
        clock_t t3 = clock();
        PPR_Q(iter->first, tempWeight);
        clock_t t4 = clock();
        timeForPPR += (t4 - t3);
        for(map<int, double>::iterator iter2 = tempWeight.begin(); iter2 != tempWeight.end(); ++iter2){
            if(totalWeight.find(iter2->first) == totalWeight.end()){
                totalWeight[iter2->first] = 0.0;
            }
            totalWeight[iter2->first] += iter2->second * iter->second * ratio_;
        }
    }
	
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

vector<RecNode> EQFG::rec_QFGP(int qid, vector<string> entities)
{
	vector<RecNode> QFG_rec;
	if(qid != -1){
		string thequery = queries_[qid];
		QFG_rec = rec_QFG(thequery);
	}
	
    vector<RecNode> ret;
	double sum = 0.0;
    map<int, double> expanded_eids = expandEntity_sp(entities);
	for(map<int, double>::iterator i = expanded_eids.begin(); i != expanded_eids.end(); ++i){
		sum += i->second;
	}
	if(sum != 0)
	for(map<int, double>::iterator i = expanded_eids.begin(); i != expanded_eids.end(); ++i){
		i->second /= sum;
	}
    map<int, double> totalWeight;
	// initialize by the QFG results
	for(int i = 0; i < QFG_rec.size(); ++i){
		totalWeight[QFG_rec[i].qid_] = QFG_rec[i].w_;
	}
	// revised by the YAGO result
	for(map<int, double>::iterator iter = expanded_eids.begin(); iter != expanded_eids.end(); ++iter){
        map<int, double> tempWeight;
        PPR_Q(iter->first, tempWeight);
        for(map<int, double>::iterator iter2 = tempWeight.begin(); iter2 != tempWeight.end(); ++iter2){
            if(totalWeight.find(iter2->first) == totalWeight.end()){
                totalWeight[iter2->first] = 0.0;
            }
            totalWeight[iter2->first] += iter2->second * iter->second * ratio_;
        }
    }
	
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

vector<RecNode> EQFG::rec_DQG(int qid, vector<string> entities)
{
	vector<int> eids;
	for(int i = 0; i < entities.size(); ++i){
		string fullEntity = "http://en.wikipedia.org/wiki/";
		fullEntity += entities[i].substr(1, entities[i].size() - 2);
		//cerr << fullEntity << endl;
		if(entity2id_.find(fullEntity) != entity2id_.end())
			eids.push_back(entity2id_[fullEntity]);
	}
	vector<RecNode> ret;
	
	map<int, double> entity2weight;
	for(int i = 0; i < eids.size(); ++i){
		if (entity2docPro_.find(eids[i]) == entity2docPro_.end()){
			//cerr << entities_[eids[i]] << " not found!" << endl;
			continue;
		}
		for(map<int, double>::iterator j = entity2docPro_[eids[i]].begin(); j != entity2docPro_[eids[i]].end(); ++j){
			int docid = j->first;
			if (doc2entityPro_.find(docid) == doc2entityPro_.end())
				continue;
			for(map<int, double>::iterator k = doc2entityPro_[docid].begin(); k != doc2entityPro_[docid].end(); ++k){
				entity2weight[k->first] = j->second * k->second / eids.size();
			}
		}
	}
	
	cerr << "The DQG related entities for ";
	for(map<int, double>::iterator i = entity2weight.begin(); i != entity2weight.end(); ++i){
		cerr << '\t' << entities_[i->first];
	}
	cerr << endl;
	
	map<int, double> totalWeight;
	for(map<int, double>::iterator iter = entity2weight.begin(); iter != entity2weight.end(); ++iter){
        map<int, double> tempWeight;
        PPR_Q(iter->first, tempWeight);
        for(map<int, double>::iterator iter2 = tempWeight.begin(); iter2 != tempWeight.end(); ++iter2){
            if(totalWeight.find(iter2->first) == totalWeight.end()){
                totalWeight[iter2->first] = 0.0;
            }
            totalWeight[iter2->first] += iter2->second * iter->second;
        }
    }
	vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

vector<RecNode> EQFG::rec_EQFGP(int qid, vector<string> entities)
{
	vector<int> eids;
	for(int i = 0; i < entities.size(); ++i){
		string fullEntity = "http://en.wikipedia.org/wiki/";
		fullEntity += entities[i].substr(1, entities[i].size() - 2);
		//cerr << fullEntity << endl;
		if(entity2id_.find(fullEntity) != entity2id_.end())
			eids.push_back(entity2id_[fullEntity]);
	}
	vector<RecNode> EQFG_rec = rec_EQFG(qid, eids);
    vector<RecNode> ret;
	
    map<int, double> expanded_eids = expandEntity_sp(entities);

    map<int, double> totalWeight;
	// initialize by the EQFG results
	for(int i = 0; i < EQFG_rec.size(); ++i){
		totalWeight[EQFG_rec[i].qid_] = EQFG_rec[i].w_;
	}
	// revised by the YAGO result
	for(map<int, double>::iterator iter = expanded_eids.begin(); iter != expanded_eids.end(); ++iter){
        map<int, double> tempWeight;
        PPR_Q(iter->first, tempWeight);
        for(map<int, double>::iterator iter2 = tempWeight.begin(); iter2 != tempWeight.end(); ++iter2){
            if(totalWeight.find(iter2->first) == totalWeight.end()){
                totalWeight[iter2->first] = 0.0;
            }
            totalWeight[iter2->first] += iter2->second * iter->second * ratio_;
        }
    }
	
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

vector<RecNode> EQFG::rec_QFGPEQFG(int qid, vector<string> entities)
{
	cerr << 1<<endl;
	vector<RecNode> QFG_rec;
	if(qid != -1){
		string thequery = queries_[qid];
		QFG_rec = rec_QFG(thequery);
	}
	cerr << 2<<endl;
	vector<int> eids;
	for(int i = 0; i < entities.size(); ++i){
		string fullEntity = "http://en.wikipedia.org/wiki/";
		fullEntity += entities[i].substr(1, entities[i].size() - 2);
		//cerr << fullEntity << endl;
		if(entity2id_.find(fullEntity) != entity2id_.end())
			eids.push_back(entity2id_[fullEntity]);
	}
	cerr << 3<<endl;
	vector<RecNode> EQFG_rec = rec_EQFG(qid, eids);
	cerr << 4<<endl;
    vector<RecNode> ret;
	

    map<int, double> totalWeight;
	// initialize by the EQFG results

	for(int i = 0; i < QFG_rec.size(); ++i){
		totalWeight[QFG_rec[i].qid_] = QFG_rec[i].w_;
	}
	for(int i = 0; i < EQFG_rec.size(); ++i){
		if(totalWeight.find(EQFG_rec[i].qid_) == totalWeight.end())
			totalWeight[EQFG_rec[i].qid_] = 0.0;
		totalWeight[EQFG_rec[i].qid_] += EQFG_rec[i].w_ * ratio_;
	}
	cerr << 5<<endl;
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
	cerr << 6<<endl;
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
	cerr << 7<<endl;
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
	cerr << 8<<endl;
    return ret;
}

vector<RecNode> EQFG::rec_P(int qid, vector<string> entities)
{
    vector<RecNode> ret;
    map<int, double> expanded_eids = expandEntity_sp(entities);
	// for(map<int, double>::iterator i = expanded_eids.begin(); i != expanded_eids.end(); ++i){
	// 	cerr << i->first << '\t' << i->second << '\t';
	// }
	// cerr << endl;
    map<int, double> totalWeight;
	for(map<int, double>::iterator iter = expanded_eids.begin(); iter != expanded_eids.end(); ++iter){
        map<int, double> tempWeight;
        PPR_Q(iter->first, tempWeight);
        for(map<int, double>::iterator iter2 = tempWeight.begin(); iter2 != tempWeight.end(); ++iter2){
            if(totalWeight.find(iter2->first) == totalWeight.end()){
                totalWeight[iter2->first] = 0.0;
            }
            totalWeight[iter2->first] += iter2->second * iter->second;
        }
    }
    vector<pair<int, double> > results;
    for(map<int, double>::iterator iter = totalWeight.begin(); iter != totalWeight.end(); ++iter){
        results.push_back(*iter);
    }
    sort(results.begin(), results.end(), pairIntDoubleLargerCMP);
    for(int i = 0; i < results.size(); ++i){
        ret.push_back(RecNode(results[i].first, results[i].second));
    }
    return ret;
}

void EQFG::rec_P_fromFile(string queryEntityPath, string cachePath, string metapathPath)
{
    cachePath_EQFG_ = cachePath;
    loadMetaPaths(metapathPath);
    
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
    	while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			//cerr << tempS << endl;
			entities.push_back(tempS);
        }
        vector<RecNode> ret = rec_P(qid, entities);
        showRecNodes(strs[1], ret);
    }
}

void EQFG::rec_DQG_fromFile(string queryEntityPath, string cachePath)
{
	cachePath_EQFG_ = cachePath;
	ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
	while(getline(queryEntityIn, line)){
		cerr << line << endl;
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			//cerr << tempS << endl;
			entities.push_back(tempS);
        }
		vector<RecNode> ret = rec_DQG(qid, entities);
		showRecNodes(strs[1], ret);
	}	
}

void EQFG::rec_QTFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, string QTFGPath, double ratio)
{
    timeForPPR = 0.0;
    timeForExp = 0.0;
    clock_t t1, t2, t3, t4;
    t1 = clock();
	cachePath_EQFG_ = cachePath;
    loadMetaPaths(metapathPath);
	ratio_ = ratio;
    t2 = clock();
    cerr << (double)(t2-t1)/CLOCKS_PER_SEC << " s to load metapath" << endl;
	
    fstream QTFGin(QTFGPath.c_str(), ios::in);
	ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
	string QTFGline;
	map<string, map<string, double> > QTFGresult;
	while(getline(QTFGin, QTFGline)){
		while(QTFGline[QTFGline.size() - 1] == '\r')
        	QTFGline = QTFGline.substr(0, QTFGline.size() - 1);
		vector<string> strs = split(QTFGline, "\t");
		string query = strs[0];
		map<string, double> tempRec;
		for(int i = 1; i < strs.size(); i += 2){
			tempRec[strs[i]] = atof(strs[i + 1].c_str());
		}
		QTFGresult[query] = tempRec;
	}
	QTFGin.close();
	t3 = clock();
    cerr << (double)(t3-t2)/CLOCKS_PER_SEC << " s to load TQGraph results." << endl;

	while(getline(queryEntityIn, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			//cerr << tempS << endl;
			entities.push_back(tempS);
        }
		vector<RecNode> ret;
		if(QTFGresult.find(strs[1]) != QTFGresult.end())
			ret = rec_QTFGP(qid, entities, QTFGresult[strs[1]]);
        else{
			map<string, double> empty;
			ret = rec_QTFGP(qid, entities, empty);
		}
		showRecNodes(strs[1], ret);
	}
    t4 = clock();
	cerr << (double)(t4-t3)/CLOCKS_PER_SEC << " s to recommend." << endl;
    cerr << timeForPPR/CLOCKS_PER_SEC << " s to PPR." << endl;
    cerr << timeForExp/CLOCKS_PER_SEC << " s to expand entity." << endl;
}


void EQFG::rec_QTFG_DQGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, string QTFGPath, string DQGPath, double ratio)
{
    timeForPPR = 0.0;
    timeForExp = 0.0;
    clock_t t1, t2, t3, t4;
    t1 = clock();
	cachePath_EQFG_ = cachePath;
    loadMetaPaths(metapathPath);
	ratio_ = ratio;
    t2 = clock();
    cerr << (double)(t2-t1)/CLOCKS_PER_SEC << " s to load metapath" << endl;
	
    ifstream QTFGin(QTFGPath.c_str(), ios::in);
	ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
	string QTFGline;
	map<string, map<string, double> > QTFGresult;
	while(getline(QTFGin, QTFGline)){
		while(QTFGline[QTFGline.size() - 1] == '\r')
        	QTFGline = QTFGline.substr(0, QTFGline.size() - 1);
		vector<string> strs = split(QTFGline, "\t");
		string query = strs[0];
		map<string, double> tempRec;
		for(int i = 1; i < strs.size(); i += 2){
			tempRec[strs[i]] = atof(strs[i + 1].c_str());
		}
		QTFGresult[query] = tempRec;
	}
	QTFGin.close();
	t3 = clock();
    cerr << (double)(t3-t2)/CLOCKS_PER_SEC << " s to load TQGraph results." << endl;

	ifstream DQGin(DQGPath.c_str(), ios::in);
	string DQGline;
	while(getline(DQGin, DQGline)){
		while(DQGline[DQGline.size() - 1] == '\r')
        	DQGline = DQGline.substr(0, DQGline.size() - 1);
		vector<string> strs = split(DQGline, "\t");
		string query = strs[0];
		map<string, double> tempRec;
		for(int i = 1; i < strs.size(); i += 2){
			tempRec[strs[i]] = atof(strs[i + 1].c_str());
		}
		if(QTFGresult.find(query) == QTFGresult.end()){
			QTFGresult[query] = map<string, double>();
		}
		for(map<string, double>::iterator ii = tempRec.begin(); ii != tempRec.end(); ++ii){
			if(QTFGresult[query].find(ii->first) == QTFGresult[query].end()){
				QTFGresult[query][ii->first] = 0.0;
			}
			QTFGresult[query][ii->first] += ii->second;
		}
	}
	DQGin.close();
	t3 = clock();
    cerr << (double)(t3-t2)/CLOCKS_PER_SEC << " s to load DQG results." << endl;
	
	
	while(getline(queryEntityIn, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			//cerr << tempS << endl;
			entities.push_back(tempS);
        }
		vector<RecNode> ret;
		if(QTFGresult.find(strs[1]) != QTFGresult.end())
			ret = rec_QTFGP(qid, entities, QTFGresult[strs[1]]);
        else{
			map<string, double> empty;
			ret = rec_QTFGP(qid, entities, empty);
		}
		showRecNodes(strs[1], ret);
	}
    t4 = clock();
	cerr << (double)(t4-t3)/CLOCKS_PER_SEC << " s to recommend." << endl;
    cerr << timeForPPR/CLOCKS_PER_SEC << " s to PPR." << endl;
    cerr << timeForExp/CLOCKS_PER_SEC << " s to expand entity." << endl;
}

void EQFG::combineTwoResult(string testpath, string resultpath1, string resultpath2, double alpha)
{
	ifstream in1(resultpath1.c_str(), ios::in);
    string line;
	map<string, map<string, double> > result;
	while(getline(in1, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
		vector<string> strs = split(line, "\t");
		string query = strs[0];
		map<string, double> tempRec;
		for(int i = 1; i < strs.size(); i += 2){
			tempRec[strs[i]] = alpha * atof(strs[i + 1].c_str());
		}
		result[query] = tempRec;
	}
	in1.close();
	
	ifstream in2(resultpath2.c_str(), ios::in);
	while(getline(in2, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
		vector<string> strs = split(line, "\t");
		string query = strs[0];
		if (result.find(query) == result.end()){
			result[query] = map<string, double>();
		}
		for(int i = 1; i < strs.size(); i += 2){
			if (result[query].find(strs[i]) == result[query].end()){
				result[query][strs[i]] = 0.0;
			}
			result[query][strs[i]] += (1 - alpha) * atof(strs[i + 1].c_str());
		}
	}
	in2.close();
	
	ifstream queryEntityIn(testpath.c_str(), ios::in);
	while(getline(queryEntityIn, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        string q = strs[1];
		
		vector<pair<string, double> > results;
		
		for(map<string, double>::iterator iter = result[q].begin(); iter != result[q].end(); ++iter){
			results.push_back(*iter);
		}
		sort(results.begin(), results.end(), pairStringDoubleLargerCMP);
		
		cout << q;
		for(int i = 0; i < results.size() && i < 5 ; ++i){
			cout << '\t' << results[i].first << '\t' << results[i].second;
		}
		cout << endl;
	}
}

void EQFG::joinTwoResult(string testpath, string resultpath1, string resultpath2)
{
    ifstream in1(resultpath1.c_str(), ios::in);
    string line;
    map<string, vector<string> > result;
    while(getline(in1, line)){
        while(line[line.size() - 1] == '\r')
            line = line.substr(0, line.size() - 1);
        vector<string> strs = split(line, "\t");
        string query = strs[0];
        vector<string> tempRec;
        for(int i = 1; i < strs.size(); i += 2){
            tempRec.push_back(strs[i]);// [strs[i]] = alpha * atof(strs[i + 1].c_str());
        }
        result[query] = tempRec;
    }
    in1.close();
    
    ifstream in2(resultpath2.c_str(), ios::in);
    while(getline(in2, line)){
        while(line[line.size() - 1] == '\r')
            line = line.substr(0, line.size() - 1);
        vector<string> strs = split(line, "\t");
        string query = strs[0];
        if (result.find(query) == result.end()){
            result[query] = vector<string>();
        }
        for(int i = 1; result[query].size() < 5 && i < strs.size(); i += 2){
            bool found = false;
            string newRec = strs[i];
            for (int j = 0; j < result[query].size(); ++j){
                if(result[query][j] == newRec){
                    found = true;
                    break;
                }
            }
            if (!found){
                result[query].push_back(newRec);
            }
        }
    }
    in2.close();
    
    ifstream queryEntityIn(testpath.c_str(), ios::in);
    while(getline(queryEntityIn, line)){
        while(line[line.size() - 1] == '\r')
            line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        string q = strs[1];
        
        vector<pair<string, double> > results;
        
        for(vector<string>::iterator iter = result[q].begin(); iter != result[q].end(); ++iter){
            results.push_back(make_pair(*iter, 1.0));
        }
        cout << q;
        for(int i = 0; i < results.size() && i < 5 ; ++i){
            cout << '\t' << results[i].first << '\t' << results[i].second;
        }
        cout << endl;
    }
}

void EQFG::combineThreeResult(string testpath, string resultpath1, string resultpath2, string resultpath3, double alpha, double beta)
{
	ifstream in1(resultpath1.c_str(), ios::in);
    string line;
	map<string, map<string, double> > result;
	while(getline(in1, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
		vector<string> strs = split(line, "\t");
		string query = strs[0];
		map<string, double> tempRec;
		for(int i = 1; i < strs.size(); i += 2){
			tempRec[strs[i]] = alpha * atof(strs[i + 1].c_str());
		}
		result[query] = tempRec;
	}
	in1.close();
	
	ifstream in2(resultpath2.c_str(), ios::in);
	while(getline(in2, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
		vector<string> strs = split(line, "\t");
		string query = strs[0];
		if (result.find(query) == result.end()){
			result[query] = map<string, double>();
		}
		for(int i = 1; i < strs.size(); i += 2){
			if (result[query].find(strs[i]) == result[query].end()){
				result[query][strs[i]] = 0.0;
			}
			result[query][strs[i]] += beta * atof(strs[i + 1].c_str());
		}
	}
	in2.close();
	
	ifstream in3(resultpath3.c_str(), ios::in);
	while(getline(in3, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
		vector<string> strs = split(line, "\t");
		string query = strs[0];
		if (result.find(query) == result.end()){
			result[query] = map<string, double>();
		}
		for(int i = 1; i < strs.size(); i += 2){
			if (result[query].find(strs[i]) == result[query].end()){
				result[query][strs[i]] = 0.0;
			}
			result[query][strs[i]] += (1.0 - alpha - beta) * atof(strs[i + 1].c_str());
		}
	}
	in3.close();
	
	ifstream queryEntityIn(testpath.c_str(), ios::in);
	while(getline(queryEntityIn, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        string q = strs[1];
		
		vector<pair<string, double> > results;
		
		for(map<string, double>::iterator iter = result[q].begin(); iter != result[q].end(); ++iter){
			results.push_back(*iter);
		}
		sort(results.begin(), results.end(), pairStringDoubleLargerCMP);
		
		cout << q;
		for(int i = 0; i < results.size() && i < 5 ; ++i){
			cout << '\t' << results[i].first << '\t' << results[i].second;
		}
		cout << endl;
	}
}


void EQFG::rec_QFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, double ratio)
{
	ratio_ = ratio;
    cachePath_EQFG_ = cachePath;
    loadMetaPaths(metapathPath);
    
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
    	while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			cerr << tempS << endl;
			entities.push_back(tempS);
        }
        vector<RecNode> ret = rec_QFGP(qid, entities);
        showRecNodes(strs[1], ret);
    }
}

void EQFG::rec_EQFGP_fromFile(string queryEntityPath, string cachePath, string metapathPath, double ratio)
{
	ratio_ = ratio;
    cachePath_EQFG_ = cachePath;
    loadMetaPaths(metapathPath);
    
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
    	while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			cerr << tempS << endl;
			entities.push_back(tempS);
        }
        vector<RecNode> ret = rec_EQFGP(qid, entities);
        showRecNodes(strs[1], ret);
    }
}

void EQFG::rec_QFGPEQFG_fromFile(string queryEntityPath, string cachePath, double ratio)
{
	ratio_ = ratio;
    cachePath_EQFG_ = cachePath;
    
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
    	while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			//cerr << tempS << endl;
			entities.push_back(tempS);
        }
        vector<RecNode> ret = rec_QFGPEQFG(qid, entities);
        showRecNodes(strs[1], ret);
    }
}

void EQFG::rec_EQFG_fromFile(string queryEntityPath, string cachePath)
{
    cachePath_EQFG_ = cachePath;
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
	while(line[line.size() - 1] == '\r')
                line = line.substr(0, line.size() - 1);

        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<int> eids;
        for(int i = 3; i < strs.size(); i += 2){
            if(entity2id_.find(strs[i]) != entity2id_.end())
                eids.push_back(entity2id_[strs[i]]);
        }
        vector<RecNode> ret = rec_EQFG(qid, eids);
        //vector<RecNode> ret = rec_EQFG_t(qid, eids);
        showRecNodes(strs[1], ret);
    }
}



void EQFG::rec_EQFG_t_fromFile(string queryEntityPath, string cachePath)
{
    cachePath_EQFG_ = cachePath;
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    string line;
    while(getline(queryEntityIn, line)){
	while(line[line.size() - 1] == '\r')
                line = line.substr(0, line.size() - 1); 
       	vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<int> eids;
        for(int i = 3; i < strs.size(); i += 2){
            if(entity2id_.find(strs[i]) != entity2id_.end())
                eids.push_back(entity2id_[strs[i]]);
        }
        //vector<RecNode> ret = rec_EQFG(qid, eids);
        vector<RecNode> ret = rec_EQFG_t(qid, eids);
        showRecNodes(strs[1], ret);
    }
}
void EQFG::rec_QFG_fromFile(string queryPath)
{
    ifstream queryIn(queryPath.c_str(), ios::in);
    string line;
    while(getline(queryIn, line)){
     	while(line[line.size() - 1] == '\r')
                line = line.substr(0, line.size() - 1);
	vector<RecNode> ret = rec_QFG(line);
        showRecNodes(line, ret);
    }
}

string EQFG::metapath2string(const vector<int> &metapath)
{
    string rtn;
    stringstream ss;
    for(int i = 0; i < metapath.size(); ++i){
        ss << metapath[i] << '\t';
    }
    rtn = ss.str();
    return rtn;
}
vector<int> EQFG::string2metapath(string s){
    stringstream ss(s);
    int id;
    vector<int> rtn;
    while(ss >> id)
        rtn.push_back(id);
    return rtn;
}

void EQFG::loadMetaPaths(string path)
{
    ifstream in(path.c_str(), ios::in);
    int numofmetapath;
    in >> numofmetapath;
    for(int i = 0; i < numofmetapath; ++i){
        int size;
        in >> size;
        vector<int> path;
        for(int j = 0; j < size; ++j){
            int tid;
            in >> tid;
            path.push_back(tid);
        }
        metapaths_.push_back(path);
    }
    int sizeofmetapathList;
    in >> sizeofmetapathList;
    for(int i = 0; i < sizeofmetapathList; ++i){
        int eid, size;
        in >> eid >> size;
        map<int, double> tempmap;
		double sum = 0.0;
        for(int j = 0; j < size; ++j){
            int pid;
            double w;
            in >> pid >> w;
            tempmap[pid] = w;
			sum += w;
        }
		for(map<int, double>::iterator j = tempmap.begin(); j != tempmap.end(); ++j){
			j->second /= sum;
		}
        MetaPathsLists_[eid] = tempmap;
    }
    in.close();
}

void EQFG::saveMetaPaths(string path)
{
    ofstream out(path, ios::out);
    /*
    
    for(map<int, Meta_Paths>::iterator i = MetaPathsLists_.begin(); i != MetaPathsLists_.end(); ++i){
        int len = i->second.weights_.size();
        out << i->first << len << '\t' << i->second.nodeTypes_.size() << endl;
        for(int j = 0; j < len; ++j){
            out << i->second.weights_[j] << '\t';
        }
        out << endl;
        for(int j = 0; j < len; ++j){
            for(int k = 0; k < i->second.nodeTypes_.size(); ++k){
                out << i->second.nodeTypes_[j][k] << '\t';
            }
            out << endl;
        }
        for(int j = 0; j < len; ++j){
            for(int k = 0; k < i->second.nodeTypes_.size() - 1; ++k){
                out << i->second.linkTypes_[j][k] << '\t';
            }
            out << endl;
        }
    }
     */
    out << metapaths_.size() << endl;
    for(int i = 0; i < metapaths_.size(); ++i){
        out << metapaths_[i].size();
        for(int j = 0; j < metapaths_[i].size(); ++j){
            out << '\t' << metapaths_[i][j];
        }
        out << endl;
    }
    out << MetaPathsLists_.size() << endl;
    for(map<int, map<int, double> >::iterator i = MetaPathsLists_.begin(); i != MetaPathsLists_.end(); ++i){
        out << i->first << '\t' << i->second.size();
        for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
            out << '\t' << j->first << '\t' << j->second;
        }
        out << endl;
    }
    out.close();
}

void EQFG::saveMetaPathsList(string path, map<int, map<int, double> > & PathList)
{
    ofstream out(path, ios::out);
    
    out << metapaths_.size() << endl;
    for(int i = 0; i < metapaths_.size(); ++i){
        out << metapaths_[i].size();
        for(int j = 0; j < metapaths_[i].size(); ++j){
            out << '\t' << metapaths_[i][j];
        }
        out << endl;
    }
    out << PathList.size() << endl;
    for(map<int, map<int, double> >::iterator i = PathList.begin(); i != PathList.end(); ++i){
        out << i->first << '\t' << i->second.size();
        for(map<int, double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
            out << '\t' << j->first << '\t' << j->second;
        }
        out << endl;
    }
    out.close();
}

void EQFG::getMetaPath(const vector<int> &pathNodes, vector<int> &nodeTypes, vector<int> &edgeTypes)
{
    nodeTypes.clear();
    edgeTypes.clear();
    for(int i = 0; i < pathNodes.size(); ++i){
        if(hin_.nodes_[pathNodes[i]].types_id_.size() > 0)
            nodeTypes.push_back(hin_.nodes_[pathNodes[i]].types_id_[0]);
        else
            nodeTypes.push_back(-1);
    }
    for(int i = 0; i < pathNodes.size() - 1; ++i){
        int sid = pathNodes[i];
        int eid = pathNodes[i + 1];
        int edgetype = -1;
        for(int j = 0; j < hin_.edges_src_[sid].size(); ++j){
            if(hin_.edges_src_[sid][j].dst_ == eid){
                edgetype = hin_.edges_src_[sid][j].edge_type_;
                break;
            }
        }
        edgeTypes.push_back(edgetype);
    }
}

bool EQFG::entityAllShortestPaths(int id1, int id2, int dist, vector< vector< vector<int> > > &pathNodes)
{
    pathNodes.clear();

    vector< vector<int> > tempPaths;
    tempPaths.push_back(vector<int>());
    tempPaths[0].push_back(id1);

    if (id1 == id2){
        pathNodes.push_back(tempPaths);
    }else{
        pathNodes.push_back(vector< vector<int>>());
    }

    for (int i = 1; i <= dist; ++i){
        vector< vector<int> > newPaths;
        for(int j = 0; j < tempPaths.size(); ++j){
            int endid = tempPaths[j][tempPaths[j].size() - 1];
            for(int i = 0; i < hin_.edges_dst_[endid].size(); ++i){
                int nextid = hin_.edges_dst_[endid][i].src_;
                newPaths.push_back(tempPaths[j]);
                newPaths[newPaths.size() - 1].push_back(nextid);
            }
        }
        pathNodes.push_back( vector < vector <int> >());
        for (int j = 0; j < newPaths.size(); ++j){
            if (newPaths[j][newPaths[j].size() - 1] == id2){
                pathNodes[pathNodes.size() - 1].push_back(newPaths[j]);
            }
        }
        tempPaths = newPaths;
    }
    return true;
}



bool EQFG::entityShortestPath(int id1, int id2, int dist, vector<int> &pathNodes)
{
    if(id1 == id2){
        pathNodes.push_back(id1);
        return true;
    }
    map<int, int> id_dis, id_pre;
    id_dis[id1] = 0;
    
    queue<int> q;
    q.push(id1);
    bool found = false;
    while(!q.empty() && !found){
        int node = q.front();
        q.pop();
        int nowdis = id_dis[node];
        if(nowdis >= dist)
            continue;
        for(int i = 0; i < hin_.edges_src_[node].size(); ++i){
            int nextid = hin_.edges_src_[node][i].dst_;
            if(id_dis.find(nextid) != id_dis.end())
                continue;
            id_pre[nextid] = node;
            id_dis[nextid] = nowdis + 1;
            if(nextid == id2){
                found = true;
                break;
            }
            q.push(nextid);
        }
        if(found)
            break;
        for(int i = 0; i < hin_.edges_dst_[node].size(); ++i){
            int nextid = hin_.edges_dst_[node][i].src_;
            if(id_dis.find(nextid) != id_dis.end())
                continue;
            id_pre[nextid] = node;
            id_dis[nextid] = nowdis + 1;
            if(nextid == id2){
                found = true;
                break;
            }
            q.push(nextid);
        }
    }
    if(!found)
        return false;
    vector<int> tempv;
    int t = id2;
    while(true){
        tempv.push_back(t);
        if(t == id1)
            break;
        t = id_pre[t];
    }
    for(int i = tempv.size() - 1; i >= 0; --i){
        pathNodes.push_back(tempv[i]);
    }
    return true;
}

void EQFG::insertMetaPath(int YAGOeid, const vector<int> &nodeTypes, const vector<int> &edgeTypes, double w)
{
    string s_path = metapath2string(edgeTypes);
    if(metapath_2id_.find(s_path) == metapath_2id_.end()){
        metapath_2id_[s_path] = metapaths_.size();
        metapaths_.push_back(edgeTypes);
    }
    int metapathid = metapath_2id_[s_path];
    for(int i = 0 ; i < hin_.nodes_[YAGOeid].types_id_.size(); ++i){
        int typeID = hin_.nodes_[YAGOeid].types_id_[i];
        if(MetaPathsLists_.find(typeID) == MetaPathsLists_.end())
            MetaPathsLists_[typeID] = map<int, double>();
        if(MetaPathsLists_[typeID].find(metapathid) == MetaPathsLists_[typeID].end())
            MetaPathsLists_[typeID][metapathid] = 0.0;
        MetaPathsLists_[typeID][metapathid] += w;
    }
    /*
    for(int i = 0 ; i < hin_.nodes_[YAGOeid].types_id_.size(); ++i){
        int typeID = hin_.nodes_[YAGOeid].types_id_[i];
        map<int, Meta_Paths>::iterator iter;
        if((iter = MetaPathsLists_.find(typeID)) == MetaPathsLists_.end()){
            //Meta_Paths tempPaths(hin_);
            iter = MetaPathsLists_.insert(make_pair(typeID, Meta_Paths(hin_))).first;
        }
        iter->second.nodeTypes_.push_back(nodeTypes);
        iter->second.linkTypes_.push_back(edgeTypes);
        iter->second.weights_.push_back(w);
        //MetaPathsLists_[typeID].nodeTypes_.push_back(nodeTypes);
        //MetaPathsLists_[typeID].weights_.push_back(w);
        //MetaPathsLists_[typeID].linkTypes_.push_back(edgeTypes);
    }
     */
}
void EQFG::insertMetaPathToOneList(int YAGOeid, const vector<int> &nodeTypes, const vector<int> &edgeTypes, double w, map<int, map<int, double> > & PathList)
{
    string s_path = metapath2string(edgeTypes);
    if(metapath_2id_.find(s_path) == metapath_2id_.end()){
        metapath_2id_[s_path] = metapaths_.size();
        metapaths_.push_back(edgeTypes);
    }
    int metapathid = metapath_2id_[s_path];
    for(int i = 0 ; i < hin_.nodes_[YAGOeid].types_id_.size(); ++i){
        int typeID = hin_.nodes_[YAGOeid].types_id_[i];
        if(PathList.find(typeID) == PathList.end())
            PathList[typeID] = map<int, double>();
        if(PathList[typeID].find(metapathid) == PathList[typeID].end())
            PathList[typeID][metapathid] = 0.0;
        PathList[typeID][metapathid] += w;
    }
}
void EQFG::trainMetaPaths_all(string outDir)
{
    loadYagoByDefault();
    pll_.LoadIndex("./data/YAGO/yagoGraphTwoColumn.index");

    map<int, map<int, double> > PathList_randomShortestOne;
    map<int, map<int, double> > PathList_allShortestOnes;
    map<int, map<int, double> > PathList_length_0;
    map<int, map<int, double> > PathList_length_1;
    map<int, map<int, double> > PathList_length_2;
    map<int, map<int, double> > PathList_length_3;



    for(int i = 0; i < ENodes_.size(); ++i){
        cerr << "running the " << i<< " entity" << '\t' << ENodes_[i].toEntityEdges_.size() << endl;
        int YAGOeid1 = YAGOEID[i];
        if(YAGOeid1 == -1)
            continue;
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            int YAGOeid2 = YAGOEID[ENodes_[i].toEntityEdges_[j].eid_];
            if(YAGOeid2 == -1)
                continue;
            double w = ENodes_[i].toEntityEdges_[j].w_;
            int dis = pll_.QueryDistance(YAGOeid1, YAGOeid2);
            //cerr << ENodes_[i].toEntityEdges_[j].w_ << endl;
	       if(ENodes_[i].toEntityEdges_[j].w_ < 0.1)
                continue;
            if(dis > 3)
                continue;
	       vector< vector< vector< int> > > length2paths;
            //vector<int> pathNodes;
            //entityShortestPath(YAGOeid1, YAGOeid2, MAXDISTANCE_T, pathNodes);
            entityAllShortestPaths(YAGOeid1, YAGOeid2, 3, length2paths);

            for(int k = 0; k < length2paths[0].size(); ++k){
                vector<int> nodetypes, edgetypes;
                vector<int> pathNodes = length2paths[0][k];
                getMetaPath(pathNodes, nodetypes, edgetypes);
                double newW = w / (hin_.nodes_[YAGOeid1].types_id_.size() * length2paths[0].size());
                insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_length_0);
            }
            for(int k = 0; k < length2paths[1].size(); ++k){
                vector<int> nodetypes, edgetypes;
                vector<int> pathNodes = length2paths[1][k];
                getMetaPath(pathNodes, nodetypes, edgetypes);
                double newW = w / (hin_.nodes_[YAGOeid1].types_id_.size() * length2paths[1].size());
                insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_length_1);
            }
            for(int k = 0; k < length2paths[2].size(); ++k){
                vector<int> nodetypes, edgetypes;
                vector<int> pathNodes = length2paths[2][k];
                getMetaPath(pathNodes, nodetypes, edgetypes);
                double newW = w / (hin_.nodes_[YAGOeid1].types_id_.size() * length2paths[2].size());
                insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_length_2);
            }
            for(int k = 0; k < length2paths[3].size(); ++k){
                vector<int> nodetypes, edgetypes;
                vector<int> pathNodes = length2paths[3][k];
                getMetaPath(pathNodes, nodetypes, edgetypes);
                double newW = w / (hin_.nodes_[YAGOeid1].types_id_.size() * length2paths[3].size());
                insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_length_3);
            }
            
            int shortestLen = -1;
            for(int k = 0; k <= 3; ++k){
                if (length2paths[k].size() > 0){
                    shortestLen = k;
                    break;
                }
            }
            if (shortestLen != -1){
                vector<int> nodetypes, edgetypes;
                vector<int> pathNodes = length2paths[shortestLen][0];
                getMetaPath(pathNodes, nodetypes, edgetypes);
                double newW = w / hin_.nodes_[YAGOeid1].types_id_.size();
                insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_randomShortestOne);

                for(int k = 0; k < length2paths[shortestLen].size(); ++k){
                    vector<int> nodetypes, edgetypes;
                    vector<int> pathNodes = length2paths[shortestLen][k];
                    getMetaPath(pathNodes, nodetypes, edgetypes);
                    double newW = w / (hin_.nodes_[YAGOeid1].types_id_.size() * length2paths[shortestLen].size());
                    insertMetaPathToOneList(YAGOeid1, nodetypes, edgetypes, newW, PathList_allShortestOnes);
                }
            }
        }
    }

    saveMetaPathsList(outDir + "len_0.metapath", PathList_length_0);
    saveMetaPathsList(outDir + "len_1.metapath", PathList_length_1);
    saveMetaPathsList(outDir + "len_2.metapath", PathList_length_2);
    saveMetaPathsList(outDir + "len_3.metapath", PathList_length_3);

    saveMetaPathsList(outDir + "random_one.metapath", PathList_randomShortestOne);
    saveMetaPathsList(outDir + "all_shortest.metapath", PathList_allShortestOnes);

}


void EQFG::trainMetaPaths()
{
    loadYagoByDefault();
    
    pll_.LoadIndex("./data/YAGO/yagoGraphTwoColumn.index");
    
    for(int i = 0; i < ENodes_.size(); ++i){
        cerr << "running the " << i<< " entity" << '\t' << ENodes_[i].toEntityEdges_.size() << endl;
        int YAGOeid1 = YAGOEID[i];
        if(YAGOeid1 == -1)
            continue;
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            
            int YAGOeid2 = YAGOEID[ENodes_[i].toEntityEdges_[j].eid_];
            if(YAGOeid2 == -1)
                continue;
            double w = ENodes_[i].toEntityEdges_[j].w_;
            int dis = pll_.QueryDistance(YAGOeid1, YAGOeid2);
            if(ENodes_[i].toEntityEdges_[j].w_ < 0.01)
                continue;
            if(dis > MAXDISTANCE_T)
                continue;
            vector<int> pathNodes;
            entityShortestPath(YAGOeid1, YAGOeid2, MAXDISTANCE_T, pathNodes);
            vector<int> nodetypes, edgetypes;
            getMetaPath(pathNodes, nodetypes, edgetypes);
            double newW = w / hin_.nodes_[YAGOeid1].types_id_.size();
            insertMetaPath(YAGOeid1, nodetypes, edgetypes, newW);
        }
    }
}

void EQFG::selectTwoTypeEntity(int t1, int t2)
{
    //17924	3475257
    for(int i = 0; i < ENodes_.size(); ++i){
        int YAGOeid1 = YAGOEID[i];
        if(YAGOeid1 == -1)
            continue;
        if(find(hin_.nodes_[YAGOeid1].types_id_.begin(), hin_.nodes_[YAGOeid1].types_id_.end(), t1) == hin_.nodes_[YAGOeid1].types_id_.end())
            continue;
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            int YAGOeid2 = YAGOEID[ENodes_[i].toEntityEdges_[j].eid_];
            if(YAGOeid2 == -1)
                continue;
            if(find(hin_.nodes_[YAGOeid2].types_id_.begin(), hin_.nodes_[YAGOeid2].types_id_.end(), t2) == hin_.nodes_[YAGOeid2].types_id_.end())
                continue;
            cout << YAGOeid1 << '\t' << YAGOeid2 << endl;
        }
    }
}

void EQFG::loadYagoByDefault()
{
    hin_.buildYAGOGraphbyDefault();
    // do the entity matching
    int non_match = 0;
    for(int i = 0; i < entities_.size(); ++i){
        int tempid = entities_[i].rfind("/wiki/");
        string name = "<" + entities_[i].substr(tempid + 6) + ">";
        if(hin_.key2id_.find(name) == hin_.key2id_.end()){
            non_match ++;
            YAGOEID.push_back(-1);
        }else{
            YAGOEID.push_back(hin_.key2id_[name]);
        }
    }
    cerr << "num of non-match entities are: " << non_match << endl;
    //trainMetaPaths();
}


void EQFG::expand_exp(string queryEntityPath, string metapathPath)
{
	loadMetaPaths(metapathPath);
	string line;
	ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
	
	clock_t total = 0;
	int times = 0;
	
	while(getline(queryEntityIn, line)){
		while(line[line.size() - 1] == '\r')
        	line = line.substr(0, line.size() - 1);
        vector<string> strs= split(line, "\t");
        int qid = -1;
        if(query2id_.find(strs[1]) != query2id_.end())
            qid = query2id_[strs[1]];
        vector<string> entities;
        for(int i = 3; i < strs.size(); i += 2){
			string prefix = "http://en.wikipedia.org/wiki/";
			string tempS = "<";
			tempS += strs[i].substr(prefix.size());
			tempS += ">";
			cerr << tempS << endl;
			entities.push_back(tempS);
        }
		clock_t t1 = clock();
		map<int, double> expanded_eids = expandEntity_sp(entities);
		clock_t t2 = clock();
		total += t2 - t1;
		times ++;
	}
	double aver = (total + 0.0) / (CLOCKS_PER_SEC * times);
	cout << "average time for each expansion is \t" << aver << endl;
	
}

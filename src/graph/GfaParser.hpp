

#ifndef MDBG_METAG_GFAPARSER
#define MDBG_METAG_GFAPARSER

using namespace std;
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "Graph.hpp"
#include "../Commons.hpp"
//#include "GraphSimplify.hpp"

struct UnitigLength{
	//u_int32_t _index;
    u_int32_t _length;
    //u_int32_t _nodeIndex;
    
    float _abundance;
    u_int32_t _startNodeIndex;
};




class GfaParser{

public:

    static AdjGraph* createGraph_lol(const string& filename){

        AdjGraph* graph = new AdjGraph();

        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();
        u_int64_t unitig_id = 0;

        u_int32_t nbNodes = 0;

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            

            if((*fields)[0] == "S"){
                nbNodes += 1;
            }

            
        }

        for(size_t i=0; i<nbNodes; i++){
            graph->addNode(i);
        }

        u_int32_t nbEdges = 0;
        infile.clear();                 // clear fail and eof bits
        infile.seekg(0, std::ios::beg);

        unordered_set<DbgEdge, hash_pair> isEdgeVisited;

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << line << endl;
            //cout << (*fields)[0] << endl;

           if((*fields)[0] == "L"){
                string& from = (*fields)[1];
                bool directionFrom = (*fields)[2] == "+";
                string& to = (*fields)[3];
                bool directionTo = (*fields)[4] == "+";
                //u_int64_t overlap = std::stoull((*fields)[5]);
                //cout << line << endl;

                u_int32_t from_id = std::stoull(from);
                u_int32_t to_id = std::stoull(to);

                graph->addNode(from_id);
                graph->addNode(to_id);

                DbgEdge edge = {from_id, to_id};
				edge = edge.normalize();

				if(isEdgeVisited.find(edge) != isEdgeVisited.end()){
					continue;
				}

				isEdgeVisited.insert(edge);

                //if(from_id == 1651){
                //    cout << to_id << " " << (*fields)[2] << " " << (*fields)[4] << " " << isSuccessor << endl;
                //}

                bool isSuccessor = false;
                if(directionFrom && directionTo){
                    graph->addEdge_checkDuplicate(from_id, to_id, 0, true);
                    graph->addEdge_checkDuplicate(to_id, from_id, 0, false);
                }
                else if(!directionFrom && !directionTo){
                    graph->addEdge_checkDuplicate(from_id, to_id, 0, false);
                    graph->addEdge_checkDuplicate(to_id, from_id, 0, true);
                }
                else if(directionFrom && !directionTo){
                    
                    graph->addEdge_checkDuplicate(from_id, to_id, 0, false);
                    graph->addEdge_checkDuplicate(to_id, from_id, 0, true);
                }
                else if(!directionFrom && directionTo){
                    
                    graph->addEdge_checkDuplicate(from_id, to_id, 0, true);
                    graph->addEdge_checkDuplicate(to_id, from_id, 0, false);
                }


                

                //graph->addEdge_checkDuplicate(to_id, from_id, 0, !isSuccessor);

                nbEdges += 1;
                //graphInfo->addNode(from);
                //graphInfo->addNode(to);

                //graphEdges->push_back({graphInfo->name_to_id(from), graphInfo->name_to_id(to), overlap, directionFrom, directionTo, false});
                //graphEdges->push_back({graph->name_to_id(to), graph->name_to_id(from), overlap, directionTo, directionFrom, true});
            }

            //cout << line.substr(start, end - start) << endl;
            
        }

        cout << "hiiii: " << nbEdges << endl;
        delete fields;
        delete fields_optional;



        return graph;
    }

    static BiGraph* createBiGraph_lol(const string& filename, bool indexSuccessors, u_int32_t nbNodes){

        unordered_set<DbgEdge, hash_pair> seenEdges;

        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        if(nbNodes == 0){
            while (std::getline(infile, line)){
                if(line[0] == 'S') nbNodes += 1;
            }
        }


        //nbNodes = 5000000;
        //cout << "to remove" << endl;
        BiGraph* graph = new BiGraph(nbNodes);

        infile.clear();
        infile.seekg(0, std::ios::beg);


        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

           if((*fields)[0] == "L"){
                string& from = (*fields)[1];
                bool fromOrient = (*fields)[2] == "+";
                string& to = (*fields)[3];
                bool toOrient = (*fields)[4] == "+";
                u_int16_t overlap = std::stoull((*fields)[5]);

                u_int32_t from_id = std::stoull(from);
                u_int32_t to_id = std::stoull(to);
                
                //DbgEdge edge = {from_id, to_id};
				//edge = edge.normalize();

                //e_in, e_out = line.split()[1], line.split()[3]
                //e_tuple = tuple(sorted([e_in, e_out]))
                //to_remove = e_tuple in seen_edges
                //if e_in == e_out: # self loop
                //    to_remove = True
                //seen_edges.add(e_tuple)
                //if not to_remove:
                //    print(line.strip())

                //bool toRemove = seenEdges.find(edge) != seenEdges.end() || from_id == to_id;
                //seenEdges.insert(edge);
                //if(toRemove) continue;
                if(from_id == to_id) continue; //self loop

                if(indexSuccessors){
                    graph->addEdge(from_id, fromOrient, to_id, toOrient, overlap);
                }
                else{
                    graph->addEdge(to_id, toOrient, from_id, fromOrient, overlap);
                }
            }

            
        }

        delete fields;
        delete fields_optional;

        //cout << "allo ?? : " << graph->_nbNodes << endl;
        graph->_nodeAbundances.resize(graph->_nbNodes/2, 0);
        graph->_nodeLengths.resize(graph->_nbNodes/2, 0);
        //cout << graph->_nodeAbundances.size() << endl;
        GfaParser::getNodeData(filename, graph->_nodeAbundances, graph->_nodeLengths);
        //cout << graph->_nodeAbundances.size() << endl;

        return graph;
    }


    static BiGraph* createBiGraph_binary(const string& filename, bool indexSuccessors, u_int32_t nbNodes, float kminmerLengthMean, float kminmerOverlapMean){

        u_int8_t isS = 0;
        u_int8_t isL = 1;
        u_int8_t edgePlus = 0;
        u_int8_t edgeMinus = 1;

        
        ifstream infile(filename);
        /*


        if(nbNodes == 0){

            while(true){
                u_int8_t type;
                infile.read((char*)&type, sizeof(type));

                if(infile.eof())break;


                if(type == isS){

                    u_int32_t nodeName;
                    u_int32_t abundance;

                    infile.read((char*)&nodeName, sizeof(nodeName));
                    infile.read((char*)&abundance, sizeof(abundance));

                    nbNodes += 1;
                }
                else{

                    u_int32_t nodeName1;
                    u_int32_t nodeName2;
                    u_int8_t ori1;
                    u_int8_t ori2;

                    infile.read((char*)&nodeName1, sizeof(nodeName1));
                    infile.read((char*)&ori1, sizeof(ori1));
                    infile.read((char*)&nodeName2, sizeof(nodeName2));
                    infile.read((char*)&ori2, sizeof(ori2));
                }

            }
            
        }
        */

        infile.clear();
        infile.seekg(0, std::ios::beg);

        infile.read((char*)&nbNodes, sizeof(nbNodes));
        //nbNodes = 5000000;
        //cout << "to remove" << endl;
        BiGraph* graph = new BiGraph(nbNodes);


        graph->_nodeAbundances.resize(graph->_nbNodes/2, 0);
        graph->_nodeLengths.resize(graph->_nbNodes/2, 0);

        while(true){
            //u_int8_t type;
            //infile.read((char*)&type, sizeof(type));



            //if(type == isS){

            //    u_int32_t nodeName;
            //    u_int32_t abundance;

            //    infile.read((char*)&nodeName, sizeof(nodeName));
            //    infile.read((char*)&abundance, sizeof(abundance));
                //graph->_nodeAbundances[nodeName] = abundance;
                //graph->_nodeLengths[nodeName] = kminmerLengthMean;

            //    nbNodes += 1;
            //}
            //else{

                u_int32_t nodeName1;
                u_int32_t nodeName2;
                u_int8_t ori1;
                u_int8_t ori2;

                infile.read((char*)&nodeName1, sizeof(nodeName1));

                if(infile.eof())break;

                infile.read((char*)&ori1, sizeof(ori1));
                infile.read((char*)&nodeName2, sizeof(nodeName2));
                infile.read((char*)&ori2, sizeof(ori2));

                bool fromOrient = ori1 == edgePlus;
                bool toOrient = ori2 == edgePlus;
                /*
                string& from = (*fields)[1];
                bool fromOrient = (*fields)[2] == "+";
                string& to = (*fields)[3];
                bool toOrient = (*fields)[4] == "+";
                u_int16_t overlap = std::stoull((*fields)[5]);

                u_int32_t from_id = std::stoull(from);
                u_int32_t to_id = std::stoull(to);
                */
                
                if(nodeName1 == nodeName2) continue; //self loop

                graph->addEdge(nodeName1, fromOrient, nodeName2, toOrient, kminmerOverlapMean);

            //}

        }





        //cout << "allo ?? : " << graph->_nbNodes << endl;
        //cout << graph->_nodeAbundances.size() << endl;
        //GfaParser::getNodeData(filename, graph->_nodeAbundances, graph->_nodeLengths);
        //cout << graph->_nodeAbundances.size() << endl;

        return graph;
    }

    /*
    UnitigGraph* createGraph(const string& filename, vector<int32_t>& node_to_unitig, vector<u_int32_t>& unitigLengths){

        //unordered_map<u_int64_t, u_int64_t> counts;
        //GraphInfo* graphInfo = new GraphInfo();

        //vector<GraphEdge>* graphEdges = new vector<GraphEdge>();
        //vector<GraphNode>* graphNodes = new vector<GraphNode>();
        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();
        //u_int64_t unitig_id = 0;

        u_int32_t nbNodes = 0;

        while (std::getline(infile, line)){
            if(line[0] == 'S') nbNodes += 1;
        }

        UnitigGraph* graph = new UnitigGraph(nbNodes);

        infile.clear();
        infile.seekg(0, std::ios::beg);


        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                
                size_t i=0;
                for(string& field : (*fields)){
                    if(i> 2){
                        //cout << field << endl;
                        tokenize(field, fields_optional, ':');
                        if((*fields_optional)[0] == "LN"){
                            string& unitig_name = (*fields)[1];
                            u_int64_t unitig_length = std::stoull((*fields_optional)[2]);

                            //graph->addNode(unitigName_to_id(unitig_name));
                            unitigLengths.push_back(unitig_length);

                            //cout << unitig_name << " " <<unitigName_to_id(unitig_name) << endl;
                            //graphInfo->addNode(unitig_name);
                            //cout << graph->name_to_id(unitig_name) << " " << unitig_length << endl;
                            //graphNodes->push_back({graphInfo->name_to_id(unitig_name), unitig_length});
                            //cout << unitig_length << endl;
                        }
                    }
                    i += 1;
                }

                //unitig_id += 1;
            }
            else if((*fields)[0] == "L"){

                string& from = (*fields)[1];
                bool fromOrient = (*fields)[2] == "+";
                string& to = (*fields)[3];
                bool toOrient = (*fields)[4] == "+";
                u_int64_t overlap = std::stoull((*fields)[5]);

                u_int32_t from_id = unitigName_to_id(from);
                u_int32_t to_id = unitigName_to_id(to);

                graph->addEdge(from_id, fromOrient, to_id, toOrient, 0);

                
            }
            else if((*fields)[0] == "A"){

                
                //cout << line << endl;

                string unitig_name = (*fields)[1];
                string node = (*fields)[4];

                //if(151064 == std::stoull(node.c_str())){
                //    cout << "lala" << endl;
                //    cout << unitig_name << " " << unitigName_to_id(unitig_name) << endl;
                //}
                //if(unitigName_to_id(unitig_name) == 0){
                //    cout << unitig_name << " " << node << endl;
                //}
                //cout << unitig_name << endl;

                //cout << unitig_name << endl;
                //cout << node << " " << unitig_name << endl;
                node_to_unitig[std::stoull(node.c_str())] = unitigName_to_id(unitig_name);

                //counts[unitigName_to_id(unitig_name)] += 1;
                //cout << counts[unitigName_to_id(unitig_name)] << endl;
                //cout << line << endl;
                //cout << std::stoull(node.c_str()) << " " << std::stoull(unitig_name) << endl;
            }

            //cout << line.substr(start, end - start) << endl;
            
        }


        delete fields;
        delete fields_optional;



        return graph;
    }*/

    static u_int64_t getNodeToUnitig(const string& filename, vector<int32_t>& node_to_unitig){

        u_int64_t nbUnitigs = 0;

        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            if(line[0] == 'S'){
                nbUnitigs += 1;
            }
            else if(line[0] == 'A'){
                tokenize(line, fields, '\t');

                string unitig_name = (*fields)[1];
                string node = (*fields)[4];

                node_to_unitig[std::stoull(node.c_str())] = unitigName_to_id(unitig_name);

            }

        }


        delete fields;
        delete fields_optional;

        return nbUnitigs;
    }

    static u_int64_t getUnitigLengths(const string& filename, vector<UnitigLength>& unitigLengths){

        u_int64_t nbUnitigs = 0;

        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            if(line[0] == 'S'){

                tokenize(line, fields, '\t');

                size_t i=0;
                for(string& field : (*fields)){
                    if(i> 2){
                        //cout << field << endl;
                        tokenize(field, fields_optional, ':');
                        if((*fields_optional)[0] == "LN"){
                            string& unitig_name = (*fields)[1];
                            u_int64_t unitig_length = std::stoull((*fields_optional)[2]);
                            //cout << unitig_length << endl;
                            //graph->addNode(unitigName_to_id(unitig_name));
                            unitigLengths.push_back({unitigName_to_id(unitig_name), unitig_length});

                            //cout << unitig_name << " " << unitigName_to_id(unitig_name) << " " << unitig_length << endl;
                            //graphInfo->addNode(unitig_name);
                            //cout << graph->name_to_id(unitig_name) << " " << unitig_length << endl;
                            //graphNodes->push_back({graphInfo->name_to_id(unitig_name), unitig_length});
                            //cout << unitig_length << endl;
                        }
                    }
                    i += 1;
                }

                //unitigLengths.push_back();
                nbUnitigs += 1;
            }

        }


        delete fields;
        delete fields_optional;

        return nbUnitigs;
    }

    static void getNodeData(const string& filename, vector<u_int32_t>& nodeAbundances, vector<u_int32_t>& nodeLengths){ //unordered_map<u_int32_t, u_int32_t>& nodeAbundances, unordered_map<u_int32_t, u_int32_t>& nodeLengths){

        //nodeAbundances.clear();
        //nodeLengths.clear();

        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            if(line[0] == 'S'){

                tokenize(line, fields, '\t');

                size_t i=0;
                for(string& field : (*fields)){
                    if(i> 2){
                        
                        tokenize(field, fields_optional, ':');
                        if((*fields_optional)[0] == "dp"){
                            u_int32_t index = std::stoull((*fields)[1]);
                            u_int32_t abundance = std::stoull((*fields_optional)[2]);
                            nodeAbundances[index] = abundance;

                        }
                        else if((*fields_optional)[0] == "LN"){
                            u_int32_t index = std::stoull((*fields)[1]);
                            u_int32_t length = std::stoull((*fields_optional)[2]);
                            nodeLengths[index] = length;

                        }
                    }
                    i += 1;
                }

            }

        }


        delete fields;
        delete fields_optional;
    }

    BiGraph* createBiGraph(const string& filename, vector<int32_t>& node_to_unitig, vector<u_int32_t>& unitigLengths){


        ifstream infile(filename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        u_int32_t nbNodes = 0;
        while (std::getline(infile, line)){
            if(line[0] == 'S') nbNodes += 1;
        }
        //cout << "Gfa nb nodes: " << nbNodes << endl;

        BiGraph* graph = new BiGraph(nbNodes);


        infile.clear();
        infile.seekg(0, std::ios::beg);

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                
                size_t i=0;
                for(string& field : (*fields)){
                    if(i> 2){
                        //cout << field << endl;
                        tokenize(field, fields_optional, ':');
                        if((*fields_optional)[0] == "LN"){
                            string& unitig_name = (*fields)[1];
                            u_int64_t unitig_length = std::stoull((*fields_optional)[2]);

                            //graph->addNode(unitigName_to_id(unitig_name));
                            unitigLengths.push_back(unitig_length);

                            //cout << unitig_name << " " <<unitigName_to_id(unitig_name) << endl;
                            //graphInfo->addNode(unitig_name);
                            //cout << graph->name_to_id(unitig_name) << " " << unitig_length << endl;
                            //graphNodes->push_back({graphInfo->name_to_id(unitig_name), unitig_length});
                            //cout << unitig_length << endl;
                        }
                    }
                    i += 1;
                }

            }
            else if((*fields)[0] == "L"){
                string& from = (*fields)[1];
                bool fromOrient = (*fields)[2] == "+";
                string& to = (*fields)[3];
                bool toOrient = (*fields)[4] == "+";
                u_int64_t overlap = std::stoull((*fields)[5]);

                u_int32_t from_id = unitigName_to_id(from);
                u_int32_t to_id = unitigName_to_id(to);

                graph->addEdge(from_id, fromOrient, to_id, toOrient, 0);

                //cout << "Add edge" << endl;
            }
            else if((*fields)[0] == "A"){

                string unitig_name = (*fields)[1];
                string node = (*fields)[4];

                node_to_unitig[std::stoull(node.c_str())] = unitigName_to_id(unitig_name);

            }

            
        }


        delete fields;
        delete fields_optional;


        return graph;
    }


    static void tokenize(const string& line, vector<string>* tokens, char delim){

        tokens->clear();

        int start = 0;
        int end = line.find(delim);
        
        while (end != -1) {
            tokens->push_back(line.substr(start, end - start));
            //cout << line.substr(start, end - start) << endl;
            start = end + 1; //del.size();
            end = line.find(delim, start);
        }

        tokens->push_back(line.substr(start, end - start));

    }

    static u_int32_t unitigName_to_id(string unitig_name){
        unitig_name.erase(unitig_name.begin());
        unitig_name.erase(unitig_name.begin());
        unitig_name.erase(unitig_name.begin());
        unitig_name.erase(unitig_name.end()-1);
        return std::stoull(unitig_name)-1;
    }
    
    static string unitigIndex_to_unitigName(u_int32_t unitigIndex){

        string unitigName = "utg";
		string unitig_name_id = to_string(unitigIndex+1);
		size_t nbZeros = 7 - unitig_name_id.size();
		//cout << unitigIndex << " " << nbZeros << endl;
		for(size_t i=0; i<nbZeros; i++){
			unitigName += "0";
		}
		unitigName += unitig_name_id + "l";

        return unitigName;
    }

    static void rewriteGfa_withoutEdges(const string& filename, const string& outputFilename, const unordered_set<DbgEdge, hash_pair>& removedEdges){

        
        int nbEdges = 0;
        ifstream infile(filename);
	    ofstream outputFile(outputFilename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                outputFile << line << endl;
            }
            else if((*fields)[0] == "L"){
                string& from = (*fields)[1];
                string& to = (*fields)[3];

                u_int32_t utg_from = stoull(from); //unitigName_to_id(from);
                u_int32_t utg_to = stoull(to);; //unitigName_to_id(to);

                if(removedEdges.find({utg_from, utg_to}) != removedEdges.end()) continue;
                if(removedEdges.find({utg_to, utg_from}) != removedEdges.end()) continue;

                outputFile << line << endl;
                nbEdges += 1;
            }
            else {
                outputFile << line << endl;
            }
            
        }

        cout << "lalala: " << nbEdges << endl;

        delete fields;
        delete fields_optional;
        
        

    }

    static void rewriteGfa_withNodes(const string& filename, const string& outputFilename, const unordered_set<u_int32_t>& nodes){

        
        cout << filename << endl;
        cout << outputFilename << endl;

        ifstream infile(filename);
	    ofstream outputFile(outputFilename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                u_int32_t node = std::stoull((*fields)[1]);
                if(nodes.find(node) == nodes.end()) continue;
                if(node == 55479) cout << line << endl;
                outputFile << line << endl;
            }
            else if((*fields)[0] == "L"){
                u_int32_t from = std::stoull((*fields)[1]);
                u_int32_t to = std::stoull((*fields)[3]);

                if(nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) continue;

                outputFile << line << endl;
            }
            else {
                outputFile << line << endl;
            }
            
        }


        delete fields;
        delete fields_optional;
        
    }

    static void rewriteGfa_withoutNodes(const string& filename, const string& outputFilename, const unordered_set<u_int32_t>& nodes, const unordered_set<DbgEdge, hash_pair>& edges, BiGraph* graph){

        
        cout << filename << endl;
        cout << outputFilename << endl;

        ifstream infile(filename);
	    ofstream outputFile(outputFilename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                u_int32_t node = std::stoull((*fields)[1]);
                //if(node == 214){
                //    cout << "FOUNDED" << endl;
                //    cout << (nodes.find(node) != nodes.end()) << endl;
                //}
                if(nodes.find(node) == nodes.end()) continue;
                outputFile << line << endl;
            }
            else if((*fields)[0] == "L"){
                u_int32_t from = std::stoull((*fields)[1]);
                bool from_orient = (*fields)[2] == "+";
                u_int32_t to = std::stoull((*fields)[3]);
                bool to_orient = (*fields)[4] == "+";

                if(nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) continue;

                u_int32_t nodeIndex_from = graph->nodeName_to_nodeIndex(from, from_orient);
                u_int32_t nodeIndex_to = graph->nodeName_to_nodeIndex(to, to_orient);
                
                if(!graph->edgeExists(nodeIndex_from, nodeIndex_to)) continue;
                //DbgEdge edge = {nodeIndex_from, nodeIndex_to};
                //edge = edge.normalize();
                //if(edges.find(edge) != edges.end()) continue;

                outputFile << line << endl;
            }
            else {
                outputFile << line << endl;
            }
            
        }


        delete fields;
        delete fields_optional;
        
    }



    void rewriteGfa_withUnitigs(const string& filename, const string& outputFilename, const unordered_set<u_int32_t>& nodes, const vector<UnitigData>& unitigDatas){

        
        cout << filename << endl;
        cout << outputFilename << endl;

        ifstream infile(filename);
	    ofstream outputFile(outputFilename);

        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            
            tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            if((*fields)[0] == "S"){
                u_int32_t node = unitigName_to_id((*fields)[1]);
                //u_int32_t node = std::stoull((*fields)[1]);
                if(nodes.find(node) == nodes.end()) continue;
                outputFile << line << endl; //<< "\tdp:f:" << unitigDatas[node]._meanAbundance << endl;
            }
            else if((*fields)[0] == "L"){
                u_int32_t from = unitigName_to_id((*fields)[1]);//std::stoull((*fields)[1]);
                u_int32_t to = unitigName_to_id((*fields)[3]);//std::stoull((*fields)[3]);

                if(nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) continue;

                outputFile << line << endl;
            }
            else {
                //outputFile << line << endl;
            }
            
        }


        delete fields;
        delete fields_optional;
        
    }

    /*
    void writeGfa(AdjGraph* graph, const string& outputFilename, const vector<u_int32_t>& unitigsLength){

        cout << "Dumping graph: " << outputFilename << endl;

        //cout << graph->_nbNodes << endl;
	    //vector<u_int64_t> neighbors;

	    ofstream outputFile(outputFilename);

        for(size_t n=0; n<graph->_nbNodes; n++){
            //if(graph->_nodes[n]->isBidirection) continue;

            string utg_name = unitigIndex_to_unitigName(n);
            //cout << n << " " << graph->_graphInfo->_unitigs_length.size() << endl;
            outputFile << "S" << "\t" << utg_name << "\t" << "*" << "\t" << "LN:i:" << unitigsLength[n] << endl;


            adjNode* node = graph->_nodes[n];
            while (node != nullptr) {

                string from = "";
                if(node->from_direction) from = "+"; else from = "-";

                string to = "";
                if(node->to_direction) to = "+"; else to = "-";

                outputFile << "L" << "\t" << utg_name << "\t" << from << "\t" << unitigIndex_to_unitigName(node->val) << "\t" << to << "\t" << "0M" << endl;

                node = node->next;

            }


        }

        outputFile.close();
    }*/

};



#endif
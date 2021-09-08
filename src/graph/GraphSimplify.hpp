

//#include "Graph.hpp"
//TODO:
// remove bubble, tips: consider length

//#define PRINT_DEBUG_SIMPLIFICATION

using namespace std;


struct Unitig{
    u_int32_t _index;
    u_int32_t _startNode;
    u_int32_t _endNode;
    float _abundance;
    u_int32_t _length;
};

class GraphSimplify{

public:

    BiGraph* _graphSuccessors;
    BiGraph* _graphPredecessors;
    string _outputDir;
    string _inputGfaFilename;
    vector<Unitig> _unitigs;
    vector<u_int32_t> _nodeToUnitig;
    vector<u_int32_t> _nodeAbundances;
    vector<u_int32_t> _nodeLengths;
    //unordered_map<u_int32_t, u_int32_t> _nodeAbundances;
    //unordered_map<u_int32_t, u_int32_t> _nodeLengths;
    vector<bool> _isNodeRemoved;

    GraphSimplify(const string& inputGfaFilename, const string& outputDir){
        _inputGfaFilename = inputGfaFilename;
        _outputDir = outputDir;

        _graphSuccessors = GfaParser::createBiGraph_lol(inputGfaFilename, true);
	    _graphPredecessors = GfaParser::createBiGraph_lol(inputGfaFilename, false);
        _nodeAbundances.resize(_graphSuccessors->_nbNodes/2, false);
        _nodeLengths.resize(_graphSuccessors->_nbNodes/2, false);
        GfaParser::getNodeData(inputGfaFilename, _nodeAbundances, _nodeLengths);
        _isNodeRemoved.resize(_graphSuccessors->_nbNodes, false);
    }


    void execute(){


        bool dummy;
        unordered_set<u_int32_t> removedNodes;

        
        while(true){


            while(true){
                compact();
                u_int64_t nbBubblesRemoved = bubble(10000);
                //cout << nbBubblesRemoved << endl;
                if(nbBubblesRemoved == 0) break;
            }

            u_int64_t nbTipsRemoved = tip();
            //cout << nbTipsRemoved << endl;
            if(nbTipsRemoved == 0) break;

        }

        


        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            if(!_isNodeRemoved[nodeIndex]) continue;
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
            removedNodes.insert(nodeName);
        }

        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_tmp.gfa", removedNodes);


    }

    u_int64_t tip(){

        u_int64_t nbRemoved = 0;

        vector<u_int32_t> neighbors;
        vector<u_int32_t> unitigNodes;

        vector<bool> isVisited(_graphSuccessors->_nbNodes, false);
        bool dummy = false;

        for(Unitig& unitig : _unitigs){

            if(isVisited[unitig._startNode]) continue;
            if(isVisited[unitig._endNode]) continue;

            getPredecessors(unitig._startNode, 0, neighbors);
            if(neighbors.size() == 0) continue;

            getSuccessors(unitig._endNode, 0, neighbors);
            if(neighbors.size() > 0) continue;

            isVisited[unitig._startNode] = true;
            isVisited[unitig._endNode] = true;

            getUnitigNodes(unitig, unitigNodes);
            for(u_int32_t node : unitigNodes){
                _isNodeRemoved[node] = true;
            }

            #ifdef PRINT_DEBUG_SIMPLIFICATION
                cout << "\tTip: " << _graphSuccessors->nodeIndex_to_nodeName(unitig._endNode, dummy) << endl;
            #endif 

            nbRemoved += 1;
        }

        return nbRemoved;
    }

    u_int64_t bubble(float maxLength){

        u_int64_t nbBubblesRemoved = 0;
        //vector<u_int32_t> removedUnitigs;

        vector<bool> isVisited(_graphSuccessors->_nbNodes, false);
        vector<u_int32_t> unitigNodes;
        bool dummy = false;

        for(Unitig& unitig : _unitigs){
            //cout << "-------------" << endl;



            //u_int32_t startNodeName = _graphSuccessors->nodeIndex_to_nodeName(unitig._startNode, dummy);
            //u_int32_t visitedNodeIndex = unitig._startNode;
            if(isVisited[unitig._startNode]) continue;
            

            //cout << "--------------" << endl;
            //if(isVisited[unitig._startNode]) continue;
            //if(isVisited[unitig._endNode]) continue;

            //isVisited[unitig._startNode] = true;
            //isVisited[unitig._endNode] = true;


            Unitig& utg_1 = _unitigs[unitig._index];

            //getNeighbors(utg_1._endNode, orientation, neighbors);
            vector<u_int32_t> neighbors_utg1;
            getSuccessors(utg_1._endNode, 0, neighbors_utg1);
            

            //cout << _graphSuccessors->nodeIndex_to_nodeName(utg_1._endNode, dummy) << " " << neighbors.size() << endl;

            //else if(_graphSuccessors->nodeIndex_to_nodeName(unitig._startNode, dummy) == 6734){
            //    cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
                //cout << _graphSuccessors->nodeIndex_to_nodeName(utg_1._endNode, dummy) << " " << neighbors.size() << endl;
            //}
            
            if(neighbors_utg1.size() <= 1) continue;
                
            bool isBubble = false;

            for(size_t i=0; i<neighbors_utg1.size(); i++) {
                

                Unitig& utg_2 = _unitigs[_nodeToUnitig[neighbors_utg1[i]]];
                //visitedNodeIndex = _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy);
                if(isVisited[utg_2._startNode]) continue;
                //isVisited[startNodeName] = true;



                for(size_t j=i+1; j<neighbors_utg1.size(); j++){
                    Unitig& utg_3 = _unitigs[_nodeToUnitig[neighbors_utg1[j]]];
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy);
                    if(isVisited[utg_3._startNode]) continue;
                    //isVisited[startNodeName] = true;


                    vector<u_int32_t> neighbors_utg2;
                    getPredecessors(utg_2._startNode, 0, neighbors_utg2);
                    if(neighbors_utg2.size() != 1) continue;
                    getSuccessors(utg_2._endNode, 0, neighbors_utg2);
                    if(neighbors_utg2.size() != 1) continue;
                    
                    Unitig& utg_4 = _unitigs[_nodeToUnitig[neighbors_utg2[0]]];

                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy);
                    if(isVisited[utg_4._startNode]) continue;
                    //isVisited[startNodeName] = true;


                    vector<u_int32_t> neighbors_utg3;
                    getPredecessors(utg_3._startNode, 0, neighbors_utg3);
                    if(neighbors_utg3.size() != 1) continue;
                    getSuccessors(utg_3._endNode, 0, neighbors_utg3);
                    if(neighbors_utg3.size() != 1) continue;

                    if(_unitigs[_nodeToUnitig[neighbors_utg3[0]]]._index != utg_4._index) continue;
                    if(_graphSuccessors->nodeIndex_to_nodeName(utg_1._endNode, dummy) == _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy)) continue; //Repeated unitig with a cycle on on side
                    if(utg_2._length > maxLength || utg_3._length > maxLength) continue;

                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "\tBubble: " << _graphSuccessors->nodeIndex_to_nodeName(utg_1._endNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy) << " " << utg_2._length << " " << utg_3._length << endl;
                    #endif

                    if(utg_2._abundance > utg_3._abundance){
                        getUnitigNodes(utg_3, unitigNodes);

                        isVisited[utg_3._startNode] = true;
                        isVisited[utg_3._endNode] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy);
                        //isVisited[startNodeName] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._endNode, dummy);
                        //isVisited[startNodeName] = true;
                    }
                    else{
                        getUnitigNodes(utg_2, unitigNodes);

                        isVisited[utg_2._startNode] = true;
                        isVisited[utg_2._endNode] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy);
                        //isVisited[startNodeName] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_2._endNode, dummy);
                        //isVisited[startNodeName] = true;
                    }

                    for(u_int32_t node : unitigNodes){
                        _isNodeRemoved[node] = true;
                        //removedNodes.insert(node);
                    }
                    
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_1._startNode, dummy);
                    //isVisited[startNodeName] = true;
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy);
                    //isVisited[startNodeName] = true;
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy);
                    //isVisited[startNodeName] = true;
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy);
                    //isVisited[startNodeName] = true;

                    isBubble = true;
                    nbBubblesRemoved += 1;

                    break;

                }

                if(isBubble) break;
            }

        }


        return nbBubblesRemoved;
    }



    void compact(){

        _unitigs.clear();
        _nodeToUnitig.clear();

        _nodeToUnitig.resize(_graphSuccessors->_nbNodes, -1);
        //ofstream outfile(_outputDir + "/" + "tmp_compacted_graph.gfa");

        //vector<u_int32_t> nodeToUnitig(_graphSuccessors->_nbNodes, -1);
        //vector<Unitig> unitigs;
        //unordered_set<u_int32_t> nodes;
        //unordered_map<u_int32_t, u_int32_t> startNodes;
        //unordered_map<u_int32_t, u_int32_t> endNodes;

        u_int64_t unitigIndex = 0;
        
        vector<u_int32_t> neighbors;
        bool dummy = false;
        vector<bool> isVisited(_graphSuccessors->_nbNodes, false);


        

        //size_t nodeIndex = _graphSuccessors->nodeName_to_nodeIndex(720, true);

        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            //if(nodeIndex % 2 == 1) continue;

            if(_isNodeRemoved[nodeIndex]) continue;
            if(isVisited[nodeIndex]) continue;

            isVisited[nodeIndex] = true;
            //u_int32_t nodeName = _graphSuccessors->nodeName_to_nodeIndex(7701, true);

            u_int32_t startNode = nodeIndex;
            u_int32_t endNode = nodeIndex;

            //forward
            while(true){
                getSuccessors(endNode, 0, neighbors);
                if(neighbors.size() != 1) break;
                u_int32_t successor = neighbors[0];
                getPredecessors(successor, 0, neighbors);
                if(neighbors.size() != 1) break;
                endNode = successor;
            }

            
            //backward
            while(true){
                getPredecessors(startNode, 0, neighbors);
                if(neighbors.size() != 1) break;
                u_int32_t predecessor = neighbors[0];
                getSuccessors(predecessor, 0, neighbors);
                if(neighbors.size() != 1) break;
                startNode = predecessor;
            }
            

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy);

            //startNodes[startNode] = unitigIndex;
            //endNodes[endNode] = unitigIndex;
            //nodes.insert(endNode);

            /*
            //if(unitigIndex == 1575){
                cout << "-----------------" << endl;
                cout << "Unitig index: " << unitigIndex << endl;
                cout << "Original node: " << _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy) << endl;
                cout << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
                cout << startNode << " " << endNode << endl;
                
            //}
            */

            isVisited[startNode] = true;
            isVisited[endNode] = true;

            //cout << unitigIndex << " " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << endl;
            //int lala = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
            //if(lala == 4452 && dummy){
            //    cout << "lal " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
            //}
            //lala = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
            //if(lala == 4452 && !dummy){
            //    cout << "loul " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
            //}
            _nodeToUnitig[startNode] = unitigIndex;
            _nodeToUnitig[endNode] = unitigIndex;

            u_int32_t abundance_sum = 0;
            u_int32_t nbNodes = 0;
            u_int32_t abundance_max = 0;

            u_int32_t length = 0;
            u_int32_t node = startNode;
            //cout << "---------------" << endl;
            while(true){

                //cout << _graphSuccessors->nodeIndex_to_nodeName(node, dummy) << endl;
                u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(node, dummy);

                if(_nodeAbundances[nodeName] > abundance_max){
                    abundance_max = _nodeAbundances[nodeName];
                }

                length += _nodeLengths[nodeName];
                //abundance_sum += _nodeAbundances[nodeName];
                //nbNodes += 1;

                getSuccessors(node, 0, neighbors);
                isVisited[node] = true;
                _nodeToUnitig[node] = unitigIndex;

                if(node == endNode) break;
                node = neighbors[0];

            }

            //float abundance = ((float) abundance_sum) / nbNodes;
            //cout << abundance << endl;
            //if(_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy) == 4453){
            //    cout << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
            //}

            //if(nodeIndex % 2 == 0) continue;

            _unitigs.push_back({unitigIndex, startNode, endNode, abundance_max, length});
            //outfile << "S" << "\t" << unitigIndex << "\t" << "*" << endl;
            unitigIndex += 1;

            //cout << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << endl;
            //cout << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << endl;
            //cout << _graphSuccessors->nodeIndex_to_nodeName(neighbors[0], dummy) << endl;
            //getPredecessors(nodeName, neighbors);
            //cout << _graphSuccessors->nodeIndex_to_nodeName(neighbors[0], dummy) << endl;

        }

        #ifdef PRINT_DEBUG_SIMPLIFICATION
            cout << "Nb unitigs: " << unitigIndex << endl;
        #endif

        /*
        for(Unitig& unitig : unitigs){

            int lala = _graphSuccessors->nodeIndex_to_nodeName(unitig._startNode, dummy);
            if(lala == 4452 && dummy){
                cout << "4452+ " << nodeToUnitig[unitig._startNode] << endl;
                getPredecessors(unitig._startNode, neighbors);
                for(u_int32_t nn : neighbors){
                    bool orientation;
                    int loulou = _graphSuccessors->nodeIndex_to_nodeName(nn, orientation);
                    cout << loulou << " " << orientation << "    " << nodeToUnitig[nn] << endl;
                }
            }
            else if(lala == 4453 && dummy){
                cout << "4453+ " << nodeToUnitig[unitig._startNode] << endl;
                getPredecessors(unitig._startNode, neighbors);
                for(u_int32_t nn : neighbors){
                    bool orientation;
                    int loulou = _graphSuccessors->nodeIndex_to_nodeName(nn, orientation);
                    cout << loulou << " " << orientation << "    " << nodeToUnitig[nn] << endl;
                }
            }
            
        }
        */

        /*
        ifstream infile(_inputGfaFilename);


        std::string line;
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

        while (std::getline(infile, line)){
            
            GfaParser::tokenize(line, fields, '\t');
            
            //cout << (*fields)[0] << endl;

            //if((*fields)[0] == "S"){
            //    outfile << line << endl;
            //}
            //else 
            
            if((*fields)[0] == "L"){
                string& from = (*fields)[1];
                string& fromOrient = (*fields)[2];
                string& to = (*fields)[3];
                string& toOrient = (*fields)[4];

                u_int32_t fromName = stoull(from); //unitigName_to_id(from);
                u_int32_t toName = stoull(to);; //unitigName_to_id(to);

                u_int32_t fromNodeIndex = _graphSuccessors->nodeName_to_nodeIndex(fromName, fromOrient == "+");
                u_int32_t toNodeIndex = _graphSuccessors->nodeName_to_nodeIndex(toName, toOrient == "+");

               
                //if(nodes.find(fromNodeIndex) == nodes.end()) continue;
                //if(nodes.find(toNodeIndex) == nodes.end()) continue;

                //if(nodes[fromNodeIndex] == nodes[toNodeIndex]) continue;
                //if(fromNodeIndex == 5432 || toNodeIndex == 5432){
                  //  cout << nodes.find(fromNodeIndex) << " " << nodes.find(toNodeIndex) << endl;
                //}
                //if(nodes.find(fromNodeIndex) == nodes.end()) continue;
                //if(nodes.find(toNodeIndex) == nodes.end()) continue;
                //u_int32_t unitigIndex_from = nodeToUnitigs[fromNodeIndex];
                //u_int32_t unitigIndex_to = nodeToUnitigs[toNodeIndex];

                //if(unitigIndex_from == -1 || unitigIndex_to == -1) continue;
                //cout << fromName << " " << unitigIndex_from << " " << toName << " " << unitigIndex_to << endl;
            
                //outfile << "L" << "\t" << nodes[fromNodeIndex] << "\t" <<  fromOrient << "\t" << nodes[toNodeIndex]  << "\t" << toOrient << "\t" << "0M" << endl;
            }
            //else {
            //    outfile << line << endl;
            //}
            
        }


        delete fields;
        delete fields_optional;
        */


        //outfile.close();
    }

    void getSuccessors(u_int32_t n, float abundanceCutoff_min, vector<u_int32_t>& successors){

        bool dummy;
        successors.clear();
        adjNode* node = _graphSuccessors->_nodes[n];

        while (node != nullptr) {

			u_int32_t nn = node->val;
            if(_isNodeRemoved[nn]){
                node = node->next;
                continue;
            }

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
			if(abundanceCutoff_min != 0 && _unitigs[_nodeToUnitig[nn]]._abundance < abundanceCutoff_min){ //Abundance min cutoff
				node = node->next;
				continue;
			}

            successors.push_back(nn);
            node = node->next;
        }
    }

    void getPredecessors(u_int32_t n, float abundanceCutoff_min, vector<u_int32_t>& predecessors){

        bool dummy;
        predecessors.clear();
        adjNode* node = _graphPredecessors->_nodes[n];
        
        while (node != nullptr) {

			u_int32_t nn = node->val;
            if(_isNodeRemoved[nn]){
                node = node->next;
                continue;
            }

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
			if(abundanceCutoff_min != 0 && _unitigs[_nodeToUnitig[nn]]._abundance < abundanceCutoff_min){ //Abundance min cutoff
				node = node->next;
				continue;
			}

            predecessors.push_back(nn);
            node = node->next;
        }
    }

    void getUnitigNodes(Unitig& unitig, vector<u_int32_t>& nodes){

        nodes.clear();
        bool dummy;
        u_int32_t node = unitig._startNode;
        vector<u_int32_t> neighbors;

        while(true){

            //if(_isNodeRemoved[node]) continue;

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(node, dummy);
            nodes.push_back(node);

            getSuccessors(node, 0, neighbors);

            if(node == unitig._endNode) break;
            node = neighbors[0];

        }

    }

};
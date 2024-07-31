


#ifndef MDBG_METAG_GRAPH
#define MDBG_METAG_GRAPH


#include <iostream>
#include <unordered_map>
#include <vector>
#include <queue>

using namespace std;
typedef pair<int, int> iPair;




struct adjNode {
    u_int32_t val;
    float weight;
    adjNode* next;
};

struct GraphEdge {
    u_int64_t start_ver;
    u_int64_t end_ver;
    float weight;
};

struct GraphNode {
    u_int64_t _nodeID, _length;
};



class BiGraph{




public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;
    //u_int32_t _nodeIndex;

    vector<vector<u_int32_t>> _nodes;

    
    BiGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes * 2;
        _nodes.resize(_nbNodes);


        _nbEdges = 0;
    }

    ~BiGraph() {
    }


    static u_int32_t nodeName_to_nodeIndex(u_int32_t nodeName, bool isOrientationForward){
        if(isOrientationForward){
            return nodeName * 2;
        }
        else{
            return nodeName * 2 + 1;
        }
    }

    static u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex, bool& isOrientationForward){
        isOrientationForward = (nodeIndex % 2) == 0;
        return nodeIndex / 2;
    }

    static u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex){
        //isOrientationForward = (nodeIndex % 2) == 0;
        return nodeIndex / 2;
    }

    bool addEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient, u_int16_t overlap){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, fromOrient);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, toOrient);

        _nodes[nodeIndex_from].push_back(nodeIndex_to);
        _nbEdges += 1;

        return true;
    }


    bool edgeExists(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){

        vector<u_int32_t>& successors = _nodes[nodeIndex_from];
        // Traversing through the first vector list
        // and removing the second element from it
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i] == nodeIndex_to) {
                return true;
            }
        }
    
        return false;
    }

    static string nodeToString(u_int32_t nodeIndex){
        bool orient;
        u_int32_t nodeName = nodeIndex_to_nodeName(nodeIndex, orient);
        string orientStr = "";
        if(orient){
            orientStr = "+";
        }
        else{
            orientStr = "-";
        }
        //string unitigName = unitigIndex_to_unitigName(nodeName);
        string unitigName = to_string(nodeName); //unitigIndex_to_unitigName(nodeName);
        return unitigName + orientStr;
    }



};

class UnitigGraph{

public:

    class Node{
        public:

        u_int32_t _unitigIndex;
        vector<u_int32_t> _nodes;
        //u_int32_t _startNode;
        //u_int32_t _endNode;
        float _abundance;
        u_int32_t _length;
        //u_int32_t _nbNodes;
        //vector<u_int32_t> _nodes;
        vector<u_int32_t> _successors;
        //vector<Node*> _predecessors;
        vector<float> _abundances;
        //vector<u_int32_t> _readIndexes;
        //vector<vector<u_int32_t>> _successorUnitigPaths;
        //u_int32_t _sortingIndex;

        //u_int32_t _nodeIndexStart_subStart;
        //u_int32_t _nodeIndexStart_subEnd;
        //u_int32_t _nodeIndexEnd_subStart;
        //u_int32_t _nodeIndexEnd_subEnd;

        Node(UnitigGraph* graph, u_int32_t unitigIndex, const vector<u_int32_t>& nodes, const vector<u_int32_t>& nodeDatas){
            _unitigIndex = unitigIndex;
            _nodes = nodes;

            _length = graph->_kminmerLength + ((_nodes.size()-1)*graph->_kminmerLengthNonOverlap);
            

            for(u_int32_t nodeIndex : nodes){
                _abundances.push_back(nodeDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]);
            }

            _abundance = computeMedianAbundance(_abundances);

            //_sortingIndex = startNode();
        }


        vector<Node*> getSuccessors(UnitigGraph* unitigGraph) const{

            vector<Node*> successors;

            for(u_int32_t unitigIndex : _successors){
                successors.push_back(unitigGraph->_nodes[unitigIndex]);
            }

            return successors;

        }

        vector<Node*> getPredecessors(UnitigGraph* unitigGraph) const{

            const Node* unitigRC = unitigGraph->unitigIndex_toReverseDirection(this);
            vector<Node*> predecessors;
            
            for(u_int32_t unitigIndex : unitigRC->_successors){

                predecessors.push_back(unitigGraph->_nodes[unitigGraph->unitigIndex_toReverseDirection(unitigIndex)]);
            }

            return predecessors;
        }

        u_int64_t nodeSum(){
            u_int64_t sum = 0;
            for(u_int32_t nodeIndex : _nodes){
                sum += nodeIndex;
            }
            return sum;
        }

        u_int64_t successorSum(UnitigGraph* unitigGraph){
            
            u_int64_t sum = 0;
            for(Node* nn : getSuccessors(unitigGraph)){
                sum += nn->startNode() + nn->endNode();
            }
            return sum;
            /*
            u_int64_t sum = 0;
            for(Node* nn : _successors){
                for(u_int32_t nodeIndex :  nn->_nodes){
                    sum += nodeIndex;
                }
            }
            return sum;
            */
        }

        u_int64_t predecessorSum(UnitigGraph* unitigGraph){

            u_int64_t sum = 0;
            for(Node* nn : getPredecessors(unitigGraph)){
                sum += nn->startNode() + nn->endNode();
            }
            return sum;
            /*
            u_int64_t sum = 0;
            for(Node* nn : _predecessors){
                for(u_int32_t nodeIndex :  nn->_nodes){
                    sum += nodeIndex;
                }
            }
            return sum;
            */
        }

        u_int32_t startNode(){
            return _nodes[0];
        }

        u_int32_t endNode(){
            return _nodes[_nodes.size()-1];
        }

        u_int8_t isCircular(UnitigGraph* unitigGraph){
            const vector<Node*>& successors = getSuccessors(unitigGraph);
            const vector<Node*>& predecessors = getPredecessors(unitigGraph);
            return _nodes.size() > 1 && (successors.size() == 1) && (predecessors.size() == 1) && (successors[0] == this) && (predecessors[0] == this);  //_nodes.size() > 1 && (startNode() == endNode());
        }
        void mergeWith(Node* node2, u_int32_t kminmerLength, u_int32_t kminmerLengthNonOverlap){


            //_length += (node2->_nodes.size()*kminmerLengthNonOverlap);
            for(float abundance : node2->_abundances){
                _abundances.push_back(abundance);
            }



            _abundance = computeMedianAbundance(_abundances);
            
            _nodes.insert( _nodes.end(), node2->_nodes.begin(), node2->_nodes.end());

            _length = kminmerLength + ((_nodes.size()-1)*kminmerLengthNonOverlap);

            //unordered_set<u_int32_t> readIndexes;
            //for(u_int32_t readIndex : _readIndexes){
            //    readIndexes.insert(readIndex);
            //}
            //for(u_int32_t readIndex : node2->_readIndexes){
            //    readIndexes.insert(readIndex);
            //}

            //_readIndexes.clear();
            //for(u_int32_t readIndex : readIndexes){
            //    _readIndexes.push_back(readIndex);
            //}

            //if(_abundances.size() != _nodes.size()){
            //    cout << "omg" << endl;
            //    exit(1);
            //}

            /*
            for(u_int32_t nodeIndex : _nodes){
                if(nodeIndex == 585233){
                    cout << "haha " << _length << endl;
                    getchar();
                }
                if(nodeIndex == 11822266){
                    cout << "houhou " << _length << endl;
                    getchar();
                }
            }*/

        }

        double compute_median_float(vector<float>& scores){
            size_t size = scores.size();

            if (size == 0){
                return 0;  // Undefined, really.
            }
            else{
                sort(scores.begin(), scores.end());
                if (size % 2 == 0){
                    return (scores[size / 2 - 1] + scores[size / 2]) / 2.0; //scores[size / 2 - 1];
                }
                else {
                    return scores[size / 2];
                }
            }
        }

        double computeMedianAbundance(vector<float> abundances){

            return compute_median_float(abundances);

            /*
            size_t n = abundances.size();

            if (n % 2 == 0) {

                nth_element(abundances.begin(), abundances.begin() + n / 2, abundances.end());
                nth_element(abundances.begin(), abundances.begin() + (n - 1) / 2, abundances.end());

                return (double)(abundances[(n - 1) / 2]+ abundances[n / 2])/ 2.0;
            }

            else {

                nth_element(abundances.begin(), abundances.begin() + n / 2, abundances.end());

                return (double)abundances[n / 2];
            }
            */
        }


    };

	static bool NodeComparatorByLength(UnitigGraph::Node*a, UnitigGraph::Node*b){

        if(a->_nodes.size() == b->_nodes.size()){
            return a->startNode() > b->startNode();
        }
        return a->_nodes.size() > b->_nodes.size();
        //if(a._startNodeIndex == b._startNodeIndex){
        //    return a._length < b._length;
        //}
        //return a._startNodeIndex < b._startNodeIndex;
        //if(a._length == b._length){
        //    return a._startNodeIndex < b._startNodeIndex;
        //}
		//return a._length < b._length;
	}

    
	static bool NodeComparatorByLengthRev(UnitigGraph::Node*a, UnitigGraph::Node*b){

        if(a->_nodes.size() == b->_nodes.size()){
            return a->startNode() > b->startNode();
        }
        return a->_nodes.size() < b->_nodes.size();
        
	}

	static bool NodeComparator(UnitigGraph::Node*a, UnitigGraph::Node*b){
        return a->startNode() < b->startNode();
	}

    //class Edge{
    //    u_int32_t _indexFrom;
    //    u_int32_t _indexTo;
    //};


    u_int32_t _currentUnitigIndex;
    vector<Node*> _nodes;
    //vector<u_int32_t> _nodeIndex_to_unitigIndex;
	float _kminmerLength;
	float _kminmerOverlapMean;
	float _kminmerLengthNonOverlap;
    //const vector<u_int32_t>& _nodeDatas;

    class GraphChangedIndex{
        public:
        unordered_set<u_int32_t> _superBubbleIndex2_notSuperbubble;
        unordered_map<u_int32_t, vector<u_int32_t>> _superBubbleIndex2_notSuperbubbleNodes;
    };

    GraphChangedIndex _graphChangedIndex;

    UnitigGraph(float kminmerLength, float kminmerOverlapMean, float kminmerLengthNonOverlap){//}, const vector<u_int32_t>& nodeDatas) : _nodeDatas(nodeDatas){
        _currentUnitigIndex = 0;
        _kminmerLength = kminmerLength;
        _kminmerOverlapMean = kminmerOverlapMean;
        _kminmerLengthNonOverlap = kminmerLengthNonOverlap;
        //_nodeIndex_to_unitigIndex.resize(nbNodes, -1);
    }

    ~UnitigGraph(){
        for(size_t i=0; i<_nodes.size(); i++){
            delete _nodes[i];
        }
    }

    void load(const string& nodeFilename, const string& edgeSuccessorFilename, const string& edgePredecessorFilename, vector<u_int32_t>& nodeDatas){

        cout << "Loading unitig graph nodes" << endl;

		ifstream nodeFile(nodeFilename);

		while(true){

			u_int32_t size;
			nodeFile.read((char*)&size, sizeof(size));
			

			if(nodeFile.eof()) break;

			vector<u_int32_t> nodePath;
			nodePath.resize(size);
			nodeFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));


			u_int32_t unitigIndex;
			nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

            //cout << size << " " << unitigIndex << endl;

            addNode(nodePath, nodeDatas);

            vector<u_int32_t> nodePathRC = nodePath;
            std::reverse(nodePathRC.begin(), nodePathRC.end());
            for(size_t i=0; i<nodePathRC.size(); i++){
                nodePathRC[i] = nodeIndex_toReverseDirection(nodePathRC[i]);
            }

            addNode(nodePathRC, nodeDatas);

        }


        cout << "Nb unitigs: " << _nodes.size() << endl;


        cout << "Clearing node datas" << endl;
        nodeDatas.clear();
        
        cout << "Loading unitig graph edges successors" << endl;

        u_int64_t nbEdgeSuccessors = 0;
		ifstream edgeSuccessorFile(edgeSuccessorFilename);

		while(true){

			u_int32_t unitigIndexFrom;
			edgeSuccessorFile.read((char*)&unitigIndexFrom, sizeof(unitigIndexFrom));
			
			if(edgeSuccessorFile.eof()) break;

			u_int32_t unitigIndexTo;
			edgeSuccessorFile.read((char*)&unitigIndexTo, sizeof(unitigIndexTo));

            addSuccessor(unitigIndexFrom, unitigIndexTo);
            //cout << unitigIndexFrom << " " << unitigIndexTo << endl;
            nbEdgeSuccessors += 1;
        }

        cout << "Nb edges: " << nbEdgeSuccessors<< endl;

        /*
        cout << "Loading unitig graph edges predecessors" << endl;

        int nbEdgePredecessors = 0;

		ifstream edgePredecessorFile(edgePredecessorFilename);

		while(true){

			u_int32_t unitigIndexFrom;
			edgePredecessorFile.read((char*)&unitigIndexFrom, sizeof(unitigIndexFrom));
			
			if(edgePredecessorFile.eof()) break;

			u_int32_t unitigIndexTo;
			edgePredecessorFile.read((char*)&unitigIndexTo, sizeof(unitigIndexTo));

            addPredecessors(unitigIndexFrom, unitigIndexTo);
            nbEdgePredecessors += 1;
            //cout << unitigIndexFrom << " " << unitigIndexTo << endl;
        }

        cout << nbEdgePredecessors << endl;
        */
        cout << "Loading unitig graph sorting" << endl;

        for(size_t i=0; i<_nodes.size(); i++){
            
            UnitigGraph::Node* node = _nodes[i];

            vector<UnitigGraph::Node*> successors = node->getSuccessors(this);


            std::sort(successors.begin(), successors.end(), UnitigGraph::NodeComparator);

            for(size_t i=0; i<successors.size(); i++){
                node->_successors[i] = successors[i]->_unitigIndex;
            }
            //std::sort(node->_predecessors.begin(), node->_predecessors.end(), UnitigGraph::NodeComparator);

        }
        
        cout << "Loading unitig graph done" << endl;
    }


    void addNode(const vector<u_int32_t>& nodes, const vector<u_int32_t>& nodeDatas){

        u_int32_t unitigIndex = _currentUnitigIndex;

        Node* node = new Node(this, unitigIndex, nodes, nodeDatas);

        //for(u_int32_t nodeIndex : nodes){
        //    _nodeIndex_to_unitigIndex[nodeIndex] = unitigIndex;
        //}

        _nodes.push_back(node);
        _currentUnitigIndex += 1;
    }

    void addSuccessor(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex){
        _nodes[fromUnitigIndex]->_successors.push_back(toUnitigIndex);
        //Edge edge = new Edge(fromUnitigIndex, toUnitigIndex);
    }

    //void addPredecessors(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex){
    //    _nodes[fromUnitigIndex]->_predecessors.push_back(_nodes[toUnitigIndex]);
        //Edge edge = new Edge(fromUnitigIndex, toUnitigIndex);
    //}

    void removePredecessor(Node* fromNode, Node* toNode){
        


        /*
        for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){
            if(node->startNode() == 2713422){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->startNode() == 2713423){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->endNode() == 2713422){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->endNode() == 2713423){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
        }
        */

        //fromNode->_predecessors.erase(std::remove(fromNode->_predecessors.begin(), fromNode->_predecessors.end(), toNode), fromNode->_predecessors.end());
        toNode->_successors.erase(std::remove(toNode->_successors.begin(), toNode->_successors.end(), fromNode->_unitigIndex), toNode->_successors.end());

        Node* fromNode_rc = unitigIndex_toReverseDirection(fromNode);
        Node* toNode_rc = unitigIndex_toReverseDirection(toNode);

        fromNode_rc->_successors.erase(std::remove(fromNode_rc->_successors.begin(), fromNode_rc->_successors.end(), toNode_rc->_unitigIndex), fromNode_rc->_successors.end());
        //toNode_rc->_predecessors.erase(std::remove(toNode_rc->_predecessors.begin(), toNode_rc->_predecessors.end(), fromNode_rc), toNode_rc->_predecessors.end());
    
        /*
        if(fromNode->startNode() == 2713422 || fromNode->startNode() == 2713423 || fromNode->endNode() == 2713422 || fromNode->endNode() == 2713423){

            cout << "lala" << endl;
        }
        
        if(toNode->startNode() == 2713422 || toNode->startNode() == 2713423 || toNode->endNode() == 2713422 || toNode->endNode() == 2713423){

            cout << "lala" << endl;
        }
        if(fromNode_rc->startNode() == 2713422 || fromNode_rc->startNode() == 2713423 || fromNode_rc->endNode() == 2713422 || fromNode_rc->endNode() == 2713423){

            cout << "lala" << endl;
        }
        if(toNode_rc->startNode() == 2713422 || toNode_rc->startNode() == 2713423 || toNode_rc->endNode() == 2713422 || toNode_rc->endNode() == 2713423){

            cout << "lala" << endl;
        }
        */
        nodeChanged(fromNode);
        nodeChanged(toNode);
        //nodeChanged(fromNode_rc);
        //nodeChanged(toNode_rc);
    }

    /*
    void removePredecessor(Node* fromNode, Node* toNode, ofstream& debugFile){

        //debugFile << "RP" << fromNode->startNode() << " " << fromNode->endNode() << " " << toNode->startNode() << " " << toNode->endNode() << endl;
        fromNode->_predecessors.erase(std::remove(fromNode->_predecessors.begin(), fromNode->_predecessors.end(), toNode), fromNode->_predecessors.end());
        toNode->_successors.erase(std::remove(toNode->_successors.begin(), toNode->_successors.end(), fromNode), toNode->_successors.end());

        Node* fromNode_rc = unitigIndex_toReverseDirection(fromNode);
        Node* toNode_rc = unitigIndex_toReverseDirection(toNode);

        //debugFile << "RP" << fromNode_rc->startNode() << " " << fromNode_rc->endNode() << " " << toNode_rc->startNode() << " " << toNode_rc->endNode() << endl;

        fromNode_rc->_successors.erase(std::remove(fromNode_rc->_successors.begin(), fromNode_rc->_successors.end(), toNode_rc), fromNode_rc->_successors.end());
        toNode_rc->_predecessors.erase(std::remove(toNode_rc->_predecessors.begin(), toNode_rc->_predecessors.end(), fromNode_rc), toNode_rc->_predecessors.end());
    }
    */
    
    void nodeChanged(Node* node){

        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex].clear();


        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex].clear();
    }

    void removeNode(Node* node){



        //if(node->_unitigIndex > 10000000){
        //    cout << "rlala" << endl;
        //    getchar();
        //}
        
        for(Node* nn : node->getSuccessors(this)){
            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            nodeChanged(nn);
        }
        for(Node* nn : node->getPredecessors(this)){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node->_unitigIndex), nn->_successors.end());
            nodeChanged(nn);
        }

        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(Node* nn : node_rc->getSuccessors(this)){
            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node_rc), nn->_predecessors.end());
            nodeChanged(nn);
        }
        for(Node* nn : node_rc->getPredecessors(this)){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node_rc->_unitigIndex), nn->_successors.end());
            nodeChanged(nn);
        }


        nodeChanged(node);
        //nodeChanged(node_rc);

        node->_unitigIndex = -1;
        node_rc->_unitigIndex = -1;


        //_nodes[node->_unitigIndex] = nullptr;
        //_nodes[node_rc->_unitigIndex] = nullptr;
        //delete node;
        //delete node_rc;
        //node = nullptr;
        //node_rc = nullptr;
    }

    void recompact(Node* node, ofstream& debugFile){

        //debugFile << "recompact: " << node->startNode() << endl;

        while(true){

            bool isChanged = false;

            if(node->_successors.size() == 1){

                Node* nodeSuccessor = node->getSuccessors(this)[0];

                if(nodeSuccessor->getPredecessors(this).size() == 1){

                    if(node->getSuccessors(this)[0] != nodeSuccessor->getPredecessors(this)[0]){
                        isChanged = true;
                        mergeNode(node, nodeSuccessor);
                        //debugFile << "recompact: " << node->startNode() << " " << nodeSuccessor->startNode() << endl;
                    }

                }
            }

            /*
            if(node->_predecessors.size() == 1){

                Node* nodePredecessor = node->_predecessors[0];

                isChanged = true;
                mergeNode(node, nodeSuccessor);
                node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), nodeSuccessor), node1->_successors.end());
                node1->_successors = nodeSuccessor->_successors;

                delete successor;
            }
            */

            if(!isChanged) break;
        }

        //debugFile << "recompact: " << node->startNode() << " " << node->_nodes.size() << " " << node->_successors.size() << " " << node->_predecessors.size() << endl;

        /*
        cout << "a" << endl;
        for(Node* nn : node->_successors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_predecessors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        for(Node* nn : node->_predecessors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_successors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        cout << "b" << endl;
        */

    }

    void recompact(Node* node){

        while(true){

            bool isChanged = false;

            //if(node->_predecessors.size() == 1 && node->_predecessors[0]->_successors.size() == 1) break;

            if(node->_successors.size() == 1){

                Node* nodeSuccessor = node->getSuccessors(this)[0];

                if(nodeSuccessor->getPredecessors(this).size() == 1){

                    if(node->getSuccessors(this)[0] != nodeSuccessor->getPredecessors(this)[0]){
                        isChanged = true;
                        mergeNode(node, nodeSuccessor);
                    }

                }
            }

            /*
            if(node->_predecessors.size() == 1){

                Node* nodePredecessor = node->_predecessors[0];

                isChanged = true;
                mergeNode(node, nodeSuccessor);
                node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), nodeSuccessor), node1->_successors.end());
                node1->_successors = nodeSuccessor->_successors;

                delete successor;
            }
            */

            if(!isChanged) break;
        }

        /*
        cout << "a" << endl;
        for(Node* nn : node->_successors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_predecessors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        for(Node* nn : node->_predecessors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_successors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        cout << "b" << endl;
        */

    }

    
    void mergeNode(Node* node1, Node* node2){

        /*
        //if(node2->_unitigIndex == 236531 || node2->_unitigIndex == 236532){
        if(node1->_abundance < 100 && node2->_abundance > 1000){
            cout << "merged" << endl;
            cout << node1->_unitigIndex << " " << node1->_abundance << endl;
            cout << node2->_unitigIndex << " " << node2->_abundance << endl;
            cout << node1->_successors.size() << " " << node1->_predecessors.size() << endl;
            for(Node* nn : node1->_successors){
                cout << nn->_unitigIndex << endl;
            }
            for(Node* nn : node1->_predecessors){
                cout << nn->_unitigIndex << endl;
            }
            getchar();
        }
        */
       
        Node* node1_rc = unitigIndex_toReverseDirection(node1);
        Node* node2_rc = unitigIndex_toReverseDirection(node2);

        nodeChanged(node1);
        //nodeChanged(node1_rc);
        nodeChanged(node2);
        //nodeChanged(node2_rc);

        /*
        for(Node* nn : node1->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1_rc->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1_rc->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2_rc->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2_rc->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        */

        //cout << "\tMerging " << node1->_unitigIndex << " with " << node2->_unitigIndex << endl;
        node1->mergeWith(node2, _kminmerLength, _kminmerLengthNonOverlap);

        /*
        cout << "----" << endl;
        for(Node* successor : node1->_successors){
            cout << successor->_unitigIndex << endl;
        }
        cout << "----" << endl;
        for(Node* successor : node1_rc->_predecessors){
            cout << successor->_unitigIndex << endl;
        }
        */
        
        //node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), node2), node1->_successors.end());
        node1->_successors = node2->_successors;

        for(Node* nn : node2->getSuccessors(this)){


            //for(Node* nnP : nn->_predecessors){
            //    cout << "haha: " << nnP->_unitigIndex << endl;
            //}

            //for(size_t i=0; i<nn->_predecessors.size(); i++){
            //    if(nn->_predecessors[i] == node2){
            //        nn->_predecessors[i] = node1;
            //        break;
            //    }
            //}

            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node2), nn->_predecessors.end());
            //nn->_predecessors.push_back(node1);
            /*
            for(Node* nnP : nn->_predecessors){
                if(nnP == node2){
                    cout << "lala1" << endl;
                    nnP = node1;
                    break;
                }
            }
            */

            //for(Node* nnP : nn->_predecessors){
            //    cout << "houhou: " << nnP->_unitigIndex << endl;
            //}

        }
        
        //node1_rc->_readIndexes = node1->_readIndexes;
        node1_rc->_length = node1->_length;
        node1_rc->_abundance = node1->_abundance;
        node1_rc->_abundances = node1->_abundances;
        node1_rc->_nodes = node1->_nodes;
        std::reverse(node1_rc->_nodes.begin(), node1_rc->_nodes.end());
        for(size_t i=0; i<node1_rc->_nodes.size(); i++){
            node1_rc->_nodes[i] = nodeIndex_toReverseDirection(node1_rc->_nodes[i]);
        }


        //node1_rc->_predecessors.erase(std::remove(node1_rc->_predecessors.begin(), node1_rc->_predecessors.end(), node2_rc), node1_rc->_predecessors.end());
        //node1_rc->_predecessors = node2_rc->_predecessors;

        for(Node* nn : node2_rc->getPredecessors(this)){
            for(size_t i=0; i<nn->_successors.size(); i++){
                if(nn->_successors[i] == node2_rc->_unitigIndex){
                    nn->_successors[i] = node1_rc->_unitigIndex;
                    break;
                }
            }
        }


        //cout << node1->_unitigIndex << " " << node1_rc->_unitigIndex << " " << node2->_unitigIndex << " " << node2_rc->_unitigIndex << endl;

        //cout << "\tResult: " << node1->_unitigIndex << "  " << node2->_unitigIndex << endl;
        //cout << "\tResult: " << node1_rc->_unitigIndex << "  " << node2_rc->_unitigIndex << endl;

        //_nodes[node2->_unitigIndex] = nullptr;
        //_nodes[node2_rc->_unitigIndex] = nullptr;




        node2_rc->_unitigIndex = -1;
        node2->_unitigIndex = -1;
        node2->_nodes.clear();
        node2->_abundances.clear();
        node2_rc->_nodes.clear();
        node2_rc->_abundances.clear();



        //cout << "\tResult: " << node1->_unitigIndex << "  " << node2->_unitigIndex << endl;
        //cout << "\tResult: " << node1_rc->_unitigIndex << "  " << node2_rc->_unitigIndex << endl;
        //delete node2;
        //delete node2_rc;

        /*
        node1->_nodes.insert( node1->_nodes.end(), node2->_nodes.begin(), node2->_nodes.end());
        for(u_int32_t nodeIndex : node1._nodes){

        }
        */
        /*
        cout << "----" << endl;
        for(Node* successor : node1->_successors){
            cout << "S: " << successor->_unitigIndex << endl;
        }
        for(Node* successor : node1->_predecessors){
            cout << "P: " << successor->_unitigIndex << endl;
        }
        cout << "----" << endl;
        for(Node* successor : node1_rc->_successors){
            cout << "S: " << successor->_unitigIndex << endl;
        }
        for(Node* successor : node1_rc->_predecessors){
            cout << "P: " << successor->_unitigIndex << endl;
        }
        */

    }
    //void removeEdgeTip(Node* fromNode, Node* toNode){
        
    //}

    static u_int32_t nodeIndex_toReverseDirection(u_int32_t nodeIndex){
        bool orient;
        u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orient);

        u_int32_t nodeIndex_prev;
        if(orient){
            nodeIndex_prev = BiGraph::nodeName_to_nodeIndex(nodeName, false);
        }
        else{
            nodeIndex_prev = BiGraph::nodeName_to_nodeIndex(nodeName, true);
        }

        return nodeIndex_prev;
    }

    u_int32_t unitigIndex_toReverseDirection(u_int32_t unitigIndex) const{
        if(unitigIndex % 2 == 0){
            return unitigIndex+1;
        }
        else{
            return unitigIndex-1;
        }
    }

    Node* unitigIndex_toReverseDirection(const Node* node) const{
        if(node->_unitigIndex % 2 == 0){
            return _nodes[node->_unitigIndex+1];
        }
        else{
            return _nodes[node->_unitigIndex-1];
        }
    }

    u_int64_t computeChecksum_successor(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            //u_int64_t checksum = 0;
            //for(u_int32_t nodeIndex : node->_nodes){
            //    checksum += nodeIndex;
            //}
            //checksum *= node->_nodes.size();

            //checksumTotal += checksum;

            checksumTotal += node->startNode() + node->endNode();

            for(Node* nn : node->getSuccessors(this)){
                checksumTotal += (nn->startNode() + nn->endNode()) * node->_successors.size();
            }
            for(Node* nn : node->getPredecessors(this)){
                checksumTotal += (nn->startNode() + nn->endNode()) * node->getPredecessors(this).size();
            }
        }

        return checksumTotal;
    }

    u_int64_t computeChecksum_nodeIndex(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            u_int64_t checksum = 0;
            for(u_int32_t nodeIndex : node->_nodes){
                checksum += nodeIndex;
            }
            checksum *= node->_nodes.size();

            checksumTotal += checksum;
        }

        return checksumTotal;
    }

    u_int64_t computeChecksum_abundance(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            u_int64_t checksum = 0;
            for(u_int32_t ab : node->_abundances){
                checksum += ab;
            }
            checksum *= node->_nodes.size();

            checksumTotal += checksum;
        }

        return checksumTotal;
    }

    u_int64_t nbUnitigs(){
        u_int64_t nbUnitigs = 0;
        for(size_t i=0; i<_nodes.size(); i++){

            UnitigGraph::Node* node = _nodes[i];
            if(node->_unitigIndex == -1) continue;
            nbUnitigs += 1;
        }

        return nbUnitigs;
    }

    void save(const string& outputFilename, float referenceAbundance){

        //_logFile << "Saving unitig graph: " << outputFilename << endl;

		ofstream outputContigFile(outputFilename + ".unitigs.nodepath");
		//ofstream outputContigFileToUnitigIndex(outputFilename + ".unitigs.index");

		unordered_set<u_int32_t> writtenUnitigs;

        unordered_set<u_int32_t> selectedUnitigIndex;

        //ofstream colorFile(outputFilename + "_color.csv");
        //colorFile << "Name,Color" << endl;


        //ofstream file_nodeNameToUnitigIndex(outputFilename + ".nodeToUnitig");

        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;
            //cout << node->_unitigIndex << endl;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->startNode())) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->endNode())) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->startNode()));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->endNode()));

            selectedUnitigIndex.insert(node->_unitigIndex);

            //for(u_int32_t nodeIndex : node->_nodes){
            //    u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		    //    file_nodeNameToUnitigIndex.write((const char*)&nodeName, sizeof(nodeName));
		    //    file_nodeNameToUnitigIndex.write((const char*)&node->_unitigIndex, sizeof(node->_unitigIndex));
            //}

            /*
            if(referenceAbundance != 0){
                if(node->_abundance / referenceAbundance < 1.3){
                    colorFile << node->_unitigIndex << "," << "blue" << endl;
                }
                else{
                    colorFile << node->_unitigIndex << "," << "red" << endl;
                }
            }
            */
        }

        //file_nodeNameToUnitigIndex.close();

        
        ofstream outputFile(outputFilename);


        unordered_set<u_int32_t> isVisited;



        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;


            if(selectedUnitigIndex.find(node->_unitigIndex) == selectedUnitigIndex.end()) continue;
            if(isVisited.find(node->_unitigIndex) != isVisited.end()) continue;

            queue<u_int32_t> queue;
            queue.push(node->_unitigIndex);

            while(queue.size() > 0){
                u_int32_t unitigIndex = queue.front();
                queue.pop();

                //cout << unitigIndex << endl;
                //cout << unitigIndex << endl;
                //cout << (_nodes[unitigIndex] == nullptr) << endl;
                if(isVisited.find(unitigIndex) != isVisited.end()) continue;
                isVisited.insert(unitigIndex);

                string ori = "+";
                if(selectedUnitigIndex.find(unitigIndex) == selectedUnitigIndex.end()){
                    unitigIndex = unitigIndex_toReverseDirection(unitigIndex);
                    ori = "-";
                }

                //isVisited.insert(unitigIndex_toReverseDirection(unitigIndex));
                //linkedUnitigIndex.insert(unitigIndex);

                outputFile << "S" << "\tutg" << unitigIndex << "\t" << "*" << "\t" << "LN:i:" << _nodes[unitigIndex]->_length << "\t" << "dp:i:" << _nodes[unitigIndex]->_abundance << endl;

                UnitigGraph::Node* nodeCurrent = _nodes[unitigIndex];

                u_int8_t isCircular = nodeCurrent->isCircular(this);
				vector<u_int32_t> nodePath = nodeCurrent->_nodes;
				if(isCircular && nodePath.size() > 1){
					nodePath.push_back(nodePath[0]);
				}

				u_int32_t size = nodePath.size();
				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&nodePath[0], size * sizeof(u_int32_t));

                //outputContigFileToUnitigIndex.write((const char*)&unitigIndex, sizeof(unitigIndex));
                //vector<u_int32_t> successors;
                //getSuccessors_unitig(unitigIndex, 0, successors);

                for(const Node* nodeSuccessor : _nodes[unitigIndex]->getSuccessors(this)){
                    

                    u_int32_t unitigIndexN = nodeSuccessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigIndex << "\t" << ori << "\tutg" << unitigIndexN << "\t" << ori2 << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);

                    //if(unitigIndexN > 100000){
                        //cout << "wtf succ: " << unitigIndex << " -> " << unitigIndexN << endl;
                        //getchar();
                    //}
                }

                
                //vector<u_int32_t> predecessors;
                //getPredecessors_unitig(unitigIndex, 0, predecessors);
                //for(u_int32_t unitigIndexN : predecessors){

                for(Node* nodePredecessor : _nodes[unitigIndex]->getPredecessors(this)){


                    u_int32_t unitigIndexN = nodePredecessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigIndexN << "\t" << ori2 << "\tutg" << unitigIndex << "\t" << ori << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);


                    //if(unitigIndexN > 100000){
                    //    cout << "wtf pred: " << unitigIndex << " -> " << unitigIndexN << endl;
                    //    getchar();
                    //}
                }
                

            }
        }



        outputFile.close();
        //colorFile.close();
        outputContigFile.close();
        //outputContigFileToUnitigIndex.close();

        //_logFile << "\tdone" << endl;
    }

	/*
    void saveBin(const string& outputFilename){


        ofstream outputFile(outputFilename);

        u_int32_t nbNodes = _nodes.size();

	    outputFile.write((const char*)&nbNodes, sizeof(nbNodes));
        for(Node* node : _nodes){
            if(node->_unitigIndex == -1) continue;


	        outputFile.write((const char*)&node->_unitigIndex, sizeof(node->_unitigIndex));
	        outputFile.write((const char*)&node->_abundance, sizeof(node->_abundance));

            u_int32_t unitigNbNodes = node->_nodes.size();
	        outputFile.write((const char*)&unitigNbNodes, sizeof(unitigNbNodes));

            for(u_int32_t nodeIndex : node->_nodes){
	            outputFile.write((const char*)&nodeIndex, sizeof(nodeIndex));
            }

            u_int32_t nbSuccessors = node->_successors.size();
            outputFile.write((const char*)&nbSuccessors, sizeof(nbSuccessors));

            for(Node* node_nn : node->_successors){
	            outputFile.write((const char*)&node_nn->_unitigIndex, sizeof(node_nn->_unitigIndex));
            }


            u_int32_t nbPredecessors = node->_predecessors.size();
            outputFile.write((const char*)&nbPredecessors, sizeof(nbPredecessors));

            for(Node* node_nn : node->_predecessors){
	            outputFile.write((const char*)&node_nn->_unitigIndex, sizeof(node_nn->_unitigIndex));
            }
        }

		unordered_set<u_int32_t> writtenUnitigs;

        unordered_set<u_int32_t> selectedUnitigIndex;

        //ofstream colorFile(outputFilename + "_color.csv");
        //colorFile << "Name,Color" << endl;


        //ofstream file_nodeNameToUnitigIndex(outputFilename + ".nodeToUnitig");

        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;
            //cout << node->_unitigIndex << endl;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->startNode())) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->endNode())) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->startNode()));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->endNode()));

            selectedUnitigIndex.insert(node->_unitigIndex);

            //for(u_int32_t nodeIndex : node->_nodes){
            //    u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		    //    file_nodeNameToUnitigIndex.write((const char*)&nodeName, sizeof(nodeName));
		    //    file_nodeNameToUnitigIndex.write((const char*)&node->_unitigIndex, sizeof(node->_unitigIndex));
            //}

        }

        //file_nodeNameToUnitigIndex.close();

        
        ofstream outputFile(outputFilename);


        unordered_set<u_int32_t> isVisited;



        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;


            if(selectedUnitigIndex.find(node->_unitigIndex) == selectedUnitigIndex.end()) continue;
            if(isVisited.find(node->_unitigIndex) != isVisited.end()) continue;

            queue<u_int32_t> queue;
            queue.push(node->_unitigIndex);

            while(queue.size() > 0){
                u_int32_t unitigIndex = queue.front();
                queue.pop();

                //cout << unitigIndex << endl;
                //cout << unitigIndex << endl;
                //cout << (_nodes[unitigIndex] == nullptr) << endl;
                if(isVisited.find(unitigIndex) != isVisited.end()) continue;
                isVisited.insert(unitigIndex);

                string ori = "+";
                if(selectedUnitigIndex.find(unitigIndex) == selectedUnitigIndex.end()){
                    unitigIndex = unitigIndex_toReverseDirection(unitigIndex);
                    ori = "-";
                }

                //isVisited.insert(unitigIndex_toReverseDirection(unitigIndex));
                //linkedUnitigIndex.insert(unitigIndex);

                outputFile << "S" << "\tutg" << unitigIndex << "\t" << "*" << "\t" << "LN:i:" << _nodes[unitigIndex]->_length << "\t" << "dp:i:" << _nodes[unitigIndex]->_abundance << endl;

                UnitigGraph::Node* nodeCurrent = _nodes[unitigIndex];

                u_int8_t isCircular = nodeCurrent->isCircular();
				vector<u_int32_t> nodePath = nodeCurrent->_nodes;
				if(isCircular && nodePath.size() > 1){
					nodePath.push_back(nodePath[0]);
				}

				u_int32_t size = nodePath.size();
				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&nodePath[0], size * sizeof(u_int32_t));

                //outputContigFileToUnitigIndex.write((const char*)&unitigIndex, sizeof(unitigIndex));
                //vector<u_int32_t> successors;
                //getSuccessors_unitig(unitigIndex, 0, successors);

                for(Node* nodeSuccessor : _nodes[unitigIndex]->_successors){
                    

                    u_int32_t unitigIndexN = nodeSuccessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigIndex << "\t" << ori << "\tutg" << unitigIndexN << "\t" << ori2 << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);

                    //if(unitigIndexN > 100000){
                        //cout << "wtf succ: " << unitigIndex << " -> " << unitigIndexN << endl;
                        //getchar();
                    //}
                }

                
                //vector<u_int32_t> predecessors;
                //getPredecessors_unitig(unitigIndex, 0, predecessors);
                //for(u_int32_t unitigIndexN : predecessors){

                for(Node* nodePredecessor : _nodes[unitigIndex]->_predecessors){


                    u_int32_t unitigIndexN = nodePredecessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigIndexN << "\t" << ori2 << "\tutg" << unitigIndex << "\t" << ori << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);


                    //if(unitigIndexN > 100000){
                    //    cout << "wtf pred: " << unitigIndex << " -> " << unitigIndexN << endl;
                    //    getchar();
                    //}
                }
                

            }
        }



        outputFile.close();
        //colorFile.close();
        outputContigFile.close();
        //outputContigFileToUnitigIndex.close();

        //_logFile << "\tdone" << endl;
    }
    */

};

/*
class Graph2{

public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;

    vector<vector<AdjNode>> _successors;
    vector<vector<AdjNode>> _predecessors;
    
    Graph2(){
        //_nbNodes = nbNodes;
        //_successors.resize(_nbNodes);
        //_predecessors.resize(_nbNodes);

        _successors.clear();
        _predecessors.clear();
        _nbNodes = 0;
        _nbEdges = 0;
    }

    unordered_map<u_int32_t, u_int32_t> _nodeName_to_nodeIndex;
    vector<u_int32_t> _nodeIndex_to_nodeName;


    //vector<u_int64_t> _unitigs_length;

    void addNode(u_int32_t nodeName){
        if (_nodeName_to_nodeIndex.find(nodeName) != _nodeName_to_nodeIndex.end()) return;

        _nodeName_to_nodeIndex[nodeName] = _nodeIndex_to_nodeName.size();
        _nodeIndex_to_nodeName.push_back(nodeName);
    }

    u_int32_t nodeName_to_nodeIndex(u_int32_t nodeName){
        return _nodeName_to_nodeIndex[nodeName];
    }

    u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex){
        return _nodeIndex_to_nodeName[nodeIndex];
    }

    void addEdge(u_int32_t from, u_int32_t to){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to);
        _successors[nodeIndex_from].push_back({nodeIndex_to, 0}); //= getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _predecessors[nodeIndex_to].push_back({nodeIndex_from, 0}); //= getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _nbEdges += 1;
    }

    void removeEdge(u_int32_t from, u_int32_t to){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to);

        vector<AdjNode>& successors = _successors[nodeIndex_from];
        
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                successors.erase(successors.begin() + i);
                break;
            }
        }
    
        vector<AdjNode>& predecessors = _predecessors[nodeIndex_to];
        
        for (size_t i=0; i < predecessors.size(); i++) {
            if (predecessors[i]._index == nodeIndex_from) {
                predecessors.erase(predecessors.begin() + i);
                break;
            }
        }

    }
    
    u_int32_t in_degree(u_int32_t v){
        return _predecessors[nodeName_to_nodeIndex(v)].size();
    }
    
    u_int32_t out_degree(u_int32_t v){
        return _successors[nodeName_to_nodeIndex(v)].size();
    }

    void getSuccessors(u_int32_t v, vector<u_int32_t>& successors){
        successors.clear();
        for(AdjNode& node : _successors[nodeName_to_nodeIndex(v)]){
            successors.push_back(nodeIndex_to_nodeName(node._index));
        }
    }
    
    void getPredecessors(u_int32_t v, vector<u_int32_t>& predecessors){
        predecessors.clear();
        for(AdjNode& node : _predecessors[nodeName_to_nodeIndex(v)]){
            predecessors.push_back(nodeIndex_to_nodeName(node._index));
        }
    }

};

*/
/*
class UnitigGraph{
    
    adjNode* getAdjListNode(u_int32_t to, float weight, adjNode* head)   {
        adjNode* newNode = new adjNode;
        newNode->val = to;
        newNode->weight = weight;

        newNode->next = head;
        return newNode;
    }

public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;
    //u_int32_t _nodeIndex;

    vector<adjNode*> _nodes;
    //vector<bool> isVisited;
    //vector<u_int32_t> distance;
    //vector<u_int32_t> prev;

    UnitigGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes;
        _nodes.resize(_nbNodes, nullptr);

        //prev.resize(_nbNodes, 0);
        //isVisited.resize(_nbNodes, 0);
        //distance.resize(_nbNodes, 0);

        _nbEdges = 0;
    }

    ~UnitigGraph() {
        for (int i = 0; i < _nbNodes; i++){
            delete[] _nodes[i];
        }
    }

    bool addEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient, float weight){

        //if(from == 798 || to == 798){
        //  cout << unitigIndex_to_unitigName(from) << " " << fromOrient << " " << unitigIndex_to_unitigName(to) << " " << toOrient << endl;
        //}
        if((fromOrient && toOrient)){
            //if(from == 798 || to == 798){
            //    cout << "1:     " << unitigIndex_to_unitigName(from) << " -> " << unitigIndex_to_unitigName(to) << endl;
            //}
            _nodes[from] = getAdjListNode(to, weight, _nodes[from]);
            _nbEdges += 1;
            return true; 
        }
        else if((fromOrient && !toOrient)){
            //if(from == 798 || to == 798){
            //    cout << "2:    " << unitigIndex_to_unitigName(to) << " -> " << unitigIndex_to_unitigName(from) << endl;
            //}
            _nodes[to] = getAdjListNode(from, weight, _nodes[to]);
            _nbEdges += 1;
            return true; 
        }


        return false;
    }


};
*/



#endif

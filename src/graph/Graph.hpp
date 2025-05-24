


#ifndef MDBG_METAG_GRAPH
#define MDBG_METAG_GRAPH

#include "../Commons.hpp"

//#include "../graph/GfaParser.hpp"



class BiGraph{

public:

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

};


/*

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
*/



class UnitigGraph2{

public:

    struct UnitigNode{

    public:

        UnitigType _unitigName;
        vector<UnitigType> _unitigMerge;
        //vector<UnitigType> _unitigs;
        vector<UnitigType> _successors_forward;
        vector<UnitigType> _successors_reverse;
        u_int32_t _nbMinimizers;
        float _abundance;
        vector<float> _abundances;
        bool _isReversed;
        //u_int32_t _startNode_forward;
        //u_int32_t _startNode_reverse;

        vector<MinimizerType> _nodes;

        UnitigNode(){

        }

        UnitigNode(const UnitigType& unitigName, const vector<MinimizerType>& minimizers){
            _unitigName = unitigName;
            _nbMinimizers = minimizers.size();
            //_nodes = minimizers;
            _isReversed = false;
        }

        u_int32_t startNode(){
            return _unitigName;
            /*
            if(_isReversed){
                return _startNode_reverse;
            }
            else{
                return _startNode_forward;
            }
            */
        }

        u_int64_t getLength(UnitigGraph2* unitigGraph) const{
            return unitigGraph->_kminmerLength + ((_nbMinimizers - unitigGraph->_kminmerSize + 1 - 1)* unitigGraph->_kminmerLengthNonOverlap);
        }

        double computeMedianAbundance(){
            
            long double sum = 0;
            long double n = 0;

            for(const float& ab : _abundances){
                sum += ab;
                n += 1;
            }

            return sum / n;
            /*
            size_t size = _abundances.size();

            if (size == 0){
                return 0;  // Undefined, really.
            }
            else{
                //sort(_abundances.begin(), _abundances.end());
                if (size % 2 == 0){
                    return (_abundances[size / 2 - 1] + _abundances[size / 2]) / 2.0; //scores[size / 2 - 1];
                }
                else {
                    return _abundances[size / 2];
                }
            }
            */
        }

        UnitigType startNode2() const{
            return _unitigName;
        }

        void mergeWith(UnitigNode* node2, UnitigGraph2* graph){

            //cout << "Merge: " << _unitigName*2 << " " << node2->_unitigName*2 << " " << _abundance << endl;
            //for(float ab : _abundances){
            //    cout << "\t" << ab << endl;
            //}
            //cout << endl;
            //for(float ab : node2->_abundances){
            //    cout << "\t" << ab << endl;
            //}
            //cout << endl;

            //for(float abundance : node2->_abundances){
            //    _abundances.push_back(abundance);
            //}
            //mergeAbundances(node2);
            vector<float> abundances;
            std::merge(_abundances.begin(), _abundances.end(), node2->_abundances.begin(), node2->_abundances.end(), std::back_inserter(abundances));

            _abundances = abundances;

            _abundance = computeMedianAbundance();
        
            _nbMinimizers += (node2->_nbMinimizers - graph->_kminmerSize + 1);
            

            //cout << endl;
            //for(float ab : _abundances){
            //    cout << "\t" << ab << endl;
            //}
            
            //_nodes.insert( _nodes.end(), node2->_nodes.begin()+graph->_kminmerSize-1, node2->_nodes.end());

            //_length = computeUnitigLength(_nodes.size(), graph->_kminmerSize, graph->_kminmerLength, graph->_kminmerLengthNonOverlap);


        }

        /*
        void mergeAbundances(UnitigNode* node2){


            const vector<float> abundance1 = _abundances;
            const vector<float>& abundance2 = node2->_abundances;

            _abundances.clear();
            
            size_t i=0;
            size_t j=0;
            cout << abundance1.size() << " " << abundance2.size() << endl;

            while(true){

                if(i == abundance1.size() && j == abundance2.size()) break;

                cout << i << " " << j << " " << abundance1[i] << " " << abundance2[j] << endl;

                if(abundance1[i] == abundance2[j]){
                    _abundances.push_back(abundance1[i]);
                    _abundances.push_back(abundance2[j]);
                    i += 1;
                    j += 1;
                }
                else if(abundance1[i] < abundance2[j]){
                    _abundances.push_back(abundance1[i]);
                    i += 1;
                }
                else{
                    _abundances.push_back(abundance2[j]);
                    j += 1;
                }

                getchar();
            }

            cout << _abundances.size() << endl;
            //_abundances = abundances;
            //return nbShared;

        }
        */
    };


    vector<UnitigNode*> _unitigs;
    //vector<u_int32_t> _nodeIndex_to_unitigIndex;
    size_t _kminmerSize;
	float _kminmerLength;
	float _kminmerOverlapMean;
	float _kminmerLengthNonOverlap;
	unordered_map<KmerVec, u_int32_t>& _kmerVec_to_nodeName;
    int _nbCores;
    //const vector<u_int32_t>& _nodeDatas;

    //class GraphChangedIndex{
    //    public:
    //    unordered_set<u_int32_t> _superBubbleIndex2_notSuperbubble;
    //    unordered_map<u_int32_t, vector<u_int32_t>> _superBubbleIndex2_notSuperbubbleNodes;
    //};

    //GraphChangedIndex _graphChangedIndex;

    UnitigGraph2(size_t kminmerSize, float kminmerLength, float kminmerOverlapMean, float kminmerLengthNonOverlap, unordered_map<KmerVec, u_int32_t>& kmerVec_to_nodeName, const int nbCores) : _kmerVec_to_nodeName(kmerVec_to_nodeName){//}, const vector<u_int32_t>& nodeDatas) : _nodeDatas(nodeDatas){
        //_currentUnitigIndex = 0;
        _kminmerSize = kminmerSize;
        _kminmerLength = kminmerLength;
        _kminmerOverlapMean = kminmerOverlapMean;
        _kminmerLengthNonOverlap = kminmerLengthNonOverlap;
        _nbCores = nbCores;
        //_nodeIndex_to_unitigIndex.resize(nbNodes, -1);
    }

    ~UnitigGraph2(){
        for(size_t i=0; i<_unitigs.size(); i++){
            delete _unitigs[i];
        }
    }


    void load(const string& inputDir){

        //const string& nodeFilename, const string& abundanceFilename, const string& edgeSuccessorFilename
        //    + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph.nodes.abundances.bin", _inputDir + "/unitigGraph.edges.successors.bin"

        //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala2.txt");

        Logger::get().debug() << "Load unitigs";

        u_int64_t nbUnitigs = 0;
        u_int64_t nbMinimizers = 0;

		ifstream nodeFile(inputDir + "/unitigGraph.nodes.bin");

		while(true){

			u_int32_t size;
			nodeFile.read((char*)&size, sizeof(size));
			

			if(nodeFile.eof()) break;

			vector<MinimizerType> minimizers;
			minimizers.resize(size);
			nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));


			UnitigType unitigIndex;
			nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

            //u_int32_t nodeName;
			//nodeFile.read((char*)&nodeName, sizeof(nodeName));



            //cout << "Add node: " << unitigIndex << " " << minimizers.size() << endl;
            addNode(unitigIndex/2, minimizers);

            //vector<MinimizerType> minimizersRC = minimizers;
            //std::reverse(minimizersRC.begin(), minimizersRC.end());
            
            //addNode(unitigIndex+1, minimizersRC);

            //debugFile << unitigIndex << " " << minimizers.size() << endl;
            nbUnitigs += 1;
            nbMinimizers += minimizers.size();
        }

        nodeFile.close();
        //return;

        //debugFile.close();
        Logger::get().debug() << "Nb unitigs: " << nbUnitigs;
        Logger::get().debug() << "Nb minimizers: " << nbMinimizers;
        //getchar();

        /*
        cout << "Add start end node" << endl;
		ifstream startEndFile(inputDir + "/unitigGraph.startEnd.bin");

		while(true){


			UnitigType unitigIndex;
			startEndFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			u_int32_t startNode;
			startEndFile.read((char*)&startNode, sizeof(startNode));
			

			if(startEndFile.eof()) break;

            bool isReversed = false;
            const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

            if(isReversed){
                _unitigs[unitigName]->_startNode_reverse = startNode;
            }
            else{
                _unitigs[unitigName]->_startNode_forward = startNode;
            }
            
            //_unitigs[unitigName]->_abundance = _unitigs[unitigName]->computeMedianAbundance();

            //u_int32_t unitigIndexRC = unitigIndex_toReverseDirection(unitigIndex);
            //_nodes[unitigIndexRC]->_abundances = abundances;
            //_nodes[unitigIndexRC]->_abundance = _nodes[unitigIndexRC]->computeMedianAbundance();

        }

        startEndFile.close();
        */



        Logger::get().debug() << "Add unitig abundances";


		ifstream abundanceFile(inputDir + "/unitigGraph.nodes.abundances.bin");

		while(true){


			u_int32_t unitigIndex;
			abundanceFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			u_int32_t nbAbundances;
			abundanceFile.read((char*)&nbAbundances, sizeof(nbAbundances));
			

			if(abundanceFile.eof()) break;

			vector<float> abundances;
			abundances.resize(nbAbundances);
			abundanceFile.read((char*)&abundances[0], nbAbundances * sizeof(float));

            const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex);

            _unitigs[unitigName]->_abundances = abundances;
            //_unitigs[unitigName]->_abundance = _unitigs[unitigName]->computeMedianAbundance();

            //u_int32_t unitigIndexRC = unitigIndex_toReverseDirection(unitigIndex);
            //_nodes[unitigIndexRC]->_abundances = abundances;
            //_nodes[unitigIndexRC]->_abundance = _nodes[unitigIndexRC]->computeMedianAbundance();

        }

        abundanceFile.close();
        

        Logger::get().debug() << "Sorting abundances";
        processUnitigsParallel(AbundanceSortFunctor(this), _nbCores);

        //cout << "Clearing node datas" << endl;
        
        Logger::get().debug() << "Add unitig edges";
        //cout << "Loading unitig graph edges successors" << endl;

        u_int64_t nbEdges = 0;
		ifstream edgeSuccessorFile(inputDir + "/unitigGraph.edges.successors.bin");

		while(true){

			u_int32_t unitigIndexFrom;
			edgeSuccessorFile.read((char*)&unitigIndexFrom, sizeof(unitigIndexFrom));
			
			if(edgeSuccessorFile.eof()) break;

			u_int32_t unitigIndexTo;
			edgeSuccessorFile.read((char*)&unitigIndexTo, sizeof(unitigIndexTo));

            addSuccessor(unitigIndexFrom, unitigIndexTo);
            //cout << unitigIndexFrom << " " << unitigIndexTo << endl;
            nbEdges += 1;
        }

        edgeSuccessorFile.close();

        Logger::get().debug() << "Nb unitig edges: " << nbEdges;
       
        //cout << "Loading unitig graph sorting" << endl;
        
        Logger::get().debug() << "Sorting unitig edges";

        /*
        vector<UnitigType> successors;
        getSuccessors(50091, successors);


        vector<UnitigType> preds;
        getPredecessors(50091, preds);

        cout << "Succs:" << endl;
        for(UnitigType unitigIndex : successors){
            cout << "\tutg" << unitigIndex << " " << _unitigs[unitigIndex_to_unitigName(unitigIndex)]->_startNode_forward << endl;
        }
        cout << "Preds:" << endl;
        for(UnitigType unitigIndex : preds){
            cout << "\tutg" << unitigIndex << " " << _unitigs[unitigIndex_to_unitigName(unitigIndex)]->startNode() << endl;
        }
        */

        //vector<SuccessorComparatorData> sorters;

        for(UnitigNode* unitig : _unitigs){
            
            /*
            if(unitig->_unitigName == 25045){
                cout << unitig->startNode() << " " << unitig->_startNode_forward << " " << unitig->_startNode_reverse << endl;

                cout << "Succs:" << endl;
                for(const UnitigType& unitigIndex : unitig->_successors_forward){
                    cout << "\tutg" << unitigIndex << " " << _unitigs[unitigIndex_to_unitigName(unitigIndex)]->startNode() << endl;
                }

                cout << "Preds:" << endl;

                for(const UnitigType& unitigIndex : unitig->_successors_reverse){
                    cout << "\tutg" << unitigIndex << " " << _unitigs[unitigIndex_to_unitigName(unitigIndex)]->_startNode_forward << " " << _unitigs[unitigIndex_to_unitigName(unitigIndex)]->_startNode_reverse << endl;
                }

            }
            */
            /*
            sorters.clear();

            for(const UnitigType& unitigIndex : unitig->_successors_forward){

                sorters.push_back({unitigIndex, getStartNode(unitigIndex)});
            }
            //vector<UnitigGraph::Node*> successors = node->_successors;


            std::sort(sorters.begin(), sorters.end(), UnitigGraph2::SuccessorComparator);

            unitig->_successors_forward.clear();

            for(const SuccessorComparatorData& data : sorters){
                unitig->_successors_forward.push_back(data._unitigName);
            }
            */
            //for(size_t i=0; i<successors.size(); i++){
            //    node->_successors[i] = successors[i]->_unitigIndex;
            //}
            std::sort(unitig->_successors_forward.begin(), unitig->_successors_forward.end(), UnitigGraph2::NodeComparator);
            std::sort(unitig->_successors_reverse.begin(), unitig->_successors_reverse.end(), UnitigGraph2::NodeComparator);

            /*
            if(unitig->_unitigName == 50091){

                cout << "succs:" << endl;
                for(const SuccessorComparatorData& data : sorters){
                    cout << "\tutg" << data._unitigName/2 << " " << data._nodeName << endl;
                }

            }
            */


            /*
            sorters.clear();

            for(const UnitigType& unitigIndex : unitig->_successors_reverse){
                sorters.push_back({unitigIndex, getStartNode(unitigIndex)});
            }
            //vector<UnitigGraph::Node*> successors = node->_successors;


            std::sort(sorters.begin(), sorters.end(), UnitigGraph2::SuccessorComparator);

            unitig->_successors_reverse.clear();

            for(const SuccessorComparatorData& data : sorters){
                unitig->_successors_reverse.push_back(data._unitigName);
            }
            */
            /*
            if(unitig->_unitigName == 50091){

                cout << "preds:" << endl;
                for(const SuccessorComparatorData& data : sorters){
                    cout << "\tutg" << data._unitigName/2 << " " << data._nodeName << endl;
                }

                //getchar();

            }
            */
        }

        
        

        Logger::get().debug() << "Done loading unitig graph";
    }

    /*
    struct SuccessorComparatorData{
        UnitigType _unitigName;
        u_int32_t _nodeName;
    };

	static bool SuccessorComparator(const SuccessorComparatorData& unitigName1, const SuccessorComparatorData& unitigName2){
        return unitigName1._nodeName < unitigName2._nodeName;
	}
    */

	static bool NodeComparator(const UnitigType& unitigName1, const UnitigType& unitigName2){
        return unitigName1 < unitigName2;
	}

	static bool NodeComparatorByLengthRev(UnitigGraph2::UnitigNode*a, UnitigGraph2::UnitigNode*b){

        if(a->_nbMinimizers == b->_nbMinimizers){
            return a->startNode2() > b->startNode2();
        }
        return a->_nbMinimizers < b->_nbMinimizers;
        
	}
    
    unordered_map<UnitigType, vector<MinimizerType>> _unitigName_to_minimizers;

    void addNode(const UnitigType unitigIndex, const vector<MinimizerType>& minimizers){

        while(_unitigs.size() <= unitigIndex){
            _unitigs.push_back(nullptr);
        }
        //cout << nodes.size() << " " << _kminmerSize << endl;
        //vector<u_int32_t> nodesDummy(nodes.size()-_kminmerSize+1, 0);

        //u_int32_t unitigIndex = _currentUnitigIndex;

        UnitigNode* unitig = new UnitigNode(unitigIndex, minimizers);

        //cout << _nodes.size() << " " << unitigIndex << endl;
        _unitigs[unitigIndex] = unitig;
        //for(u_int32_t nodeIndex : nodes){
        //    _nodeIndex_to_unitigIndex[nodeIndex] = unitigIndex;
        //}
        //cout << "\tAdd node: " << _currentUnitigIndex << endl;
        //_nodes.push_back(node);
        //_currentUnitigIndex += 1;

        //_unitigName_to_minimizers[unitigIndex] = minimizers;

    }

    void addSuccessor(UnitigType fromUnitigIndex, UnitigType toUnitigIndex){


        //if(fromUnitigIndex == 1916 || toUnitigIndex == 1916 || fromUnitigIndex == 1917 || toUnitigIndex == 1917){
        //    cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
        //}

        //cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;

        bool isReversedFrom;
        UnitigType fromUnitigName = unitigIndex_to_unitigName(fromUnitigIndex, isReversedFrom);

        bool isReversedTo;
        UnitigType toUnitigName = unitigIndex_to_unitigName(toUnitigIndex, isReversedTo);

        //if(fromUnitigName == 17 || toUnitigName == 17){
        //    cout << fromUnitigIndex << " -> " << toUnitigIndex << endl;
        //    getchar();
        //}


        if(isReversedFrom){
            //cout << "\tAdd edge reverse: " << fromUnitigName << " " << toUnitigIndex << endl;
            _unitigs[fromUnitigName]->_successors_reverse.push_back(toUnitigIndex);
            //_unitigs[fromUnitigName]
        }
        else{
            //cout << "\tAdd edge forward: " << fromUnitigName << " " << toUnitigIndex << endl;
            _unitigs[fromUnitigName]->_successors_forward.push_back(toUnitigIndex);
        }

        //cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
        //_nodes[fromUnitigIndex]->_successors.push_back(_nodes[toUnitigIndex]);
        //_nodes[toUnitigIndex]->_predecessors.push_back(_nodes[fromUnitigIndex]);

    }

    static UnitigType unitigName_to_unitigIndex(const UnitigType& unitigName, const bool& isReversed){
        if(isReversed){
            return unitigName * 2 + 1;
        }
        else{
            return unitigName * 2;
        }
    }

    static UnitigType unitigIndex_to_unitigName(const UnitigType& unitigIndex, bool& isReversed){
        isReversed = (unitigIndex % 2) == 1;
        return unitigIndex / 2;
    }

    static UnitigType unitigIndex_to_unitigName(const UnitigType& unitigIndex){
        //isOrientationForward = (nodeIndex % 2) == 0;
        return unitigIndex / 2;
    }

    static UnitigType unitigIndex_to_reverseDirection(const UnitigType& unitigIndex){
        if(unitigIndex % 2 == 0){
            return unitigIndex+1;
        }
        else{
            return unitigIndex-1;
        }
    }

    static string unitigIndexToString(const UnitigType& unitigIndex){

        bool isReversed;
        const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

        if(isReversed){
            return to_string(unitigName) + "-";
        }
        else{
            return to_string(unitigName) + "+";
        }
    }

    UnitigType getStartNode(const UnitigType& unitigIndex){

        return unitigIndex;
        /*
        bool isReversed;
        const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

        UnitigGraph2::UnitigNode* unitig = _unitigs[unitigName];

        if(isReversed){
            return unitig->_startNode_reverse;
        }
        else{
            return unitig->_startNode_forward;
        }
        */
    }

    /*
    Node* unitigIndex_toReverseDirection(const Node* node) const{
        if(node->_unitigIndex % 2 == 0){
            return _nodes[node->_unitigIndex+1];
        }
        else{
            return _nodes[node->_unitigIndex-1];
        }
    }
    */
    
    void nodeChanged(UnitigNode* node){
        /*
        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex].clear();


        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex].clear();
        */
    }



    void removeNode(UnitigNode* node){

        const UnitigType unitigName = node->_unitigName;

        //cout << "Removing node: utg" << node->_unitigName << endl;



        //for(const UnitigType& unitigIndexSuccessor : node->_successors_forward){
        //    cout << "\tSuccessor forward: utg" << unitigIndex_to_unitigName(unitigIndexSuccessor) << endl;
        //}
        //for(const UnitigType& unitigIndexSuccessor : node->_successors_reverse){
        //    cout << "\tSuccessor reverse: utg" << unitigIndex_to_unitigName(unitigIndexSuccessor) << endl;
        //}

        removeEdges(node, false);
        removeEdges(node, true);
        nodeChanged(node);

        //if(unitigName == 778){
        //    cout << "plop" << endl;
        //    exit(1);
        //}

        delete node;
        _unitigs[unitigName] = nullptr;
        

        /*
        for(UnitigNode* nn : node->_successors){
            nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            nodeChanged(nn);
        }
        for(UnitigNode* nn : node->_predecessors){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node), nn->_successors.end());
            nodeChanged(nn);
        }

        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(UnitigNode* nn : node_rc->_successors){
            nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node_rc), nn->_predecessors.end());
            //nodeChanged(nn);
        }
        for(UnitigNode* nn : node_rc->_predecessors){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node_rc), nn->_successors.end());
            //nodeChanged(nn);
        }
        */


        //delete node;
        //delete node_rc;
        //nodeChanged(node_rc);

        //node->clear();
        //node_rc->clear();

    }

    void removeEdges(UnitigNode* node, bool isReversed){

        //cout << "\tRemove edges: " << isReversed << endl;
        UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, isReversed);
        

        vector<UnitigType> successors;
        getSuccessors(unitigIndex, successors);
        
        u_int32_t unitigIndexToRemove = unitigIndex_to_reverseDirection(unitigIndex);
        //vector<UnitigType> predecessors;
        //getPredecessors(unitigIndex, predecessors);

        for(const UnitigType& unitigIndexSuccessor : successors){


            bool isReversedSuccessor;
            const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);

            //cout << "\tEdge: " << node->_unitigName << "-" << isReversed << " -> " << unitigNameSuccessor << "-" << isReversedSuccessor << endl;

            if(isReversedSuccessor){

                //if(_unitigs[unitigNameSuccessor] == nullptr){
                //    cout << "omg" << endl;
                //    getchar();
                //}

                vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

                int nbSuccessors = successorsFoward.size();
                //cout << "\t\tNb successors: " << successorsFoward.size() << endl;
                successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexToRemove), successorsFoward.end());
                //cout << "\t\tNb successors: " << successorsFoward.size() << endl;

                //if(nbSuccessors == successorsFoward.size()){
                //    cout << "issue" << endl;
                //    getchar();
                //}

            }
            else{


                //if(_unitigs[unitigNameSuccessor] == nullptr){
                //    cout << "omg" << endl;
                //    getchar();
                //}

                vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

                int nbSuccessors = successorsReverse.size();
                //cout << "\t\tNb successors: " << successorsReverse.size() << endl;
                successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexToRemove), successorsReverse.end());
                //cout << "\t\tNb successors: " << successorsReverse.size() << endl;

                //if(nbSuccessors == successorsReverse.size()){
                //    cout << "issue" << endl;
                //    getchar();
                //}


            }
            /*
            

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_forward){

                cout << "\t\tEdge neighbor forward: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\t\tFound forward" << endl;
                }
            }

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_reverse){

                cout << "\t\tEdge neighbor reverse: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\tFound reverse" << endl;
                }
            }
            */
            //const UnitigType& 
            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            //nodeChanged(nn);
        }
        //for(UnitigNode* nn : node->_predecessors){
        //    nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node), nn->_successors.end());
        //    nodeChanged(nn);
        //}


    }


    void recompact(UnitigNode* node){

        //cout << "Recompact: " << node->_unitigName << endl;

        recompact(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true));
        recompact(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));
        /*
        while(true){

            bool isChanged = false;

            //if(node->_predecessors.size() == 1 && node->_predecessors[0]->_successors.size() == 1) break;

            if(node->_successors.size() == 1){

                Node* nodeSuccessor = node->_successors[0];

                if(nodeSuccessor->_predecessors.size() == 1){

                    if(node->_successors[0] != nodeSuccessor->_predecessors[0]){
                        isChanged = true;
                        mergeNode(node, nodeSuccessor);
                    }

                }
            }



            if(!isChanged) break;
        }
        */
    }

    void recompact(const UnitigType& unitigIndex){



        //cout << "\tRecompact: " << isReversed << endl;
        

        while(true){

            bool isChanged = false;

            vector<UnitigType> successors;
            getSuccessors(unitigIndex, successors);
            
            
            //cout << "\t\tNb successors: " << successors.size() << endl;

            if(successors.size() == 1){

                UnitigNode* nodeSuccessor = _unitigs[successors[0]];
                //if(nodeSuccessor == nullptr){
                //    cout << "wtf" << endl;
                //    getchar();
                //    continue;
                //}

                //if(node != nullptr && node->_unitigName == 3153 && nodeSuccessor != nullptr){
                    //cout << "\t\t" << nodeSuccessor->_successors_forward.size() << " " << nodeSuccessor->_successors_reverse.size() << endl;
                //}

                vector<UnitigType> predecessors;
                getPredecessors(successors[0], predecessors);

                //cout << "\t\tNb predecessors: " << predecessors.size() << endl;

                if(predecessors.size() == 1){

                    if(successors[0] != predecessors[0]){ //Do not compact a circular unitig
                        isChanged = true;
                        mergeNode(unitigIndex, successors[0]);
                    }
                }


            }

            if(!isChanged) break;
        }

        /*
        for(const UnitigType& unitigIndexSuccessor : successors){

            bool isReversedSuccessor;
            const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);

            cout << "\t\tEdge: " << node->_unitigName << "-" << isReversed << " -> " << unitigNameSuccessor << "-" << isReversedSuccessor << endl;
        }
        u_int32_t unitigIndexToRemove = unitigIndex_to_reverseDirection(unitigIndex);
        //vector<UnitigType> predecessors;
        //getPredecessors(unitigIndex, predecessors);

        for(const UnitigType& unitigIndexSuccessor : successors){


            bool isReversedSuccessor;
            const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);

            cout << "\tEdge: " << node->_unitigName << "-" << isReversed << " -> " << unitigNameSuccessor << "-" << isReversedSuccessor << endl;

            if(isReversedSuccessor){

                vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

                cout << "\t\tNb successors: " << successorsFoward.size() << endl;
                successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexToRemove), successorsFoward.end());
                cout << "\t\tNb successors: " << successorsFoward.size() << endl;


            }
            else{

                vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

                cout << "\t\tNb successors: " << successorsReverse.size() << endl;
                successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexToRemove), successorsReverse.end());
                cout << "\t\tNb successors: " << successorsReverse.size() << endl;



            }
            */
            /*
            

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_forward){

                cout << "\t\tEdge neighbor forward: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\t\tFound forward" << endl;
                }
            }

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_reverse){

                cout << "\t\tEdge neighbor reverse: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\tFound reverse" << endl;
                }
            }
            //const UnitigType& 
            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            //nodeChanged(nn);
        }
        //for(UnitigNode* nn : node->_predecessors){
        //    nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node), nn->_successors.end());
        //    nodeChanged(nn);
        //}

        */

    }

    void removePredecessors(const UnitigType& unitigIndex, ofstream& debugFile){
        
        //cout << "Remove predecessors" << endl;

        const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex);

        vector<UnitigType> predecessors;
        getPredecessors(unitigIndex, predecessors);
        
        u_int32_t unitigIndexToRemove = unitigIndex; //unitigIndex_to_reverseDirection(unitigIndex);
        //vector<UnitigType> predecessors;
        //getPredecessors(unitigIndex, predecessors);

        for(const UnitigType& unitigIndexSuccessor : predecessors){

            debugFile << "Tip removePred: " << unitigIndex << " " << unitigIndexSuccessor << endl;

            //cout << "\t" << unitigIndexToString(unitigIndexSuccessor) << endl;

            bool isReversedSuccessor;
            const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);

            if(isReversedSuccessor){


                vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

                //for(UnitigType unitigIndex : successorsReverse){
                //    cout << "\t\t" << unitigIndex << endl;
                //}

                //int nbSuccessors = successorsReverse.size();
                //cout << successorsReverse.size() << endl;
                successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexToRemove), successorsReverse.end());
                //cout << successorsReverse.size() << endl;
                


            }
            else{


                vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

                //for(UnitigType unitigName : successorsFoward){
                //    cout << "\t\t" << unitigName << endl;
                //}

                //int nbSuccessors = successorsFoward.size();
                //cout << successorsFoward.size() << endl;
                successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexToRemove), successorsFoward.end());
                //cout << successorsFoward.size() << endl;


            }

            debugFile << "Tip recompacting: " << unitigIndexSuccessor << endl;
            recompact(unitigIndexSuccessor);



        }


        _unitigs[unitigName]->_successors_forward.clear();
        _unitigs[unitigName]->_successors_reverse.clear();
        

        //fromNode->_predecessors.erase(std::remove(fromNode->_predecessors.begin(), fromNode->_predecessors.end(), toNode), fromNode->_predecessors.end());
        //toNode->_successors.erase(std::remove(toNode->_successors.begin(), toNode->_successors.end(), fromNode), toNode->_successors.end());

        //Node* fromNode_rc = unitigIndex_toReverseDirection(fromNode);
        //Node* toNode_rc = unitigIndex_toReverseDirection(toNode);

        //fromNode_rc->_successors.erase(std::remove(fromNode_rc->_successors.begin(), fromNode_rc->_successors.end(), toNode_rc), fromNode_rc->_successors.end());
        //toNode_rc->_predecessors.erase(std::remove(toNode_rc->_predecessors.begin(), toNode_rc->_predecessors.end(), fromNode_rc), toNode_rc->_predecessors.end());
    

        //nodeChanged(fromNode);
        //nodeChanged(toNode);
    }


    static void reverseComplementUnitigs(const vector<UnitigType>& unitigs, vector<UnitigType>& unitigsRC){

        unitigsRC.resize(unitigs.size(), 0);
        

        for(size_t i=0; i<unitigs.size(); i++){
            unitigsRC[i] = unitigIndex_to_reverseDirection(unitigs[unitigs.size()-i-1]);
        }

    }

    vector<MinimizerType> unitigSequenceToMinimizerSequence(const vector<UnitigType>& unitigs){

        vector<MinimizerType> minimizersSequence;

        for(size_t i=0; i<unitigs.size(); i++){

            UnitigType unitigIndex = unitigs[i];

            bool isReversed;
            const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

            vector<MinimizerType> minimizers = _unitigName_to_minimizers[unitigName];

            if(isReversed){
                std::reverse(minimizers.begin(), minimizers.end());
            }

            if(i==0){

                minimizersSequence = minimizers;
            }
            else{

                minimizersSequence.insert( minimizersSequence.end(), minimizers.begin()+_kminmerSize-1, minimizers.end());
                //for(size_t i=_kminmerSize-1; i<minimizers.size(); i++){
                //    minimizersSequence.push_back(minimizers[i]);
                //}
            }
        }

        return minimizersSequence;
    }

    void mergeNode(const UnitigType& unitigIndex1, const UnitigType& unitigIndex2){

        //cout << "Merge unitigs: " << unitigIndexToString(unitigIndex1) << " -> " << unitigIndexToString(unitigIndex2) << endl;


        bool isReversed1;
        const UnitigType& unitigName1 = unitigIndex_to_unitigName(unitigIndex1, isReversed1);

        bool isReversed2;
        const UnitigType& unitigName2 = unitigIndex_to_unitigName(unitigIndex2, isReversed2);

        UnitigNode* unitig1 = _unitigs[unitigName1];
        UnitigNode* unitig2 = _unitigs[unitigName2];

        /*

        bool print = false;
        if(unitigName1 == 1353 || unitigName2 == 1353){
            cout << "\t" << unitigIndexToString(unitigIndex1) << endl;
            for(const MinimizerType& m : unitig1->_nodes){
                cout << "\t\t" << m << endl;
            }
            cout << "\t" << unitigIndexToString(unitigIndex2) << endl;
            for(const MinimizerType& m : unitig2->_nodes){
                cout << "\t\t" << m << endl;
            }

            //print = true;
        }
        */

        if(unitig1->_unitigMerge.size() == 0){
            unitig1->_isReversed = isReversed1;
            unitig1->_unitigMerge.push_back(unitigIndex1);
        }

        if(unitig1->_isReversed != isReversed1){
            unitig1->_isReversed = isReversed1;

            vector<UnitigType> unitigsRC;
            reverseComplementUnitigs(unitig1->_unitigMerge, unitigsRC);

            unitig1->_unitigMerge = unitigsRC;

            //u_int32_t memo = unitig1->_startNode_forward;
            //unitig1->_startNode_forward = unitig1->_startNode_reverse;
            //unitig1->_startNode_reverse = memo;
        }

        if(unitig2->_unitigMerge.size() == 0){
            unitig1->_unitigMerge.push_back(unitigIndex2);
        }
        else{
            if(unitig2->_isReversed != isReversed2){

                vector<UnitigType> unitigsRC;
                reverseComplementUnitigs(unitig2->_unitigMerge, unitigsRC);
                
                unitig1->_unitigMerge.insert( unitig1->_unitigMerge.end(), unitigsRC.begin(), unitigsRC.end());

                //unitig1->_startNode_reverse = unitig2->_startNode_forward;
            }
            else{

                unitig1->_unitigMerge.insert( unitig1->_unitigMerge.end(), unitig2->_unitigMerge.begin(), unitig2->_unitigMerge.end());
                
                //unitig1->_startNode_reverse = unitig2->_startNode_reverse;
            }
        }

        /*
        if(!isReversed1 && !isReversed2){

            //if(print) cout << "1" << endl;

            vector<MinimizerType> minimizers = unitig2->_nodes;


            unitig1->_nodes.insert( unitig1->_nodes.end(), minimizers.begin()+_kminmerSize-1, minimizers.end());



        }
        else if(!isReversed1 && isReversed2){

            //if(print) cout << "2" << endl;

            vector<MinimizerType> minimizers = unitig2->_nodes;
            std::reverse(minimizers.begin(), minimizers.end());


            unitig1->_nodes.insert( unitig1->_nodes.end(), minimizers.begin()+_kminmerSize-1, minimizers.end());



            //getchar();

        }
        else if(isReversed1 && !isReversed2){
            //if(print) cout << "3" << endl;

            vector<MinimizerType> minimizers = unitig2->_nodes;
            std::reverse(minimizers.begin(), minimizers.end());


            unitig1->_nodes.insert( unitig1->_nodes.begin(), minimizers.begin(), minimizers.end()-_kminmerSize+1);


        }
        else if(isReversed1 && isReversed2){
            
            //if(print) cout << "4" << endl;

            vector<MinimizerType> minimizers = unitig2->_nodes;
            //std::reverse(minimizers.begin(), minimizers.end());


            unitig1->_nodes.insert( unitig1->_nodes.begin(), minimizers.begin(), minimizers.end()-_kminmerSize+1);


        }
        */

        /*
        vector<UnitigType> unitigsRC;
        reverseComplementUnitigs(unitig1->_unitigMerge, unitigsRC);

        const vector<MinimizerType>& minimizersRC = unitigSequenceToMinimizerSequence(unitigsRC);
        const vector<MinimizerType>& minimizers = unitigSequenceToMinimizerSequence(unitig1->_unitigMerge);

        bool print = true;
        
        if(minimizers == unitig1->_nodes || minimizersRC == unitig1->_nodes){
            print = false;
        }

        


        
        if(print){

            cout << "issue " << endl;

            cout << "\tresult:" << endl;
            for(auto m : unitig1->_nodes){
                cout << "\t\t" << m << endl;
            }

            cout << "\tresult2:" << endl;

            for(const UnitigType& unitigIndex : unitig1->_unitigMerge){
                cout << "\t\t" << unitigIndexToString(unitigIndex) << endl;
            }

            vector<UnitigType> unitigsRC;
            reverseComplementUnitigs(unitig1->_unitigMerge, unitigsRC);

            cout << "\tresult3:" << endl;

            for(const UnitigType& unitigIndex : unitigsRC){
                cout << "\t\t" << unitigIndexToString(unitigIndex) << endl;
            }

            for(MinimizerType m : unitigSequenceToMinimizerSequence(unitigsRC)){
                cout << "\t\t" << m << endl;
            }
            getchar();
        }
        */

        nodeChanged(unitig1);
        nodeChanged(unitig2);
        
        unitig1->mergeWith(unitig2, this);

        


        vector<UnitigType> successors;
        getSuccessors(unitigIndex2, successors);
        
        u_int32_t unitigIndexToReplace = unitigIndex_to_reverseDirection(unitigIndex2);
        u_int32_t unitigIndex1_rev = unitigIndex_to_reverseDirection(unitigIndex1);
        //vector<UnitigType> predecessors;
        //getPredecessors(unitigIndex, predecessors);

        for(const UnitigType& unitigIndexSuccessor : successors){

            bool isReversedSuccessor;
            const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);

            //cout << "\tEdge: " << unitigIndexToString(unitigIndex2) << " -> " << unitigIndexToString(unitigIndexSuccessor) << endl;

            if(isReversedSuccessor){

                vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

                std::replace (successorsFoward.begin(), successorsFoward.end(), unitigIndexToReplace, unitigIndex1_rev);
                //for(const UnitigTyp)
            }
            else{

                vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

                std::replace (successorsReverse.begin(), successorsReverse.end(), unitigIndexToReplace, unitigIndex1_rev);
            }
            /*

            if(isReversedSuccessor){

                vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

                cout << "\t\tNb successors: " << successorsFoward.size() << endl;
                successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexToRemove), successorsFoward.end());
                cout << "\t\tNb successors: " << successorsFoward.size() << endl;


            }
            else{

                vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

                cout << "\t\tNb successors: " << successorsReverse.size() << endl;
                successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexToRemove), successorsReverse.end());
                cout << "\t\tNb successors: " << successorsReverse.size() << endl;



            }
            */
            /*
            

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_forward){

                cout << "\t\tEdge neighbor forward: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\t\tFound forward" << endl;
                }
            }

            for(const UnitigType& unitigIndexLol : _unitigs[unitigNameSuccessor]->_successors_reverse){

                cout << "\t\tEdge neighbor reverse: " << unitigIndexSuccessor << " -> " << unitigIndexLol << endl;

                if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
                    cout << "\t\tFound reverse" << endl;
                }
            }
            */
            //const UnitigType& 
            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            //nodeChanged(nn);
        }

        if(isReversed1){
            if(isReversed2){
                unitig1->_successors_reverse = unitig2->_successors_reverse;
            }
            else{
                unitig1->_successors_reverse = unitig2->_successors_forward;
            }
        }
        else{
            if(isReversed2){
                unitig1->_successors_forward = unitig2->_successors_reverse;
            }
            else{
                unitig1->_successors_forward = unitig2->_successors_forward;
            }
        }


        /*
        for(Node* nn : node2->_successors){


            for(size_t i=0; i<nn->_predecessors.size(); i++){
                if(nn->_predecessors[i] == node2){
                    nn->_predecessors[i] = node1;
                    break;
                }
            }

        }
        */

        //if(unitig2->_unitigName == 2025){
        //    cout << "plop" << endl;
        //    exit(1);
        //}

        delete unitig2;
        _unitigs[unitigName2] = nullptr;





        /*
        for(const UnitigType& unitigIndexLol : _unitigs[unitigName1]->_successors_forward){

            cout << "\t\tEdge neighbor forward: " << unitigName1 << " -> " << unitigIndex_to_unitigName(unitigIndexLol) << endl;

            //if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
            //    cout << "\t\t\tFound forward" << endl;
            //}
        }

        for(const UnitigType& unitigIndexLol : _unitigs[unitigName1]->_successors_reverse){

            cout << "\t\tEdge neighbor reverse: " << unitigName1 << " -> " << unitigIndex_to_unitigName(unitigIndexLol) << endl;

            //if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
            //    cout << "\t\tFound reverse" << endl;
            //}
        }

        
        for(const UnitigType& unitigIndexLol : _unitigs[unitigName2]->_successors_forward){

            cout << "\t\tEdge neighbor forward: " << unitigName2 << " -> " << unitigIndex_to_unitigName(unitigIndexLol) << endl;

            //if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
            //    cout << "\t\t\tFound forward" << endl;
            //}
        }

        for(const UnitigType& unitigIndexLol : _unitigs[unitigName2]->_successors_reverse){

            cout << "\t\tEdge neighbor reverse: " << unitigName2 << " -> " << unitigIndex_to_unitigName(unitigIndexLol) << endl;

            //if(unitigIndex_to_unitigName(unitigIndexLol) == node->_unitigName){
            //    cout << "\t\tFound reverse" << endl;
            //}
        }

        */
        /*
        if(isReversedSuccessor){

            vector<UnitigType>& successorsFoward = _unitigs[unitigNameSuccessor]->_successors_forward;

            cout << "\t\tNb successors: " << successorsFoward.size() << endl;
            successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexToRemove), successorsFoward.end());
            cout << "\t\tNb successors: " << successorsFoward.size() << endl;


        }
        else{

            vector<UnitigType>& successorsReverse = _unitigs[unitigNameSuccessor]->_successors_reverse;

            cout << "\t\tNb successors: " << successorsReverse.size() << endl;
            successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexToRemove), successorsReverse.end());
            cout << "\t\tNb successors: " << successorsReverse.size() << endl;



        }
        */


        /*
        Node* node1_rc = unitigIndex_toReverseDirection(node1);
        Node* node2_rc = unitigIndex_toReverseDirection(node2);

        nodeChanged(node1);
        nodeChanged(node2);


        node1->mergeWith(node2, this);

        node1->_successors = node2->_successors;

        for(Node* nn : node2->_successors){


            for(size_t i=0; i<nn->_predecessors.size(); i++){
                if(nn->_predecessors[i] == node2){
                    nn->_predecessors[i] = node1;
                    break;
                }
            }

        }

        node1_rc->_length = node1->_length;
        node1_rc->_abundance = node1->_abundance;
        node1_rc->_abundances = node1->_abundances;
        node1_rc->_nodes = node1->_nodes;
        std::reverse(node1_rc->_nodes.begin(), node1_rc->_nodes.end());
        //for(size_t i=0; i<node1_rc->_nodes.size(); i++){
        //    node1_rc->_nodes[i] = nodeIndex_toReverseDirection(node1_rc->_nodes[i]);
        //}


        //node1_rc->_predecessors.erase(std::remove(node1_rc->_predecessors.begin(), node1_rc->_predecessors.end(), node2_rc), node1_rc->_predecessors.end());
        node1_rc->_predecessors = node2_rc->_predecessors;

        for(Node* nn : node2_rc->_predecessors){
            for(size_t i=0; i<nn->_successors.size(); i++){
                if(nn->_successors[i] == node2_rc){
                    nn->_successors[i] = node1_rc;
                    break;
                }
            }
        }



        node2_rc->clear();
        node2->clear();
        */
    }


    u_int64_t getNbUnitigs() const{

        u_int64_t nbUnitigs = 0;

        for(const UnitigNode* unitig : _unitigs){
            if(unitig == nullptr) continue;

            nbUnitigs += 1;
        }

        return nbUnitigs;
    }

    void save(const string& outputFilename){

        u_int64_t checksum = 0;
        u_int64_t sumAbundance = 0;
        u_int64_t sumLength = 0;

        //cout << outputFilename << endl;
        //cout << "\ttodo: remettre la partie qui correspond aux sequences" << endl;
        //cout << "\tNb unitigs: " << getNbUnitigs() << endl;
        //return;


        Logger::get().debug() << "Saving unitig graph: " << outputFilename;

		ofstream outputContigFile(outputFilename + ".unitigs.nodepath");

		unordered_set<UnitigType> writtenUnitigs;

        unordered_set<UnitigType> selectedUnitigIndex;

        u_int64_t nbMinimizers = 0;

        /*
        for(const UnitigNode* unitig : _unitigs){
            //if(node == nullptr) continue;
            if(unitig->_unitigName == -1) continue;
            //cout << node->_unitigIndex << endl;


			if(writtenUnitigs.find(unitig->_unitigName) != writtenUnitigs.end()) continue;
			//if(writtenUnitigs.find(unitigIndex_toReverseDirection(node->_unitigIndex)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(unitig->_unitigName);
			//writtenUnitigs.insert(unitigIndex_toReverseDirection(node->_unitigIndex));

            UnitigType unitigIndex = unitigName_to_unitigIndex(unitig->_unitigName, false);
            selectedUnitigIndex.insert(unitigIndex);

            //for(u_int32_t nodeIndex : node->_nodes){
            //    u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		    //    file_nodeNameToUnitigIndex.write((const char*)&nodeName, sizeof(nodeName));
		    //    file_nodeNameToUnitigIndex.write((const char*)&node->_unitigIndex, sizeof(node->_unitigIndex));
            //}


        }
        */
        //file_nodeNameToUnitigIndex.close();

        
        ofstream outputFile(outputFilename);


        unordered_set<UnitigType> isVisited;

        /*
        for(UnitigType sourceUnitigIndex = 0; sourceUnitigIndex < _unitigs.size()*2; sourceUnitigIndex++){

            //bool isReversed;
            UnitigType sourceUnitigName = unitigIndex_to_reverseDirection(sourceUnitigIndex);
            const UnitigNode* unitig = _unitigs[sourceUnitigName];
            
            //if(node == nullptr) continue;
            if(unitig == nullptr || unitig->_unitigName == -1) continue;
        */

        for(const UnitigNode* unitig : _unitigs){
            //if(node == nullptr) continue;
            if(unitig == nullptr) continue;

            UnitigType sourceUnitigIndex = unitigName_to_unitigIndex(unitig->_unitigName, false);

            //if(selectedUnitigIndex.find(unitig->_unitigName) == selectedUnitigIndex.end()) continue;
            if(isVisited.find(unitig->_unitigName) != isVisited.end()) continue;

            queue<UnitigType> queue;
            queue.push(sourceUnitigIndex);

            while(queue.size() > 0){

                UnitigType unitigIndex = queue.front();
                queue.pop();


                bool isReversed;
                UnitigType unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);
                

                if(isVisited.find(unitigName) != isVisited.end()) continue;
                isVisited.insert(unitigName);

                //UnitigType unitigIndex = unitigName_to_unitigIndex(unitigName, false);
                
                string ori1 = "+";
                if(isReversed){
                    ori1 = "-";
                }

                //string ori = "+";
                //if(selectedUnitigIndex.find(unitigIndex) == selectedUnitigIndex.end()){
                //    unitigIndex = unitigIndex_to_reverseDirection(unitigIndex);
                //    ori = "-";
                //}

                //isVisited.insert(unitigIndex_toReverseDirection(unitigIndex));
                //linkedUnitigIndex.insert(unitigIndex);
                //float abundance = 1;//_nodes[unitigIndex]->_abundance;

                UnitigGraph2::UnitigNode* unitig = _unitigs[unitigName];

                outputFile << "S" << "\tutg" << unitigName << "\t" << "*" << "\t" << "LN:i:" << unitig->getLength(this) << "\t" << "dp:i:" << unitig->_abundance << endl;
                

                u_int8_t isCircular = 0; //isCircular(unitig);
				//vector<u_int32_t> nodePath = nodeCurrent->_nodes;
				//if(isCircular && nodePath.size() > 1){
				//	nodePath.push_back(nodePath[0]);
				//}

                vector<UnitigType> unitigs;
                if(unitig->_unitigMerge.size() == 0){
                    unitigs.push_back(unitigIndex);    
                }
                else{
                    unitigs = unitig->_unitigMerge;
                }

                //if(isReversed){
                //    vector<UnitigType> unitigsRC;
                //    reverseComplementUnitigs(unitigs, unitigsRC);
                //    unitigs = unitigsRC;
                //}

                u_int32_t size = unitigs.size();

				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&unitigs[0], size * sizeof(UnitigType));

                //cout << endl;
                //cout << "utg" << unitigName << endl;
                //for(float abundance : _unitigs[unitigName]->_abundances){
                //    cout << abundance << " ";
                //}
                //cout << endl;

                sumAbundance += _unitigs[unitigName]->_abundance;
                sumLength += _unitigs[unitigName]->getLength(this);
                checksum += (_unitigs[unitigName]->getLength(this)*_unitigs[unitigName]->_abundance);
                
                //cout << unitigName << " " << _unitigs[unitigName]->getLength(this) << " " << _unitigs[unitigName]->_abundance << " " << checksum << endl;
                nbMinimizers += _unitigs[unitigName]->_nbMinimizers;


                //outputContigFileToUnitigIndex.write((const char*)&unitigIndex, sizeof(unitigIndex));
                //vector<u_int32_t> successors;
                //getSuccessors_unitig(unitigIndex, 0, successors);

                vector<UnitigType> successors;
                getSuccessors(unitigIndex, successors);


                for(const UnitigType& unitigIndexSuccessor : successors){
                    
                    //if(unitigIndex == 0){
                    //    cout << "\tSucc: " << unitigIndexSuccessor << endl;
                    //}
                    bool isReversedSuccessor;
                    const UnitigType& unitigNameSuccessor = unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);
                    //u_int32_t unitigIndexN = nodeSuccessor->_unitigIndex;

                    string ori2 = "+";
                    if(isReversedSuccessor){
                        ori2 = "-";
                    }

                    //if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                    //    unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                    //    ori2 = "-";
                    //}
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigName << "\t" << ori1 << "\tutg" << unitigNameSuccessor << "\t" << ori2 << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexSuccessor);

                    //if(unitigIndexN > 100000){
                        //cout << "wtf succ: " << unitigIndex << " -> " << unitigIndexN << endl;
                        //getchar();
                    //}
                }

                
                //vector<u_int32_t> predecessors;
                //getPredecessors_unitig(unitigIndex, 0, predecessors);
                //for(u_int32_t unitigIndexN : predecessors){


                vector<UnitigType> predecessors;
                getPredecessors(unitigIndex, predecessors);

                for(const UnitigType& unitigIndexPredecessor : predecessors){


                    //if(unitigIndex == 0){
                    //    cout << "\tPred: " << unitigIndexPredecessor << endl;
                    //}

                    bool isReversedPredecessor;
                    const UnitigType& unitigNamePredecessor = unitigIndex_to_unitigName(unitigIndexPredecessor, isReversedPredecessor);
                    //u_int32_t unitigIndexN = nodeSuccessor->_unitigIndex;

                    string ori2 = "+";
                    if(isReversedPredecessor){
                        ori2 = "-";
                    }

                    //if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                    //    unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                    //    ori2 = "-";
                    //}
                    
                    u_int32_t overlap = 1;
                    outputFile << "L" << "\tutg" << unitigNamePredecessor << "\t" << ori2 << "\tutg" << unitigName << "\t" << ori1 << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexPredecessor);
                }

                

            }
        }


        //if(_kminmerSize >= 18) getchar();

        outputFile.close();
        //colorFile.close();
        outputContigFile.close();
        //outputContigFileToUnitigIndex.close();

        //_logFile << "\tdone" << endl;


        //cout << "\tNb minimizers: " << nbMinimizers << endl;
        //cout << "\tSum abundance: " << sumAbundance << endl;
        //cout << "\tSum length: " << sumLength << endl;
        //cout << "\tGraph abundance: " << getChecksumAbundance() << endl;
        //cout << "\tChecksum: " << checksum << endl;

        //getchar();


        
    }
    
    u_int64_t getChecksumAbundance(){

        u_int64_t sumAbundance = 0;

        for(UnitigNode* unitig : _unitigs){
            if(unitig == nullptr) continue;

            sumAbundance += unitig->_abundance;
            
        }

        return sumAbundance;
    }

    u_int64_t getChecksum(){

        u_int64_t sumAbundance = 0;

        for(UnitigNode* unitig : _unitigs){
            if(unitig == nullptr) continue;

            sumAbundance += (unitig->_nbMinimizers* unitig->_abundance);
            
        }

        return sumAbundance;
    }

    void getSuccessors(const UnitigType& unitigIndex, vector<UnitigType>& successors){

        successors.clear();

        bool isReversed;
        UnitigType unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

        if(isReversed){
            for(const UnitigType& unitigIndexSuccessor : _unitigs[unitigName]->_successors_reverse){
                successors.push_back(unitigIndexSuccessor);
            } 
        }
        else{
            for(const UnitigType& unitigIndexSuccessor : _unitigs[unitigName]->_successors_forward){
                successors.push_back(unitigIndexSuccessor);
            } 
        }
        
    }

    void getPredecessors(const UnitigType& unitigIndex, vector<UnitigType>& predecessors){

        predecessors.clear();

        const u_int32_t& unitigIndexRev = unitigIndex_to_reverseDirection(unitigIndex);
        
        bool isReversed;
        UnitigType unitigName = unitigIndex_to_unitigName(unitigIndexRev, isReversed);

        if(isReversed){
            for(const UnitigType& unitigIndexSuccessor : _unitigs[unitigName]->_successors_reverse){
                predecessors.push_back(unitigIndex_to_reverseDirection(unitigIndexSuccessor));
            } 
        }
        else{
            for(const UnitigType& unitigIndexSuccessor : _unitigs[unitigName]->_successors_forward){
                predecessors.push_back(unitigIndex_to_reverseDirection(unitigIndexSuccessor));
            } 
        }
        
    }

    bool hasSuccessors(const UnitigType& unitigIndex){


        bool isReversed;
        UnitigType unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

        if(isReversed){
            return _unitigs[unitigName]->_successors_reverse.size() > 0;
        }
        else{
            return _unitigs[unitigName]->_successors_forward.size() > 0;
        }
        
    }

    u_int64_t nbSuccessors(const UnitigType& unitigIndex){

        u_int64_t nbSuccessors = 0;

        bool isReversed;
        UnitigType unitigName = unitigIndex_to_unitigName(unitigIndex, isReversed);

        if(isReversed){
            return _unitigs[unitigName]->_successors_reverse.size();
        }
        else{
            return _unitigs[unitigName]->_successors_forward.size();
        }
        
    }

    bool hasNoPredecessors(const UnitigType& unitigIndex){


        const u_int32_t& unitigIndexRev = unitigIndex_to_reverseDirection(unitigIndex);
        
        bool isReversed;
        UnitigType unitigName = unitigIndex_to_unitigName(unitigIndexRev, isReversed);

        if(isReversed){
            return _unitigs[unitigName]->_successors_reverse.size() == 0;
        }
        else{
            return _unitigs[unitigName]->_successors_forward.size() == 0;
        }
        
    }

    u_int8_t isCircular(const UnitigNode* unitig){

        const UnitigType& unitigIndex = unitigName_to_unitigIndex(unitig->_unitigName, false);

        vector<UnitigType> successors;
        getSuccessors(unitigIndex, successors);

        vector<UnitigType> predecessors;
        getPredecessors(unitigIndex, predecessors);

        //const vector<Node*>& successors = _successors;
        //const vector<Node*>& predecessors = getPredecessors(unitigGraph);
        return (unitig->_nbMinimizers -_kminmerSize+1) > 1 && (successors.size() == 1) && (predecessors.size() == 1) && (successors[0] == unitigIndex) && (predecessors[0] == unitigIndex);  //_nodes.size() > 1 && (startNode() == endNode());
    }

	bool isRepeatSide(UnitigNode* unitig){

		//if(node->_length > 50000) return false;
		if((unitig->_nbMinimizers-_kminmerSize+1) > _kminmerSize*2) return false;
		if(unitig->_successors_forward.size() == 0) return false;
		if(unitig->_successors_reverse.size() == 0) return false;


        const UnitigType& unitigIndex = unitigName_to_unitigIndex(unitig->_unitigName, false);

        vector<UnitigType> successors;
        getSuccessors(unitigIndex, successors);

        vector<UnitigType> predecessors;
        getPredecessors(unitigIndex, predecessors);


		for(const UnitigType& unitigIndexSuccessor : successors){

            const UnitigType& unitigNameSuccessor = unitigIndex_to_unitigName(unitigIndexSuccessor);
			if(unitigNameSuccessor == unitig->_unitigName) continue;

		    for(const UnitigType& unitigIndexPredecessor : predecessors){

                const UnitigType& unitigNamePredecessor = unitigIndex_to_unitigName(unitigIndexPredecessor);
			    if(unitigNamePredecessor == unitig->_unitigName) continue;

				if(unitigIndexSuccessor == unitigIndexPredecessor){
					return true;
				}
			}
		}

		return false;
	}

    u_int64_t computeChecksum_successors(){

        //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala_1.txt");

        u_int64_t checksumTotal = 0;

        for(UnitigNode* node : _unitigs){

            if(node == nullptr) continue;


            checksumTotal += (node->_abundance * (node->_nbMinimizers-_kminmerSize+1) * node->_successors_forward.size() * node->_successors_reverse.size());

            //debugFile << node->_unitigName << " " << node->_abundance << " " << (node->_nbMinimizers-_kminmerSize+1) << " " << checksumTotal << endl;
        }

        //debugFile.close();

        return checksumTotal;
    }

    u_int64_t computeChecksum_abundance(){

        u_int64_t checksumTotal = 0;

        for(UnitigNode* node : _unitigs){

            if(node == nullptr) continue;

            u_int64_t checksum = 0;
            for(u_int32_t ab : node->_abundances){
                //checksum += ab;
            }

            checksum += (node->_abundance * (node->_nbMinimizers-_kminmerSize+1));

            checksumTotal += checksum;
        }

        return checksumTotal;
    }

	template<typename Functor>
	void processUnitigsParallel(const Functor& functor, int nbCores){

		#pragma omp parallel num_threads(nbCores)
		{

			Functor functorSub(functor);

            #pragma omp for
            for(size_t i=0; i < _unitigs.size(); i++){

                if(_unitigs[i] == nullptr) continue;

                functorSub(_unitigs[i]);
            }
			
		}
	}


	class AbundanceSortFunctor {

		public:

		UnitigGraph2* _unitigGraph;

		AbundanceSortFunctor(UnitigGraph2* unitigGraph) : _unitigGraph(unitigGraph){
		}

		AbundanceSortFunctor(const AbundanceSortFunctor& copy) : _unitigGraph(copy._unitigGraph){
		}

		~AbundanceSortFunctor(){
		}


		void operator () (UnitigNode* unitig) const {
        
            std::sort(unitig->_abundances.begin(), unitig->_abundances.end());
            unitig->_abundance = unitig->computeMedianAbundance();
		}
	};

	void partitionNodes(int k){
		
        int nbNodesPerPartition = getNbUnitigs() / k;
        cout << "Nb nodes per partition: " << nbNodesPerPartition << endl;
		//ofstream partitionFile(_partitionFilename);

        vector<UnitigType> currentPartition;
		//ReadPartition currentPartition;

		//vector<bool> isCorrected(_unitigs.size(), false);
		vector<int> unitigName_to_partition(_unitigs.size(), false);
		vector<bool> isVisited(_unitigs.size(), false);
		//size_t nbPartitions = 0;
        size_t partitionIndex = 0;

		for(UnitigType unitigName = 0; unitigName < _unitigs.size(); unitigName++){

            if(_unitigs[unitigName] == nullptr) continue;
            if(isVisited[unitigName]) continue;

            queue<UnitigType> q; 

            isVisited[unitigName] = true;
            q.push(unitigName);

            while (!q.empty()) {
                UnitigType curr = q.front(); 
                q.pop();
                
                currentPartition.push_back(curr);

                for(UnitigType x : _unitigs[curr]->_successors_forward){
                    x /= 2;
                    if(isVisited[x]) continue;
                    if(currentPartition.size() > nbNodesPerPartition) continue;
                    isVisited[x] = true; 
                    q.push(x); 
                }

                for(UnitigType x : _unitigs[curr]->_successors_reverse){
                    x /= 2;
                    if(isVisited[x]) continue;
                    if(currentPartition.size() > nbNodesPerPartition) continue;
                    isVisited[x] = true; 
                    q.push(x); 
                }

				
            }

            if(currentPartition.size() > nbNodesPerPartition){

                for(UnitigType unitigName : currentPartition){
                    unitigName_to_partition[unitigName] = partitionIndex;
                }

                cout << "Partiton detected: " <<    currentPartition.size() << endl;
                currentPartition.clear();
                partitionIndex += 1;
            }
        }


        if(currentPartition.size() > 0){
            for(UnitigType unitigName : currentPartition){
                unitigName_to_partition[unitigName] = partitionIndex;
            }

            cout << "Partiton detected: " <<    currentPartition.size() << endl;
            currentPartition.clear();
            partitionIndex += 1;
        }




        int nbIndependentNodes = 0;
        int nbCriticalNodes = 0;

		for(UnitigType unitigName = 0; unitigName < _unitigs.size(); unitigName++){

            if(_unitigs[unitigName] == nullptr) continue;

            int currentPartition = unitigName_to_partition[unitigName];

            bool isIndependent = true;

            for(UnitigType x : _unitigs[unitigName]->_successors_forward){
                x /= 2;

                if(unitigName_to_partition[x] != currentPartition){
                    isIndependent = false;
                }
            }

            for(UnitigType x : _unitigs[unitigName]->_successors_reverse){
                x /= 2;

                if(unitigName_to_partition[x] != currentPartition){
                    isIndependent = false;
                }
            }

            if(isIndependent){
                nbIndependentNodes += 1;
            }
            else{
                nbCriticalNodes += 1;
            }
        }

        cout << "Nb independant nodes: " << nbIndependentNodes << endl;
        cout << "Nb critical nodes:    " << nbCriticalNodes << endl;
        getchar();
	}

};


#endif

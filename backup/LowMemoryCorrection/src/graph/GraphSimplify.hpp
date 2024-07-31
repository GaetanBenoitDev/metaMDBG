

//#include "Graph.hpp"
//TODO:
// remove bubble, tips: consider length

//#define PRINT_DEBUG_SIMPLIFICATION



#ifndef MDBG_METAG_GRAPHSIMPLIFY
#define MDBG_METAG_GRAPHSIMPLIFY

#include "Commons.hpp"
#include "ProgressiveAbundanceFilter.hpp"
#include <limits>
//#include "contigFeatures/ContigFeature.hpp"
//#include <experimental/filesystem>

//using namespace std;
//namespace fs = std::experimental::filesystem;

//class GraphSimplify;

struct NodeAbundance{
    u_int32_t _nodeIndex;
    double _abundance;
};

struct UnitigTip{
    u_int32_t _unitigIndex;
    u_int32_t _startNodeIndex;
    u_int32_t _length;
};

struct BubbleSide{
    u_int32_t _unitigIndex;
    u_int32_t _nbNodes;
    u_int32_t _nodeIndexSum;
    //u_int32_t _quality;
};

struct Bubble{
    u_int32_t _startNode;
    u_int32_t _endNode;
    vector<NodeAbundance> _nodes;
};


//const float _globalCutoff = 3;

struct Unitig{
    u_int32_t _index;
    u_int32_t _startNode;
    u_int32_t _endNode;
    float _abundance;
    u_int32_t _length;
    u_int32_t _nbNodes;
    vector<u_int32_t> _nodes;
    //u_int32_t _quality;
    //u_int32_t _debug_nodeIndexSource;
    //u_int32_t _debug_nodeIndexSum;
};

struct ComplexArea{
    u_int32_t _nodeIndex_source;
    u_int32_t _nodeIndex_sink;
    vector<u_int32_t> _unitigs;
    vector<u_int32_t> _path;
    //float _abundanceCutoff_min;
};


//"todo: take superbubble dans algo d'assemblage, ça peut loop actuellement"
//"calculer la vrai length des unitigs en prenant en compte les overlaps"
//"ne plus utiliser de component connexe qunad abondance < 30"

struct SaveState{
    //float abundanceCutoff_min;
    unordered_set<u_int32_t> _isNodeValid2;
    vector<bool> _isBubble;
    unordered_set<DbgEdge, hash_pair> _isEdgeRemoved;
    vector<Bubble> _bubbles;
    unordered_map<u_int32_t, u_int32_t> _nodeToUnitig;
    vector<Unitig> _unitigs;
    unordered_set<u_int32_t> _isLongUnitig;
    unordered_map<u_int32_t, u_int32_t> _longUnitigNodeAbundance;
};

struct SaveState2{
    //float abundanceCutoff_min;
    float _abundanceCutoff_min;
    unordered_set<u_int32_t> _nodeNameRemoved_tmp;
    vector<u_int32_t> _nodeNameRemoved;
    vector<DbgEdge> _isEdgeRemoved;
    vector<u_int32_t> _isBubble;
    vector<NodeAbundance> _longUnitigNodeAbundance;
    vector<u_int64_t> _rareReads;
    //vector<Bubble> _bubbles;
    //unordered_map<u_int32_t, u_int32_t> _longUnitigNodeAbundance;
};

class GraphSimplify{

public:


    ofstream& _logFile;
    BiGraph* _graphSuccessors;
    //BiGraph* _graphPredecessors;
    string _outputDir;
    vector<Unitig> _unitigs;
    //vector<u_int32_t> _nodeToUnitig;
    //vector<u_int32_t> _nodeToUnitig;
    //vector<u_int32_t> _nodeAbundances;
    //vector<u_int32_t> _nodeLengths;
    //unordered_map<Destroy path: 2u_int32_t, u_int32_t> _nodeAbundances;
    //unordered_map<u_int32_t, u_int32_t> _nodeLengths;
    //vector<bool> _isNodeValid;
    unordered_set<u_int32_t> _isNodeErased;
    //unordered_set<u_int32_t> _isNodeValid2;
    //unordered_set<u_int32_t> _isNodeValid2_cache;
    vector<bool> _isBubble;
    vector<bool> _isNodeIndexIndexed;
    //unordered_set<DbgEdge, hash_pair> _isEdgeRemoved;
    //unordered_set<DbgEdge, hash_pair> _isEdgeUnsupported;

    u_int32_t _nbNodes;
    vector<Bubble> _bubbles;
    unordered_set<u_int32_t> _isNodeInvalid_tmp;
    unordered_set<u_int32_t> _isUnitigInvalid_tmp;
    //GraphSimplify::SuperbubbleSimplify* _superbubbleSimplify;
    //u_int32_t _lastAbundanceCutoff;
    unordered_map<u_int32_t, ComplexArea> _complexAreas_source; 
    unordered_map<u_int32_t, ComplexArea> _complexAreas_sink;
    //unordered_set<u_int32_t> _isLongUnitig;
    unordered_map<u_int32_t, u_int32_t> _longUnitigNodeAbundance;
    vector<vector<UnitigLength>> _startingUnitigs;
    vector<Unitig> _startingUnitigstest;
    unordered_set<u_int32_t> _currentUnitigNodes;
    unordered_set<u_int64_t> _removedReadIndex;
    unordered_set<u_int32_t> _isBoundUnitig;
    unordered_set<u_int32_t> _debug_groundTruthNodeNames;
    unordered_set<u_int64_t> _rareReads;
    unordered_set<u_int32_t> _repeatedNodenames;
    size_t _kminmerSize;
    unordered_set<u_int32_t> _isNodenameRoundabout;
    bool _isCopy;
    vector<Unitig> _cleanedLongTips;
    //ContigFeature* _contigFeature;
    int _nbCores;
    u_int64_t _superbubble_checksum;
    u_int64_t _tip_checksum;
    u_int64_t _abcutoff_checksum;
    vector<bool> _isNodeRemoved;

	//float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    unordered_map<u_int32_t, u_int32_t> _endNode_to_unitigIndex;
    unordered_map<u_int32_t, u_int32_t> _startNode_to_unitigIndex;

    //const vector<u_int32_t>& _nodeName_to_abundance;

    GraphSimplify(BiGraph* graphSuccessors, const string& outputDir, u_int32_t nbNodes, size_t kminmerSize, int nbCores, float kminmerLengthMean, float kminmerOverlapMean, ofstream& logFile) : _logFile(logFile){

        //_unitigGraph = nullptr;
        _kminmerLengthMean = kminmerLengthMean;
        _kminmerOverlapMean = kminmerOverlapMean;

        _superbubble_checksum = 0;
        _tip_checksum = 0;
        _abcutoff_checksum = 0;
        //_bubble_checksum = 0;

        //_contigFeature = nullptr;
        _graphSuccessors = graphSuccessors;
        _outputDir = outputDir;
        _kminmerSize = kminmerSize;
        _nbCores = nbCores;

    
        _nbNodes = _graphSuccessors->_nbNodes;
        _isCopy = false;


    }



    ~GraphSimplify(){
        //if(!_isCopy) delete _graphSuccessors;
        //delete _graphPredecessors;
    }


    void clear(float abundanceCutoff_min){

        _nextUnitigIndex = 0;
        _unitigs.clear();
        //_nodeToUnitig = vector<u_int32_t>(_graphSuccessors->_nbNodes, -1);

        _cleanedLongTips.clear();
        //_isNodenameRoundabout.clear();
        _repeatedNodenames.clear();
        _rareReads.clear();
        _isBoundUnitig.clear();
        _removedReadIndex.clear();
        _currentUnitigNodes.clear();
        //_isLongUnitig.clear();
        _longUnitigNodeAbundance.clear();
        _isNodeInvalid_tmp.clear();
        //_isNodeValid2.clear();
        _isUnitigInvalid_tmp.clear();


        _isNodeIndexIndexed = vector<bool>(_graphSuccessors->_nbNodes, false);
        //_isNodeRemoved = vector<bool>(_graphSuccessors->_nbNodes, false);
        //_isEdgeRemoved.clear();
        _bubbles.clear();
        
    }



    u_int32_t _nextUnitigIndex;


    //u_int32_t nodeIndex_to_unitigIndex(u_int32_t nodeIndex){
    //    return _unitigs[_nodeToUnitig[nodeIndex]]._index;
    //}

    //Unitig& nodeIndex_to_unitig(u_int32_t nodeIndex){
    //    return _unitigs[_nodeToUnitig[nodeIndex]];
    //}

    //float getNodeUnitigAbundance(u_int32_t nodeIndex){
    //    return _unitigs[_nodeToUnitig[nodeIndex]]._abundance;
    //}

    unordered_set<u_int32_t> _removedUnitigs;


    void addUnitigGraphEdges(){

        _logFile << "Add unitig graph edges" << endl;

        //cout << "a " << _startNode_to_unitigIndex.size() << endl;
        //_logFile << "Unitig graph nb nodes: " << _unitigGraph->_nodes.size() << endl;
        /*
        //for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){
        for(const auto& it : _startNode_to_unitigIndex){

            u_int32_t startNode = it.first;
            u_int32_t unitigIndex = it.second;

            //cout << "a \t" << startNode << " " << unitigIndex << endl;

            vector<u_int32_t> predecessors;
            getPredecessors(startNode, predecessors);
                

            for(u_int32_t nodeIndex : predecessors){
                u_int32_t toUnitigIndex = _endNode_to_unitigIndex[nodeIndex];
                dumpUnitigEdge(unitigIndex, toUnitigIndex, false);
                //_unitigGraph->addPredecessors(node->_unitigIndex, toUnitigIndex);
            }

        }
        */

        //cout << "b " << _startNode_to_unitigIndex.size() << endl;

        for(const auto& it : _endNode_to_unitigIndex){

            u_int32_t endNode = it.first;
            u_int32_t unitigIndex = it.second;

            //cout << "b \t" << endNode << " " << unitigIndex << endl;

            vector<u_int32_t> successors;
            getSuccessors(endNode, successors);
                
            //_logFile << "c" << endl;

            for(u_int32_t nodeIndex : successors){
                u_int32_t toUnitigIndex = _startNode_to_unitigIndex[nodeIndex];
                dumpUnitigEdge(unitigIndex, toUnitigIndex, true);
                //_unitigGraph->addSuccessor(node->_unitigIndex, toUnitigIndex);
            }

        }


    }

    unordered_set<u_int32_t> _unitigIndexToClean;

    void compact(bool rebuild, const vector<UnitigData>& unitigDatas){


        _logFile << "\tCompacting" << endl;
        
        if(!rebuild){
            //if(_unitigGraph != nullptr){
            //    delete _unitigGraph;
            //}

            
        }




        //vector<u_int32_t> nodes;


        _nbSingletonUnitigs = 0;
        _nextUnitigIndex = 0;
        _unitigs.clear();

        //for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
        //    if(!_isNodeRemoved[nodeIndex]){
        //        nodes.push_back(nodeIndex);
        //    }
        //    //for (u_int32_t nodeIndex : _isNodeValid2){
        //}
        //std::sort(nodes.begin(), nodes.end());




        //srand(time(NULL));



        
        UnitigFunctor functor(this, rebuild);
        size_t i=0;

        #pragma omp parallel num_threads(1)
        {

            UnitigFunctor functorSub(functor);

            while(true){

                u_int32_t nodeIndex = -1;
                bool isEof = false;

                #pragma omp critical
                {
                    
                    if(i >= _graphSuccessors->_nbNodes){
                        isEof = true;
                    }
                    else{
                        nodeIndex = i;//nodes[i];
                        //_logFile << nodeIndex << endl;
                    }

                    i += 2;
                }

                if(isEof) break;

                functorSub(nodeIndex);

            }
            

        }

        addUnitigGraphEdges();

        _logFile << "\tdone" << endl;
        
    }

    void dumpUnitigNode(u_int32_t unitigIndex, const vector<u_int32_t>& nodePath){

        u_int32_t size = nodePath.size();
        
        _unitigGraphFile_nodes.write((const char*)&size, sizeof(size));
        _unitigGraphFile_nodes.write((const char*)&nodePath[0], size * sizeof(u_int32_t));
        _unitigGraphFile_nodes.write((const char*)&unitigIndex, sizeof(unitigIndex));

    }

    void dumpUnitigEdge(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex, bool isSuccessor){

        if(isSuccessor){
            _unitigGraphFile_edges_successors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
            _unitigGraphFile_edges_successors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));
        }
        else{
            _unitigGraphFile_edges_predecessors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
            _unitigGraphFile_edges_predecessors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));
        }
    }

    u_int64_t _nbSingletonUnitigs;

	class UnitigFunctor {

		public:

		GraphSimplify* _graph;
        vector<Unitig>& _unitigs;
        bool _rebuild;

		UnitigFunctor(GraphSimplify* graph, bool rebuild) : _graph(graph), _unitigs(graph->_unitigs){
            _rebuild = rebuild;
		}

		UnitigFunctor(const UnitigFunctor& copy) : _graph(copy._graph), _unitigs(copy._graph->_unitigs){
            _rebuild = copy._rebuild;
		}

		~UnitigFunctor(){
		}

		void operator () (u_int32_t nodeIndex) {


            bool exist = false;
            #pragma omp critical(unitig)
            {
                if(_graph->_isNodeIndexIndexed[nodeIndex]) exist = true;
                //if(_nodeToUnitig.find(nodeIndex) != _nodeToUnitig.end()) exist = true;
            }



            if(exist) return;

            vector<u_int32_t> neighbors;
            bool dummy = false;


            u_int32_t startNode = nodeIndex;
            u_int32_t endNode = nodeIndex;

            u_int8_t isCircular = CONTIG_LINEAR;




            u_int32_t smallestNodeIndex = -1;

            //forward
            while(true){
                //_logFile << endNode << " " << nodeIndex << endl;

                if(endNode < smallestNodeIndex){
                    smallestNodeIndex = endNode;
                }

                _graph->getSuccessors(endNode, neighbors);


                //_logFile << "lala1: " << neighbors.size() << endl;
                if(neighbors.size() != 1) break;

                //_logFile << endNode << " " << neighbors[0] << endl;

                u_int32_t successor = neighbors[0];
                _graph->getPredecessors(successor, neighbors);
                
                if(neighbors.size() != 1) break;

                if(successor == nodeIndex){
                    endNode = successor;
                    isCircular = true;
                    //_logFile << "circular" << endl;
                    break; //Circular
                }


                endNode = successor;

            }

            if(!isCircular){

                //backward
                while(true){
                    _graph->getPredecessors(startNode, neighbors);
                    //_logFile << "loulou1: " << neighbors.size() << endl;
                    if(neighbors.size() != 1) break;


                    u_int32_t predecessor = neighbors[0];
                    _graph->getSuccessors(predecessor, neighbors);
                    //_logFile << "loulou2: " << neighbors.size() << endl;
                    if(neighbors.size() != 1) break;

                    startNode = predecessor;


                }

            }
            else{
                startNode = smallestNodeIndex;
                endNode = smallestNodeIndex;

            }


            //u_int32_t length = 0;
            u_int32_t node = startNode;


            size_t i=0;
            u_int32_t nbNodes = 0;
            u_int32_t lastNode = -1;

            vector<u_int32_t> nodes;
            vector<u_int32_t> nodesRC;
            //u_int32_t qualitySum = 0;
            //u_int32_t minQuality = -1;
            //u_int32_t maxQuality = 0;

            //vector<float> abundances;

            while(true){

                u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(node);


                //if(i == 0){
                //    length += _graph->_kminmerLengthMean;
                //}
                //else{
                //    length += _graph->_kminmerLengthMean - _graph->_kminmerOverlapMean;
                    
                //}

                //abundances.push_back(_graph->_nodeName_to_abundance[nodeName]);

                _graph->getSuccessors(node, neighbors);

                nodes.push_back(node);
                nodesRC.push_back(UnitigGraph::nodeIndex_toReverseDirection(node));
                nbNodes += 1;


                if(isCircular){
                    if(i == 0){ 
                        if(neighbors[0] == endNode) break; //Unitig with single node circular, or single node without edges
                    }
                    else{
                        if(node == endNode) break;
                    }
                }
                else{
                    if(node == endNode) break;
                }

                lastNode = node;
                node = neighbors[0];

                i += 1;


            }

            //if(isCircular){
            //    abundances.pop_back();
            //}


            std::reverse(nodesRC.begin(), nodesRC.end());
            //float median = Utils::compute_median_float(abundances);// _graph->compute_median_float(abundances);


            #pragma omp critical(unitig)
            {


                bool isValid = false; 
                
                for(u_int32_t nodeIndex : nodes){
                    if(!_graph->_isNodeIndexIndexed[nodeIndex]){
                        isValid = true;
                        break;
                    }
                }



                if(isValid){

                    u_int32_t unitigIndexRC = _graph->_nextUnitigIndex + 1;

                    for(u_int32_t nodeIndex : nodes){

                        _graph->_isNodeIndexIndexed[nodeIndex] = true;
                        _graph->_isNodeIndexIndexed[UnitigGraph::nodeIndex_toReverseDirection(nodeIndex)] = true;


                    }


                    if(isCircular && nodes.size() > 1){
                        nodes.pop_back();
                        nodesRC.erase(nodesRC.begin());
                    }

                    _graph->dumpUnitigNode(_graph->_nextUnitigIndex, nodes);
                    //_graph->dumpUnitigNode(_graph->_nextUnitigIndex, nodes);


                    u_int32_t startNode = nodes[0];
                    u_int32_t endNode = nodes[nodes.size()-1];
                    u_int32_t startNodeRC = nodesRC[0];
                    u_int32_t endNodeRC = nodesRC[nodesRC.size()-1];


                    _graph->_startNode_to_unitigIndex[startNode] = _graph->_nextUnitigIndex;
                    _graph->_endNode_to_unitigIndex[endNode] = _graph->_nextUnitigIndex;
                    _graph->_startNode_to_unitigIndex[startNodeRC] = _graph->_nextUnitigIndex+1;
                    _graph->_endNode_to_unitigIndex[endNodeRC] = _graph->_nextUnitigIndex+1;
                    
                    //_graph->_unitigGraph->addNode(nodes, median, length);
                    //_graph->_unitigGraph->addNode(nodesRC, median, length);

                    /*
                    u_int32_t endNode = nodes[nodes.size()-1];
                    vector<u_int32_t> successors;
                    _graph->getSuccessors(endNode, successors);
                        
                    u_int32_t startNode = nodes[0];
                    vector<u_int32_t> predecessors;
                    _graph->getPredecessors(startNode, predecessors);

                    if(successors.size() == 0 && predecessors.size() == 0){
                        _graph->_nbSingletonUnitigs += 1;
                        cout << _graph->_nbSingletonUnitigs << " " << _graph->_unitigGraph->_nodes.size() /2 << endl;
                    }
                    */

                    _graph->_nextUnitigIndex += 2;
                }

            }

            


        }

	};




    void getSuccessors(u_int32_t n, vector<u_int32_t>& successors){


        successors.clear();

        const vector<u_int32_t>& nodes = _graphSuccessors->_nodes[n];

        for(u_int32_t nodeIndex : nodes){

            successors.push_back(nodeIndex);
        }

    }


    void getPredecessors(u_int32_t n, vector<u_int32_t>& predecessors){

        predecessors.clear();
        

        vector<u_int32_t>& nodes = _graphSuccessors->_nodes[UnitigGraph::nodeIndex_toReverseDirection(n)];

        for(u_int32_t nodeIndex : nodes){


            predecessors.push_back(UnitigGraph::nodeIndex_toReverseDirection(nodeIndex));
        }

    }
    
    


    u_int32_t unitigIndex_toReverseDirection(u_int32_t unitigIndex){
        if(unitigIndex % 2 == 0){
            return unitigIndex+1;
        }
        else{
            return unitigIndex-1;
        }

    }



    float _prevAbudanceCutoff;

    struct LongNeighbor{
        u_int32_t _nodeIndex;
        float _abundance;
    };

    vector<u_int32_t> sccDebug;
    ofstream file_debug;
    int __loadState2_index = 0;


    unordered_set<u_int32_t> _possibleTips;
	ofstream _unitigGraphFile_nodes;
	ofstream _unitigGraphFile_edges_successors;
	ofstream _unitigGraphFile_edges_predecessors;

    void debug_writeGfaErrorfree(u_int32_t currentAbundance, float abundanceCutoff_min, u_int32_t nodeIndex_source, u_int64_t k, bool saveGfa, bool doesSaveState, bool doesLoadState, const vector<UnitigData>& unitigDatas, bool crushBubble, bool smallBubbleOnly, bool detectRoundabout, bool insertBubble, bool saveAllState, bool doesSaveUnitigGraph, MDBG* mdbg, size_t minimizerSize, size_t nbCores, bool useLocalAbundanceFilter, bool removeLongTips){

        _unitigGraphFile_nodes = ofstream(_outputDir + "/unitigGraph.nodes.bin");
        _unitigGraphFile_edges_successors = ofstream(_outputDir + "/unitigGraph.edges.successors.bin");
        _unitigGraphFile_edges_predecessors = ofstream(_outputDir + "/unitigGraph.edges.predecessors.bin");

        clear(0);

        compact(false, unitigDatas);

        _unitigGraphFile_nodes.close();
        _unitigGraphFile_edges_successors.close();
        _unitigGraphFile_edges_predecessors.close();

    }





	ofstream _file_readPath;


	unordered_map<u_int32_t, unordered_set<u_int64_t>> _unitigIndex_to_readIndexes;
    MDBG* _mdbg;


    

};


#endif

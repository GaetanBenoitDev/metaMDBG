

//#include "Graph.hpp"
//TODO:
// remove bubble, tips: consider length

#define PRINT_DEBUG_SIMPLIFICATION



#ifndef MDBG_METAG_GRAPHSIMPLIFY
#define MDBG_METAG_GRAPHSIMPLIFY

#include "Commons.hpp"
#include <limits>
//#include <experimental/filesystem>

//using namespace std;
//namespace fs = std::experimental::filesystem;

//class GraphSimplify;

struct NodeAbundance{
    u_int32_t _nodeIndex;
    double _abundance;
};

struct Bubble{
    vector<NodeAbundance> _nodes;
};


struct Unitig{
    u_int32_t _index;
    u_int32_t _startNode;
    u_int32_t _endNode;
    float _abundance;
    u_int32_t _length;
    u_int32_t _nbNodes;
    vector<u_int32_t> _nodes;
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
};

class GraphSimplify{

public:




    BiGraph* _graphSuccessors;
    //BiGraph* _graphPredecessors;
    string _outputDir;
    string _inputGfaFilename;
    vector<Unitig> _unitigs;
    //vector<u_int32_t> _nodeToUnitig;
    unordered_map<u_int32_t, u_int32_t> _nodeToUnitig;
    vector<u_int32_t> _nodeAbundances;
    vector<u_int32_t> _nodeLengths;
    //unordered_map<Destroy path: 2u_int32_t, u_int32_t> _nodeAbundances;
    //unordered_map<u_int32_t, u_int32_t> _nodeLengths;
    //vector<bool> _isNodeValid;
    unordered_set<u_int32_t> _isNodeValid2;
    vector<bool> _isBubble;
    unordered_set<DbgEdge, hash_pair> _isEdgeRemoved;
    unordered_set<DbgEdge, hash_pair> _isEdgeUnsupported;

    vector<Bubble> _bubbles;
    unordered_set<u_int32_t> _isNodeInvalid_tmp;
    //GraphSimplify::SuperbubbleSimplify* _superbubbleSimplify;
    //u_int32_t _lastAbundanceCutoff;
    unordered_map<u_int32_t, ComplexArea> _complexAreas_source; 
    unordered_map<u_int32_t, ComplexArea> _complexAreas_sink;
    unordered_set<u_int32_t> _isLongUnitig;

    GraphSimplify(const string& inputGfaFilename, const string& outputDir, u_int32_t nbNodes){
        _inputGfaFilename = inputGfaFilename;
        _outputDir = outputDir;

        _graphSuccessors = GfaParser::createBiGraph_lol(inputGfaFilename, true, nbNodes);
	    //_graphPredecessors = GfaParser::createBiGraph_lol(inputGfaFilename, false);
        _nodeAbundances.resize(_graphSuccessors->_nbNodes/2, false);
        _nodeLengths.resize(_graphSuccessors->_nbNodes/2, false);
        GfaParser::getNodeData(inputGfaFilename, _nodeAbundances, _nodeLengths);
        //clear(0);
        //_isNodeValid = vector<bool>(_graphSuccessors->_nbNodes, true);
        //_isBubble = vector<bool>(_graphSuccessors->_nbNodes, false);

        //_superbubbleSimplify = new SuperbubbleSimplify(this);

    }

    GraphSimplify(BiGraph* graph){
        _inputGfaFilename = "";
        _outputDir = "";

        _graphSuccessors = graph;
        
        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            _isNodeValid2.insert(nodeIndex);
        }

        _nodeAbundances = vector<u_int32_t>(_graphSuccessors->_nbNodes/2, 10000);
        _nodeLengths = vector<u_int32_t>(_graphSuccessors->_nbNodes/2, 10);
        


        //cout << "getSuccessors 564 debug to remove" << endl;


	    //_graphPredecessors = GfaParser::createBiGraph_lol(inputGfaFilename, false);
        //_nodeAbundances.resize(_graphSuccessors->_nbNodes/2, false);
        //_nodeLengths.resize(_graphSuccessors->_nbNodes/2, false);
        //GfaParser::getNodeData(inputGfaFilename, _nodeAbundances, _nodeLengths);
        //_isNodeValid = vector<bool>(_graphSuccessors->_nbNodes, true);
        //_isBubble = vector<bool>(_graphSuccessors->_nbNodes, false);

    }

    ~GraphSimplify(){
        delete _graphSuccessors;
        //delete _graphPredecessors;
    }

    SaveState _saveState;

    void saveState(){
        _saveState = {_isNodeValid2, _isBubble,_isEdgeRemoved, _bubbles, _nodeToUnitig, _unitigs};
    }

    void loadState(float abundanceCutoff_min){
        _isNodeValid2.clear();
        _isNodeInvalid_tmp.clear();
        //_isNodeValid2 = _saveState._isNodeValid2;
        _isBubble = _saveState._isBubble;
        _isEdgeRemoved = _saveState._isEdgeRemoved;
        _bubbles = _saveState._bubbles;
        _nodeToUnitig = _saveState._nodeToUnitig;
        _unitigs = _saveState._unitigs;

        /*
        for(u_int32_t nodeIndex : _saveState._isNodeValid2){
            if(_nodeAbundances[BiGraph::nodeIndex_to_nodeName(nodeIndex)] <= abundanceCutoff_min/3) continue;
            _isNodeValid2.insert(nodeIndex);
        }
        */
    }

    void clear(float abundanceCutoff_min){

        _isLongUnitig.clear();
        _isNodeInvalid_tmp.clear();
        _isNodeValid2.clear();

        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            if(_nodeAbundances[BiGraph::nodeIndex_to_nodeName(nodeIndex)] <= abundanceCutoff_min/3) continue;
            _isNodeValid2.insert(nodeIndex);
        }
        
        /*
        //_isNodeRemoved.resize(_graphSuccessors->_nbNodes, false);
        _isNodeValid = vector<bool>(_graphSuccessors->_nbNodes, true);
        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            if(_nodeAbundances[BiGraph::nodeIndex_to_nodeName(nodeIndex)] <= abundanceCutoff_min/3){
                _isNodeValid[nodeIndex] = false;
            }
        }*/
        
        _isBubble = vector<bool>(_graphSuccessors->_nbNodes, false);
        _isEdgeRemoved.clear();
        _bubbles.clear();
        _complexAreas_source.clear();
        _complexAreas_sink.clear();
        
        _isEdgeRemoved = _isEdgeUnsupported;
        //cout << _isEdgeRemoved.size() << endl;
        //_isNodeValid2.erase(BiGraph::nodeName_to_nodeIndex(1009227, false));
        //_isNodeValid2.erase(BiGraph::nodeName_to_nodeIndex(1009227, true));
    }

    /*
    GraphSimplify* clone(){
        GraphSimplify* clone = new GraphSimplify(_inputGfaFilename, _outputDir);

        clone->_graphSuccessors = _graphSuccessors;
        clone->_graphPredecessors = _graphPredecessors;
        clone->_nodeAbundances = _nodeAbundances;
        clone->_nodeLengths = _nodeLengths;
        clone->_isNodeRemoved = _isNodeRemoved;
        clone->_isEdgeRemoved = _isEdgeRemoved;

        return clone;
    }
    */

    void execute(float abundanceCutoff_min, u_int64_t k){

        clear(0);
        compact(false);

        bool dummy;
        unordered_set<u_int32_t> removedNodes;

        //cout << "to remove" << endl;
        //abundanceCutoff_min = 12;
        
        while(true){

            u_int64_t nbTipsRemoved = 0;
            u_int64_t nbErrorRemoved = 0;
            u_int64_t nbBubblesRemoved = 0;
            u_int64_t nbSuperbubblesRemoved = 0;
            u_int64_t nbSmallLoopRemoved = 0;

            bool isModification = false;

            cout << "1" << endl;
            while(true){
                cout << "lala" << endl;
                compact(true);
                cout << "loulou" << endl;
                nbErrorRemoved = removeLowAbundantNodes(abundanceCutoff_min);
                cout << "lili" << endl;
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb error removed: " << nbErrorRemoved << endl;
                #endif
                if(nbErrorRemoved == 0) break;
                isModification = true;
            }

            
            cout << "2" << endl;
            
            while(true){
                cout << "lala" << endl;
                compact(true);
                cout << "loulou" << endl;
                unordered_set<u_int32_t> isTips;
                nbTipsRemoved = tip(4*k, true, isTips);
                nbTipsRemoved = tip(4*k, false, isTips);
                cout << "lili" << endl;
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb tip removed: " << nbTipsRemoved << endl;
                #endif
                if(nbTipsRemoved == 0) break;
                isModification = true;
            }
            cout << "3" << endl;

            /*
            while(true){

                u_int64_t nbErrorRemoved = 0;

                compact();

                for(Unitig& unitig : _unitigs){
                    
                    cout << unitig._index << " " << unitig._nbNodes << endl;
                    if(unitig._nbNodes >= 2*k) continue;

                    vector<u_int32_t> successors;
                    getSuccessors_unitig(unitig._index, successors);

                    vector<u_int32_t> predecessors;
                    getPredecessors_unitig(unitig._index, predecessors);


                    float neighborAbundanceMin = std::numeric_limits<float>::max();
                    
                    for(u_int32_t unitigIndex : successors){
                        if(_unitigs[unitigIndex]._abundance < neighborAbundanceMin){
                            neighborAbundanceMin = _unitigs[unitigIndex]._abundance;
                        }
                    }
                    for(u_int32_t unitigIndex : predecessors){
                        if(_unitigs[unitigIndex]._abundance < neighborAbundanceMin){
                            neighborAbundanceMin = _unitigs[unitigIndex]._abundance;
                        }
                    }

                    
                    if(unitig._abundance < neighborAbundanceMin*0.5){
                        vector<u_int32_t> unitigNodes;
                        cout << "lala" << endl;
                        getUnitigNodes(unitig, unitigNodes);
                        for(u_int32_t node : unitigNodes){
                            _isNodeValid2.erase(node);
                        }
                        cout << "loulou" << endl;
                        isModification = true;
                        nbErrorRemoved += 1;
                    }

                }

                if(nbErrorRemoved == 0) break;
            }*/


            /*
            while(true){
                compact();
                nbSuperbubblesRemoved = superbubble(50000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb superbubble removed: " << nbSuperbubblesRemoved << endl;
                #endif
                if(nbSuperbubblesRemoved == 0) break;
                isModification = true;
            }
            */
            /*
            while(true){
                compact();
                nbBubblesRemoved = bubble(50000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb bubble removed: " << nbBubblesRemoved << endl;
                #endif
                if(nbBubblesRemoved == 0) break;
                isModification = true;
            }
            */



            cout << "4" << endl;
            
            while(true){
                compact(true);
                nbSmallLoopRemoved = removeSmallLoop(10000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb small loop removed: " << nbSmallLoopRemoved << endl;
                #endif
                if(nbSmallLoopRemoved == 0) break;
                isModification = true;
            }
            
            
            

            if(!isModification) break;
        }


        

        //exit(1);

        /*
        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            if(_isNodeValid2.find(nodeIndex) != _isNodeValid2.end()) continue;
            //if(_isNodeValid[nodeIndex]) continue;
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
            removedNodes.insert(nodeName);
        }

        cout << "REWRITE GFA WITHOUT NODe: faire l'inverse rewrite gfa WITH node et utiliser _isNodeValid2 direct" << endl;
        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_tmp.gfa", removedNodes, _isEdgeRemoved, _graphSuccessors);
        */

    }

    u_int64_t removeLowAbundantNodes(float abundanceCutoff_min){

        unordered_set<u_int32_t> removedUnitigs;
        
        u_int64_t nbRemoved = 0;

        if(abundanceCutoff_min == 0) return nbRemoved;
        
        unordered_set<u_int32_t> removedNodes;
        for(u_int32_t nodeIndex : _isNodeValid2){
            if(_unitigs[_nodeToUnitig[nodeIndex]]._abundance < abundanceCutoff_min){
                removedNodes.insert(nodeIndex);
                removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
            }
        }

        removeUnitigs(removedUnitigs);
        for(u_int32_t nodeIndex : removedNodes){
            _isNodeValid2.erase(nodeIndex);
            //_removedFrom[nodeIndex] = 5;
        }

        /*
        for (auto it = _isNodeValid2.begin(); it != _isNodeValid2.end(); ){

            u_int32_t nodeIndex = *it;

            //cout << (_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()) << endl;
            if(_unitigs[_nodeToUnitig[nodeIndex]]._abundance < abundanceCutoff_min){
                removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
                it = _isNodeValid2.erase(it);
                nbRemoved += 1;
            }
            else{
                ++it;
            }

        }
        */

        /*
        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            //cout << nodeIndex << " " <<_graphSuccessors->_nbNodes << endl;
            if(!_isNodeValid[nodeIndex]) continue;
            
            //cout << "lala: " << _nodeToUnitig[nodeIndex] << endl;
            //cout << _unitigs[_nodeToUnitig[nodeIndex]]._abundance << endl;

            if(_unitigs[_nodeToUnitig[nodeIndex]]._abundance < abundanceCutoff_min){
                _isNodeValid[nodeIndex] = false;
                nbRemoved += 1;
            }

        }*/
        
        //cout << "miam1" << endl;
        //removeUnitigs(removedUnitigs);
        //cout << "miam2" << endl;

        /*
        for(u_int32_t unitigIndex : removedUnitigs){
            vector<u_int32_t> unitigNodes; 
            getUnitigNodes(_unitigs[unitigIndex], unitigNodes);
            for(u_int32_t node : unitigNodes){
                _isNodeValid2.erase(node);
                _isNodeValid2.erase(nodeIndex_toReverseDirection(node));
            }
        }
        */

        return nbRemoved;
    }

    //unordered_map<u_int32_t, u_int32_t> _removedFrom;

    u_int64_t tip(float maxLength, bool indexingTips, unordered_set<u_int32_t>& isTips){

        unordered_set<u_int32_t> removedNodes;

        u_int64_t nbRemoved = 0;

        vector<u_int32_t> neighbors;

        //vector<bool> isVisited(_graphSuccessors->_nbNodes, false);
        bool dummy = false;

        for(Unitig& unitig : _unitigs){
            if(unitig._startNode == -1) continue;


            //if(unitig._length > maxLength) continue;
            if(unitig._nbNodes >= maxLength) continue;
            //if(_isNodeValid2.find(unitig._startNode) == _isNodeValid2.end()) continue; //already removed


            //if(isVisited[unitig._startNode]) continue;
            //if(isVisited[unitig._endNode]) continue;

            getPredecessors(unitig._startNode, 0, neighbors);
            if(indexingTips){
                if(neighbors.size() == 0) continue;
            }
            else{

                u_int64_t nbPredecessors = 0;
                for(u_int32_t nodeIndex : neighbors){
                    u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                    if(isTips.find(unitigIndex) == isTips.end()) nbPredecessors += 1;
                }
                
                if(nbPredecessors != 1) continue;
            }



            getSuccessors(unitig._endNode, 0, neighbors);
            if(neighbors.size() > 0) continue;


            //isVisited[unitig._startNode] = true;
            //isVisited[unitig._endNode] = true;

            /*
            getUnitigNodes(unitig, unitigNodes);
            for(u_int32_t node : unitigNodes){
                _isNodeValid2.erase(node);
                _isNodeValid2.erase(nodeIndex_toReverseDirection(node));

                //_isNodeValid[node] = false;
            }
            */

            if(indexingTips){
                isTips.insert(unitig._index);
                isTips.insert(unitigIndex_toReverseDirection(unitig._index));
                continue;
            }

            vector<u_int32_t> unitigNodes; 
            getUnitigNodes(unitig, unitigNodes);
            for(u_int32_t node : unitigNodes){
                //if(_isNodeValid2.find(node) == _isNodeValid2.end()){
                //    cout << "omg1" << endl;
                //    getchar();
                //}
                //if(_isNodeValid2.find(nodeIndex_toReverseDirection(node)) == _isNodeValid2.end()){
                //    cout << (_removedFrom.find(nodeIndex_toReverseDirection(node)) != _removedFrom.end()) << endl;
                //    cout << _removedFrom[nodeIndex_toReverseDirection(node)] << endl;
                //    cout << "omg2" << endl;
                //    getchar();
                //}
                removedNodes.insert(node);
                removedNodes.insert(nodeIndex_toReverseDirection(node));
            }
            //cout << "blabla2" << endl;

            //removedUnitigs.insert(unitig._index);
            //removedUnitigs.insert(unitigIndex_toReverseDirection(unitig._index));

            /*
            _availableUnitigIndexes.insert(nodeIndex_to_unitigIndex(node));
            _availableUnitigIndexes.insert(nodeIndex_to_unitigIndex(nodeIndex_toReverseDirection(node)));

            getSuccessors_unitig(unitig._index, neighbors);
            for(u_int32_t unitigIndex : neighbors){
                _rebuildInvalidUnitigs.insert(unitigIndex);
                _rebuildInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex));
            }
            getPredecessors_unitig(unitig._index, neighbors);
            for(u_int32_t unitigIndex : neighbors){
                _rebuildInvalidUnitigs.insert(unitigIndex);
                _rebuildInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex));
            }
            */
            

            #ifdef PRINT_DEBUG_SIMPLIFICATION
                cout << "\tTip: " << _graphSuccessors->nodeIndex_to_nodeName(unitig._endNode, dummy) << " " << unitig._length << endl;
            #endif 

            nbRemoved += 1;
        }

        if(indexingTips) return 0;

        unordered_set<u_int32_t> removedUnitigs;
        for(u_int32_t nodeIndex : removedNodes){
            //if(_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()){
            //    cout << "KOUERK" << endl;
            //}
            removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
        }
        removeUnitigs(removedUnitigs);
        for(u_int32_t nodeIndex : removedNodes){
            //_removedFrom[nodeIndex] = 1;
            _isNodeValid2.erase(nodeIndex);
        }


        return nbRemoved;
    }

    void removeUnitigs(unordered_set<u_int32_t>& unitigs){

        _removedUnitigs.clear();
        unordered_set<u_int32_t> allInvalidUnitigs;

        for(u_int32_t unitigIndex : unitigs){
            collectInvalidUnitig(unitigIndex, allInvalidUnitigs);
        }

        for(u_int32_t unitigIndex : allInvalidUnitigs){
            disableUnitig(unitigIndex);
        }


    }

    void collectInvalidUnitig(u_int32_t unitigIndex, unordered_set<u_int32_t>& allInvalidUnitigs){
        //_removedUnitigs.insert(unitigIndex);
        collectInvalidUnitig2(unitigIndex, allInvalidUnitigs);

        vector<u_int32_t> neighbors;

        getSuccessors_unitig(unitigIndex, 0, neighbors);
        
        for(u_int32_t unitigIndex_nn : neighbors){
            collectInvalidUnitig2(unitigIndex_nn, allInvalidUnitigs);
            //_rebuildInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex));
        }
        getPredecessors_unitig(unitigIndex, 0, neighbors);
        for(u_int32_t unitigIndex_nn : neighbors){
            collectInvalidUnitig2(unitigIndex_nn, allInvalidUnitigs);
            //_rebuildInvalidUnitigs.insert(unitigIndex);
            //_rebuildInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex));
        }
    }

    void collectInvalidUnitig2(u_int32_t unitigIndex, unordered_set<u_int32_t>& allInvalidUnitigs){
        allInvalidUnitigs.insert(unitigIndex);
        //allInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex));
        vector<u_int32_t> neighbors;

        getSuccessors_unitig(unitigIndex, 0, neighbors);
        for(u_int32_t unitigIndex_nn : neighbors){
            allInvalidUnitigs.insert(unitigIndex_nn);
            //allInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex_nn));
        }
        getPredecessors_unitig(unitigIndex, 0, neighbors);
        for(u_int32_t unitigIndex_nn : neighbors){
            allInvalidUnitigs.insert(unitigIndex_nn);
            //allInvalidUnitigs.insert(unitigIndex_toReverseDirection(unitigIndex_nn));
        }

    }

    void disableUnitig(u_int32_t unitigIndex){
        //cout << _unitigs[unitigIndex]._startNode << endl;
        //cout << _unitigs[unitigIndex]._endNode << endl;
        //if(_unitigs[unitigIndex]._startNode == -1) return;
        //_rebuildInvalidUnitigs.insert(unitigIndex);
        //_availableUnitigIndexes.insert(unitigIndex);


        vector<u_int32_t> unitigNodes; 
        getUnitigNodes(_unitigs[unitigIndex], unitigNodes);
        for(u_int32_t nodeIndex : unitigNodes){
            _nodeToUnitig.erase(nodeIndex);
        }


        _unitigs[unitigIndex]._startNode = -1;
        _unitigs[unitigIndex]._endNode = -1;
        _removedUnitigs.insert(unitigIndex);
        //_unitigs[unitigIndex]._nodes;
    }

    u_int64_t superbubble(u_int64_t maxLength){


        unordered_set<u_int32_t> removedNodes;
        //unordered_set<u_int32_t> removedUnitigs;

        #ifdef PRINT_DEBUG_SIMPLIFICATION
            cout << "\tStart superbubble" << endl;
        #endif 

        u_int64_t nbRemoved = 0;

        //cout << "Super bubble!" << endl;
        for(Unitig& unitig : _unitigs){
            if(unitig._startNode == -1) continue;
            //if(BiGraph::nodeIndex_to_nodeName(unitig._endNode) != 1307) continue;
            //cout << BiGraph::nodeIndex_to_nodeName(unitig._endNode) << endl;

            if(_isBubble[unitig._endNode]) continue;
            //if(_isBubble[nodeIndex_toReverseDirection(unitig._endNode)]) continue;
            
            u_int32_t unitigIndex_exit = detectSuperbubble(unitig._index, maxLength);
            if(unitigIndex_exit == -1) continue;

            #ifdef PRINT_DEBUG_SIMPLIFICATION
                cout << "\tSuperbubble: " << BiGraph::nodeToString(unitig._endNode) << " " << BiGraph::nodeToString(_unitigs[unitigIndex_exit]._startNode) << endl;
                cout << "\tSuperbubble: " << BiGraph::nodeToString(nodeIndex_toReverseDirection(unitig._endNode)) << " " << BiGraph::nodeToString(nodeIndex_toReverseDirection(_unitigs[unitigIndex_exit]._startNode)) << endl;
                //getchar();
            #endif

            //cout << "Found superbubble: " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_exit]._startNode) << endl;
            simplifySuperbubble(unitig._index, unitigIndex_exit, removedNodes);


            nbRemoved += 1;
        }

        unordered_set<u_int32_t> removedUnitigs;
        for(u_int32_t nodeIndex : removedNodes){
            removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
        }
        removeUnitigs(removedUnitigs);
        for(u_int32_t nodeIndex : removedNodes){
            //_removedFrom[nodeIndex] = 2;
            _isNodeValid2.erase(nodeIndex);
        }


        //exit(1);
        return nbRemoved;
    }

    u_int32_t _nextUnitigIndex;

    //https://arxiv.org/pdf/1307.7925.pdf
    u_int32_t detectSuperbubble(u_int32_t unitigIndex_source, u_int64_t maxLength){

        unordered_set<u_int32_t> isVisited;
        unordered_set<u_int32_t> seen;
        unordered_map<u_int32_t, u_int64_t> pathLength;
        vector<u_int32_t> queue;

        queue.push_back(unitigIndex_source);
        pathLength[unitigIndex_source] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];
            //cout << "\tVisited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._endNode) << endl;
            queue.pop_back();

            if(pathLength[v] > maxLength) continue;

            isVisited.insert(v);
            if(seen.find(v) != seen.end()){
                seen.erase(v);
            }

            vector<u_int32_t> successors;
            getSuccessors_unitig(v, 0, successors);
            if(successors.size() == 0) return -1; //abort tip

            for(u_int32_t u : successors){
                if(_isBubble[_unitigs[u]._startNode]) return -1;
                if(_isBubble[nodeIndex_toReverseDirection(_unitigs[u]._startNode)]) return -1;
                if(u == unitigIndex_source) return -1; //cycle including s

                if(isVisited.find(u) == isVisited.end()){
                    seen.insert(u);
                    pathLength[u] = pathLength[v] + _unitigs[u]._length;
                }

            }

            for(u_int32_t u : successors){

                //cout << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

                vector<u_int32_t> predecessors;
                getPredecessors_unitig(u, 0, predecessors);
                bool allPredecessorsAreVisited = true;
                for(u_int32_t p : predecessors){
                    if(isVisited.find(p) == isVisited.end()){
                        allPredecessorsAreVisited = false;
                        break;
                    }
                }

                if(allPredecessorsAreVisited){
                    //cout << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                    queue.push_back(u);
                }

                //cout << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                    u_int32_t t = queue[0];
                    vector<u_int32_t> successors_t;
                    getSuccessors_unitig(t, 0, successors_t);
                    if(std::find(successors_t.begin(), successors_t.end(), unitigIndex_source) == successors_t.end()){
                        return t;
                    }
                    else{
                        return -1; // cycle including s
                    }

                }


            }

        }

        return -1;
    }




    void simplifySuperbubble(u_int32_t source_unitigIndex, u_int32_t sink_unitigIndex, unordered_set<u_int32_t>& removedNodes){

        vector<u_int32_t> unitigIndexes;
        collectNodes_betweenSourceSink_unitig(source_unitigIndex, sink_unitigIndex, unitigIndexes, 0, {}, true);
        //if(unitigIndexes.size() == 0) return;

        unordered_set<u_int32_t> keepNodes;
        vector<u_int32_t> nodes;
        //vector<u_int32_t> removedNodes;
        Bubble bubble;
        Bubble bubbleRC;

        //Choose one path in the superbubble
        u_int32_t unitigIndex = source_unitigIndex;
        vector<u_int32_t> path = {source_unitigIndex};
        while(true){
            vector<u_int32_t> successors;
            getSuccessors_unitig(unitigIndex, 0, successors);

            float maxAbundance = 0;
            u_int32_t maxV = -1;
            for(u_int32_t v : successors){
                if(_unitigs[v]._abundance > maxAbundance){
                    maxAbundance = _unitigs[v]._abundance;
                    maxV = v;
                }
            }

            unitigIndex = maxV;
            path.push_back(unitigIndex);

            if(unitigIndex == sink_unitigIndex) break;
        }
        
		//shortestPath_unitig(source_unitigIndex, sink_unitigIndex, path, false, false);
        for(u_int32_t unitigIndex : path){
            getUnitigNodes(_unitigs[unitigIndex], nodes);
            for(u_int32_t nodeIndex : nodes){
                //cout << "\tKeep nodes: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
                keepNodes.insert(nodeIndex);
                //keepNodes.insert(nodeIndex_toReverseDirection(nodeIndex));
            }
        }
        
        
        //Mark all nodes as bubble
        for(u_int32_t unitigIndex : unitigIndexes){
            getUnitigNodes(_unitigs[unitigIndex], nodes);
            for(u_int32_t nodeIndex : nodes){
                _isBubble[nodeIndex] = true;
                _isBubble[nodeIndex_toReverseDirection(nodeIndex)] = true;
                if(keepNodes.find(nodeIndex) != keepNodes.end()) continue;
                //if(keepNodes.find(nodeIndex_toReverseDirection(nodeIndex)) != keepNodes.end()) continue;
                //removedNodes.push_back(nodeIndex);
                bubble._nodes.push_back({nodeIndex, getNodeUnitigAbundance(nodeIndex)});
                bubbleRC._nodes.push_back({nodeIndex_toReverseDirection(nodeIndex), getNodeUnitigAbundance(nodeIndex_toReverseDirection(nodeIndex)) });
                //cout << "\tRemove nodes: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
            }
        }
        
        std::reverse(bubbleRC._nodes.begin(), bubbleRC._nodes.end());

        for(NodeAbundance& nodeIndex : bubble._nodes){ //removed nodes
             //cout << "1" << endl;
            //_isNodeValid2.erase(nodeIndex);
            removedNodes.insert(nodeIndex._nodeIndex);
            removedNodes.insert(nodeIndex_toReverseDirection(nodeIndex._nodeIndex));
            //_isNodeValid2.erase(nodeIndex_toReverseDirection(nodeIndex));
             //cout << "2" << endl;
        }

        if(bubble._nodes.size() == 0){ 
            cout << "1" << endl;
            cout << source_unitigIndex << endl;
            cout << BiGraph::nodeIndex_to_nodeName(_unitigs[source_unitigIndex]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[sink_unitigIndex]._startNode) << endl;
            getchar();
        }
        //if(bubbleRC._nodes.size() == 0){ 
        //    cout << "2" << endl;
        //    getchar();
        //}
        _bubbles.push_back(bubble);
        _bubbles.push_back(bubbleRC);
        /*
        Bubble bubble;
        for(u_int32_t node : unitigNodes){
            _isNodeValid2.erase(node);
            //_isNodeValid[node] = false;
            bubble._nodes.push_back(node);
            //removedNodes.insert(node);
        }
        _bubbles.push_back(bubble);
        */
	}

/*
push s into S
2: repeat
3: pick out an arbitrary v ∈ S
4: label v as visited
5: if v does not have a child then
    6: abort // tip
7: for u in v’s children do
    8: if u = s then
    9: abort // cycle including s
    10: label u as seen
    11: if all of u’s parents are visited then
    12:     push u into S
    13: if only one vertex t is left in S and no other vertex is seen then
    14:     if edge (t, s) does not exist then
    15:         return t
    16:     else
    17:         abort // cycle including s
    18: until |S| = 0
    }
    */

    u_int64_t bubble(float maxLength){

        unordered_set<u_int32_t> removedNodes;

        u_int64_t nbBubblesRemoved = 0;
        //vector<u_int32_t> removedUnitigs;

        //vector<bool> isVisited(_graphSuccessors->_nbNodes, false);
        vector<u_int32_t> unitigNodes;
        bool dummy = false;

        for(Unitig& unitig : _unitigs){
            if(unitig._startNode == -1) continue;
            //cout << "-------------" << endl;



            //u_int32_t startNodeName = _graphSuccessors->nodeIndex_to_nodeName(unitig._startNode, dummy);
            //u_int32_t visitedNodeIndex = unitig._startNode;
            //if(isVisited[unitig._startNode]) continue;
            if(_isBubble[unitig._startNode]) continue;
            if(_isBubble[nodeIndex_toReverseDirection(unitig._startNode)]) continue;
            

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
                if(_isBubble[utg_2._startNode]) continue;
                if(_isBubble[nodeIndex_toReverseDirection(utg_2._startNode)]) continue;
                //if(isVisited[utg_2._startNode]) continue;
                //isVisited[startNodeName] = true;



                for(size_t j=i+1; j<neighbors_utg1.size(); j++){
                    Unitig& utg_3 = _unitigs[_nodeToUnitig[neighbors_utg1[j]]];
                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy);
                    //if(isVisited[utg_3._startNode]) continue;
                    if(_isBubble[utg_3._startNode]) continue;
                    if(_isBubble[nodeIndex_toReverseDirection(utg_3._startNode)]) continue;
                    //isVisited[startNodeName] = true;


                    vector<u_int32_t> neighbors_utg2;
                    getPredecessors(utg_2._startNode, 0, neighbors_utg2);
                    if(neighbors_utg2.size() != 1) continue;
                    getSuccessors(utg_2._endNode, 0, neighbors_utg2);
                    if(neighbors_utg2.size() != 1) continue;
                    
                    Unitig& utg_4 = _unitigs[_nodeToUnitig[neighbors_utg2[0]]];

                    //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy);
                    //if(isVisited[utg_4._startNode]) continue;
                    if(_isBubble[utg_4._startNode]) continue;
                    if(_isBubble[nodeIndex_toReverseDirection(utg_4._startNode)]) continue;
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
                        cout << "\tBubble: " << _graphSuccessors->nodeIndex_to_nodeName(utg_1._endNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(utg_4._startNode, dummy) << " " << utg_2._length << " " << utg_3._length << endl;
                    #endif

                    
                    //Memorise bubble unitigs
                    getUnitigNodes(utg_3, unitigNodes);
                    for(u_int32_t node : unitigNodes){
                        _isBubble[node] = true;
                        _isBubble[nodeIndex_toReverseDirection(node)] = true;
                    }
                    getUnitigNodes(utg_2, unitigNodes);
                    for(u_int32_t node : unitigNodes){
                        _isBubble[node] = true;
                        _isBubble[nodeIndex_toReverseDirection(node)] = true;
                    }
                    
                    
                    //remove bubble
                    if(utg_2._abundance > utg_3._abundance){
                        getUnitigNodes(utg_3, unitigNodes);

                        //isVisited[utg_3._startNode] = true;
                        //isVisited[utg_3._endNode] = true;
                        //isVisited[nodeIndex_toReverseDirection(utg_3._startNode)] = true;
                        //isVisited[nodeIndex_toReverseDirection(utg_3._endNode)] = true;


                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._startNode, dummy);
                        //isVisited[startNodeName] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_3._endNode, dummy);
                        //isVisited[startNodeName] = true;
                    }
                    else{
                        getUnitigNodes(utg_2, unitigNodes);

                        //isVisited[utg_2._startNode] = true;
                        //isVisited[utg_2._endNode] = true;
                        //isVisited[nodeIndex_toReverseDirection(utg_2._startNode)] = true;
                        //isVisited[nodeIndex_toReverseDirection(utg_2._endNode)] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_2._startNode, dummy);
                        //isVisited[startNodeName] = true;
                        //startNodeName = _graphSuccessors->nodeIndex_to_nodeName(utg_2._endNode, dummy);
                        //isVisited[startNodeName] = true;
                    }
                    

                    Bubble bubble;
                    Bubble bubbleRC;
                    for(u_int32_t node : unitigNodes){
                        removedNodes.insert(node);
                        removedNodes.insert(nodeIndex_toReverseDirection(node));
                        //_isNodeValid2.erase(node);
                        //_isNodeValid2.erase(nodeIndex_toReverseDirection(node));
                        //_isNodeValid[node] = false;
                        bubble._nodes.push_back({node, getNodeUnitigAbundance(node)});
                        bubbleRC._nodes.push_back({nodeIndex_toReverseDirection(node), getNodeUnitigAbundance(nodeIndex_toReverseDirection(node))});


                        //removedUnitigs.insert(nodeIndex_to_unitigIndex(node));
                        //removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex_toReverseDirection(node)));
                        //removedNodes.insert(node);
                    }
                    std::reverse(bubbleRC._nodes.begin(), bubbleRC._nodes.end());
                    _bubbles.push_back(bubble);
                    _bubbles.push_back(bubbleRC);

                    
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

        unordered_set<u_int32_t> removedUnitigs;
        for(u_int32_t nodeIndex : removedNodes){
            removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
        }
        removeUnitigs(removedUnitigs);
        for(u_int32_t nodeIndex : removedNodes){
            _isNodeValid2.erase(nodeIndex);
            //_removedFrom[nodeIndex] = 3;
        }


        return nbBubblesRemoved;
    }


    u_int64_t removeSmallLoop(float maxLength){

        unordered_set<u_int32_t> removedUnitigs;

        //cout << "removing small loops" << endl;
        u_int64_t nbRemoved = 0;

        //vector<u_int32_t> unitigNodes;

        //vector<bool> isVisited(_graphSuccessors->_nbNodes, false);
        bool dummy = false;

        for(Unitig& unitig_source : _unitigs){
            if(unitig_source._startNode == -1) continue;


            //if(isVisited[unitig_source._startNode]) continue;
            //if(isVisited[unitig_source._endNode]) continue;


            vector<u_int32_t> neighbors;
            getSuccessors(unitig_source._endNode, 0, neighbors);

            /*
            if(unitig_source._endNode == _graphSuccessors->nodeName_to_nodeIndex(1997, dummy)){
                cout << "lala " << neighbors.size() << endl;
                for(u_int32_t nn : neighbors){
                    cout << _graphSuccessors->nodeToString(nn) << endl;
                }
            }
            */

            if(neighbors.size() < 2) continue;

            bool foundSmallLoop = false;

            for(size_t i=0; i<neighbors.size(); i++){

                Unitig& unitig_dest = nodeIndex_to_unitig(neighbors[i]);

                for(size_t j=0; j<neighbors.size(); j++){
                    if(i == j) continue;

                    Unitig& unitig_loop = nodeIndex_to_unitig(neighbors[j]);
                    if(unitig_loop._length > maxLength) continue;



                    vector<u_int32_t> predecessors;
                    getPredecessors(unitig_loop._startNode, 0, predecessors);
                    if(predecessors.size() != 2) continue;
                    vector<u_int32_t> successors;
                    getSuccessors(unitig_loop._endNode, 0, successors);
                    if(successors.size() != 2) continue;

                    /*
                    if(unitig_source._endNode == _graphSuccessors->nodeName_to_nodeIndex(1997, dummy)){
                        cout << "Pred: " << predecessors.size() << endl;
                        for(u_int32_t nn : predecessors){
                            cout << _graphSuccessors->nodeToString(nn) << endl;
                        }
                    }
                    
                    if(unitig_source._endNode == _graphSuccessors->nodeName_to_nodeIndex(1997, dummy)){
                        cout << "Succ: " << successors.size() << endl;
                        for(u_int32_t nn : successors){
                            cout << _graphSuccessors->nodeToString(nn) << endl;
                        }
                    }*/

                    if((nodeIndex_to_unitigIndex(predecessors[0]) == unitig_loop._index && nodeIndex_to_unitigIndex(predecessors[1]) == unitig_source._index) || (nodeIndex_to_unitigIndex(predecessors[0]) == unitig_source._index && nodeIndex_to_unitigIndex(predecessors[1]) == unitig_loop._index)){
                        if((nodeIndex_to_unitigIndex(successors[0]) == unitig_loop._index && nodeIndex_to_unitigIndex(successors[1]) == unitig_dest._index) || (nodeIndex_to_unitigIndex(successors[0]) == unitig_dest._index && nodeIndex_to_unitigIndex(successors[1]) == unitig_loop._index)){
                            
                            #ifdef PRINT_DEBUG_SIMPLIFICATION
                                cout << "Small loop: " <<  _graphSuccessors->nodeToString(unitig_loop._startNode) << " " << _graphSuccessors->nodeToString(unitig_loop._endNode)  << endl;
                            #endif
                            
                            removedUnitigs.insert(unitig_source._index);
                            removedUnitigs.insert(unitig_dest._index);
                            removedUnitigs.insert(unitig_loop._index);
                            removedUnitigs.insert(unitigIndex_toReverseDirection(unitig_source._index));
                            removedUnitigs.insert(unitigIndex_toReverseDirection(unitig_dest._index));
                            removedUnitigs.insert(unitigIndex_toReverseDirection(unitig_loop._index));

                            DbgEdge edge = {unitig_source._endNode, unitig_dest._startNode};
                            edge = edge.normalize();
                            _isEdgeRemoved.insert(edge);
                            edge = {unitig_loop._endNode, unitig_loop._startNode};
                            edge = edge.normalize();
                            _isEdgeRemoved.insert(edge);

                            
                            edge = {nodeIndex_toReverseDirection(unitig_dest._startNode), nodeIndex_toReverseDirection(unitig_source._endNode)};
                            edge = edge.normalize();
                            _isEdgeRemoved.insert(edge);
                            edge = {nodeIndex_toReverseDirection(unitig_loop._startNode), nodeIndex_toReverseDirection(unitig_loop._endNode)};
                            edge = edge.normalize();
                            _isEdgeRemoved.insert(edge);
                            

                            nbRemoved += 1;
                            foundSmallLoop = true;
                            break;
                            //isVisited[unitig_source._startNode] = true;
                            //isVisited[unitig_source._destNode] = true;
                        }
                    }

                }

                if(foundSmallLoop) break;
            }

            /*
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
            */

            //#ifdef PRINT_DEBUG_SIMPLIFICATION
            //    cout << "\tSmall loop: " << _graphSuccessors->nodeIndex_to_nodeName(unitig._endNode, dummy) << endl;
            //#endif 

            //nbRemoved += 1;
        }

        removeUnitigs(removedUnitigs);

        return nbRemoved;
    }

    u_int32_t nodeIndex_to_unitigIndex(u_int32_t nodeIndex){
        return _unitigs[_nodeToUnitig[nodeIndex]]._index;
    }

    Unitig& nodeIndex_to_unitig(u_int32_t nodeIndex){
        return _unitigs[_nodeToUnitig[nodeIndex]];
    }

    float getNodeUnitigAbundance(u_int32_t nodeIndex){
        return _unitigs[_nodeToUnitig[nodeIndex]]._abundance;
    }

    unordered_set<u_int32_t> _removedUnitigs;

    void compact(bool rebuild){

        //if(rebuild && _rebuildInvalidUnitigs.size() == 0) return;

        #ifdef PRINT_DEBUG_SIMPLIFICATION
            cout << "\tCompacting" << endl;
        #endif

        if(rebuild){

            for(u_int32_t unitigIndex : _removedUnitigs){
                //cout << unitigIndex << endl;
                vector<u_int32_t> unitigNodes; 
                getUnitigNodes(_unitigs[unitigIndex], unitigNodes);
                for(u_int32_t nodeIndex : unitigNodes){
                    if(_isNodeValid2.find(nodeIndex) == _isNodeValid2.end()) continue;
                    computeUnitig(nodeIndex);
                }
            }
            //for(u_int32_t invalidUnitigIndex : _rebuildInvalidUnitigs){
            //    Unitig& unitig = _unitigs[invalidUnitigIndex];
            //    unitig._startNode = -1;
            //}
        }
        else{
            _nextUnitigIndex = 0;
            _unitigs.clear();
            _nodeToUnitig.clear();

            for (u_int32_t nodeIndex : _isNodeValid2){
                computeUnitig(nodeIndex);
            }
        }


        //_nodeToUnitig.resize(_graphSuccessors->_nbNodes, -1);
        //ofstream outfile(_outputDir + "/" + "tmp_compacted_graph.gfa");

        //vector<u_int32_t> nodeToUnitig(_graphSuccessors->_nbNodes, -1);
        //vector<Unitig> unitigs;
        //unordered_set<u_int32_t> nodes;
        //unordered_map<u_int32_t, u_int32_t> startNodes;
        //unordered_map<u_int32_t, u_int32_t> endNodes;

        //u_int64_t unitigIndex = 0;
        
        //vector<bool> isVisited(_graphSuccessors->_nbNodes, false);




        //size_t nodeIndex = _graphSuccessors->nodeName_to_nodeIndex(720, true);

        //for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){


        u_int32_t unitigIndex = 0;
        #ifdef PRINT_DEBUG_SIMPLIFICATION
            cout << "Nb unitigs: " << unitigIndex << endl;
        #endif

        cout << "Nb unitigs: " << _unitigs.size() << endl;
        //cout << "Unitigs: " << _unitigs.size() << " " << _nodeToUnitig.size() << endl;
        //getchar();
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

    void computeUnitig(u_int32_t nodeIndex){
        if(_nodeToUnitig.find(nodeIndex) != _nodeToUnitig.end()) return;
        
        vector<u_int32_t> neighbors;
        bool dummy = false;

        //cout << nodeIndex << " " << _isNodeValid2.size() << endl;
        //u_int32_t nodeIndex = *it;
        //cout << nodeIndex << " " << _isNodeValid2.size() << endl;
        //if(nodeIndex % 10000 == 0){
        //    cout << nodeIndex << " " << _graphSuccessors->_nbNodes << endl;
        //}
        //if(nodeIndex % 2 == 1) continue;

        //if(!_isNodeValid[nodeIndex]) continue;
        //if(isVisited[nodeIndex]) continue;

        //isVisited[nodeIndex] = true;
        //u_int32_t nodeName = _graphSuccessors->nodeName_to_nodeIndex(7701, true);

        //cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
        u_int32_t startNode = nodeIndex;
        u_int32_t endNode = nodeIndex;

        bool isCircular = false;

        int lala1 = 0;
        int lala2 = 0;

        //forward
        while(true){
            getSuccessors(endNode, 0, neighbors);
            //cout << "lala1: " << neighbors.size() << endl;
            if(neighbors.size() != 1) break;
            if(neighbors[0] == nodeIndex){
                endNode = neighbors[0];
                isCircular = true;
                //cout << "circular" << endl;
                break; //Circular
            }
            //cout << endNode << " " << neighbors[0] << endl;

            u_int32_t successor = neighbors[0];
            getPredecessors(successor, 0, neighbors);
            //cout << "lala2: " << neighbors.size() << endl;
            if(neighbors.size() != 1) break;

            //cout << successor << " " << neighbors[0] << endl;

            //cout << "1: " << BiGraph::nodeIndex_to_nodeName(successor) << endl;
            endNode = successor;
            lala1 += 1;


        }

        if(!isCircular){

            //backward
            while(true){
                getPredecessors(startNode, 0, neighbors);
                //cout << "loulou1: " << neighbors.size() << endl;
                if(neighbors.size() != 1) break;
                //if(neighbors[0] == nodeIndex){
                //    startNode = neighbors[0];
                //    break; //Circular, todo: don't need to backward if circular
                //}

                //cout << startNode << " " << neighbors[0] << endl;

                u_int32_t predecessor = neighbors[0];
                getSuccessors(predecessor, 0, neighbors);
                //cout << "loulou2: " << neighbors.size() << endl;
                if(neighbors.size() != 1) break;

                //cout << predecessor << " " << neighbors[0] << endl;

                startNode = predecessor;

                lala2 += 1;
            }

        }

        
        //if(isCircular){
        //cout << startNode << " " << endNode << endl;
        //}

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

        //isVisited[startNode] = true;
        //isVisited[endNode] = true;

        //cout << unitigIndex << " " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << endl;
        //int lala = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
        //if(lala == 4452 && dummy){
        //    cout << "lal " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
        //}
        //lala = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy);
        //if(lala == 4452 && !dummy){
        //    cout << "loul " << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
        //}
        _nodeToUnitig[startNode] = _nextUnitigIndex;
        _nodeToUnitig[endNode] = _nextUnitigIndex;
        
        u_int32_t unitigIndexRC = _nextUnitigIndex + 1;
        _nodeToUnitig[nodeIndex_toReverseDirection(startNode)] = unitigIndexRC;
        _nodeToUnitig[nodeIndex_toReverseDirection(endNode)] = unitigIndexRC;

        u_int32_t abundance_sum = 0;
        u_int32_t abundance_max = 0;

        u_int32_t length = 0;
        u_int32_t node = startNode;
        //cout << "---------------" << endl;

        //cout << BiGraph::nodeIndex_to_nodeName(startNode) << " " << BiGraph::nodeIndex_to_nodeName(endNode) << endl;
        
        size_t i=0;
        u_int32_t nbNodes = 0;
        u_int32_t lastNode = -1;

        vector<u_int32_t> nodes;
        vector<u_int32_t> nodesRC;
        //if(unitigIndex == 114)
        //cout << "----" << endl;

        while(true){

            //if(unitigIndex == 114){
            //    cout << "HIIII: " << BiGraph::nodeIndex_to_nodeName(node) << " " << length << endl;
            //}

            //cout << _graphSuccessors->nodeIndex_to_nodeName(node, dummy) << endl;
            u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(node);

            if(_nodeAbundances[nodeName] > abundance_max){
                abundance_max = _nodeAbundances[nodeName];
            }

            if(i == 0){
                length += _nodeLengths[nodeName];
            }
            else{
                u_int16_t overlapLength = _graphSuccessors->getOverlap(lastNode, node);
                length += _nodeLengths[nodeName] - overlapLength;
            }

            //if(BiGraph::nodeIndex_to_nodeName(startNode) == 5304){
            //    cout << "\t" << length << endl;
            //}
            abundance_sum += _nodeAbundances[nodeName];
            //nbNodes += 1;

            getSuccessors(node, 0, neighbors);

            //if(unitigIndex == 114)
            //cout << neighbors.size() << endl;

            //cout << "lili: " << node << " " << neighbors.size() << endl;

            //isVisited[node] = true;
            _nodeToUnitig[node] = _nextUnitigIndex;
            nodes.push_back(node);
            nodesRC.push_back(nodeIndex_toReverseDirection(node));
            _nodeToUnitig[nodeIndex_toReverseDirection(node)] = unitigIndexRC;
            nbNodes += 1;

            //if(isCircular){
                //cout << "\tnode: " << BiGraph::nodeIndex_to_nodeName(node) << endl;
            //}

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

        /*
        if(startNode == 162558){
            cout << "-------------" << endl;
            cout << startNode << " " << endNode << endl;
            cout << lala1 << " " << lala2 << endl;
            getSuccessors(endNode, 0, neighbors);
            cout << "nb successors: " << neighbors.size() << endl;
            getPredecessors(startNode, 0, neighbors);
            cout << "nb predecessors: " << neighbors.size() << endl;
            //exit(1);
        }
        */

        //float abundance = ((float) abundance_sum) / nbNodes;
        //cout << abundance << endl;
        //if(_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, dummy) == 4453){
        //    cout << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << " " << dummy << " " << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << " " << dummy << endl;
        //}

        //if(nodeIndex % 2 == 0) continue;

        //cout << "Unitig: " << BiGraph::nodeIndex_to_nodeName(startNode) << " " << length << endl;
        //cout << BiGraph::nodeIndex_to_nodeName(startNode) << " " << BiGraph::nodeIndex_to_nodeName(endNode) << endl;


        std::reverse(nodesRC.begin(), nodesRC.end());
        _unitigs.push_back({_nextUnitigIndex, startNode, endNode, ((float)abundance_sum) / ((float)nbNodes), length, nbNodes, nodes});
        _unitigs.push_back({unitigIndexRC, nodeIndex_toReverseDirection(endNode), nodeIndex_toReverseDirection(startNode), ((float)abundance_sum) / ((float)nbNodes), length, nbNodes, nodesRC});




        //outfile << "S" << "\t" << unitigIndex << "\t" << "*" << endl;





        _nextUnitigIndex += 2;
        //cout << _graphSuccessors->nodeIndex_to_nodeName(startNode, dummy) << endl;
        //cout << _graphSuccessors->nodeIndex_to_nodeName(endNode, dummy) << endl;
        //cout << _graphSuccessors->nodeIndex_to_nodeName(neighbors[0], dummy) << endl;
        //getPredecessors(nodeName, neighbors);
        //cout << _graphSuccessors->nodeIndex_to_nodeName(neighbors[0], dummy) << endl;
    }

    void getSuccessors(u_int32_t n, float abundanceCutoff_min, vector<u_int32_t>& successors){

        bool dummy;
        successors.clear();

        vector<AdjNode>& nodes = _graphSuccessors->_nodes[n];
        //adjNode* node = _graphSuccessors->_nodes[n];

        //while (node != nullptr) {
        for(AdjNode& node : nodes){

			u_int32_t nn = node._index;
            //if(!_isNodeValid[nn]){
            if(_isNodeValid2.find(nn) == _isNodeValid2.end()){
                //node = node->next;
                continue;
            }

            //if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()){
            //    continue;
            //}

            if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()) continue;
            if(_allowedNodeIndex.size() > 0 && _allowedNodeIndex.find(BiGraph::nodeIndex_to_nodeName(nn)) == _allowedNodeIndex.end()) continue;

            if(isEdgeRemoved(n, nn)){
                //node = node->next;
                continue;
            }
            
			//if(abundanceCutoff_min != 0 && getNodeUnitigAbundance(nn) < abundanceCutoff_min) continue;

            
            /*
            if(BiGraph::nodeName_to_nodeIndex(2068, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(2068, true) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, true) == nn){
                continue;
            }
            */
            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
                
                
                //Abundance min cutoff
			//	node = node->next;
			//	continue;
			//}

            successors.push_back(nn);
            //node = node->next;
        }
    }

    void getSuccessors_overlap(u_int32_t n, float abundanceCutoff_min, vector<AdjNode>& successors){

        //cout << "succ overlap " << endl;

        bool dummy;
        successors.clear();

        vector<AdjNode>& nodes = _graphSuccessors->_nodes[n];
        //adjNode* node = _graphSuccessors->_nodes[n];

        //while (node != nullptr) {
        for(AdjNode& node : nodes){

			u_int32_t nn = node._index;
            //if(!_isNodeValid[nn]){
            if(_isNodeValid2.find(nn) == _isNodeValid2.end()){
                //node = node->next;
                continue;
            }

            //if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()){
            //    continue;
            //}


            if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()) continue;
            if(_allowedNodeIndex.size() > 0 && _allowedNodeIndex.find(BiGraph::nodeIndex_to_nodeName(nn)) == _allowedNodeIndex.end()) continue;

            if(isEdgeRemoved(n, nn)){
                //node = node->next;
                continue;
            }
            

			//if(abundanceCutoff_min != 0 && getNodeUnitigAbundance(nn) < abundanceCutoff_min) continue;
            
            /*
            if(BiGraph::nodeName_to_nodeIndex(2068, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(2068, true) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, true) == nn){
                continue;
            }
            */
            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
			//if(abundanceCutoff_min != 0 && _unitigs[_nodeToUnitig[nn]]._abundance < abundanceCutoff_min){ //Abundance min cutoff
			//	node = node->next;
			//	continue;
			//}

            successors.push_back(node);
            //node = node->next;
        }
    }

    void getPredecessors_overlap(u_int32_t n, float abundanceCutoff_min, vector<AdjNode>& successors){

        //cout << "pred overlap " << endl;
        bool dummy;
        successors.clear();

        vector<AdjNode>& nodes = _graphSuccessors->_nodes[nodeIndex_toReverseDirection(n)];
        //adjNode* node = _graphSuccessors->_nodes[n];

        //while (node != nullptr) {
        for(AdjNode& node : nodes){

			u_int32_t nn = nodeIndex_toReverseDirection(node._index);
            AdjNode node_nn = {nn, node._overlap};

            //cout << "YEY " << BiGraph::nodeIndex_to_nodeName(nn) << endl;
            //if(!_isNodeValid[nn]){
            if(_isNodeValid2.find(nn) == _isNodeValid2.end()){
                //node = node->next;
                continue;
            }

            if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()) continue;
            if(_allowedNodeIndex.size() > 0 && _allowedNodeIndex.find(BiGraph::nodeIndex_to_nodeName(nn)) == _allowedNodeIndex.end()) continue;

            //if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()){
            //    continue;
            //}


            if(isEdgeRemoved(n, nn)){
                //node = node->next;
                continue;
            }
            

			//if(abundanceCutoff_min != 0 && getNodeUnitigAbundance(nn) < abundanceCutoff_min) continue;
            
            /*
            if(BiGraph::nodeName_to_nodeIndex(2068, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(2068, true) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, true) == nn){
                continue;
            }
            */
            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
			//if(abundanceCutoff_min != 0 && _unitigs[_nodeToUnitig[nn]]._abundance < abundanceCutoff_min){ //Abundance min cutoff
			//	node = node->next;
			//	continue;
			//}

            successors.push_back(node_nn);
            //node = node->next;
        }
    }

    void getPredecessors(u_int32_t n, float abundanceCutoff_min, vector<u_int32_t>& predecessors){

        bool dummy;
        predecessors.clear();
        



        vector<AdjNode>& nodes = _graphSuccessors->_nodes[nodeIndex_toReverseDirection(n)];
        //adjNode* node = _graphSuccessors->_nodes[nodeIndex_toReverseDirection(n)];

        for(AdjNode& node : nodes){
        //while (node != nullptr) {

			u_int32_t nn = nodeIndex_toReverseDirection(node._index);
            //cout << "Pred: " << BiGraph::nodeIndex_to_nodeName(nn) << endl;

            //if(!_isNodeValid[nn]){
            if(_isNodeValid2.find(nn) == _isNodeValid2.end()){
                //node = node->next;
                continue;
            }

            if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()) continue;

            if(_allowedNodeIndex.size() > 0 && _allowedNodeIndex.find(BiGraph::nodeIndex_to_nodeName(nn)) == _allowedNodeIndex.end()) continue;
            //if(_isNodeInvalid_tmp.find(nn) != _isNodeInvalid_tmp.end()){
            //    continue;
            //}

            if(isEdgeRemoved(n, nn)){
                //node = node->next;
                continue;
            }


			//if(abundanceCutoff_min != 0 && getNodeUnitigAbundance(nn) < abundanceCutoff_min) continue;

            /*
            if(BiGraph::nodeName_to_nodeIndex(2068, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(2068, true) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, false) == nn){
                continue;
            }
            if(BiGraph::nodeName_to_nodeIndex(5645, true) == nn){
                continue;
            }*/

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nn, dummy);
			//if(abundanceCutoff_min != 0 && _unitigs[_nodeToUnitig[nn]]._abundance < abundanceCutoff_min){ //Abundance min cutoff
			//	node = node->next;
			//	continue;
			//}

            predecessors.push_back(nn);
            //node = node->next;
        }
    }

    
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

    u_int32_t unitigIndex_toReverseDirection(u_int32_t unitigIndex){
        return nodeIndex_to_unitigIndex(nodeIndex_toReverseDirection(_unitigs[unitigIndex]._startNode));

    }

    void getNeighbors(u_int32_t n, float abundanceCutoff_min, vector<u_int32_t>& neighbors){

        neighbors.clear();

        vector<u_int32_t> successors;
        getSuccessors(n, abundanceCutoff_min, successors);

        vector<u_int32_t> predecessors;
        getPredecessors(n, abundanceCutoff_min, predecessors);

        neighbors.insert(neighbors.end(), successors.begin(), successors.end());
        neighbors.insert(neighbors.end(), predecessors.begin(), predecessors.end());
    }

    
    void getSuccessors_unitig(u_int32_t unitigIndex, float abundanceCutoff, vector<u_int32_t>& successors){

        //cout << unitigIndex << " " << _unitigs.size() << endl;
        successors.clear();

        Unitig& unitig = _unitigs[unitigIndex];

        //cout << (unitig._endNode == -1)<< endl;
        //cout << BiGraph::nodeIndex_to_nodeName(unitig._endNode) << endl;
        vector<u_int32_t> successors_nodeIndex;
        getSuccessors(unitig._endNode, 0, successors_nodeIndex);
            

        for(u_int32_t nodeIndex : successors_nodeIndex){

            if(_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()) continue; //reinserting bubble

            //if(nodeIndex_to_unitigIndex(nodeIndex) == nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(1302, true))) cout << "hhhhhhha" << endl;
            //if(nodeIndex_to_unitigIndex(nodeIndex) == nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(1302, false))) cout << "hhhhhhha" << endl;
            //cout << (_isNodeInvalid_tmp.find(nodeIndex_to_unitigIndex(nodeIndex)) != _isNodeInvalid_tmp.end()) << endl;

            //if(_isNodeInvalid_tmp.find(nodeIndex_to_unitigIndex(nodeIndex)) != _isNodeInvalid_tmp.end()) continue;

            successors.push_back(nodeIndex_to_unitigIndex(nodeIndex));
        }
    }

    void getPredecessors_unitig(u_int32_t unitigIndex, float abundanceCutoff, vector<u_int32_t>& predecessors){

        predecessors.clear();

        Unitig& unitig = _unitigs[unitigIndex];

        vector<u_int32_t> predecessors_nodeIndex;
        getPredecessors(unitig._startNode, 0, predecessors_nodeIndex);
            
        for(u_int32_t nodeIndex : predecessors_nodeIndex){
            //if(_isNodeInvalid_tmp.find(nodeIndex_to_unitigIndex(nodeIndex)) != _isNodeInvalid_tmp.end()) continue;

            if(_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()) continue; //reinserting bubble

            predecessors.push_back(nodeIndex_to_unitigIndex(nodeIndex));
        }
    }

    /*
    u_int32_t in_degree(u_int32_t v){
        vector<u_int32_t> predecessors;
        getPredecessors(v, 0, predecessors);
        return predecessors.size();
    }
    
    u_int32_t out_degree(u_int32_t v){
        vector<u_int32_t> successors;
        getSuccessors(v, 0, successors);
        return successors.size();
    }
    */

    /*
    u_int32_t in_degree(u_int32_t unitigIndex){
        vector<u_int32_t> predecessors;
        getPredecessors_unitig(unitigIndex, predecessors);
        return predecessors.size();
    }
    
    u_int32_t out_degree(u_int32_t unitigIndex){
        vector<u_int32_t> successors;
        getSuccessors_unitig(unitigIndex, successors);
        return successors.size();
    }*/

    bool isEdgeRemoved(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){
        DbgEdge edge = {nodeIndex_from, nodeIndex_to};
		edge = edge.normalize();
        return _isEdgeRemoved.find(edge) != _isEdgeRemoved.end();
    }


    void getUnitigNodes_list(const vector<u_int32_t>& unitigs, unordered_set<u_int32_t>& nodeIndexes){
        
        nodeIndexes.clear();

        for(u_int32_t unitigIndex : unitigs){
            vector<u_int32_t> nodes;
            getUnitigNodes(_unitigs[unitigIndex], nodes);
            for(u_int32_t nodeIndex : nodes){
                //cout << "Scc: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
                nodeIndexes.insert(nodeIndex);
            }
        }

    }

    void getUnitigNodes(const Unitig& unitig, vector<u_int32_t>& nodes){


        nodes.clear();
        for(u_int32_t nodeIndex : unitig._nodes){
            nodes.push_back(nodeIndex);
        }

        return; // unitig._nodes;

        bool isCircular = false;

        nodes.clear();
        u_int32_t node = unitig._startNode;
        vector<u_int32_t> neighbors;


        while(true){
            //cout << node << endl;
            nodes.push_back(node);

            //cout << nodes.size() << " " << unitig._nbNodes << endl;
            if(nodes.size() == unitig._nbNodes) break;

            getSuccessors(node, 0, neighbors);
            //cout << neighbors.size() << endl;
            node = neighbors[0];
        }
        /*
            if(unitig._startNode == 162558){
                cout << "-------------" << endl;
                cout << unitig._startNode << " " << unitig._endNode << endl;
                getSuccessors(unitig._endNode, 0, neighbors);
                cout << "nb successors: " << neighbors.size() << endl;
                getPredecessors(unitig._startNode, 0, neighbors);
                cout << "nb predecessors: " << neighbors.size() << endl;
                exit(1);
            }

        getSuccessors(node, 0, neighbors);

        if(unitig._startNode == unitig._endNode && neighbors.size() == 1){
            isCircular = true;
        }

        size_t i = 0;
        while(true){

            //if(_isNodeRemoved[node]) continue;

            //u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(node, dummy);
            nodes.push_back(node);

            getSuccessors(node, 0, neighbors);
            //if(neighbors.size() != 1) break;

            if(isCircular){
                if(i == 0){ 
                    if(neighbors[0] == unitig._endNode) break; //Unitig with single node circular, or single node without edges
                }
                else{
                    if(node == unitig._endNode) break;
                }
            }
            else if(unitig._startNode == unitig._endNode){ //Single node circular and with several neighbors
                break;
            }
            else{
                if(node == unitig._endNode) break;
            }

            node = neighbors[0];

            i += 1;
        }

        */
    }

	double compute_median_float(vector<float> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0){
				return scores[size / 2 - 1]; //(scores[size / 2 - 1] + scores[size / 2]) / 2;
			}
			else {
				return scores[size / 2];
			}
		}
	}

    float _prevAbudanceCutoff;

    struct LongNeighbor{
        u_int32_t _nodeIndex;
        float _abundance;
    };

    void debug_writeGfaErrorfree(u_int32_t currentAbundance, float abundanceCutoff_min, u_int32_t nodeIndex_source, u_int64_t k){


        unordered_map<u_int32_t, vector<LongNeighbor>> nodeLongNeighbors;
        bool useConnectComponent = false;
        /*
        if(abundanceCutoff_min < 6){
            useConnectComponent = false;
            if(((u_int64_t)_prevAbudanceCutoff) == ((u_int64_t)abundanceCutoff_min)) return;
            //useConnectComponent = true;
            //abundanceCutoff_min = 1;
        }*/
        
        //if(currentAbundance == _lastAbundanceCutoff) return;
        //_lastAbundanceCutoff = currentAbundance;
        
        //bool needCleaning = true;


        //loadState(abundanceCutoff_min);
        clear(abundanceCutoff_min);
        compact(false);
        //vector<float> lala = {40, 800};
        //cout << compute_median_float(lala) << endl;
        //exit(1);

  
        //clear(0);

        _prevAbudanceCutoff = abundanceCutoff_min;

        while(true){

            //compact(false);

            u_int64_t nbTipsRemoved = 0;
            u_int64_t nbErrorRemoved = 0;
            u_int64_t nbBubblesRemoved = 0;
            u_int64_t nbSmallLoopRemoved = 0;

            bool isModification = false;

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;
            
            if(useConnectComponent){
                //compact(true);
                removeUnconnectedNodes(nodeIndex_source);
            }
            /*
            //u_int32_t loulou = BiGraph::nodeName_to_nodeIndex(214, false);

            //remove very low abundant nodes
            while(true){
                compact();
                nbErrorRemoved = removeLowAbundantNodes(abundanceCutoff_min/2.0);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb error removed 1: " << nbErrorRemoved << endl;
                #endif
                if(nbErrorRemoved == 0) break;
                isModification = true;
            }

            //Remove abundant nodes progressively (should save variant bubbles with lower abundance)
            for(float i=0.1; i<abundanceCutoff_min; i++){
                while(true){
                    compact();
                    nbErrorRemoved = removeLowAbundantNodes(i);
                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb error removed 2: " << nbErrorRemoved << endl;
                    #endif
                    if(nbErrorRemoved == 0) break;
                    isModification = true;
                }
            }

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;
            */
            if(useConnectComponent){
                removeUnconnectedNodes(nodeIndex_source);
            }


            while(true){

                bool isModSub = false;

                while(true){
                    compact(true);

                    unordered_set<u_int32_t> isTips;
                    nbTipsRemoved = tip(4*k, true, isTips);
                    nbTipsRemoved = tip(4*k, false, isTips);

                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb tip removed: " << nbTipsRemoved << endl;
                    #endif
                    if(nbTipsRemoved == 0) break;
                    isModification = true;
                    isModSub = true;
                }

                while(true){
                    compact(true);
                    u_int64_t nbSuperbubblesRemoved = superbubble(100000);
                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb superbubble removed: " << nbSuperbubblesRemoved << endl;
                    #endif
                    if(nbSuperbubblesRemoved == 0) break;
                    isModification = true;
                    isModSub = true;
                }
                

                //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

                while(true){
                    cout << "bubbulu 1" << endl;
                    compact(true);
                    cout << "bubbulu 2" << endl;
                    nbBubblesRemoved = bubble(100000);
                    cout << "bubbulu 3" << endl;
                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb bubble removed: " << nbBubblesRemoved << endl;
                    #endif
                    if(nbBubblesRemoved == 0) break;
                    isModification = true;
                    isModSub = true;
                }

                while(true){
                    compact(true);
                    nbSmallLoopRemoved = removeSmallLoop(10000);
                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb small loop removed: " << nbSmallLoopRemoved << endl;
                    #endif
                    if(nbSmallLoopRemoved == 0) break;
                    isModification = true;
                }

                if(!isModSub) break;
            }
            
            //unordered_set<u_int32_t> 
            u_int64_t nbErrorsRemoved = 0;
            float localCutoffRate = 0.5;
            float cutoffGlobal = 1;
            float aplha = 0.1;
            double t=2.5;
            u_int32_t prevCutoff = -1;
            while(t < 500){ //500 = max unitig abundance over graph
                
                unordered_set<u_int32_t> removedNodes;
                //if(((u_int32_t)t) != prevCutoff){
                //    prevCutoff = t;

                cout << t << endl;
                //cout << "Nb node valid: " << _isNodeValid2.size() << endl;

                compact(true);
                    
                //t = t * (1+aplha);
                //continue;
                    
                /*
                for(Bubble& bubble : _bubbles){
                    u_int32_t startNode = bubble._nodes[0];
                    u_int32_t endNode = bubble._nodes[bubble._nodes.size()-1];
                    if(_isNodeValid2.find(startNode) != _isNodeValid2.end()) continue;
                    //if(_isNodeValid[startNode]) continue; //&& !_isNodeRemoved[endNode]

                    getPredecessors(startNode, 0, neighbors);
                    if(neighbors.size() == 0) continue;

                    bool isPathValid = false;
                    for(u_int32_t n : neighbors){
                        //if(_isNodeValid[n]){
                        if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                            isPathValid = true;
                            break;
                        }
                    }

                    if(!isPathValid) continue;

                    getSuccessors(endNode, 0, neighbors);
                    if(neighbors.size() == 0) continue;

                    isPathValid = false;
                    for(u_int32_t n : neighbors){
                        if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                        //if(_isNodeValid[n]){
                            isPathValid = true;
                            break;
                        }
                    }

                    if(!isPathValid) continue;

                    for(u_int32_t nodeIndex : bubble._nodes){
                        _isNodeValid2.insert(nodeIndex);
                        //_isNodeValid2.insert(nodeIndex_toReverseDirection(nodeIndex));
                        //cout << "Bubble saved: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
                        //_isNodeValid[nodeIndex] = true;
                    }

                    bubbleAdded += 1;
                }
                */

                for(Unitig& unitig : _unitigs){
                    if(unitig._startNode == -1) continue;
                    //cout << unitig._nbNodes << " " << unitig._abundance << endl;
                    if(unitig._nbNodes >= 4*k) continue;

                    //cout << "mioum: " << _unitigs[0]._startNode << endl;
                    //cout << unitig._index << " " << _unitigs.size() << endl;
                    /*
                    vector<u_int32_t> successors;
                    getSuccessors_unitig(unitig._index, successors);

                    vector<u_int32_t> predecessors;
                    getPredecessors_unitig(unitig._index, predecessors);

                    vector<float> values;// = {40, 800};

                    double neighborAbundanceMin = 1000000;
                    double neighborAbundanceMax = 0;
                    double neighborAbundance = 0;
                    for(u_int32_t unitigIndex : successors){
                        if(_unitigs[unitigIndex]._abundance < neighborAbundanceMin){
                            neighborAbundanceMin = _unitigs[unitigIndex]._abundance;
                        }
                        if(_unitigs[unitigIndex]._abundance > neighborAbundanceMax){
                            neighborAbundanceMax = _unitigs[unitigIndex]._abundance;
                        }
                        neighborAbundance += _unitigs[unitigIndex]._abundance;
                        values.push_back(_unitigs[unitigIndex]._abundance);
                    }
                    for(u_int32_t unitigIndex : predecessors){
                        if(_unitigs[unitigIndex]._abundance < neighborAbundanceMin){
                            neighborAbundanceMin = _unitigs[unitigIndex]._abundance;
                        }
                        if(_unitigs[unitigIndex]._abundance > neighborAbundanceMax){
                            neighborAbundanceMax = _unitigs[unitigIndex]._abundance;
                        }
                        neighborAbundance += _unitigs[unitigIndex]._abundance;
                        values.push_back(_unitigs[unitigIndex]._abundance);
                    }
                    neighborAbundance /= (successors.size() + predecessors.size());
                    
                    //double nbRepeatedTimes = neighborAbundanceMax / currentAbundance;
                    //double minCutoff = (currentAbundance / 4);
                    //double maxCutoff = (currentAbundance / 2);
                    //double cutoff = neighborAbundanceMax / 18;
                    //float median = compute_median_float(values);
                    //cout << unitig._abundance << " " << neighborAbundance << endl;
                    //float cutoff = localCutoffRate * neighborAbundanceMin;
                    //if(cutoff < t) cutoff = t;
                    //cutoff = 12;
                    //float cutoff = min(cutoffGlobal, cutoffLocal);
                    */

                    //cout << "1" << endl;
                    double localabundance = computeLocalAbundance(unitig._index, k, nodeLongNeighbors);
                    double cutoff = min(t, localabundance*0.5);

                    //cout << "2" << endl;
                    if(unitig._abundance < cutoff){
                        vector<u_int32_t> unitigNodes;
                        getUnitigNodes(unitig, unitigNodes);
                        for(u_int32_t node : unitigNodes){
                            removedNodes.insert(node);
                            removedNodes.insert(nodeIndex_toReverseDirection(node));
                            //_isNodeValid2.erase(node);
                            nbErrorsRemoved += 1;
                        }
                        isModification = true;
                    }
                    //cout << "3" << endl;
                    //cout << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._abundance << " " << computeLocalAbundance(unitig._index, k) << endl;
                    //if(localabundance < 30) getchar();
                    //cout << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._abundance << " " << computeLocalAbundance(unitig._index, k) << endl;
                    //if(computeLocalAbundance(unitig._index, k) > 50){
                    //    getchar();
                    //}
                    /*
                    if(unitig._abundance < localCutoffRate * neighborAbundanceMin){
                        vector<u_int32_t> unitigNodes;
                        getUnitigNodes(unitig, unitigNodes);
                        for(u_int32_t node : unitigNodes){
                            _isNodeValid2.erase(node);
                            nbErrorsRemoved += 1;
                        }
                        isModification = true;
                        
                    }*/

                    //if(20545 == BiGraph::nodeIndex_to_nodeName(unitig._startNode) || 20545 == BiGraph::nodeIndex_to_nodeName(unitig._endNode)){
                    //    cout << unitig._abundance << " " << localabundance << " " << cutoff << endl;
                    //    getchar();
                    //}

                }
                //}
               
                unordered_set<u_int32_t> removedUnitigs;
                for(u_int32_t nodeIndex : removedNodes){
                    //if(_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()){
                    //    cout << "KOUERK" << endl;
                    //}
                    //_unitigs[_nodeToUnitig[nodeIndex]]._index

                    removedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex));
                }
                removeUnitigs(removedUnitigs);
                for(u_int32_t nodeIndex : removedNodes){
                    _isNodeValid2.erase(nodeIndex);
                    //_removedFrom[nodeIndex] = 4;
                }
                

                t = t * (1+aplha);
            }

            //cout << "Error: " << _isNodeRemoved[loulou] << endl;
            //compact();
            



            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            //exit(1);






            /*
            unordered_set<u_int32_t> validNodes;
            for (auto& nodeIndex : _isNodeValid2){
                u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
                validNodes.insert(nodeName);
            }
            GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_errorFree.gfa", validNodes, _isEdgeRemoved, _graphSuccessors);
            */

            





            if(!isModification) break;
        }

        cout << "inserting bubbles" << endl;

        vector<u_int32_t> isLongUnitigNodes;
        for(u_int32_t nodeIndex : _isNodeValid2){
            //if(_nodeToUnitig.find(nodeIndex) == _nodeToUnitig.end()) continue;

            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
            if(_unitigs[unitigIndex]._length > 10000 && !_isBubble[nodeIndex]){
                isLongUnitigNodes.push_back(nodeIndex);
            }
        }

        //cout << "Nb bubbles: " << _bubbles.size() << endl;
        //for(Bubble& bubble : _bubbles){
        //    for(u_int32_t nodeIndex : bubble._nodes){
        //        cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
        //    }
        //}

        /*
        cout << "inserting bubbles" << endl;
        while(true){
            vector<u_int32_t> neighbors;
            u_int64_t bubbleAdded = 0;

            u_int32_t startUnitigIndex = -1;
            
            for(Bubble& bubble : _bubbles){

                u_int32_t startNode = bubble._nodes[0]._nodeIndex;
                u_int32_t endNode = bubble._nodes[bubble._nodes.size()-1]._nodeIndex;
                if(_isNodeValid2.find(startNode) != _isNodeValid2.end()) continue;
                //if(_isNodeValid2.find(nodeIndex_toReverseDirection(startNode)) != _isNodeValid2.end()) continue;
                //if(_isNodeValid[startNode]) continue; //&& !_isNodeRemoved[endNode]

                getPredecessors(startNode, 0, neighbors);
                if(neighbors.size() == 0) continue;


                bool isPathValid = false;
                for(u_int32_t n : neighbors){
                    //if(_isNodeValid[n]){
                    if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                        isPathValid = true;
                        if(_nodeToUnitig.find(n) != _nodeToUnitig.end()) startUnitigIndex = nodeIndex_to_unitigIndex(n);
                        break;
                    }
                }

                if(!isPathValid) continue;

                getSuccessors(endNode, 0, neighbors);
                if(neighbors.size() == 0) continue;


                isPathValid = false;
                for(u_int32_t n : neighbors){
                    if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                    //if(_isNodeValid[n]){
                        isPathValid = true;
                        if(_nodeToUnitig.find(n) != _nodeToUnitig.end()) startUnitigIndex = nodeIndex_to_unitigIndex(n);
                        break;
                    }
                }

                if(!isPathValid) continue;

                //cout << bubble._nodes[0]._abundance << " " << computeLocalAbundance(startUnitigIndex, k, nodeLongNeighbors) << endl;

                for(NodeAbundance& node : bubble._nodes){

                    //u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                    double localabundance = computeLocalAbundance(startUnitigIndex, k, nodeLongNeighbors);
                    //cout << localabundance << endl;
                    double cutoff = localabundance*0.5;

                    
                    //cout << "2" << endl;
                    if(node._abundance >= cutoff){
                        if(_isNodeValid2.find(node._nodeIndex) == _isNodeValid2.end()){
                            bubbleAdded += 1;
                            _isNodeValid2.insert(node._nodeIndex);
                            //if(node._abundance < 5){
                                //cout << node._abundance << " " << localabundance << endl;
                            //    getchar();
                            //}
                        }

                        //cout << node._abundance << endl;
                        //cout << "allo: " <<  startUnitigIndex << endl;
                    }

                    //_isNodeValid2.insert(nodeIndex_toReverseDirection(nodeIndex));
                    //cout << "Bubble saved: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
                    //_isNodeValid[nodeIndex] = true;
                }

            }

            //cout << bubbleAdded << endl;
            if(bubbleAdded == 0) break;
        }
        cout << "done" << endl;
        */

        compact(false);


        //cout << (_isNodeValid2.find(BiGraph::nodeName_to_nodeIndex(1840, false)) != _isNodeValid2.end())<< endl;
        cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;
        //getchar();


        
        unordered_set<u_int32_t> validNodes;
        for (auto& nodeIndex : _isNodeValid2){
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            validNodes.insert(nodeName);
        }
        string outputFilename = _outputDir + "/minimizer_graph_cleaned.gfa";
        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, outputFilename, validNodes, _isEdgeRemoved, _graphSuccessors);
        

        //compact();
        //superbubble(50000);
        
        for(u_int32_t nodeIndex : isLongUnitigNodes){
            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
            _isLongUnitig.insert(unitigIndex);
        }

        /*
        //if(currentAbundance == _lastAbundanceCutoff) return;
        //_lastAbundanceCutoff = currentAbundance;
        
        //bool needCleaning = true;
        bool useConnectComponent = true;
        if(abundanceCutoff_min < 6){
            useConnectComponent = false;
            if(((u_int64_t)_prevAbudanceCutoff) == ((u_int64_t)abundanceCutoff_min)) return;
            //useConnectComponent = true;
            //abundanceCutoff_min = 1;
        }

        clear(abundanceCutoff_min);
		//loadState();



        _prevAbudanceCutoff = abundanceCutoff_min;

        while(true){

            u_int64_t nbTipsRemoved = 0;
            u_int64_t nbErrorRemoved = 0;
            u_int64_t nbBubblesRemoved = 0;
            u_int64_t nbSmallLoopRemoved = 0;

            bool isModification = false;

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            if(useConnectComponent){
                compact();
                removeUnconnectedNodes(nodeIndex_source);
            }

            //u_int32_t loulou = BiGraph::nodeName_to_nodeIndex(214, false);

            //remove very low abundant nodes
            while(true){
                compact();
                nbErrorRemoved = removeLowAbundantNodes(abundanceCutoff_min/2.0);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb error removed 1: " << nbErrorRemoved << endl;
                #endif
                if(nbErrorRemoved == 0) break;
                isModification = true;
            }

            //Remove abundant nodes progressively (should save variant bubbles with lower abundance)
            for(float i=0.1; i<abundanceCutoff_min; i++){
                while(true){
                    compact();
                    nbErrorRemoved = removeLowAbundantNodes(i);
                    #ifdef PRINT_DEBUG_SIMPLIFICATION
                        cout << "Nb error removed 2: " << nbErrorRemoved << endl;
                    #endif
                    if(nbErrorRemoved == 0) break;
                    isModification = true;
                }
            }

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            if(useConnectComponent){
                removeUnconnectedNodes(nodeIndex_source);
            }
getStronglyConnectedComponent_node
            //cout << "Error: " << _isNodeRemoved[loulou] << endl;
            //compact();

            while(true){
                compact();
                nbTipsRemoved = tip(50000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb tip removed: " << nbTipsRemoved << endl;
                #endif
                if(nbTipsRemoved == 0) break;
                isModification = true;
            }

            //cout << "tip: " << _isNodeRemoved[loulou] << endl;

            

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;
            
            
            while(true){
                compact();
                u_int64_t nbSuperbubblesRemoved = superbubble(50000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb superbubble removed: " << nbSuperbubblesRemoved << endl;
                #endif
                if(nbSuperbubblesRemoved == 0) break;
                isModification = true;
            }
            nodeIndex_neighbor
            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            while(true){
                compact();
                nbBubblesRemoved = bubble(50000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb bubble removed: " << nbBubblesRemoved << endl;
                #endif
                if(nbBubblesRemoved == 0) break;
                isModification = true;
            }

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            //exit(1);
            while(true){
                compact();
                nbSmallLoopRemoved = removeSmallLoop(10000);
                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    cout << "Nb small loop removed: " << nbSmallLoopRemoved << endl;
                #endif
                if(nbSmallLoopRemoved == 0) break;
                isModification = true;
            }

            //cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;

            //cout << "Small bubble: " << _isNodeRemoved[loulou] << endl;

            if(!isModification) break;
        }


        //cout << "Nb bubbles: " << _bubbles.size() << endl;
        //for(Bubble& bubble : _bubbles){
        //    for(u_int32_t nodeIndex : bubble._nodes){
        //        cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
        //    }
        //}

        while(true){
            vector<u_int32_t> neighbors;
            u_int64_t bubbleAdded = 0;

            for(Bubble& bubble : _bubbles){
                u_int32_t startNode = bubble._nodes[0];
                u_int32_t endNode = bubble._nodes[bubble._nodes.size()-1];
                if(_isNodeValid2.find(startNode) != _isNodeValid2.end()) continue;
                //if(_isNodeValid[startNode]) continue; //&& !_isNodeRemoved[endNode]

                getPredecessors(startNode, 0, neighbors);
                if(neighbors.size() == 0) continue;

                bool isPathValid = false;
                for(u_int32_t n : neighbors){
                    //if(_isNodeValid[n]){
                    if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                        isPathValid = true;
                        break;
                    }
                }

                if(!isPathValid) continue;

                getSuccessors(endNode, 0, neighbors);
                if(neighbors.size() == 0) continue;

                isPathValid = false;
                for(u_int32_t n : neighbors){
                    if(_isNodeValid2.find(n) != _isNodeValid2.end()){
                    //if(_isNodeValid[n]){
                        isPathValid = true;
                        break;
                    }
                }

                if(!isPathValid) continue;

                for(u_int32_t nodeIndex : bubble._nodes){
                    _isNodeValid2.insert(nodeIndex);
                    _isNodeValid2.insert(nodeIndex_toReverseDirection(nodeIndex));
                    //cout << "Bubble saved: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
                    //_isNodeValid[nodeIndex] = true;
                }

                bubbleAdded += 1;
            }

            if(bubbleAdded == 0) break;
        }

        compact();

        //cout << (_isNodeValid2.find(BiGraph::nodeName_to_nodeIndex(1840, false)) != _isNodeValid2.end())<< endl;
        cout << "Nb nodes valid: " << _isNodeValid2.size() << endl;
        //compact();
        //superbubble(50000);
        */
        
        /*
        //Debug gfa
        unordered_set<u_int32_t> validNodes;
        for (auto& nodeIndex : _isNodeValid2){
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            validNodes.insert(nodeName);
        }
        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_errorFree.gfa", validNodes, _isEdgeRemoved, _graphSuccessors);
        */

        //exit(1);
        /*
        //cout << "orlolo" << endl;
        for(size_t nodeIndex=0; nodeIndex<_graphSuccessors->_nbNodes; nodeIndex++){
            //cout << nodeIndex << " " << BiGraph::nodeName_to_nodeIndex(214, true) << endl;
            //if(nodeIndex == BiGraph::nodeName_to_nodeIndex(214, true)){
            //    cout << "haha" << endl;
            //    cout << _isNodeRemoved[nodeIndex] << endl;
            //}
            //if(nodeIndex == BiGraph::nodeName_to_nodeIndex(214, false)){
            //    cout << "haha" << endl;
            //    cout << _isNodeRemoved[nodeIndex] << endl;
            //}

            //if(_isNodeValid[nodeIndex]) continue;
            if(_isNodeValid2.find(nodeIndex) != _isNodeValid2.end()){
                cout << _graphSuccessors->nodeIndex_to_nodeName(nodeIndex) << endl;
                cout << "lala" << _isNodeValid2.size() << endl;
                continue;
            }
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            removedNodes.insert(nodeName);
        }
        */



        //cout << "lalalal" << endl;
        //cout << _isNodeRemoved[BiGraph::nodeName_to_nodeIndex(214, true)] << endl;

        //exit(1);
    }

    void saveCurrentGfa(u_int32_t nodeIndex_source, u_int32_t abundance, u_int32_t binIndex){



        cout << "Saving current gfa" << endl;
        //Bizarre mais peut etre detruit dans les region pourri composé de 3 unitig relié a leur extreminité (donc 3 tips) par exemple
        if(_isNodeValid2.find(nodeIndex_source) == _isNodeValid2.end()){
            return;
        }


        unordered_set<u_int32_t> memoNodes = _isNodeValid2;
        for(u_int32_t nodeIndex : memoNodes){
            if(getNodeUnitigAbundance(nodeIndex) < abundance/4){
                _isNodeValid2.erase(nodeIndex);
            }
        }

        unordered_set<u_int32_t> component;
        getConnectedComponent(nodeIndex_source, component);

        _isNodeValid2 = memoNodes;
        /*
        cout << "Compute strongly connected component" << endl;
        vector<u_int32_t> component;
        getStronglyConnectedComponent(nodeIndex_source, component); //ATTENTION UTILISE MAINTENANT UNITIG_INDEX EN ENTREE et pas nodeIndex_source
        cout << "done" << endl;
        //exit(1);
        */

        unordered_set<u_int32_t> isNodeValid2;

        //cout << component.size() << endl;
        //u_int32_t componentSize = 0;

        vector<u_int32_t> nodes;
        for(u_int32_t unitigIndex : component){
            getUnitigNodes(_unitigs[unitigIndex], nodes);
            for(u_int32_t nodeIndex : nodes){
                isNodeValid2.insert(nodeIndex);
                //_isNodeValid[nodeIndex] = true;
                //componentSize += 1;
            }
        }


        fs::path path(_outputDir + "/gfa");
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 
        //System::file().mkdir(_outputDir + "/gfa", -1);

        //fs::path p(_outputDir + "/gfa");
        //if(!fs::exists(p)){
        //    fs::create_directory(_outputDir + "/gfa");
        //}

        cout << "Component size: " << isNodeValid2.size() << endl;
        
        unordered_set<u_int32_t> validNodes;
        for (auto& nodeIndex : isNodeValid2){
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            validNodes.insert(nodeName);
        }
        string outputFilename = _outputDir + "/gfa/" + to_string(binIndex) + "_" + to_string(abundance) + "_" + to_string(BiGraph::nodeIndex_to_nodeName(nodeIndex_source)) + ".gfa";
        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, outputFilename, validNodes, _isEdgeRemoved, _graphSuccessors);
        //getchar();

    }
    
    void getConnectedComponent(u_int32_t nodeIndex_source, unordered_set<u_int32_t>& component){

        component.clear();

        vector<u_int32_t> neighbors;
        //unordered_set<u_int32_t> isVisited;

        //for(size_t n=0; n<_nbNodes; n++){
        //    if(isVisited[n]) continue;

        queue <u_int32_t> queue;

        //cout << "----------" << endl;
        //cout << nodeIndex_source << endl;
        //cout << _unitigs.size() << endl;
        //cout << _nodeToUnitig[nodeIndex_source] << endl;


        queue.push(nodeIndex_to_unitigIndex(nodeIndex_source));
        component.insert(nodeIndex_to_unitigIndex(nodeIndex_source));
        //component.push_back(nodeIndex_source);

        while (!queue.empty()){

            u_int64_t unitigIndex = queue.front();
            //cout << "1: " << unitigIndex << endl;

            queue.pop();

            vector<u_int32_t> neighbors1;
            getNeighbors(_unitigs[unitigIndex]._startNode, 0, neighbors1);
            vector<u_int32_t> neighbors2;
            getNeighbors(_unitigs[unitigIndex]._endNode, 0, neighbors2);


            for(u_int32_t nodeIndex_neighbor : neighbors1){
                //cout << "2: " << _nodeToUnitig[nodeIndex_neighbor] << endl;
                u_int32_t unitigIndex_neighbor = nodeIndex_to_unitigIndex(nodeIndex_neighbor);
                if (component.find(unitigIndex_neighbor) != component.end()) continue;

                queue.push(unitigIndex_neighbor);

                component.insert(unitigIndex_neighbor);
                //component.push_back(nodeIndex_neighbor);
            }

            for(u_int32_t nodeIndex_neighbor : neighbors2){
                //cout << "3: " << _nodeToUnitig[nodeIndex_neighbor] << endl;
                u_int32_t unitigIndex_neighbor = nodeIndex_to_unitigIndex(nodeIndex_neighbor);
                if (component.find(unitigIndex_neighbor) != component.end()) continue;

                queue.push(unitigIndex_neighbor);

                component.insert(unitigIndex_neighbor);
                //component.push_back(nodeIndex_neighbor);
            }
        }

    }

    void removeUnconnectedNodes(u_int32_t nodeIndex_source){
        
        //Bizarre mais peut etre detruit dans les region pourri composé de 3 unitig relié a leur extreminité (donc 3 tips) par exemple
        if(_isNodeValid2.find(nodeIndex_source) == _isNodeValid2.end()){
            return;
        }

        cout << "Compute connected component" << endl;
        cout << _isNodeValid2.size() << endl;
        unordered_set<u_int32_t> component;
        getConnectedComponent(nodeIndex_source, component);
        cout << "done" << endl;

        /*
        cout << "Compute strongly connected component" << endl;
        vector<u_int32_t> component;
        getStronglyConnectedComponent(nodeIndex_source, component); //ATTENTION UTILISE MAINTENANT UNITIG_INDEX EN ENTREE et pas nodeIndex_source
        cout << "done" << endl;
        //exit(1);
        */

        unordered_set<u_int32_t> isNodeValid2;

        //cout << component.size() << endl;
        //u_int32_t componentSize = 0;

        vector<u_int32_t> nodes;
        for(u_int32_t unitigIndex : component){
            getUnitigNodes(_unitigs[unitigIndex], nodes);
            for(u_int32_t nodeIndex : nodes){
                isNodeValid2.insert(nodeIndex);
                //_isNodeValid[nodeIndex] = true;
                //componentSize += 1;
            }
        }

        cout << "Component size: " << isNodeValid2.size() << endl;
        _isNodeValid2 = isNodeValid2;

        /*
        unordered_set<u_int32_t> validNodes;
        for (auto& nodeIndex : _isNodeValid2){
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            validNodes.insert(nodeName);
        }
        //Debug gfa
        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_errorFree.gfa", validNodes, _isEdgeRemoved, _graphSuccessors);
        exit(1);
        */
        /*
        if(_isNodeValid2.size() == 90){


            unordered_set<u_int32_t> validNodes;
            
            for (auto& nodeIndex : _isNodeValid2){
                u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
                cout << BiGraph::nodeToString(nodeName) << endl;
                validNodes.insert(nodeName);
            }

            GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_component.gfa", validNodes, _isEdgeRemoved, _graphSuccessors);
        }
        */

    }

    double computeLocalAbundance(u_int32_t unitigIndex_source, u_int64_t k, unordered_map<u_int32_t, vector<LongNeighbor>>& nodeLongNeighbors){

        //cout << "---------------" << endl;

        float abundanceMinIndexed = std::numeric_limits<float>::max();

        bool needIndexing = true;
        vector<u_int32_t> unitigNodes;
        getUnitigNodes(_unitigs[unitigIndex_source], unitigNodes);


        bool hasChanged = false;
        for(u_int32_t nodeIndex : unitigNodes){
            u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

            const auto& it = nodeLongNeighbors.find(nodeName);
            if(it != nodeLongNeighbors.end()){
                needIndexing = false;
                for(const LongNeighbor& ln : it->second){

                    if(_isNodeValid2.find(ln._nodeIndex) == _isNodeValid2.end()){
                        hasChanged = true;
                        break;
                    }

                    if(_nodeToUnitig.find(ln._nodeIndex) == _nodeToUnitig.end()){
                        hasChanged = true;
                        break;
                    }

                    u_int32_t unitigIndex = nodeIndex_to_unitigIndex(ln._nodeIndex);
                    float abundance = ln._abundance;
                    if(abundance < abundanceMinIndexed){
                        abundanceMinIndexed = abundance;
                    }
                    //if(nodeName == 580)
                    //cout << "Existing: " << BiGraph::nodeIndex_to_nodeName(ln._nodeIndex) << " " << abundance << endl;

                    float currentAbundance = _unitigs[unitigIndex]._abundance;
                    if(currentAbundance != abundance){
                        hasChanged = true;
                        break;
                    }
                }
                //return 50;
            }
            if(hasChanged) break;
        }


        //cout << hasChanged << endl;
        if(hasChanged){
            needIndexing = true;
            for(u_int32_t nodeIndex : unitigNodes){
                u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
                if(nodeLongNeighbors.find(nodeName) != nodeLongNeighbors.end()){
                    //cout << "cleaned" << endl;
                    nodeLongNeighbors.erase(nodeName);
                }
            }
        }

        if(!hasChanged && !needIndexing){
            return abundanceMinIndexed;
        }
        //cout << "Indexed: " << abundanceMinIndexed << endl;

        static float maxVal = std::numeric_limits<float>::max();

        //cout << "recomp" << endl;
        int lala = 0;

        u_int64_t totalLongUnitigLength = 0;
        u_int64_t nbLongUnitigs = 0;
        float abundanceMin = maxVal;

        //unordered_map<u_int32_t, u_int32_t> distance;
        //neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);

        unordered_set<u_int32_t> isVisited;
    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        //distance[unitigIndex_source] = 0;
        queue.push(unitigIndex_source);
        isVisited.insert(unitigIndex_source);
        //cout << "------" << endl;

        while(!queue.empty()){
            
            if(nbLongUnitigs >= 2 && totalLongUnitigLength > 1500) break;

            u_int32_t unitigIndex = queue.front();
            queue.pop();

            vector<u_int32_t> successors;
            getSuccessors_unitig(unitigIndex, 0, successors);

            vector<u_int32_t> predecessors;
            getPredecessors_unitig(unitigIndex, 0, predecessors);

            lala += 1;
            //if(lala > 30000) break;
            //cout << unitigIndex << " " << nbLongUnitigs << " " << totalLongUnitigLength << endl;

            for(u_int32_t unitigIndex_nn : successors){

                if(isVisited.find(unitigIndex_nn) != isVisited.end()) continue;                    

                if(_unitigs[unitigIndex_nn]._nbNodes >= 3*k){
                //if(_unitigs[unitigIndex_nn]._nod > 1500){
                    nbLongUnitigs += 1;
                    totalLongUnitigLength += _nodeLengths[BiGraph::nodeIndex_to_nodeName(unitigIndex_nn)];
                    if(_unitigs[unitigIndex_nn]._abundance < abundanceMin){
                        abundanceMin = _unitigs[unitigIndex_nn]._abundance;
                    }

                    if(needIndexing){
                        nodeLongNeighbors[BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_source]._startNode)].push_back({_unitigs[unitigIndex_nn]._startNode, _unitigs[unitigIndex_nn]._abundance});
                    }    
                    //if(BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_source]._startNode) == 580)
                    //cout << "Truth: " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_nn]._startNode) << " " << _unitigs[unitigIndex_nn]._abundance << endl;
                }

                isVisited.insert(unitigIndex_nn);
                queue.push(unitigIndex_nn);
            }

            for(u_int32_t unitigIndex_nn : predecessors){

                if(isVisited.find(unitigIndex_nn) != isVisited.end()) continue;

                if(_unitigs[unitigIndex_nn]._nbNodes >= 4*k){
                //if(_unitigs[unitigIndex_nn]._length > 10000){
                    nbLongUnitigs += 1;
                    totalLongUnitigLength += _nodeLengths[BiGraph::nodeIndex_to_nodeName(unitigIndex_nn)];
                    if(_unitigs[unitigIndex_nn]._abundance < abundanceMin){
                        abundanceMin = _unitigs[unitigIndex_nn]._abundance;
                    }
                    
                    if(needIndexing){
                        nodeLongNeighbors[BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_source]._startNode)].push_back({_unitigs[unitigIndex_nn]._startNode, _unitigs[unitigIndex_nn]._abundance});
                    }  
                    //if(BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_source]._startNode) == 580)
                    //cout << "Truth: " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_nn]._startNode) << " " << _unitigs[unitigIndex_nn]._abundance << endl;
                }

                isVisited.insert(unitigIndex_nn);
                queue.push(unitigIndex_nn);

            }


        }

        //if(!hasChanged && !needIndexing){
            //getchar();
        //}
        //if(abundanceMin == maxVal) return 1;
        //cout << "Truth: " << abundanceMin << " " << hasChanged << " " << needIndexing << endl;
        //cout << "Nb nodes: " << lala << endl;
        //cout << abundanceMin << endl;
        return abundanceMin;

    }

    void collectNeighbors(u_int32_t nodeIndex_source, u_int32_t maxDistance, u_int32_t maxNeighbors, unordered_set<u_int32_t>& neighbors){

        unordered_map<u_int32_t, u_int32_t> distance;
        neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);


    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        distance[nodeIndex_source] = 0;
        queue.push(nodeIndex_source);
        neighbors.insert(nodeIndex_source);
    
        while(!queue.empty()){
            
            u_int32_t n = queue.front();
            queue.pop();

            vector<u_int32_t> nodeNeighbors;
            getNeighbors(n, 0, nodeNeighbors);

            for(u_int32_t nodeIndex_neighbor : nodeNeighbors){


                if(neighbors.find(nodeIndex_neighbor) != neighbors.end()){
                    continue;
                }

                //if(visitedNodes.find(nn) != visitedNodes.end()){
                //    node = node->next;
                //    continue;
                //}

                if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                    break;
                }

                neighbors.insert(nodeIndex_neighbor);
                //isVisited[nn] = true;
                distance[nodeIndex_neighbor] = distance[n] + 1;
                //neighbors.push_back(nn);


                if(distance[nodeIndex_neighbor] >= maxDistance){
                    continue;
                }

                //cout << "Push: " << nn << " " << distance[nn] << endl;
                queue.push(nodeIndex_neighbor);

            }

            if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                break;
            }

        }



    }

    void extractComponent(u_int32_t nodeName){
        
        u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);

        vector<u_int32_t> unitigNodes;
        getUnitigNodes(nodeIndex_to_unitig(nodeIndex), unitigNodes);
        for(u_int32_t nIndex : unitigNodes){
            u_int32_t nName = BiGraph::nodeIndex_to_nodeName(nIndex);
            cout << nName << " " << _nodeAbundances[nName] << endl;
        }


        cout << nodeName << " " << nodeIndex_to_unitig(nodeIndex)._abundance << " " << nodeIndex_to_unitig(nodeIndex)._length << endl;
        //cout << nodeIndex_to_unitig(nodeIndex)._abundance << endl;
        
        unordered_set<u_int32_t> nodes;
		collectNeighbors(nodeIndex, 200, 10000, nodes);

        unordered_set<u_int32_t> validNodes;

        for (auto& nodeIndex : nodes){
            u_int32_t nodeName = _graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
            validNodes.insert(nodeName);
        }

        GfaParser::rewriteGfa_withoutNodes(_inputGfaFilename, _inputGfaFilename + "_component.gfa", validNodes, _isEdgeRemoved, _graphSuccessors);
		//exit(1);
    }
    

    void getStronglyConnectedComponent(u_int32_t unitigIndex_source, bool forward, vector<u_int32_t>& component){ //, vector<u_int32_t>& scc
        //scc.clear();
    
        component.clear();
        //vector<vector<u_int32_t>> sccs;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        //u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(nodeIndex_source);

        queue.clear();
        queue.push_back(unitigIndex_source);

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

            if(preorder.find(v) == preorder.end()){
                i += 1;
                preorder[v] = i;
            }

            bool done = true;

            vector<u_int32_t> v_nbrs;
            if(forward){
                getSuccessors_unitig(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            else{
                getPredecessors_unitig(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            
            //cout << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << v_nbrs.size() << endl;

            for(u_int32_t w : v_nbrs){
                if(preorder.find(w) == preorder.end()){
                    queue.push_back(w);
                    done = false;
                    break;
                }
            }

            if(done){
                lowlink[v] = preorder[v];
                for(u_int32_t w : v_nbrs){
                    if(scc_found.find(w) == scc_found.end()){
                        if(preorder[w] > preorder[v]){
                            lowlink[v] = min(lowlink[v], lowlink[w]);
                        }
                        else{
                            lowlink[v] = min(lowlink[v], preorder[w]);
                        }
                    }
                }

                queue.pop_back();

                if(lowlink[v] == preorder[v]){
                    //cout << "----" << endl;
                    scc_found.insert(v);
                    vector<u_int32_t> scc = {v};

                    bool isSourceScc = false;

                    if(v == unitigIndex_source) isSourceScc = true;

                    while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                        u_int32_t k = scc_queue[scc_queue.size()-1];
                        scc_queue.pop_back();
                        scc_found.insert(k);
                        scc.push_back(k);

                        //cout << "Lala: " << unitigIndex_source << " " << k << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[k]._startNode) << endl;

                        if(k == unitigIndex_source) isSourceScc = true;
                    }

                    if(isSourceScc){
                        for(u_int32_t unitigIndex : scc){
                            component.push_back(unitigIndex);
                        }
                        //cout << "found" << endl;
                        //component = scc;
                        return;
                    }
                    //sccs.push_back(scc);
                }
                else{
                    scc_queue.push_back(v);
                }
            }
        }
        
        /*
		for(size_t i=0; i<sccs.size(); i++){
            
            vector<u_int32_t>& scc = sccs[i];

                //cout << scc.size() << endl;
            if(std::find(scc.begin(), scc.end(), unitigIndex_source) != scc.end()){
                //cout << "miam" << endl;
                component = scc;
                return;
            }

        }
        */


        /*
        u_int32_t largestComponent = 0;
        size_t largestComponentIndex = -1;
		for(size_t i=0; i<sccs.size(); i++){
            
            vector<u_int32_t>& scc = sccs[i];
            
            if(scc.size() > largestComponent){
                largestComponent = scc.size();
                largestComponentIndex = i;
            }
        }

        component = sccs[largestComponentIndex];
        */

        //for(u_int32_t nodeIndex : sccs[largestComponentIndex]){
        //    component.insert(nodeIndex);
        //}

        /*
        cout << "Nb sccs: " << sccs.size() << endl;
		for(vector<u_int32_t>& scc : sccs){
            //cout << scc.size() << endl;
			cout << "----" << endl;
			for(u_int32_t unitigIndex : scc){
				vector<u_int32_t> nodes;
				getUnitigNodes(_unitigs[unitigIndex], nodes);
				for(u_int32_t nodeIndex : nodes){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				}
			}
		}*/


    }

    void getStronglyConnectedComponents(u_int32_t unitigIndex_source, bool forward, vector<vector<u_int32_t>>& components){ //, vector<u_int32_t>& scc
        //scc.clear();
    
        components.clear();
        //vector<vector<u_int32_t>> sccs;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        //u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(nodeIndex_source);

        queue.clear();
        queue.push_back(unitigIndex_source);

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

            if(preorder.find(v) == preorder.end()){
                i += 1;
                preorder[v] = i;
            }

            bool done = true;

            vector<u_int32_t> v_nbrs;
            if(forward){
                getSuccessors_unitig(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            else{
                getPredecessors_unitig(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            
            //cout << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << v_nbrs.size() << endl;

            for(u_int32_t w : v_nbrs){
                if(preorder.find(w) == preorder.end()){
                    queue.push_back(w);
                    done = false;
                    break;
                }
            }

            if(done){
                lowlink[v] = preorder[v];
                for(u_int32_t w : v_nbrs){
                    if(scc_found.find(w) == scc_found.end()){
                        if(preorder[w] > preorder[v]){
                            lowlink[v] = min(lowlink[v], lowlink[w]);
                        }
                        else{
                            lowlink[v] = min(lowlink[v], preorder[w]);
                        }
                    }
                }

                queue.pop_back();

                if(lowlink[v] == preorder[v]){
                    //cout << "----" << endl;
                    scc_found.insert(v);
                    vector<u_int32_t> scc = {v};

                    bool isSourceScc = false;

                    if(v == unitigIndex_source) isSourceScc = true;

                    while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                        u_int32_t k = scc_queue[scc_queue.size()-1];
                        scc_queue.pop_back();
                        scc_found.insert(k);
                        scc.push_back(k);

                        //cout << "Lala: " << unitigIndex_source << " " << k << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[k]._startNode) << endl;

                        if(k == unitigIndex_source) isSourceScc = true;
                    }

                    if(isSourceScc){
                        //cout << "found" << endl;
                        //component = scc;
                        //return;
                    }
                    components.push_back(scc);
                }
                else{
                    scc_queue.push_back(v);
                }
            }
        }
        
   

    }

    void getStronglyConnectedComponent_node(u_int32_t nodeIndex_source, bool forward, vector<u_int32_t>& component){ //, vector<u_int32_t>& scc
        //scc.clear();
    
        component.clear();
        //vector<vector<u_int32_t>> sccs;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(nodeIndex_source);

        queue.clear();
        queue.push_back(nodeIndex_source);

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

            if(preorder.find(v) == preorder.end()){
                i += 1;
                preorder[v] = i;
            }

            bool done = true;

            vector<u_int32_t> v_nbrs;
            if(forward){
                getSuccessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            else{
                getPredecessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            
            //cout << BiGraph::nodeIndex_to_nodeName(v) << " " << v_nbrs.size() << endl;
            //getchar();

            for(u_int32_t w : v_nbrs){
                if(preorder.find(w) == preorder.end()){
                    queue.push_back(w);
                    done = false;
                    break;
                }
            }

            if(done){
                lowlink[v] = preorder[v];
                for(u_int32_t w : v_nbrs){
                    if(scc_found.find(w) == scc_found.end()){
                        if(preorder[w] > preorder[v]){
                            lowlink[v] = min(lowlink[v], lowlink[w]);
                        }
                        else{
                            lowlink[v] = min(lowlink[v], preorder[w]);
                        }
                    }
                }

                queue.pop_back();

                if(lowlink[v] == preorder[v]){
                    //cout << "----" << endl;
                    scc_found.insert(v);
                    vector<u_int32_t> scc = {v};

                    bool isSourceScc = false;

                    if(nodeIndex_to_unitigIndex(v) == unitigIndex_source) isSourceScc = true;

                    while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                        u_int32_t k = scc_queue[scc_queue.size()-1];
                        scc_queue.pop_back();
                        scc_found.insert(k);
                        scc.push_back(k);

                        //cout << "Lala: " << unitigIndex_source << " " << k << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[k]._startNode) << endl;

                        if(nodeIndex_to_unitigIndex(k) == unitigIndex_source) isSourceScc = true;
                    }

                    if(isSourceScc){
                        for(u_int32_t nodeIndex : scc){
                            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                            if(std::find(component.begin(), component.end(), unitigIndex) != component.end()) continue;
                            component.push_back(unitigIndex);
                        }
                        //cout << "found" << endl;
                        //component = scc;
                        return;
                    }
                    
                }
                else{
                    scc_queue.push_back(v);
                }
            }
        }
        
        /*
        //unordered_map<u_int32_t, u_int32_t> unitigTo;
        unordered_map<u_int32_t, vector<u_int32_t>> unitigToNodes;

        for(vector<u_int32_t>& scc : components_tmp){
            for(u_int32_t nodeIndex : scc){
                
                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                unitigToNodes[unitigIndex].push_back(nodeIndex);
            }
        }


        //unordered_map<u_int32_t, u_int32_t> unitigToScc;
        unordered_set<u_int32_t> nodeVisited;


        for(vector<u_int32_t>& scc : components_tmp){

            vector<u_int32_t> sccNew;
            
            for(u_int32_t nodeIndex : scc){
                if(nodeVisited.find(nodeIndex) != nodeVisited.end()) continue;

                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                for(u_int32_t n : unitigToNodes[unitigIndex]){
                    sccNew.push_back(n);
                    nodeVisited.insert(n);
                }
            }

            if(sccNew.size() == 0) continue;

            components.push_back(sccNew);

        }*/
    }

    void getStronglyConnectedComponents_node_2(u_int32_t nodeIndex_source, bool forward, vector<vector<u_int32_t>>& components){ //, vector<u_int32_t>& scc
        //scc.clear();
    
        components.clear();
        //vector<vector<u_int32_t>> sccs;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(nodeIndex_source);

        queue.clear();
        queue.push_back(nodeIndex_source);

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

            if(preorder.find(v) == preorder.end()){
                i += 1;
                preorder[v] = i;
            }

            bool done = true;

            vector<u_int32_t> v_nbrs;
            if(forward){
                getSuccessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            else{
                getPredecessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            
            //cout << BiGraph::nodeIndex_to_nodeName(v) << " " << v_nbrs.size() << endl;
            //getchar();

            for(u_int32_t w : v_nbrs){
                if(preorder.find(w) == preorder.end()){
                    queue.push_back(w);
                    done = false;
                    break;
                }
            }

            if(done){
                lowlink[v] = preorder[v];
                for(u_int32_t w : v_nbrs){
                    if(scc_found.find(w) == scc_found.end()){
                        if(preorder[w] > preorder[v]){
                            lowlink[v] = min(lowlink[v], lowlink[w]);
                        }
                        else{
                            lowlink[v] = min(lowlink[v], preorder[w]);
                        }
                    }
                }

                queue.pop_back();

                if(lowlink[v] == preorder[v]){
                    //cout << "----" << endl;
                    scc_found.insert(v);
                    vector<u_int32_t> scc = {v};

                    bool isSourceScc = false;

                    if(nodeIndex_to_unitigIndex(v) == unitigIndex_source) isSourceScc = true;

                    while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                        u_int32_t k = scc_queue[scc_queue.size()-1];
                        scc_queue.pop_back();
                        scc_found.insert(k);
                        scc.push_back(k);

                        //cout << "Lala: " << unitigIndex_source << " " << k << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[k]._startNode) << endl;

                        if(nodeIndex_to_unitigIndex(k) == unitigIndex_source) isSourceScc = true;
                    }

                    vector<u_int32_t> component;
                    //if(isSourceScc){
                        for(u_int32_t nodeIndex : scc){
                            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                            if(std::find(component.begin(), component.end(), unitigIndex) != component.end()) continue;
                            component.push_back(unitigIndex);
                        }
                        //cout << "found" << endl;
                        //component = scc;
                        //return;
                    //}
                    components.push_back(component);
                }
                else{
                    scc_queue.push_back(v);
                }
            }
        }
        
        /*
        //unordered_map<u_int32_t, u_int32_t> unitigTo;
        unordered_map<u_int32_t, vector<u_int32_t>> unitigToNodes;

        for(vector<u_int32_t>& scc : components_tmp){
            for(u_int32_t nodeIndex : scc){
                
                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                unitigToNodes[unitigIndex].push_back(nodeIndex);
            }
        }


        //unordered_map<u_int32_t, u_int32_t> unitigToScc;
        unordered_set<u_int32_t> nodeVisited;


        for(vector<u_int32_t>& scc : components_tmp){

            vector<u_int32_t> sccNew;
            
            for(u_int32_t nodeIndex : scc){
                if(nodeVisited.find(nodeIndex) != nodeVisited.end()) continue;

                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                for(u_int32_t n : unitigToNodes[unitigIndex]){
                    sccNew.push_back(n);
                    nodeVisited.insert(n);
                }
            }

            if(sccNew.size() == 0) continue;

            components.push_back(sccNew);

        }*/
    }

    void getStronglyConnectedComponents_node(u_int32_t nodeIndex_source, bool forward, vector<vector<u_int32_t>>& components){ //, vector<u_int32_t>& scc
        //scc.clear();
    
        vector<vector<u_int32_t>> components_tmp;
        components.clear();
        //vector<vector<u_int32_t>> sccs;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        //u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(nodeIndex_source);

        queue.clear();
        queue.push_back(nodeIndex_source);

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

            if(preorder.find(v) == preorder.end()){
                i += 1;
                preorder[v] = i;
            }

            bool done = true;

            vector<u_int32_t> v_nbrs;
            if(forward){
                getSuccessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            else{
                getPredecessors(v, 0, v_nbrs); //_unitigs[v]._endNode
            }
            
            //cout << BiGraph::nodeIndex_to_nodeName(v) << " " << v_nbrs.size() << endl;
            //getchar();

            for(u_int32_t w : v_nbrs){
                if(preorder.find(w) == preorder.end()){
                    queue.push_back(w);
                    done = false;
                    break;
                }
            }

            if(done){
                lowlink[v] = preorder[v];
                for(u_int32_t w : v_nbrs){
                    if(scc_found.find(w) == scc_found.end()){
                        if(preorder[w] > preorder[v]){
                            lowlink[v] = min(lowlink[v], lowlink[w]);
                        }
                        else{
                            lowlink[v] = min(lowlink[v], preorder[w]);
                        }
                    }
                }

                queue.pop_back();

                if(lowlink[v] == preorder[v]){
                    //cout << "----" << endl;
                    scc_found.insert(v);
                    vector<u_int32_t> scc = {v};

                    bool isSourceScc = false;

                    if(v == nodeIndex_source) isSourceScc = true;

                    while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                        u_int32_t k = scc_queue[scc_queue.size()-1];
                        scc_queue.pop_back();
                        scc_found.insert(k);
                        scc.push_back(k);

                        //cout << "Lala: " << unitigIndex_source << " " << k << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[k]._startNode) << endl;

                        if(k == nodeIndex_source) isSourceScc = true;
                    }

                    if(isSourceScc){
                        //cout << "found" << endl;
                        //component = scc;
                        //return;
                    }
                    components_tmp.push_back(scc);
                }
                else{
                    scc_queue.push_back(v);
                }
            }
        }
        
        //unordered_map<u_int32_t, u_int32_t> unitigTo;
        unordered_map<u_int32_t, vector<u_int32_t>> unitigToNodes;

        for(vector<u_int32_t>& scc : components_tmp){
            for(u_int32_t nodeIndex : scc){
                
                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                unitigToNodes[unitigIndex].push_back(nodeIndex);
            }
        }


        //unordered_map<u_int32_t, u_int32_t> unitigToScc;
        unordered_set<u_int32_t> nodeVisited;


        for(vector<u_int32_t>& scc : components_tmp){

            vector<u_int32_t> sccNew;
            
            for(u_int32_t nodeIndex : scc){
                if(nodeVisited.find(nodeIndex) != nodeVisited.end()) continue;

                u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
                for(u_int32_t n : unitigToNodes[unitigIndex]){
                    sccNew.push_back(n);
                    nodeVisited.insert(n);
                }
            }

            if(sccNew.size() == 0) continue;

            components.push_back(sccNew);

        }
        /*
        for(auto& it : preorder){
            u_int32_t nodeIndex = it.first;
            cout << BiGraph::nodeIndex_to_nodeName(nodeIndex);

            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
            if(unitigToScc.find(unitigIndex) == unitigToScc.end()){
                unitigToScc[unitigIndex] = it.second;
            }


        }*/
    }


    u_int32_t shortestPath_unitig(u_int32_t source_unitigIndex, u_int32_t sink_unitigIndex, vector<u_int32_t>& path, bool includeSource, bool includeSink){

        path.clear();
		unordered_set<u_int32_t> isVisited;
		unordered_map<u_int32_t, u_int32_t> distance;
		unordered_map<u_int32_t, u_int32_t> prev;
        queue<u_int32_t> queue;


        distance[source_unitigIndex] = 0;
        isVisited.insert(source_unitigIndex);
		prev[source_unitigIndex] = -1;

        queue.push(source_unitigIndex);
        bool found = false;

        while (!queue.empty() && !found){

            u_int32_t unitigIndex_current = queue.front();
            queue.pop();

			vector<u_int32_t> successors;
			getSuccessors_unitig(unitigIndex_current, 0, successors);


			for(u_int32_t unitigIndex_successor : successors){

                if (isVisited.find(unitigIndex_successor) != isVisited.end()) continue;

                distance[unitigIndex_successor] = distance[unitigIndex_current] + 1;
                queue.push(unitigIndex_successor);
                isVisited.insert(unitigIndex_successor);
                prev[unitigIndex_successor] = unitigIndex_current;

                if(unitigIndex_successor == sink_unitigIndex){
                    found = true;
                    break;
                }

            }



        }

        if(found){

            u_int32_t n = sink_unitigIndex;
            while(n != source_unitigIndex){
                if(n == sink_unitigIndex){
                    if(includeSink) path.push_back(n);
                }
                else{
                    path.push_back(n);
                }

                n = prev[n];
            }
            if(includeSource) path.push_back(source_unitigIndex);

            return distance[sink_unitigIndex];
        }

        return -1;
    }
    

    void collectNodes_betweenSourceSink_unitig(u_int32_t source_unitigIndex, u_int32_t sink_unitigIndex, vector<u_int32_t>& nodes, float abundanceCutoff_min, const vector<UnitigData>& _unitigDatas, bool forward){


        u_int32_t source_nodeName = BiGraph::nodeIndex_to_nodeName(_unitigs[source_unitigIndex]._startNode);

        nodes.clear();

		unordered_set<u_int32_t> isVisited;
        queue<u_int32_t> queue;

        isVisited.insert(source_unitigIndex);
        isVisited.insert(sink_unitigIndex);
        queue.push(source_unitigIndex);

        while (!queue.empty()){

            u_int32_t unitigIndex_current = queue.front();
            queue.pop();

			vector<u_int32_t> successors;
            if(forward){
			    getSuccessors_unitig(unitigIndex_current, 0, successors);
            }
            else{
			    getPredecessors_unitig(unitigIndex_current, 0, successors);
            }

			for(u_int32_t unitigIndex_successor : successors){

                if (isVisited.find(unitigIndex_successor) != isVisited.end()) continue;

                if(abundanceCutoff_min > 0){
                    u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_successor]._startNode);
                    u_int64_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[source_nodeName], _unitigDatas[nodeName]);
                    //cout << nodeName_neighbor << " " << nbSharedReads << endl;
                    if(nbSharedReads ==0 || nbSharedReads < abundanceCutoff_min) continue;
                }

                queue.push(unitigIndex_successor);
                isVisited.insert(unitigIndex_successor);
                nodes.push_back(unitigIndex_successor);
                //cout << "\t\tCollected node: " <<  BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex_successor]._startNode) << endl;

            }
        }

    }

    unordered_set<u_int32_t> _allowedNodeIndex;
    void disconnectSubGraph(u_int32_t source_nodeIndex, u_int64_t maxLength, vector<UnitigData>& _unitigDatas, float abundanceCutoff_min, bool forward){

        _allowedNodeIndex.clear();
        u_int32_t source_nodeName = BiGraph::nodeIndex_to_nodeName(source_nodeIndex);

        _isNodeInvalid_tmp.clear();

		vector<u_int32_t> predecessors;
        if(forward){
            getPredecessors(source_nodeIndex, 0, predecessors);
		    //getPredecessors_unitig(nodeIndex_to_unitigIndex(source_nodeIndex), predecessors);
        }
        else{
            getSuccessors(source_nodeIndex, 0, predecessors);
		    //getSuccessors_unitig(nodeIndex_to_unitigIndex(source_nodeIndex), predecessors);
        }

		for(u_int32_t unitigIndex : predecessors){
            disconnectUnitig(unitigIndex);
		}

        /*
        if(forward){
            source_nodeIndex = _unitigs[nodeIndex_to_unitigIndex(source_nodeIndex)]._endNode;
        }
        else{
            source_nodeIndex = _unitigs[nodeIndex_to_unitigIndex(source_nodeIndex)]._startNode; 
        }*/
		ofstream file_scc("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		file_scc << "Name,Colour" << endl;

        unordered_set<u_int32_t> isVisited;
        //unordered_set<u_int32_t> visitedUnitigs;
        queue <u_int32_t> queue;
        unordered_map<u_int32_t, u_int32_t> distance;

        queue.push(source_nodeIndex);
        isVisited.insert(source_nodeIndex);
        distance[source_nodeIndex] = 0;

        while (!queue.empty()){

            u_int64_t nodeIndex = queue.front();
            file_scc << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
            queue.pop();

            vector<AdjNode> neighbors;
            if(forward){
                getSuccessors_overlap(nodeIndex, 0, neighbors);
            }
            else{
                getPredecessors_overlap(nodeIndex, 0, neighbors);
            }

            //cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << neighbors.size() << endl;
            
            for(AdjNode& node_neighbor : neighbors){

                u_int32_t nodeIndex_neighbor = node_neighbor._index;
                u_int32_t nodeName_neighbor = BiGraph::nodeIndex_to_nodeName(nodeIndex_neighbor);

                //cout << nodeName_neighbor << " " << (isVisited.find(nodeIndex_neighbor) != isVisited.end()) << " " << distance[nodeIndex_neighbor] << endl;


                //cout << "2: " << _nodeToUnitig[nodeIndex_neighbor] << endl;
                //u_int32_t unitigIndex_neighbor = nodeIndex_to_unitigIndex(nodeIndex_neighbor);
                if (isVisited.find(nodeIndex_neighbor) != isVisited.end()) continue;

                distance[nodeIndex_neighbor] = distance[nodeIndex] + (_nodeLengths[BiGraph::nodeIndex_to_nodeName(nodeIndex_neighbor)] - node_neighbor._overlap);

                //cout << nodeName_neighbor << " " << distance[nodeIndex_neighbor] << endl;
                if(distance[nodeIndex_neighbor] > maxLength) continue;
                
                if(abundanceCutoff_min > 0){
				     u_int64_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[source_nodeName], _unitigDatas[nodeName_neighbor]);
                    //cout << nodeName_neighbor << " " << nbSharedReads << endl;
                    if(nbSharedReads == 0 || nbSharedReads < abundanceCutoff_min) continue;
                }

                queue.push(nodeIndex_neighbor);
                //visitedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex_neighbor));
                //visitedUnitigs.insert(nodeIndex_to_unitigIndex(nodeIndex_toReverseDirection(nodeIndex_neighbor)));
                isVisited.insert(nodeIndex_neighbor);

                //cout << "Visited: " << BiGraph::nodeToString(nodeIndex_neighbor) << " " << nodeIndex_to_unitigIndex(nodeIndex_neighbor) << "             " << distance[nodeIndex_neighbor] << "    " << _nodeLengths[BiGraph::nodeIndex_to_nodeName(nodeIndex_neighbor)] << " " << node_neighbor._overlap << endl;
                //component.push_back(nodeIndex_neighbor);
            }

        }

        for(u_int32_t nodeIndex : isVisited){

            _allowedNodeIndex.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
            /*
            //cout << "miam: " << unitigIndex << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._startNode) << endl;
            //cout << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._endNode) << endl;

            vector<u_int32_t> successors;
            if(forward){
                getSuccessors(nodeIndex, 0, successors);
            }
            else{
                getPredecessors(nodeIndex, 0, successors);
            }

            for(u_int32_t successor_nodeIndex : successors){
                //u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
                if(isVisited.find(successor_nodeIndex) != isVisited.end()) continue; 
                if(isVisited.find(nodeIndex_toReverseDirection(successor_nodeIndex)) != isVisited.end()) continue; 

                //cout << unitigIndex << " " << (visitedUnitigs.find(unitigIndex) != visitedUnitigs.end()) << endl;
                //if(visitedUnitigs.find(nodeIndex) != visitedUnitigs.end()) continue; 
                //if(isVisited.find(_unitigs[unitigIndex]._startNode) != isVisited.end() ||  isVisited.find(_unitigs[unitigIndex]._endNode) != isVisited.end()) continue;
                //if(isVisited.find(nodeIndex_toReverseDirection(_unitigs[unitigIndex]._startNode)) != isVisited.end() ||  isVisited.find(nodeIndex_toReverseDirection(_unitigs[unitigIndex]._endNode)) != isVisited.end()) continue;
                //cout << BiGraph::nodeToString(_unitigs[unitigIndex]._startNode) << " " << BiGraph::nodeToString(_unitigs[unitigIndex]._endNode) << endl;
                disconnectUnitig(successor_nodeIndex);
            }*/
        }

        //_allowedNodeIndex = isVisited;
        file_scc.close();
        /*
        unordered_set<u_int32_t> unitigs;


		ofstream file_scc("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		file_scc << "Name,Colour" << endl;

        for(u_int32_t nodeIndex : isVisited){
            file_scc << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;

            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
            unitigs.insert(unitigIndex);
        }
        file_scc.close();

        for(u_int32_t unitigIndex : unitigs){



            //cout << "miam: " << unitigIndex << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._startNode) << endl;
            //cout << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._endNode) << endl;

            vector<u_int32_t> neighbors;
            if(forward){
                getSuccessors_unitig(unitigIndex, neighbors);
            }
            else{
                getPredecessors_unitig(unitigIndex, neighbors);
            }

            for(u_int32_t unitigIndex : neighbors){
                cout << unitigIndex << " " << (visitedUnitigs.find(unitigIndex) != visitedUnitigs.end()) << endl;
                if(visitedUnitigs.find(unitigIndex) != visitedUnitigs.end()) continue; 
                //if(isVisited.find(_unitigs[unitigIndex]._startNode) != isVisited.end() ||  isVisited.find(_unitigs[unitigIndex]._endNode) != isVisited.end()) continue;
                //if(isVisited.find(nodeIndex_toReverseDirection(_unitigs[unitigIndex]._startNode)) != isVisited.end() ||  isVisited.find(nodeIndex_toReverseDirection(_unitigs[unitigIndex]._endNode)) != isVisited.end()) continue;
                //cout << BiGraph::nodeToString(_unitigs[unitigIndex]._startNode) << " " << BiGraph::nodeToString(_unitigs[unitigIndex]._endNode) << endl;
                disconnectUnitig(unitigIndex);
            }
        }
        */
        //exit(1);
    }

    void disconnectUnitig(u_int32_t unitigIndex){
        /*
        u_int32_t startNode = _unitigs[unitigIndex]._startNode;
        _isNodeInvalid_tmp.insert(nodeIndex_to_unitigIndex(startNode));
        _isNodeInvalid_tmp.insert(nodeIndex_to_unitigIndex(nodeIndex_toReverseDirection(startNode)));
        
        cout << "Disconnected: " << BiGraph::nodeIndex_to_nodeName(startNode) << endl;
        */
        _isNodeInvalid_tmp.insert(unitigIndex);
        _isNodeInvalid_tmp.insert(nodeIndex_toReverseDirection(unitigIndex));
        //cout << "Disconnected: " << BiGraph::nodeIndex_to_nodeName(unitigIndex) << endl;
    }

    void extractSubGraph(u_int32_t source_nodeIndex, vector<u_int32_t>& successor_nodeIndex, u_int64_t maxLength, bool forward){

        //cout << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		ofstream file_scc("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		file_scc << "Name,Colour" << endl;

        /*
        for(u_int32_t nodeIndex : sourceNodes){
            vector<u_int32_t> scc;
            getStronglyConnectedComponent(nodeIndex, scc);

        }
        */

        //u_int32_t source_nodeName = BiGraph::nodeIndex_to_nodeName(source_nodeIndex);

        _isNodeInvalid_tmp.clear();

        unordered_set<u_int32_t> isVisited;
        queue <u_int32_t> queue;
        unordered_map<u_int32_t, u_int32_t> distance;

        queue.push(source_nodeIndex);
        isVisited.insert(source_nodeIndex);
        distance[source_nodeIndex] = 0;

        while (!queue.empty()){

            u_int64_t nodeIndex = queue.front();

            queue.pop();

            vector<AdjNode> neighbors;
            //vector<AdjNode> successors;
            //vector<AdjNode> predecessors;
            if(forward){
                getSuccessors_overlap(nodeIndex, 0, neighbors);
            }
            else{
                getPredecessors_overlap(nodeIndex, 0, neighbors);
            }

            //for(AdjNode& node : successors) neighbors.push_back(node);
            //for(AdjNode& node : predecessors) neighbors.push_back(node);

            
            for(AdjNode& node_neighbor : neighbors){

                u_int32_t nodeIndex_neighbor = node_neighbor._index;
                if (isVisited.find(nodeIndex_neighbor) != isVisited.end()) continue;

                distance[nodeIndex_neighbor] = distance[nodeIndex] + (_nodeLengths[BiGraph::nodeIndex_to_nodeName(nodeIndex_neighbor)] - node_neighbor._overlap);

                if(distance[nodeIndex_neighbor] > maxLength) continue;

                queue.push(nodeIndex_neighbor);
                isVisited.insert(nodeIndex_neighbor);
            }

        }


        for(u_int32_t nodeIndex : isVisited){
            file_scc << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
        }
        /*
        unordered_set<u_int32_t> unitigs;

        for(u_int32_t nodeIndex : isVisited){
            u_int32_t unitigIndex = nodeIndex_to_unitigIndex(nodeIndex);
            unitigs.insert(unitigIndex);
        }

        for(u_int32_t unitigIndex : unitigs){



            //cout << unitigIndex << endl;
            //cout << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._endNode) << endl;

            vector<u_int32_t> neighbors;
            if(forward){
                getSuccessors_unitig(unitigIndex, neighbors);
            }
            else{
                getPredecessors_unitig(unitigIndex, neighbors);
            }

            for(u_int32_t unitigIndex : neighbors){
                if(isVisited.find(_unitigs[unitigIndex]._startNode) != isVisited.end() ||  isVisited.find(_unitigs[unitigIndex]._endNode) != isVisited.end()) continue;
                disconnectUnitig(unitigIndex);
            }
        }

        */
        file_scc.close();
    }

    /*
    void create_dfs_order(u_int32_t source_nodeIndex, vector<u_int32_t>& order){

        cout << endl << "-------------" << endl;

        order.clear();
        vector<u_int32_t> stack;

        u_int32_t unitigIndex_source = nodeIndex_to_unitigIndex(source_nodeIndex);
        stack.push_back(unitigIndex_source);

        unordered_map<u_int32_t, u_int8_t> nodeColors;

        static u_int8_t WHITE = 0;
        static u_int8_t GREY = 1;
        static u_int8_t BLACK = 2;

        while(!stack.empty()){
            u_int32_t u = stack[stack.size()-1];

            u_int32_t c = WHITE;
            if(nodeColors.find(u) != nodeColors.end()) c = nodeColors[u];

            cout << "\tStart node: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << c << endl;

            if(c != GREY && c != BLACK){
                nodeColors[u] = GREY;

                vector<u_int32_t> successors;
                getSuccessors_unitig(u, successors);

                cout << "\tNb succ: " << successors.size() << endl;

                for(u_int32_t w : successors){
                    
                    u_int32_t color = WHITE;
                    if(nodeColors.find(w) != nodeColors.end()) color = nodeColors[w];
                    if(color != GREY && color != BLACK){
                        stack.push_back(w);
                    }
                }
            }
            else{
                stack.pop_back();
                if(c == GREY){
                    cout << "Add order: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                    order.push_back(u);
                    nodeColors[u] = BLACK;
                }
            }
        }

        cout << order.size() << endl;
        for(u_int32_t unitigIndex : order){
            cout << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._startNode) << endl;
        }
        cout << "-----------------" << endl;
    }

    int out_parent(u_int32_t k, const vector<u_int32_t>& order, u_int32_t v, u_int32_t v2){
        u_int32_t u = order[k];
        if(u == v2){
            u = v;
        }

        cout << "\tIn degree: "<< BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << in_degree(u) << endl;

        if(in_degree(u) == 0){
            return INT_MAX;
        }
        int maximum = -1;

        vector<u_int32_t> predecessors;
        getPredecessors_unitig(u, predecessors);

        for(u_int32_t w : predecessors){
            int pos = find(order.begin(), order.end(), w) - order.begin();//order.index(w)
            if (pos <= k){
                return INT_MAX;
            }
            maximum = max(maximum, pos);
        }
        return maximum;
    }



    int out_child(u_int32_t k, const vector<u_int32_t>& order, u_int32_t v, u_int32_t v2){

        u_int32_t u = order[k];
        if (u == v2){
            u = v;
        }
        cout << "\tOut degree: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << out_degree(u) << endl;
        if (out_degree(u) == 0){
            return -2;
        }
        int minimum = INT_MAX;

        vector<u_int32_t> successors;
        getSuccessors_unitig(u, successors);

        for(u_int32_t w : successors){
            int pos = 0;
            if(w == v && v2 != -1){
                pos = find(order.begin(), order.end(), v2) - order.begin();
                //pos = order.index(v2);
            }
            else{
                pos = find(order.begin(), order.end(), w) - order.begin();
                //pos = order.index(w);
            }
            if (pos >= k){
                return -2;
            }
            minimum = min(minimum, pos);
        }
        return minimum;
    }

    void superbubble(const vector<u_int32_t>& order, u_int32_t v, u_int32_t v2){

        cout << endl << endl;

        vector<long> stack;
        vector<u_int32_t> out_parent_map;
        static long NONE = -2222;
        long t = NONE;

        for(u_int32_t k=0; k<order.size(); k++){
            cout << "\tK: " << k << endl;
            int child = out_child(k, order, v, v2);
            cout << "\tChild: " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << child << endl;
            if(child == k - 1){
                cout << "Push: " << t << endl; 
                stack.push_back(t);
                t = k -1;
            }
            else{
                while (t != NONE && t > child){
                    //cout << "lala" << endl;
                    long t2 = stack[stack.size()-1];
                    stack.pop_back();
                    if(t2 != NONE){
                        out_parent_map[t2] = max(out_parent_map[t], out_parent_map[t2]);
                    }
                    t = t2;
                }
            }
    
            if (t != NONE && out_parent_map[t] == k){
                //order[o: i + 1][::-1]
                cout << "Superbubble: " << k << " " << t << endl;
                for(size_t i=t; i<k+1; i++){
                    cout << BiGraph::nodeIndex_to_nodeName(_unitigs[order[i]]._startNode) << endl;
                }

                //report(k, t)
                long t2 = stack[stack.size()-1];
                stack.pop_back();
                cout << t2 << endl;
                if (t2 != NONE){
                    out_parent_map[t2] = max(out_parent_map[t], out_parent_map[t2]);
                }
                t = t2;
            }

            out_parent_map.push_back(out_parent(k, order, v, v2));
            if(t != NONE){
                out_parent_map[t] = max(out_parent_map[t], out_parent_map[k]);
            }
        }

        cout << endl << endl;
    }
    */
    /*
    stack = []
    out_parent_map = []
    t = None
    for k in range(len(order)):
        child = out_child(k, g, order, v, v2)
        if child == k - 1:
            stack.append(t)
            t = k - 1
        else:
            while t is not None and t > child:
                t2 = stack.pop()
                if t2 is not None:
                    out_parent_map[t2] = max(out_parent_map[t], out_parent_map[t2])
                t = t2
        if t is not None and out_parent_map[t] == k:
            report(k, t)
            t2 = stack.pop()
            if t2 is not None:
                out_parent_map[t2] = max(out_parent_map[t], out_parent_map[t2])
            t = t2
        out_parent_map.append(out_parent(k, g, order, v, v2))
        if t is not None:
            out_parent_map[t] = max(out_parent_map[t], out_parent_map[k])

    */
};


#endif



#ifndef MDBG_METAG_PROGRESSIVEABUNDANCEFILTER
#define MDBG_METAG_PROGRESSIVEABUNDANCEFILTER


#include "Commons.hpp"
#include "Graph.hpp"


class ProgressiveAbundanceFilter{

public:

    struct SaveState{
        float _abundanceCutoff_min;
        vector<u_int32_t> _nodeNameRemoved;
        vector<DbgEdge> _isEdgeRemoved;
    };

    struct BubbleSide{
        UnitigType _unitigIndex;
        u_int32_t _startNode;
        size_t _nbNodes;
        //u_int32_t _nodeIndexSum;
        //u_int32_t _quality;
    };

    
    static bool BubbleSideComparator(const BubbleSide &a, const BubbleSide &b){

        if(a._nbNodes == b._nbNodes){
            return a._startNode > b._startNode;
        }
        return a._nbNodes > b._nbNodes;

    }

    static bool BubbleSideComparatorRev(const BubbleSide &a, const BubbleSide &b){

        if(a._nbNodes == b._nbNodes){
            return a._startNode > b._startNode;
        }
        return a._nbNodes < b._nbNodes;

    }

    struct Bubble{
        UnitigType _unitigIndexSource;
        UnitigType _unitigIndexExit;
    };


    struct SuperbubbleSource{
        UnitigType _unitigIndex;
        u_int32_t _nodeName;
        u_int32_t _nbMinimizers;
    };

    static bool SuperbubbleComparator(const SuperbubbleSource& p1, const SuperbubbleSource& p2){
        if(p1._nbMinimizers == p2._nbMinimizers){
            return p1._nodeName > p2._nodeName;
        }
        return p1._nbMinimizers > p2._nbMinimizers;
    }

    
    class SuperbubbleRemoverOld{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph2* _unitigGraph2;
        u_int64_t _maxLength;
        u_int64_t _nbRemoved;
        u_int64_t _nbUnitigNameRemoved;

        SuperbubbleRemoverOld(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph2 = _progressiveAbundanceFilter->_unitigGraph2;
            _maxLength = maxLength;
        }


        bool execute(){

            bool isModification = false;

            while(true){


                auto start = high_resolution_clock::now();
                
                _nbRemoved = 0;
                _nbUnitigNameRemoved = 0;
                bool isMod = removeSuperbubble();

        
                Logger::get().debug() << "Nb superbubble removed: " << _nbRemoved << " " << _nbUnitigNameRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
                //cout << "Nb superbubble removed: " << _nbRemoved << " " << _nbUnitigNameRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
            
                //cout << "Nb superbubble removed: " << _nbRemoved << " " << _nbUnitigNameRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
                //cout << "Nb unitigs: " << _unitigGraph2->getNbUnitigs() << endl;
                break;
                //getchar();
                
                if(isMod) isModification = true;
                if(!isMod) break;
            }
            

            return isModification;
        }

        bool removeSuperbubble(){

           

            vector<SuperbubbleSource> queue;

            for(const UnitigType& unitigName : _progressiveAbundanceFilter->_validNodes2){

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[unitigName];

                if(node == nullptr) continue;

                const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(unitigName, false);

                if(_unitigGraph2->nbSuccessors(unitigIndex1) > 1){
                    queue.push_back({unitigIndex1, _unitigGraph2->getStartNode(unitigIndex1), node->_nbMinimizers});
                }

                const UnitigType& unitigIndex2 = UnitigGraph2::unitigName_to_unitigIndex(unitigName, true);

                if(_unitigGraph2->nbSuccessors(unitigIndex2) > 1){
                    queue.push_back({unitigIndex2, _unitigGraph2->getStartNode(unitigIndex2), node->_nbMinimizers});
                }

                //if(_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()){
                //    continue;
                //}

              
                //queue.push_back(node);
                
            }
            
            //std::sort(queue.begin(), queue.end(), SuperbubbleComparator);

            /*
            ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala_1.txt");

            for(size_t i=0; i<queue.size(); i++){

                const SuperbubbleSource& superbubbleSource = queue[i];
                const UnitigType& unitigIndex = superbubbleSource._unitigIndex;

                bool isReversed;
                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[unitigName];

                //debugFile << unitigIndex << " " << node->_startNode_forward << " " << node->_startNode_reverse << endl;
            }

            debugFile.close();
            cout << "done" << endl;
            //getchar();
            */

            //auto rng = std::default_random_engine {};
            //std::shuffle(std::begin(queue), std::end(queue), rng);


            //cout << queue.size() << endl;

            //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala.txt");
            //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala.txt");


            bool isModification = false;
            unordered_set<UnitigType> isUnitigBubble;
            unordered_set<UnitigType> isUnitigIndexBubble;

            int nbProcessed = queue.size();
            vector<Bubble> bubbles;

            #pragma omp parallel for num_threads(_progressiveAbundanceFilter->_nbCores)
            for(size_t i=0; i<queue.size(); i++){
            //for(const SuperbubbleSource& superbubbleSource : queue){

                const SuperbubbleSource& superbubbleSource = queue[i];
                const UnitigType& unitigIndex = superbubbleSource._unitigIndex;

                bool isReversed;
                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[unitigName];

                //#pragma omp atomic
                //nbProcessed += 1;

                if(node == nullptr) continue;
                if(_unitigGraph2->nbSuccessors(unitigIndex) <= 1) continue;

                //bool exists = false;
                //#pragma omp critical(superbubble)
                //{
                    //cout << i << " " << queue.size() << endl;
                //    if(isUnitigBubble.find(unitigName) != isUnitigBubble.end()) exists = true;
                //}
                //if(exists) continue;

                vector<UnitigType> seenNodes;
                UnitigType unitigIndexExit = -1;
                bool foundExit = isSuperbubble(unitigIndex, unitigIndexExit, isUnitigBubble, seenNodes);




                if(!foundExit){
                    //indexNotSuperbubble(node, seenNodes);
                    continue;
                }
                
                if(unitigIndexExit == UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex)){
                    //indexNotSuperbubble(node, seenNodes);
                    continue; //loop side of an inverse repeat
                }
                
                if(unitigIndexExit == unitigIndex){
                    if(node->getLength(_unitigGraph2) < 400000){
                        //indexNotSuperbubble(node, seenNodes);
                        continue;
                    }
                }

                //if(isUnitigBubble.find(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexExit)) != isUnitigBubble.end()){
                    //indexNotSuperbubble(node, seenNodes);
                //   continue;
                //}

                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(nodeExit)) != isUnitigBubble.end()){
                //    indexNotSuperbubble(node, seenNodes);
                //    continue;
                //}


                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                if(std::find(successors.begin(), successors.end(), unitigIndexExit) != successors.end()){
                    continue;
                }

                //bool isInvalid = false;
                //for(const UnitigType& unitigIndexSuccessor : successors){
                //    if(unitigIndexSuccessor == unitigIndexExit){
                //        isInvalid = true;
                //        break;
                //    }
                //}
                //if(isInvalid){
                //    //indexNotSuperbubble(node, seenNodes);
                //    continue;
                //}



                //debugFile << "Superbubble: " << unitigIndex << " " << unitigIndexExit << " " << isUnitigIndexBubble.size() << endl;
                //getchar();

                #pragma omp critical(superbubble)
                {


                    if(unitigIndex < unitigIndexExit){
                        vector<UnitigType> unitigIndexToRemove;
                        collapseSuperbubble2(unitigIndex, unitigIndexExit, unitigIndexToRemove);
                        for(const UnitigType& unitigIndex : unitigIndexToRemove){
                            isUnitigBubble.insert(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex));
                        }
                        //collapseSuperbubble(unitigIndex, unitigIndexExit, isUnitigBubble, isUnitigIndexBubble);
                        bubbles.push_back({unitigIndex, unitigIndexExit});
                    }
                    else{
                        vector<UnitigType> unitigIndexToRemove;
                        collapseSuperbubble2(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexExit), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex), unitigIndexToRemove);
                        for(const UnitigType& unitigIndex : unitigIndexToRemove){
                            isUnitigBubble.insert(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex));
                        }
                        bubbles.push_back({UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexExit), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex)});
                        //collapseSuperbubble(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexExit), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex), isUnitigBubble, isUnitigIndexBubble);
                    }
                
                    isModification = true;
                    _nbRemoved += 1;
                }

            }


            unordered_set<UnitigType> allUnitigIndexToRemove;

            #pragma omp parallel for num_threads(_progressiveAbundanceFilter->_nbCores)
            for(size_t i=0; i<bubbles.size(); i++){

                const Bubble& bubble = bubbles[i];
                const UnitigType& unitigIndexSource = bubble._unitigIndexSource;
                const UnitigType& unitigIndexExit = bubble._unitigIndexExit;


                const UnitigType& unitigNameSource = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSource);
                const UnitigType& unitigNameExit = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexExit);

                //if(_unitigGraph2->_unitigs[unitigNameSource] == nullptr) continue;
                //if(_unitigGraph2->_unitigs[unitigNameExit] == nullptr) continue;

                if(isUnitigBubble.find(unitigNameSource) != isUnitigIndexBubble.end()) continue;
                if(isUnitigBubble.find(unitigNameExit) != isUnitigIndexBubble.end()) continue;

                //unordered_set<UnitigType> keepNodes_forward;
                //getSuperbubblePath(unitigIndexSource, unitigIndexExit, keepNodes_forward);

                //unordered_set<UnitigType> keepNodes_reverse;
                //getSuperbubblePath(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexExit), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSource), keepNodes_reverse);

                vector<UnitigType> unitigIndexToRemove;

                //if(keepNodes_forward.size() > keepNodes_reverse.size()){
                collapseSuperbubble2(unitigIndexSource, unitigIndexExit, unitigIndexToRemove);
                //}
                //else{
                  //  collapseSuperbubble2(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexExit), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSource), unitigIndexToRemove);
                //}

                #pragma omp critical(superbubble)
                {
                    for(const UnitigType& unitigIndex : unitigIndexToRemove){
                        allUnitigIndexToRemove.insert(unitigIndex);
                    }
                }
            }


            unordered_set<UnitigType> recompactNodes;

    
            for(const UnitigType& unitigIndex : allUnitigIndexToRemove){
                
                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
                
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                //vector<UnitigType> succs = node->_successors;
                //vector<UnitigType> preds = node->_predecessors;


                vector<UnitigType> succs;
                _unitigGraph2->getSuccessors(unitigIndex, succs);

                vector<UnitigType> preds;
                _unitigGraph2->getPredecessors(unitigIndex, preds);

                //cout << "Removing node: " << unitigIndex << " " << unitigName << endl;
                _unitigGraph2->removeNode(_unitigGraph2->_unitigs[unitigName]);

                for(const UnitigType& predecessor : preds){
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessor)] == nullptr) continue;
                    recompactNodes.insert(predecessor);
                }

                for(const UnitigType& successor : succs){
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successor)] == nullptr) continue;
                    recompactNodes.insert(UnitigGraph2::unitigIndex_to_reverseDirection(successor));
                }

            }

            /*
            unordered_set<UnitigType> recompactNodes;

            for(size_t i=0; i<queue.size(); i++){

                const SuperbubbleSource& superbubbleSource = queue[i];
                const UnitigType& unitigIndex = superbubbleSource._unitigIndex;


                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
                
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;
                if(isUnitigIndexBubble.find(unitigIndex) != isUnitigIndexBubble.end()) continue;



                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                int nbBubbleSuccessors = 0;

                for(const UnitigType& unitigIndexSuccessor : successors){
                    if(isUnitigIndexBubble.find(unitigIndexSuccessor) != isUnitigIndexBubble.end()){
                        nbBubbleSuccessors += 1;
                    }
                }
                
                if(nbBubbleSuccessors < 2 || nbBubbleSuccessors != successors.size()) continue;

                //cout << "Collapsing: " << i << " " << unitigIndex << " " << unitigName << endl;

                vector<UnitigType> unitigIndexToRemove;
                collapseSuperbubble2(unitigIndex, isUnitigIndexBubble, unitigIndexToRemove);

                for(const UnitigType& unitigIndex : unitigIndexToRemove){
                    
                    const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
                    
                    if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                    //vector<UnitigType> succs = node->_successors;
                    //vector<UnitigType> preds = node->_predecessors;


                    vector<UnitigType> succs;
                    _unitigGraph2->getSuccessors(unitigIndex, succs);

                    vector<UnitigType> preds;
                    _unitigGraph2->getPredecessors(unitigIndex, preds);

                    //cout << "Removing node: " << unitigIndex << " " << unitigName << endl;
                    _unitigGraph2->removeNode(_unitigGraph2->_unitigs[unitigName]);

                    for(const UnitigType& predecessor : preds){
                        if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessor)] == nullptr) continue;
                        recompactNodes.insert(predecessor);
                    }

                    for(const UnitigType& successor : succs){
                        if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successor)] == nullptr) continue;
                        recompactNodes.insert(UnitigGraph2::unitigIndex_to_reverseDirection(successor));
                    }

                }


            }

            */



            
            //cout << isUnitigIndexBubble.size() << endl;
            
            /*
            //cout << "lala 1" << endl;
            for(const UnitigType& unitigIndex : isUnitigIndexBubble){
                
                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
                
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                //vector<UnitigType> succs = node->_successors;
                //vector<UnitigType> preds = node->_predecessors;


                vector<UnitigType> succs;
                _unitigGraph2->getSuccessors(unitigIndex, succs);

                vector<UnitigType> preds;
                _unitigGraph2->getPredecessors(unitigIndex, preds);

                
                _unitigGraph2->removeNode(_unitigGraph2->_unitigs[unitigName]);

                for(const UnitigType& predecessor : preds){
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessor)] == nullptr) continue;
                    recompactNodes.insert(predecessor);
                }

                for(const UnitigType& successor : succs){
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successor)] == nullptr) continue;
                    recompactNodes.insert(UnitigGraph2::unitigIndex_to_reverseDirection(successor));
                }

            }
            */
            //debugFile.close();
            //cout << _nbUnitigNameRemoved << endl;
            //getchar();

            //for(const UnitigType& unitigName : recompactUnitigNames){
                
            //    if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

            //    _unitigGraph2->recompact(_unitigGraph2->_unitigs[unitigName]);
            //}

            
            //cout << recompactNodes.size() << endl;

            vector<BubbleSide> recompactNodesVec;
            for(const UnitigType& unitigIndex : recompactNodes){

                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);

                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                recompactNodesVec.push_back({unitigIndex, _unitigGraph2->getStartNode(unitigIndex), _unitigGraph2->_unitigs[unitigName]->_nbMinimizers});
            }

            std::sort(recompactNodesVec.begin(), recompactNodesVec.end(), BubbleSideComparatorRev);


            for(const BubbleSide& bubbleSide : recompactNodesVec){

                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(bubbleSide._unitigIndex);

                //const UnitigType& unitigName = bubbleSide._unitigIndex;
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;


                //debugFile << "Superbubble compacting: " << bubbleSide._unitigIndex << endl;

                _unitigGraph2->recompact(bubbleSide._unitigIndex);
            }

            //debugFile.close();

            //cout << "done" << endl;
            //getchar();
            //cout << "Remove superbubble end: " << _unitigGraph2->computeChecksum_abundance() << endl;
            /*
            vector<UnitigGraph2::UnitigNode*> recompactNodesVec;
            for(UnitigGraph2::UnitigNode* node : recompactNodes){
                if(node->_unitigIndex == -1) continue;
                recompactNodesVec.push_back(node);
            }

            std::sort(recompactNodesVec.begin(), recompactNodesVec.end(), UnitigGraph::NodeComparatorByLengthRev);

            for(UnitigGraph2::UnitigNode* node : recompactNodesVec){
                if(node->_unitigIndex == -1) continue;
                _unitigGraph->recompact(node);
            }
            */

            return isModification;
        }


        void indexNotSuperbubble(UnitigGraph2::UnitigNode* nodeSource, const vector<UnitigType>& seenNodes){
            //for(u_int32_t unitigIndex : seenNodes){
            //    _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[unitigIndex].push_back(nodeSource->_unitigIndex);
            //}
            //_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.insert(nodeSource->_unitigIndex);
        }

        bool isSuperbubble(const UnitigType& unitigIndexSource, UnitigType& unitigIndexExit, unordered_set<UnitigType>& isUnitigBubble, vector<UnitigType>& seenNodes){

            
            //cout << "\tIs superbubble: " << unitigIndexSource << endl;

            bool isReversed;
            const UnitigType& unitigNameSource = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSource, isReversed);

            UnitigGraph2::UnitigNode* nodeSource = _unitigGraph2->_unitigs[unitigNameSource];

            seenNodes.clear();

            //vector<ProgressiveAbundanceFilter::SuperbubbleNode> nodes;
            unordered_set<UnitigType> isVisited;
            unordered_set<UnitigType> seen;
            unordered_map<UnitigType, u_int64_t> pathLength;
            vector<UnitigType> queue;

            queue.push_back(unitigIndexSource);
            pathLength[unitigIndexSource] = 0;
            seenNodes.push_back(unitigIndexSource);

            while(queue.size() > 0){

                UnitigType v = queue[queue.size()-1];

                //cout << "\t\t" << v << endl;

                vector<UnitigType> vSuccessors;
                _unitigGraph2->getSuccessors(v, vSuccessors);
                //u_int32_t v = vNode->_unitigIndex;

                if(_progressiveAbundanceFilter->_cutoffIndex == 0 && vSuccessors.size() > 5) return false;

                queue.pop_back();

                if(pathLength[v] > _maxLength){
                    return false;
                }

                isVisited.insert(v);
                if(seen.find(v) != seen.end()){
                    seen.erase(v);
                }


                if(vSuccessors.size() == 0){
                    return false; //abort tip
                }


                //cout << "\t\taa" << endl;

                for(const UnitigType& u : vSuccessors){

                    //cout << "\t\t\tL1" << u << endl;

                    //UnitigGraph2::UnitigNode* uNode;
                    seenNodes.push_back(u);
                    const UnitigType& uName = UnitigGraph2::unitigIndex_to_unitigName(u);

                    //if(isUnitigBubble.find(uName) != isUnitigBubble.end()){
                    //    return false;
                    //}


                    //u_int32_t u = uNode->_unitigIndex;
                    

                    if(u == unitigIndexSource){
                        return false; //cycle including s
                    }

                    if(isVisited.find(u) == isVisited.end()){
                        seen.insert(u);


                        if(v == unitigIndexSource){
                            pathLength[u] = _unitigGraph2->_unitigs[uName]->getLength(_unitigGraph2); // uNode->_length;
                        }
                        else{
                            pathLength[u] = pathLength[v] + ((_unitigGraph2->_unitigs[uName]->_nbMinimizers-_unitigGraph2->_kminmerSize+1) * _unitigGraph2->_kminmerLengthNonOverlap);
                        }
                        
                    }
                    else{

                        return false; //Cycle within superbubble
                    }

                }

                //cout << "\t\tbb" << endl;

                for(const UnitigType& u : vSuccessors){

                    vector<UnitigType> uPredecessors;
                    _unitigGraph2->getPredecessors(u, uPredecessors);
                    
                    bool allPredecessorsAreVisited = true;
                    
                    for(const UnitigType& p : uPredecessors){
                        if(isVisited.find(p) == isVisited.end()){
                            allPredecessorsAreVisited = false;
                            break;
                        }
                    }

                    if(allPredecessorsAreVisited){
                        queue.push_back(u);
                    }

                    if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ 
                    //if(seen.size() == 1){ 
                        
                        UnitigType t = (*begin(seen));
                        //UnitigType t = queue[0];
                        //u_int32_t t = tNode->_unitigIndex;
                        
                        vector<UnitigType> tSuccessors;
                        _unitigGraph2->getSuccessors(t, tSuccessors);

                        if(std::find(tSuccessors.begin(), tSuccessors.end(), unitigIndexSource) == tSuccessors.end()){
                            unitigIndexExit = t;
                            //return tNode;
                            return true;
                        }
                        else{
                            return false; // cycle including s
                        }

                    }

                }

            }

            return false;

        }
        
        void getSuperbubblePath(const UnitigType& unitigIndexSource, const UnitigType& unitigIndexExit, unordered_set<UnitigType>& keepNodes){

            
            keepNodes.clear();
            //unordered_set<UnitigType> keepNodes;
            //UnitigGraph2::UnitigNode* nodeCurrent = sourceNode;
            UnitigType unitigIndex = unitigIndexSource;

            while(true){


                vector<UnitigType> maxAbNodes;
                float maxAbundance = 0;

                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                for(const UnitigType& unitigIndexSuccessor : successors){

                    const float& abundance = _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)]->_abundance;
                    

                    if(abundance == maxAbundance){
                        maxAbNodes.push_back(unitigIndexSuccessor);
                    }
                    else if(abundance > maxAbundance){
                        maxAbundance = abundance;
                        maxAbNodes.clear();
                        maxAbNodes.push_back(unitigIndexSuccessor);
                    }
                }


                if(maxAbNodes.size() == 1){
                    keepNodes.insert(maxAbNodes[0]);
                    unitigIndex = maxAbNodes[0];
                }
                else{

                    
                    vector<BubbleSide> bubbleSides;
                    for(const UnitigType& unitigIndex : maxAbNodes){
                        
                        bubbleSides.push_back({unitigIndex, _unitigGraph2->getStartNode(unitigIndex), _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]->_nbMinimizers}); 
                    }


                    std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);
                    keepNodes.insert(bubbleSides[0]._unitigIndex);

                    unitigIndex = bubbleSides[0]._unitigIndex;

                }
                

                if(unitigIndex == unitigIndexExit) break;
            }
            

        }

        void collapseSuperbubble2(const UnitigType& unitigIndexSource, const UnitigType& unitigIndexExit, vector<UnitigType>& unitigIndexToRemove){



            //unordered_set<UnitigType> keepNodes;
            //getSuperbubblePath(unitigIndexSource, unitigIndexExit, keepNodes);

            vector<UnitigType> superbubbleNodes;
            collectSuperbubbleNodes(unitigIndexSource, unitigIndexExit, superbubbleNodes);

            unordered_set<UnitigType> keepNodes;
            BellmanFord(unitigIndexSource, unitigIndexExit, keepNodes, superbubbleNodes);

            //cout << "----" << endl;
            for(const UnitigType& unitigIndex : superbubbleNodes){

                if(keepNodes.find(unitigIndex) != keepNodes.end()){
                    //cout << "Keep: " << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << " " << _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]->_abundance << endl;
                    continue;
                }
                
                //cout << "Remove: " << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << " " << _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]->_abundance << endl;
                unitigIndexToRemove.push_back(unitigIndex);
            }

            //getchar();
            /*
            //unordered_set<UnitigType> isVisited;
            queue<UnitigType> queue;

            //isVisited.insert(unitigIndexSource);
            queue.push(unitigIndexSource);

            while (!queue.empty()){

                UnitigType unitigIndex = queue.front();
                //u_int32_t unitigIndex_current = currentNode->_unitigIndex;
                queue.pop();

                //if(currentNode->_successors.size() > _progressiveAbundanceFilter->_maxNbSuccessors){
                //    _progressiveAbundanceFilter->_maxNbSuccessors = currentNode->_successors.size();
                //    Logger::get().debug() << "Max successors: " << _progressiveAbundanceFilter->_maxNbSuccessors;
                //}
                
                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                for(const UnitigType& unitigIndexSuccessor : successors){

                    //if (isVisited.find(unitigIndexSuccessor) != isVisited.end()) continue;
                    


                    if(isUnitigIndexBubble.find(unitigIndexSuccessor) == isUnitigIndexBubble.end()){
                        continue;
                    }

                    queue.push(unitigIndexSuccessor);


                    if(keepNodes.find(unitigIndexSuccessor) == keepNodes.end()){
                        unitigIndexToRemove.push_back(unitigIndexSuccessor);
                    }

                    //isVisited.insert(unitigIndexSuccessor);
                    //nodes.push_back(unitigIndexSuccessor);

                    //if(unitigIndexSource == 211512){
                    //    cout << "lol: " << unitigIndexSuccessor << endl;
                    //}

                }
            }
            */

        }

        void collapseSuperbubble(const UnitigType& unitigIndexSource, const UnitigType& unitigIndexExit, unordered_set<UnitigType>& isUnitigBubble, unordered_set<UnitigType>& isUnitigIndexBubble){


            vector<UnitigType> superbubbleNodes;
            collectSuperbubbleNodes(unitigIndexSource, unitigIndexExit, superbubbleNodes);

            /*
            unordered_set<UnitigType> keepNodes;
            //UnitigGraph2::UnitigNode* nodeCurrent = sourceNode;
            UnitigType unitigIndex = unitigIndexSource;

            while(true){


                vector<UnitigType> maxAbNodes;
                float maxAbundance = 0;
                //UnitigType maxV = -1;

                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                for(const UnitigType& unitigIndexSuccessor : successors){
                    
                    const float& abundance = _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)]->_abundance;
                    
                    //if(unitigIndexSource == 211512){
                    //    cout << "rofl: " << unitigIndexSuccessor << " " << abundance << endl;
                    //}

                    if(abundance == maxAbundance){
                        maxAbNodes.push_back(unitigIndexSuccessor);
                    }
                    else if(abundance > maxAbundance){
                        maxAbundance = abundance;
                        maxAbNodes.clear();
                        maxAbNodes.push_back(unitigIndexSuccessor);
                    }
                }

                if(maxAbNodes.size() == 1){
                    keepNodes.insert(maxAbNodes[0]);
                    unitigIndex = maxAbNodes[0];
                }
                else{

                    
                    vector<BubbleSide> bubbleSides;
                    for(const UnitigType& unitigIndex : maxAbNodes){
                        
                        bubbleSides.push_back({unitigIndex, _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]->_nbMinimizers}); 
                    }


                    std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);
                    keepNodes.insert(bubbleSides[0]._unitigIndex);

                    unitigIndex = bubbleSides[0]._unitigIndex;

                }
                

                if(unitigIndex == unitigIndexExit) break;
            }
            */

            
            //#pragma omp critical(superbubble)
            //{
            for(const UnitigType& unitigIndex : superbubbleNodes){

                //if(keepNodes.find(unitigIndex) != keepNodes.end()) continue;
                
                const UnitigType& uName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);


                isUnitigBubble.insert(uName);
                isUnitigIndexBubble.insert(unitigIndex);
                isUnitigIndexBubble.insert(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));
                //isUnitigBubble.insert(_unitigGraph->unitigIndex_toReverseDirection(node));

                //if(unitigIndexSource == 211512){
                //    cout << "lal: " << unitigIndex << endl;
                //    cout << "lal: " << UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex) << endl;
                //}

            }
            //}
            
            /*
            #pragma omp critical(superbubble)
            {

                bool isValid = true;

                if(isUnitigBubble.find(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexExit)) != isUnitigBubble.end()){
                    isValid = false;
                }
                
                for(const UnitigType& unitigIndex : superbubbleNodes){
                    if(keepNodes.find(unitigIndex) != keepNodes.end()) continue;

                    const UnitigType& uName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);

                    if(isUnitigBubble.find(uName) != isUnitigBubble.end()){
                        isValid = false;
                        break;
                    }
                }

                if(isValid){

                    for(const UnitigType& unitigIndex : superbubbleNodes){

                        if(keepNodes.find(unitigIndex) != keepNodes.end()) continue;
                        
                        const UnitigType& uName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);


                        isUnitigBubble.insert(uName);
                        isUnitigIndexBubble.insert(unitigIndex);
                        isUnitigIndexBubble.insert(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));
                        //isUnitigBubble.insert(_unitigGraph->unitigIndex_toReverseDirection(node));

                        //if(unitigIndexSource == 211512){
                        //    cout << "lal: " << unitigIndex << endl;
                        //    cout << "lal: " << UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex) << endl;
                        //}

                    }
                }
            }
            */
        }

        void collectSuperbubbleNodes(const UnitigType& unitigIndexSource, const UnitigType& unitigIndexExit, vector<UnitigType>& nodes){

            nodes.clear();

            unordered_set<UnitigType> isVisited;
            queue<UnitigType> queue;

            isVisited.insert(unitigIndexSource);
            isVisited.insert(unitigIndexExit);
            queue.push(unitigIndexSource);

            while (!queue.empty()){

                UnitigType unitigIndex = queue.front();
                //u_int32_t unitigIndex_current = currentNode->_unitigIndex;
                queue.pop();

                //if(currentNode->_successors.size() > _progressiveAbundanceFilter->_maxNbSuccessors){
                //    _progressiveAbundanceFilter->_maxNbSuccessors = currentNode->_successors.size();
                //    Logger::get().debug() << "Max successors: " << _progressiveAbundanceFilter->_maxNbSuccessors;
                //}
                
                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                for(const UnitigType& unitigIndexSuccessor : successors){

                    if (isVisited.find(unitigIndexSuccessor) != isVisited.end()) continue;

                    queue.push(unitigIndexSuccessor);
                    isVisited.insert(unitigIndexSuccessor);
                    nodes.push_back(unitigIndexSuccessor);


                    //if(unitigIndexSource == 211512){
                    //    cout << "lol: " << unitigIndexSuccessor << endl;
                    //}

                }
            }

        }
        

        void BellmanFord(const UnitigType& unitigIndexSource, const UnitigType& unitigIndexExit, unordered_set<UnitigType>& keepNodes, vector<UnitigType> superbubbleNodes){

            //cout << "----" << endl;
            //cout << unitigIndexSource << " " << unitigIndexExit << endl;

            unordered_map<UnitigType, UnitigType> parent;
            //parent[unitigIndexSource] = -1;

            superbubbleNodes.push_back(unitigIndexSource);
            //superbubbleNodes.push_back(unitigIndexExit);

            unordered_map<UnitigType, double> dist;

            for(const UnitigType& unitigIndex : superbubbleNodes){
                dist[unitigIndex] = INT_MAX;
            }


            dist[unitigIndexSource] = 0;

            for(const UnitigType& uu : superbubbleNodes){
                for(const UnitigType& u : superbubbleNodes){

                    vector<UnitigType> successors;
                    _unitigGraph2->getSuccessors(u, successors);

                    for(const UnitigType& v : successors){

                        int weight = - _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(v)]->_abundance;

                        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]){
                            parent[v] = u;
                            dist[v] = dist[u] + weight;
                        }
                    }
                }
            }

            UnitigType currentUnitigIndex = unitigIndexExit;

            while(true){

                //cout << currentUnitigIndex << endl;
                //getchar();
                keepNodes.insert(currentUnitigIndex);
                currentUnitigIndex = parent[currentUnitigIndex];
                if(currentUnitigIndex == unitigIndexSource) break;
            }


        }

    };
    
    
    class TipRemover{

        public:

        struct TipData{
            UnitigType _unitigName;
            u_int32_t _startNode;
            float _abundance;
            u_int32_t _nbMinimizers;
        };

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph2* _unitigGraph2;
        u_int64_t _maxLength;

        u_int64_t _nbTipRemoved;

        /*
        struct TipComparator {
            bool operator()(UnitigGraph2::UnitigNode* p1, UnitigGraph2::UnitigNode* p2) const{
                if(p1->_nbMinimizers == p2->_nbMinimizers){
                    if(p1->_abundance == p2->_abundance){
                        return p1->startNode2() > p2->startNode2();
                    }
                    return p1->_abundance < p2->_abundance;
                }
                return p1->_nbMinimizers < p2->_nbMinimizers;
            }
        };
        */

        struct TipComparator2 {
            bool operator()(const TipData& p1, const TipData& p2) const{
                
                if(p1._nbMinimizers == p2._nbMinimizers){
                    if(p1._abundance == p2._abundance){
                        return p1._startNode > p2._startNode;
                    }
                    return p1._abundance < p2._abundance;
                }
                return p1._nbMinimizers < p2._nbMinimizers;
            }
        };


        /*
        struct TipComparatorAbundance {
            bool operator()(UnitigGraph2::UnitigNode* p1, UnitigGraph2::UnitigNode* p2) const{
                if(p1->_abundance == p2->_abundance){
                    if(p1->_nbMinimizers == p2->_nbMinimizers){
                        return p1->startNode2() > p2->startNode2();
                    }
                    return p1->_nbMinimizers < p2->_nbMinimizers;
                }
                return p1->_abundance < p2->_abundance;
            }
        };
        */

        //vector<UnitigGraph::Node*> _queue;
        //vector<UnitigGraph::Node*> _nextQueue;
        std::set<TipData, TipComparator2> _queue;
        //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator> _queue;
        u_int64_t _nextReadIndexWriter;

        TipRemover(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph2 = _progressiveAbundanceFilter->_unitigGraph2;
            _nbTipRemoved = 0;
            _maxLength = maxLength;
	        //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator>();
        }


        u_int64_t execute(){

            UnitigType unitigIndexDummy;

            for(const UnitigType& unitigName : _progressiveAbundanceFilter->_validNodes2){

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[unitigName];

                if(isTipAny(node, unitigIndexDummy)){
                    /*
                    bool isReversed = false;
                    UnitigGraph2::unitigIndex_to_unitigName(unitigIndexDummy, isReversed);
                    
                    u_int32_t startNode = 0;

                    if(isReversed){
                        startNode = node->_startNode_reverse;
                    }
                    else{
                        startNode = node->_startNode_forward;
                    }
                    */

                    _queue.insert({node->_unitigName, _unitigGraph2->getStartNode(unitigIndexDummy), node->_abundance, node->_nbMinimizers});
                }
            }  


            //cout << _queue.size() << endl;

            auto start = high_resolution_clock::now();
            bool isModification = removeTips(false, _maxLength);
            Logger::get().debug() << "Nb Tip removed: " << _nbTipRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
            //cout << "Nb Tip removed: " << _nbTipRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
            //getchar();

            return isModification;
            
            
        }

        bool removeTips(bool useK, u_int64_t maxLength){

            //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala_1.txt");

            UnitigType unitigIndexTip;
            bool isModification = false;
            bool lala = false;

            while(!_queue.empty()){

                auto top = _queue.begin();
                TipData tipData = *top;
                _queue.erase(top);

                if(_unitigGraph2->_unitigs[tipData._unitigName] == nullptr) continue;

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[tipData._unitigName];

                if(!isTipAny(node, unitigIndexTip)){
                    //_nextQueue.push_back(node);
                    continue;
                }


                //debugFile << "Tip: " << unitigIndexTip << " " << node->_nbMinimizers << " " << node->_abundance << " " <<  node->_startNode_forward << " " << node->_startNode_reverse << endl;
                //continue;

                isModification = true;
                _nbTipRemoved += 1;


                //cout << "Remove tip: " << node->_unitigName*2 << " " << node->startNode2() << " " << node->_abundance << " " << node->getLength(_unitigGraph2) << " " << _unitigGraph2->computeChecksum_abundance() << endl;
                //debugFile << "Tip: " << unitigIndexTip << " " << node->_abundance << " " << node->getLength(_unitigGraph2) << endl;
                //getchar();
                

                vector<UnitigType> preds;
                _unitigGraph2->getPredecessors(unitigIndexTip, preds);





                std::sort(preds.begin(), preds.end(), UnitigGraph2::NodeComparator);


                //const UnitigType& unitigName = unitigIndex_to_unitigName(unitigIndex);

                //vector<UnitigType> predecessors;
                //getPredecessors(unitigIndex, predecessors);
                
                //u_int32_t unitigIndexToRemove = unitigIndex; //unitigIndex_to_reverseDirection(unitigIndex);
                //vector<UnitigType> predecessors;
                //getPredecessors(unitigIndex, predecessors);

                for(const UnitigType& unitigIndexPredecessor : preds){

                    bool isReversedSuccessor;
                    const UnitigType& unitigNamePredecessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor, isReversedSuccessor);

                    if(_unitigGraph2->_unitigs[unitigNamePredecessor] == nullptr) continue;

                    //debugFile << "Tip removePred: " << unitigIndexTip << " " << unitigIndexPredecessor << endl;

                    //cout << "\t" << unitigIndexToString(unitigIndexSuccessor) << endl;


                    UnitigGraph2::UnitigNode* predecessor = _unitigGraph2->_unitigs[unitigNamePredecessor];

                    if(isReversedSuccessor){


                        vector<UnitigType>& successorsReverse = _unitigGraph2->_unitigs[unitigNamePredecessor]->_successors_reverse;

                        //for(UnitigType unitigIndex : successorsReverse){
                        //    cout << "\t\t" << unitigIndex << endl;
                        //}

                        //int nbSuccessors = successorsReverse.size();
                        //cout << successorsReverse.size() << endl;
                        successorsReverse.erase(std::remove(successorsReverse.begin(), successorsReverse.end(), unitigIndexTip), successorsReverse.end());
                        //cout << successorsReverse.size() << endl;
                        


                    }
                    else{


                        vector<UnitigType>& successorsFoward = _unitigGraph2->_unitigs[unitigNamePredecessor]->_successors_forward;

                        //for(UnitigType unitigName : successorsFoward){
                        //    cout << "\t\t" << unitigName << endl;
                        //}

                        //int nbSuccessors = successorsFoward.size();
                        //cout << successorsFoward.size() << endl;
                        successorsFoward.erase(std::remove(successorsFoward.begin(), successorsFoward.end(), unitigIndexTip), successorsFoward.end());
                        //cout << successorsFoward.size() << endl;


                    }


                    //debugFile << "Tip recompacting: " << unitigIndexPredecessor << endl;
                    _unitigGraph2->recompact(unitigIndexPredecessor);

                    UnitigType unitigIndexTip2;
                    if(isTipAny(predecessor, unitigIndexTip2)){
                        //debugFile << "Tip add: " << unitigIndexTip2 << " " << predecessor->_abundance << " " << predecessor->getLength(_unitigGraph2) << endl;
                        _queue.insert({predecessor->_unitigName, _unitigGraph2->getStartNode(unitigIndexTip2), predecessor->_abundance, predecessor->_nbMinimizers});
                    }
                    

                }

                node->_successors_forward.clear();
                node->_successors_reverse.clear();

                //getchar();

                /*
                _unitigGraph2->removePredecessors(unitigIndexTip, debugFile);

                
                for(const UnitigType& unitigIndexPredecessor : preds){
                    
                    
                    UnitigGraph2::UnitigNode* predecessor = _unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor)];

                    if(predecessor == nullptr) continue;
                
                    //UnitigGraph2::UnitigNode* predecessor_rc = _unitigGraph->unitigIndex_toReverseDirection(predecessor);
                 
                    
                

                    //_unitigGraph->removePredecessor(node, predecessor);
                    //recompactNodes.insert(predecessor);
                    //debugFile << "Tip recompacting: " << unitigIndexPredecessor << endl;
                    //_unitigGraph2->recompact(unitigIndexPredecessor);


                    if(isTipAny(predecessor, unitigIndexTip)){
                        debugFile << "Tip add1: " << unitigIndexTip << endl;
                        _queue.insert({predecessor->_unitigName, predecessor->_abundance, predecessor->_nbMinimizers});
                    }

                    //if(isTipAny(predecessor_rc, unitigIndexTip)){
                    //    _queue.insert(predecessor_rc);
                    //}
                    


                }
                */
                
            }

            //cout << "done" << endl;
            //debugFile.close();
            //getchar();
            
            return isModification;
        }


        bool isTipAny(const UnitigGraph2::UnitigNode* node, UnitigType& unitigIndexTip){

            if(node == nullptr) return false;
            if(node->getLength(_unitigGraph2) > _maxLength) return false;
            //if(node->_abundance > 1) return false;
            
            const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
            const UnitigType& unitigIndex2 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

            if(!_unitigGraph2->hasSuccessors(unitigIndex1) && !_unitigGraph2->hasNoPredecessors(unitigIndex1)){
                unitigIndexTip = unitigIndex1;
                return true;
            }

            if(!_unitigGraph2->hasSuccessors(unitigIndex2) && !_unitigGraph2->hasNoPredecessors(unitigIndex2)){
                unitigIndexTip = unitigIndex2;
                return true;
            }

            //if(node->_successors_forward.size() > 0 && node->_successors_reverse.size() == 0) return true;
            //if(node->_successors_forward.size() == 0 && node->_successors_reverse.size() > 0) return true;

            //vector<UnitigType> successors;
            //_unitigGraph2->getSuccessors(UnitigGraph2::unitigName_to_unitigIndex(2025, false), successors);

            //vector<UnitigType> predecessors;
            //_unitigGraph2->getPredecessors(UnitigGraph2::unitigName_to_unitigIndex(2025, false), predecessors);

            //if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() == 0) return false;

       

            return false;
        }

    };

    UnitigGraph2* _unitigGraph2;
    vector<SaveState> _cachedGraphStates;
    SaveState _currentSaveState;
    size_t _kminmerSize;


    struct CutoffIndex{
        u_int32_t _cutoffIndex;
        float _cutoffValue;
    };


    
    //unordered_map<u_int32_t, vector<SuperbubbleNode>> _superBubbleIndex;
    u_int32_t _cutoffIndex;
    vector<CutoffIndex> _cutoffIndexes;
    string _tmpDir;
    //ofstream _debugOuputFile;
    bool _removeBubble;
    bool _isFirstPass;
    int _nbCores;
    bool _genGraph;

    ProgressiveAbundanceFilter(UnitigGraph2* unitigGraph2, const string& tmpDir, size_t kminmerSize, bool removeBubble, bool isFirstPass, bool genGraph, int nbCores){
        _unitigGraph2 = unitigGraph2;
        _tmpDir = tmpDir;
        _kminmerSize = kminmerSize;
        _removeBubble = removeBubble;
        _isFirstPass = isFirstPass;
        _genGraph = genGraph;
        _nbCores = nbCores;
    }

    //vector<UnitigGraph::Node*> _validNodes;
    vector<UnitigType> _validNodes2;

	template<typename Functor>
	//void parse(const Functor& functor){
    void execute(Functor& functor){

        //functor(0);
        //return;

        //cout << "Save graph init" << endl;
        //_unitigGraph2->save(_tmpDir + "/assembly_graph_init.gfa");
        //_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp//assembly_graph_init_1.gfa");
        //exit(1);
        //getchar();
        /*
        //exit(1);
        cout << _unitigGraph2->_unitigs[2025]->_successors_forward.size() << " " << _unitigGraph2->_unitigs[2025]->_successors_reverse.size() << endl;
        

        vector<UnitigType> successors;
        _unitigGraph2->getSuccessors(UnitigGraph2::unitigName_to_unitigIndex(2025, false), successors);

        vector<UnitigType> predecessors;
        _unitigGraph2->getPredecessors(UnitigGraph2::unitigName_to_unitigIndex(2025, false), predecessors);

        cout << successors.size() << " " << predecessors.size() << endl;

        for(const UnitigType& unitigIndex : successors){
            cout << "Succ: " << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
        }
        for(const UnitigType& unitigIndex : predecessors){
            cout << "Pred: " << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
        }

        exit(1);
        */
        /*
        UnitigType unitigNameToRemove = 594;

        unordered_set<UnitigType> recompactUnitigNames;
        for(const UnitigType& unitigIndex : _unitigGraph2->_unitigs[unitigNameToRemove]->_successors_forward){
            recompactUnitigNames.insert(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex));
        }
        for(const UnitigType& unitigIndex : _unitigGraph2->_unitigs[unitigNameToRemove]->_successors_reverse){
            recompactUnitigNames.insert(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex));
        }

        _unitigGraph2->removeNode(_unitigGraph2->_unitigs[unitigNameToRemove]);

        for(const UnitigType& unitigName : recompactUnitigNames){
            _unitigGraph2->recompact(_unitigGraph2->_unitigs[unitigName]);
        }

        _unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa");
        exit(1);
        */

        //for(UnitigGraph2::UnitigNode* node : _unitigGraph2->_unitigs){
        //    _unitigGraph2->removeNode(node);
        //    _unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa");
        //    getchar();
        //}
        //_unitigGraph2->removeNode(_unitigGraph2->_unitigs[3782]);
        //exit(1);
        //exit(1);
        //_unitigGraph->save(_tmpDir + "/assembly_graph_init.gfa", 0);

        _checksumAbundanceTotal = 0;
        _checksumNodeTotal = 0;
        _cutoffIndex = 0;
        _maxAbundance = 1;
        //_debugOuputFile = ofstream("/home/gats/workspace/tmp/lala_" + to_string(_kminmerSize) + ".txt");


        //_unitigGraph->save(_tmpDir + "/assembly_graph.gfa", 0);
        //_cutoffIndex = 0;
        //checkSaveState(0);
        //return;
        
        u_int32_t maxAbundance = 0;
        for(UnitigGraph2::UnitigNode* node : _unitigGraph2->_unitigs){
            if(node == nullptr) continue;
            if(node->_abundance > maxAbundance){
                maxAbundance = node->_abundance;
            }
        }
        

        Logger::get().debug() << "Max abundance: " << maxAbundance;
        Logger::get().debug() << "Nb unitigs: " << _unitigGraph2->getNbUnitigs();
        //cout << "Checksum: " << _unitigGraph2->computeChecksum_successors() << endl;
        //cout << "Checksum: " << _unitigGraph2->computeChecksum_abundance() << endl;

        simplifyProgressive(functor);

        Logger::get().debug() << "Checksum: " << _checksumNodeTotal;
        Logger::get().debug() << "Checksum abundnace: " << _checksumAbundanceTotal;
        //getchar();

        //_debugOuputFile.close();

        //if(_kminmerSize > 4){
        //string command = "diff /home/gats/workspace/tmp/bug/lala_" + to_string(_kminmerSize) + ".txt" + " /home/gats/workspace/tmp/lala_" + to_string(_kminmerSize-1) + ".txt";
        //system(command.c_str());
        //getchar();
        //}
    }

    unordered_set<float> isCutoffProcessed;
    u_int64_t _maxNbSuccessors;

	template<typename Functor>
    void simplifyProgressive(Functor& functor){

        float maxAbundance = 2000;
        float currentCutoff = 0;
        _currentSaveState = {0, {}, {}};

        while(true){

            _maxNbSuccessors = 0;
            bool isModification = false;

            bool isMod = simplify();
            if(isMod){
                isModification = true;
            }


            //cout << "Max superbubble successors: " << currentCutoff << " " << _maxNbSuccessors << endl;
            checkSaveState(currentCutoff);
            
            //cout << "Remove abundance" << endl;
            u_int64_t nbErrorRemoved = removeAbundanceNoQueue(maxAbundance, currentCutoff, functor);
            if(nbErrorRemoved > 0){
                isModification = true;
            }
            


            //cout << _unitigGraph->computeChecksum_nodeIndex() << endl;
            //cout << _unitigGraph->computeChecksum_abundance() << endl;
            //getchar();
            //compact(true, unitigDatas);
            //resizeUnitigs();    




            if(!isModification) break;
        }
    }

    
    void checkSaveState(float currentCutoff){
        
        bool saveStateExist = false;
        for(const SaveState& saveState : _cachedGraphStates){
            if(currentCutoff == saveState._abundanceCutoff_min){
                saveStateExist = true;
                break;
            }
        }


        if(!saveStateExist){

            //_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph_" + to_string(_cutoffIndex) + "_2.gfa");
            //getchar();
            //cout << "saving state: " << currentCutoff << endl;

            //cout << "attention assemblygraph.gfa save" << endl;
            if(currentCutoff == 0){
                
                /*
                int nbUnitigs = 0;
                for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                    if(node->_unitigIndex == -1) continue;
                    nbUnitigs += 1;
                }

                cout << nbUnitigs << endl;
                */
                if(_genGraph){
                    _unitigGraph2->save(_tmpDir + "/assembly_graph.gfa");
                }
                //getchar();
                //exit(1);

                
                //exit(1);
            }

            //if(currentCutoff > 2){
            //    _unitigGraph->save(_tmpDir + "/assemblyGraph_" + to_string(currentCutoff) + ".gfa", 0);
            //}
            //compact(true, unitigDatas);
            //if(doesSaveUnitigGraph && currentCutoff == 0) saveUnitigGraph(_outputDir + "/minimizer_graph_u_cleaned.gfa", mdbg, minimizerSize, nbCores, true);
            

            //for(u_int32_t nodeName : currentSaveState._nodeNameRemoved_tmp){
            //    currentSaveState._nodeNameRemoved.push_back(nodeName);
            //}
            //currentSaveState._nodeNameRemoved_tmp.clear();

            int nbUnitigs = 0;
            for(size_t i=0; i<_unitigGraph2->_unitigs.size(); i++){

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[i];
                if(node == nullptr) continue;
                nbUnitigs += 1;
            }

            dumpUnitigs(currentCutoff);

            _currentSaveState._abundanceCutoff_min = currentCutoff;
            _cachedGraphStates.push_back(_currentSaveState);
            _currentSaveState = {0, {}, {}};


        
        }

    }
    
    float _maxAbundance;

    bool simplify(){

        
        bool isModification = false;
        static u_int64_t maxTipLength = 50000;
        static u_int64_t maxBubbleLength = 50000;

        u_int64_t maxLengthKminmer = _unitigGraph2->_kminmerLength*2.25;

        maxTipLength = max(maxTipLength, maxLengthKminmer);
        maxBubbleLength = max(maxBubbleLength, maxLengthKminmer);


        while(true){

            bool isModificationSub = false;
            
            _validNodes2.clear();
            for(UnitigGraph2::UnitigNode* node : _unitigGraph2->_unitigs){

                if(node == nullptr) continue;

                _validNodes2.push_back(node->_unitigName);
            }

            /*
            while(true){
                compact(true, unitigDatas);

                u_int64_t nbSelfLoopRemoved = removeSelfLoops(currentSaveState);

                #ifdef PRINT_DEBUG_SIMPLIFICATION
                    _logFile << "Nb self loop removed: " << nbSelfLoopRemoved << endl;
                #endif
                if(nbSelfLoopRemoved == 0) break;
                isModification = true;
                isModSub = true;

                
            }
            */
            
            
            if(_removeBubble){
                while(true){
                    bool isModBubble = false;

                    
                    //cout << "Superbubble" << endl;
                    
                    Logger::get().debug() << "Superbubble";
                    SuperbubbleRemoverOld superbubbleRemover(this, maxBubbleLength);
                    bool isModificationSuperbubble = superbubbleRemover.execute();
                    if(isModificationSuperbubble){
                        isModification = true;
                        isModBubble = true;
                        isModificationSub = true;
                    }
                    
                    break;
                    //_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph_init.gfa");
                    //cout << "check" << endl;
                    //getchar();
                    
                    //cout << "Nb superbubble removed: " << superbubbleRemover._nbRemovedTotal << endl;

                    //cout << _unitigGraph->computeChecksum_nodeIndex() << endl;
                    //cout << _unitigGraph->computeChecksum_abundance() << endl;
                    //cout << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s" << endl;
                    //getchar();
                    
                    
                    
                    /*
                    Logger::get().debug() << "Bubble";
                    auto start = high_resolution_clock::now();
                    
                    BubbleRemover bubbleRemover(this, maxBubbleLength);
                    bool isModificationBubble = bubbleRemover.execute();
                    if(isModificationBubble){
                        isModification = true;
                        isModBubble = true;
                        isModificationSub = true;
                    }
                    */
                    //getchar();
                    
                    //_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph_init2.gfa");
                    //cout << "check" << endl;
                    //getchar();
                    
                    //break;
                    //cout << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s" << endl;
                    

                    if(!isModBubble) break;
                }
            }
            


            //cout << "check" << endl;
            //getchar();
            //_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph_init2.gfa");

            Logger::get().debug() << "Tip";
            //auto start = high_resolution_clock::now();
            
            
            TipRemover tipRemover(this, maxTipLength);
            u_int64_t nbRemoved = tipRemover.execute();
            if(nbRemoved){
                isModification = true;
                isModificationSub = true;
            }
            
            //getchar();

            //_logFile << "Nb unitigs: " << _unitigGraph->nbUnitigs() << endl;
            

            //cout << "Nb removed tips: " << nbRemoved << endl;
            //cout << _unitigGraph->computeChecksum_nodeIndex() << endl;
            //cout << _unitigGraph->computeChecksum_abundance() << endl;
            //cout << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s" << endl;
            
            
            //getchar();
            
            
            //if(_kminmerSize == 21){
                //_unitigGraph->save("/home/gats/workspace/run/minimizer_graph_u_tmp.gfa");
                //getchar();
            //}

            if(!isModificationSub) break;
        }

        //getchar();

        return isModification;
    }  

    struct AbundanceComparator {
        bool operator()(UnitigGraph2::UnitigNode* p1, UnitigGraph2::UnitigNode* p2) const{
            if(p1->_nbMinimizers == p2->_nbMinimizers){
                return p1->startNode2() > p2->startNode2();
            }
            return p1->_nbMinimizers < p2->_nbMinimizers;
        }
    };


	template<typename Functor>
    u_int64_t removeAbundanceNoQueue(float abundanceCutoff_min, float& currentCutoff, Functor& functor){


        //ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp//lala.txt");

        u_int64_t nbErrorsRemoved = 0;
        float localCutoffRate = 0.5;
        float cutoffGlobal = 1;
        float aplha = 0.1;
        float t=1.1;
        currentCutoff = min(t, abundanceCutoff_min);



        while(t < abundanceCutoff_min){ 
            

            auto start = high_resolution_clock::now();

            currentCutoff = t;
            _maxAbundance = currentCutoff*2;


            
            _queue = std::set<UnitigGraph2::UnitigNode*, AbundanceComparator>();

            unordered_set<UnitigType> recompactNodes;


            for(const UnitigType& unitigName : _validNodes2){

                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[unitigName];
                //if(node == nullptr) continue;

                if(node->_abundance >= t) continue;

                const UnitigType& unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
                /*
                for(const UnitigType& unitigIndex : node->_successors_forward){
                    recompactNodes.insert(unitigIndex);
                }
                for(const UnitigType& unitigIndex : node->_successors_reverse){
                    recompactNodes.insert(unitigIndex);
                }
               
                //debugFile << "Abundance: " << node->_unitigName*2 << endl;
                //cout << "Remove node: " << node->_unitigName*2 << " " << node->startNode2() << " " << _unitigGraph2->computeChecksum_abundance() << endl;
                _unitigGraph2->removeNode(node);
                */

                vector<UnitigType> predecessors;
                _unitigGraph2->getPredecessors(unitigIndex, predecessors);

                vector<UnitigType> successors;
                _unitigGraph2->getSuccessors(unitigIndex, successors);

                _unitigGraph2->removeNode(node);

                //debugFile << "Abundance: " << unitigIndex << endl;


                for(const UnitigType& predecessor : predecessors){
                    
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessor)] == nullptr) continue;
                    //if(predecessor->_abundance < t) continue;

                    recompactNodes.insert(predecessor);
                    //recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(predecessor));

                }


                for(const UnitigType& successor : successors){
                    
                    if(_unitigGraph2->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successor)] == nullptr) continue;
                    //if(predecessor->_abundance < t) continue;

                    //recompactNodes.insert(successor);
                    recompactNodes.insert(UnitigGraph2::unitigIndex_to_reverseDirection(successor));

                }

                
                nbErrorsRemoved += 1;

            }

            /*
            for(const UnitigType& unitigName : recompactUnitigNames){
                
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                debugFile << "Abundance recompacting: " << unitigName*2 << endl;

                _unitigGraph2->recompact(_unitigGraph2->_unitigs[unitigName]);
            }
            */
            vector<BubbleSide> recompactNodesVec;
            for(const UnitigType& unitigIndex : recompactNodes){

                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);

                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;

                recompactNodesVec.push_back({unitigIndex, _unitigGraph2->getStartNode(unitigIndex), _unitigGraph2->_unitigs[unitigName]->_nbMinimizers});
            }

            std::sort(recompactNodesVec.begin(), recompactNodesVec.end(), BubbleSideComparatorRev);


            for(const BubbleSide& bubbleSide : recompactNodesVec){

                const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(bubbleSide._unitigIndex);

                //const UnitigType& unitigName = bubbleSide._unitigIndex;
                if(_unitigGraph2->_unitigs[unitigName] == nullptr) continue;


                //debugFile << "Abundance compacting: " << bubbleSide._unitigIndex << endl;

                _unitigGraph2->recompact(bubbleSide._unitigIndex);
            }


            

            

            t = t * (1+aplha);

            if(nbErrorsRemoved > 0){
                Logger::get().debug() << "Nb errors removed: " << nbErrorsRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
                

            }

            if(isCutoffProcessed.find(t) == isCutoffProcessed.end()){
                isCutoffProcessed.insert(t);
                functor(t);
            }

            if(nbErrorsRemoved > 0) break;
            
        }

        
        //debugFile.close();
        //getchar();

        return nbErrorsRemoved;

    }

    std::set<UnitigGraph2::UnitigNode*, AbundanceComparator> _queue;
    /*
    //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , AbundanceComparator> _queue;
    std::set<UnitigGraph::Node*, AbundanceComparator> _queue;

    struct TipComparatorLala {
        bool operator()(UnitigGraph::Node* p1, UnitigGraph::Node* p2) const{
            return p1->startNode2() > p2->startNode2();
        }
    };

    u_int64_t removeAbundance(float abundanceCutoff_min, float& currentCutoff){




        u_int64_t nbErrorsRemoved = 0;
        float localCutoffRate = 0.5;
        float cutoffGlobal = 1;
        float aplha = 0.1;
        float t=1.1;
        currentCutoff = min(t, abundanceCutoff_min);
        unordered_set<UnitigGraph::Node*> recompactNodes;



        while(t < abundanceCutoff_min){ 
            
            //if(saveAllState) checkSaveState(currentCutoff, unitigDatas, detectRoundabout, maxBubbleLength, saveState, insertBubble, doesSaveUnitigGraph, mdbg, minimizerSize, nbCores);

            currentCutoff = t;
            //unordered_set<u_int32_t> removedUnitigs;

            //cout << t << " " << abundanceCutoff_min << endl;



            //compact(true, unitigDatas);

            
            _queue = std::set<UnitigGraph::Node*, AbundanceComparator>();
            //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , AbundanceComparator>();
            //for(UnitigGraph::Node* node : _unitigGraph->_nodes){
            //    if(node->_unitigIndex == -1) continue;
            //    if(node->)
            //    _queue.push(node);
            //}
           
            for(size_t i=0; i<_validNodes.size(); i+=2){

                UnitigGraph::Node* node1 = _validNodes[i];
                UnitigGraph::Node* node2 = _validNodes[i+1];
                UnitigGraph::Node* node = nullptr;

                if(node1->_unitigIndex == -1) continue;
                if(node2->_unitigIndex == -1) continue;

                if(node1->startNode2() < node2->startNode2()){
                    node = node1;
                }
                else{
                    node = node2;
                }

                if(node->_abundance >= t) continue;
                _queue.insert(node);
            }

            while(!_queue.empty()){

                auto top = _queue.begin();
                UnitigGraph::Node* node = *top;
                _queue.erase(top);
                //UnitigGraph::Node* node = _queue.top();
                //_queue.pop();

                //if(node == nullptr) continue;
                //if(node->_unitigIndex == -1) continue;

                //for(size_t i=0; i<_unitigGraph->_nodes.size(); i+=2){
                //const Unitig& unitig = _unitigs[i];
                //UnitigGraph::Node* node = _unitigGraph->_nodes[i];

                if(node == nullptr) continue;
                if(node->_unitigIndex == -1) continue;

                //double cutoff = t; //min(t, localabundance*0.5);

                //if(node->_abundance < cutoff){


                    vector<UnitigGraph::Node*> preds = node->_predecessors;
                    _unitigGraph->removeNode(node);


                    for(UnitigGraph::Node* predecessor : preds){
                        
                        if(predecessor->_unitigIndex == -1) continue;
                        //cout << predecessor->_unitigIndex << endl;
                        //_unitigGraph->removePredecessor(node, predecessor);
                        _unitigGraph->recompact(predecessor);

                        UnitigGraph::Node* predecessor_rc = _unitigGraph->unitigIndex_toReverseDirection(predecessor);
                        //if(predecessor->_unitigIndex == -1){
                        //    cout << "argh" << endl;
                        //    getchar();
                        //}
                        //if(predecessor){
                            //cout << predecessor->_unitigIndex << endl;
                            //cout << "push; " << predecessor->_unitigIndex << endl;

                            if(predecessor->startNode2() < predecessor_rc->startNode2()){
                                //predecessor = node1;
                            }
                            else{
                                predecessor = predecessor_rc;
                            }

                            if(predecessor->_abundance >= t){
                                _queue.erase(predecessor);
                            }
                            else{
                                _queue.erase(predecessor);
                                _queue.insert(predecessor);
                            }

                        //}

                        //UnitigGraph::Node* predecessor_rc = _unitigGraph->unitigIndex_toReverseDirection(predecessor);
                        //if(predecessor_rc){
                            //cout << predecessor_rc->_unitigIndex << endl;
                            //cout << "push; " << predecessor_rc->_unitigIndex << endl;
                            //_queue.erase(predecessor_rc);
                            //_queue.insert(predecessor_rc);

                            //if(predecessor_rc->_abundance >= t){
                            //    _queue.erase(predecessor_rc);
                            //}
                        //}
                    }
                    
                    //_unitigGraph->recompact(node);

                    //_queue.push(node);
                    
                    nbErrorsRemoved += 1;

                    
                //}


            }
            
            

            t = t * (1+aplha);

            if(nbErrorsRemoved > 0) break;
        }

        for(UnitigGraph::Node* predecessor : recompactNodes){
            
            if(predecessor->_unitigIndex == -1) continue;
            _unitigGraph->recompact(predecessor);
        }

        

        return nbErrorsRemoved;

    }
    */

    u_int64_t _checksumAbundanceTotal;
    u_int64_t _checksumNodeTotal;

    void dumpUnitigs(float currentCutoff){

        //cout << "dump unitigs: a remettre" << endl;
        //return;

        u_int64_t skipped = 0;

        //cout << "dumpUnitigs" << _cutoffIndex << " truc a remettre: isCirulcar, isRepeatSide" << endl;
        int lala = 0;

        u_int64_t checksumTotal = 0;
        u_int64_t checksumAbundanceTotal = 0;
        _cutoffIndexes.push_back({_cutoffIndex, currentCutoff});

        //cout << "dump unitigs: " << currentCutoff << endl;
        ofstream outputFile(_tmpDir + "/filter/unitigs_" + to_string(_cutoffIndex) + ".bin");

        for(size_t i=0; i<_unitigGraph2->_unitigs.size(); i++){

            UnitigGraph2::UnitigNode* node = _unitigGraph2->_unitigs[i];
            if(node == nullptr) continue;

        //for(UnitigGraph::Node* node : _unitigGraph->_nodes){
        //for(size_t i=0; i<_unitigGraph->_nodes.size(); i+=2){

            //UnitigGraph::Node* node1 = _unitigGraph->_nodes[i];
            //UnitigGraph::Node* node2 = _unitigGraph->_nodes[i+1];
            //UnitigGraph::Node* node = nullptr;

            //if(node1->startNode2() < node2->startNode2()){
            //    node = node1;
            //}
            //else{
            //    node = node2;
            //}

            //if(node == nullptr) continue;
            //if(node->_unitigIndex == -1) continue;
            //if(node->_unitigIndex % 2 == 1) continue;

            //const vector<UnitigType>& successors;
            //_unitigGraph2->getSuccessors()

            if(node->_successors_forward.size() == 0 && node->_successors_reverse.size() == 0 && node->_abundance == 1){
            //if(node->_successors.size() == 0 && node->_predecessors.size() == 0 && node->_abundance == 1) {
                //skipped += node->_unitigIndex;
                //cout << "omg: " <<  node->_unitigIndex << endl;
                continue;
            }

            vector<UnitigType> unitigs;
            if(node->_unitigMerge.size() == 0){
                unitigs.push_back(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));    
            }
            else{
                unitigs = node->_unitigMerge;
            }

            u_int32_t size = unitigs.size();
            outputFile.write((const char*)&size, sizeof(size));
            

            u_int8_t isCircular = _unitigGraph2->isCircular(node);



            //if(isCircular) continue;
            //isCircular = false;
            outputFile.write((const char*)&isCircular, sizeof(isCircular));
            
            bool isRepeatSide_ = _unitigGraph2->isRepeatSide(node);
            outputFile.write((const char*)&isRepeatSide_, sizeof(isRepeatSide_));

            //if(isRepeatSide_){
                //cout << "derp: " << node->_unitigName*2 << " " << node->getLength(_unitigGraph2) << " " << node->_abundance << endl;
                //getchar();
            //}

            //if(node->_nodes.size() > 1000){
            //    cout << isCircular << endl;
            //    getchar();
            //}
            //cout << "Write unitig: " << i << " " << node->_unitigIndex << endl;
            float abundance = node->_abundance;
            outputFile.write((const char*)&abundance, sizeof(abundance));

            u_int32_t nbMinimizers = node->_nbMinimizers;
            outputFile.write((const char*)&nbMinimizers, sizeof(nbMinimizers));

            outputFile.write((const char*)&unitigs[0], size * sizeof(UnitigType));
            
            //cout << "dump: " << unitigs.size() << endl;
            /*
            u_int64_t checksum = 0;
            u_int64_t checksumAbundance = 0;
            for(u_int32_t nodeIndex : node->_nodes){
                checksum += nodeIndex;
                checksumAbundance += node->_abundance;
            }
            checksum *= node->_nodes.size();
            checksumAbundance *= node->_nodes.size();
            checksumTotal += checksum;
            _checksumAbundanceTotal += checksumAbundance;
            _checksumNodeTotal += checksum;
            */
            //cout << "Dump unitig: " << node->_unitigIndex << " " << node->_abundance << endl;
            lala += 1;
        }

        //cout << "yo? " << lala << endl;
        //cout << "skipped: " << skipped << endl;

        outputFile.close();

        _cutoffIndex += 1;


        //getchar();
    }



    
};



#endif
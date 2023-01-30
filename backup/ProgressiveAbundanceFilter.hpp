


#ifndef MDBG_METAG_PROGRESSIVEABUNDANCEFILTER
#define MDBG_METAG_PROGRESSIVEABUNDANCEFILTER


//#include "Commons.hpp"
#include "Graph.hpp"

class X
{
    public: 
    void operator()(float cutoff) {}
};

class ProgressiveAbundanceFilter{

public:




    struct SaveState{
        float _abundanceCutoff_min;
        vector<u_int32_t> _nodeNameRemoved;
        vector<DbgEdge> _isEdgeRemoved;
    };

    struct BubbleSide{
        UnitigGraph::Node* _unitig;
        size_t _nbNodes;
        u_int32_t _nodeIndexSum;
        //u_int32_t _quality;
    };

    static bool BubbleSideComparator(const BubbleSide &a, const BubbleSide &b){

        if(a._nbNodes == b._nbNodes){
            return a._nodeIndexSum > b._nodeIndexSum;
        }
        return a._nbNodes > b._nbNodes;

    }

    static bool SuperbubbleComparator(UnitigGraph::Node* p1, UnitigGraph::Node* p2){
        if(p1->_nodes.size() == p2->_nodes.size()){
            return p1->startNode() > p2->startNode();
        }
        return p1->_nodes.size() > p2->_nodes.size();
    }

    /*
    class SuperbubbleRemover{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph* _unitigGraph;
        u_int64_t _maxLength;
        u_int64_t _nbRemovedTotal;

        std::set<UnitigGraph::Node*, SuperbubbleComparator> _queue;
        //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , SuperbubbleComparator> _queue;
        u_int64_t _nextReadIndexWriter;
        u_int64_t _nbRemoved;

        SuperbubbleRemover(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph = _progressiveAbundanceFilter->_unitigGraph;
            _maxLength = maxLength;
            _nbRemovedTotal = 0;
        }

        bool execute(){

            return removeSuperbubble();
          
        }

        bool removeSuperbubble(){

            
            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){

                if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1) continue;

                if(_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()){
                    continue;
                }
        

                //if(node->_length > _maxLength) continue;
                //if(node->_unitigIndex % 2 == 1) continue;

                //if(isTip(node)){
                _queue.insert(node);
                //}
            }
            

            bool isModification = false;

            bool lala = false;
            while(!_queue.empty()){

                auto top = _queue.begin();
                UnitigGraph::Node* node = *top;
                _queue.erase(top);

                if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1) continue;

                

                //cout << node->_nodes.size() << endl;

                //cout << "is bubble? " << node->_unitigIndex << " " << (_lol.find(node->_unitigIndex) != _lol.end()) << " " << (_lol.find(_unitigGraph->unitigIndex_toReverseDirection(node->_unitigIndex)) != _lol.end()) << endl;
                vector<u_int32_t> seenNodes;
                UnitigGraph::Node* nodeExit = isSuperbubble(node, seenNodes);

                if(nodeExit == nullptr){
                    indexNotSuperbubble(node, seenNodes);
                    continue;
                }
                
                if(nodeExit == _unitigGraph->unitigIndex_toReverseDirection(node)){
                    indexNotSuperbubble(node, seenNodes);
                    continue; //loop side of an inverse repeat
                }
                if(nodeExit == node){
                    if(node->_length < 400000){
                        indexNotSuperbubble(node, seenNodes);
                        continue;
                    }
                }

                //cout << "\tRemove superbubble: " << node->_unitigIndex << " " << nodeExit->_unitigIndex << endl;

                isModification = true;

                collapseSuperbubble(node, nodeExit);
                _nbRemoved += 1;
                _nbRemovedTotal += 1;
                
                _unitigGraph->recompact(node);
                //_unitigGraph->removeNode(bubbleSide);
                //_unitigGraph->recompact(node);

                //_queue.erase(node);
                //_queue.erase(_unitigGraph->unitigIndex_toReverseDirection(node));
                _queue.insert(node);
                //_queue.insert(_unitigGraph->unitigIndex_toReverseDirection(node));
                //_queue.insert(nodeExit);
                //_queue.insert(_unitigGraph->unitigIndex_toReverseDirection(nodeExit));
                //_queue.push(node);
                //_queue.push(_unitigGraph->unitigIndex_toReverseDirection(node));

            }

            return isModification;
        }


        void indexNotSuperbubble(UnitigGraph::Node* nodeSource, const vector<u_int32_t>& seenNodes){
            for(u_int32_t unitigIndex : seenNodes){
                _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[unitigIndex].push_back(nodeSource->_unitigIndex);
            }
            _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.insert(nodeSource->_unitigIndex);
        }

        UnitigGraph::Node* isSuperbubble(UnitigGraph::Node* nodeSource, vector<u_int32_t>& seenNodes){

            seenNodes.clear();

            //vector<ProgressiveAbundanceFilter::SuperbubbleNode> nodes;
            unordered_set<u_int32_t> isVisited;
            unordered_set<u_int32_t> seen;
            unordered_map<u_int32_t, u_int64_t> pathLength;
            vector<UnitigGraph::Node*> queue;

            queue.push_back(nodeSource);
            pathLength[nodeSource->_unitigIndex] = 0;
            seenNodes.push_back(nodeSource->_unitigIndex);

            while(queue.size() > 0){

                UnitigGraph::Node* vNode = queue[queue.size()-1];
                u_int32_t v = vNode->_unitigIndex;


                //nodes.push_back({vNode, (u_int64_t)vNode->_nodes.size(), (u_int64_t)vNode->successorSum(), (u_int64_t)vNode->predecessorSum()});

                queue.pop_back();

                if(pathLength[v] > _maxLength){
                    //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                    return nullptr;
                }

                isVisited.insert(v);
                if(seen.find(v) != seen.end()){
                    seen.erase(v);
                }

                //vector<u_int32_t> successors;
                //_graph->getSuccessors_unitig(v, 0, successors);

                if(vNode->_successors.size() == 0){
                    //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                    return nullptr; //abort tip
                }


                //std::sort(successors.begin(),successors.end(), [this](const u_int32_t& unitigIndex1, const u_int32_t& unitigIndex2) {
                //    return _graph->_unitigs[unitigIndex1]._startNode < _graph->_unitigs[unitigIndex2]._startNode;
                //});

                //for(const u_int32_t& u : successors){
                for(UnitigGraph::Node* uNode : vNode->_successors){


                    u_int32_t u = uNode->_unitigIndex;
                    seenNodes.push_back(u);
                    

               

                    //if(isBubble[_unitigs[u]._startNode]) return -1;
                    //if(isBubble[nodeIndex_toReverseDirection(_unitigs[u]._startNode)]) return -1;
                    if(uNode == nodeSource){
                        //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                        return nullptr; //cycle including s
                    }

                    if(isVisited.find(u) == isVisited.end()){
                        seen.insert(u);
                        //long length = uNode->_length - _unitigGraph->_kminmerOverlapMean;
                        //_logFile << _unitigs[u]._length << " " << node._overlap << endl;
                        //if(length < 0) length = _unitigs[u]._length;
                        //pathLength[u] = pathLength[v] + length;
                        //pathLength[u] = pathLength[v] + _unitigs[u]._nbNodes;

                        if(vNode == nodeSource){
                            pathLength[u] = uNode->_length;
                        }
                        else{
                            pathLength[u] = pathLength[v] + (uNode->_nodes.size() * _unitigGraph->_kminmerLengthNonOverlap);
                        }
                    }
                    else{

                        //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                        //_logFile << "exit 4" << endl;
                        return nullptr; //Cycle within superbubble
                    }

                }


                for(UnitigGraph::Node* uNode : vNode->_successors){

                    u_int32_t u = uNode->_unitigIndex;
                    //u_int32_t u = node._index;
                    //_logFile << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

                    //vector<u_int32_t> predecessors;
                    //_graph->getPredecessors_unitig(u, 0, predecessors);
                    bool allPredecessorsAreVisited = true;
                    //for(u_int32_t p : predecessors){
                    for(UnitigGraph::Node* pNode : uNode->_predecessors){
                        if(isVisited.find(pNode->_unitigIndex) == isVisited.end()){
                            allPredecessorsAreVisited = false;
                            break;
                        }
                    }

                    if(allPredecessorsAreVisited){
                        //_logFile << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                        queue.push_back(uNode);
                    }

                    //_logFile << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                    if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]->_unitigIndex) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                        
                        UnitigGraph::Node* tNode = queue[0];
                        u_int32_t t = tNode->_unitigIndex;
                        //vector<u_int32_t> successors_t;
                        //_graph->getSuccessors_unitig(t, 0, successors_t);
                        

                        if(std::find(tNode->_successors.begin(), tNode->_successors.end(), nodeSource) == tNode->_successors.end()){
                            return tNode;
                        }
                        else{
                            //_logFile << "exit 5" << endl;
                            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                            return nullptr; // cycle including s
                        }

                    }


                }

            }


            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
            return nullptr;

        }
        
        void collapseSuperbubble(UnitigGraph::Node* sourceNode, UnitigGraph::Node* exitNode){

            vector<UnitigGraph::Node*> superbubbleNodes;
            collectSuperbubbleNodes(sourceNode, exitNode, superbubbleNodes);

            unordered_set<u_int32_t> keepNodes;
            //vector<u_int32_t> nodes;

            //Choose one path in the superbubble
            //u_int32_t unitigIndex = source_unitigIndex;
            UnitigGraph::Node* nodeCurrent = sourceNode;
            //vector<UnitigGraph::Node*> path = {nodeCurrent};


            while(true){

                //vector<u_int32_t> successors;
                //_graph->getSuccessors_unitig(unitigIndex, 0, successors);


                vector<UnitigGraph::Node*> maxAbNodes;
                float maxAbundance = 0;
                u_int32_t maxV = -1;

                for(UnitigGraph::Node* nn : nodeCurrent->_successors){
                    
                    if(nn->_abundance == maxAbundance){
                        maxAbNodes.push_back(nn);
                    }
                    else if(nn->_abundance > maxAbundance){
                        maxAbundance = nn->_abundance;
                        maxAbNodes.clear();
                        maxAbNodes.push_back(nn);
                    }
                }

                if(maxAbNodes.size() == 1){
                    keepNodes.insert(maxAbNodes[0]->_unitigIndex);
                    //path.push_back(maxAbNodes[0]);
                    nodeCurrent = maxAbNodes[0];
                }
                else{

                    
                    vector<BubbleSide> bubbleSides;
                    for(UnitigGraph::Node* nn : maxAbNodes){

                        
                        //u_int64_t sum = 0;
                        //for(u_int32_t nodeIndex : _unitigs[unitigIndex]._nodes){
                        //    sum += nodeIndex;
                        //}
                        
                        bubbleSides.push_back({nn, nn->_nodes.size(), nn->startNode()}); //_unitigs[unitigIndex]._quality
                    }


                    std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);
                    //path.push_back(bubbleSides[0]._unitig);
                    keepNodes.insert(bubbleSides[0]._unitig->_unitigIndex);
                    //_logFile << "Keep: " << bubbleSides[0]._unitigIndex << endl;

                    nodeCurrent = bubbleSides[0]._unitig;

                    //_logFile << "Keep node: " << _unitigs[unitigIndex]._startNode << endl;
                    //_logFile << _unitigs[unitigIndex]._nbNodes;
                    //getchar();
                }
                

                if(nodeCurrent == exitNode) break;
            }

            for(UnitigGraph::Node* node : superbubbleNodes){
                if(keepNodes.find(node->_unitigIndex) != keepNodes.end()) continue;
                
                _unitigGraph->removeNode(node);
            }
            
        }

        void collectSuperbubbleNodes(UnitigGraph::Node* sourceNode, UnitigGraph::Node* exitNode, vector<UnitigGraph::Node*>& nodes){

            nodes.clear();

            unordered_set<u_int32_t> isVisited;
            queue<UnitigGraph::Node*> queue;

            isVisited.insert(sourceNode->_unitigIndex);
            isVisited.insert(exitNode->_unitigIndex);
            queue.push(sourceNode);

            while (!queue.empty()){

                UnitigGraph::Node* currentNode = queue.front();
                u_int32_t unitigIndex_current = currentNode->_unitigIndex;
                queue.pop();

                for(UnitigGraph::Node* nn : currentNode->_successors){

                    if (isVisited.find(nn->_unitigIndex) != isVisited.end()) continue;

                    queue.push(nn);
                    isVisited.insert(nn->_unitigIndex);
                    nodes.push_back(nn);

                }
            }

        }

    };
    */
    
    struct SuperbubbleWriter{
        u_int64_t _readIndex;
        UnitigGraph::Node* _unitigSource;
        UnitigGraph::Node* _unitigExit;
        vector<UnitigGraph::Node*> _unitigRemoved;
        vector<u_int32_t> _seenNodes;
        //u_int32_t _prevNodeIndex;
    };

    struct SuperbubbleWriter_Comparator {
        bool operator()(SuperbubbleWriter const& p1, SuperbubbleWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

    priority_queue<SuperbubbleWriter, vector<SuperbubbleWriter> , SuperbubbleWriter_Comparator> _readWriterQueue;
    u_int64_t _nextReadIndexWriter;
    u_int64_t _nbBubbleRemoved;

    class SuperbubbleRemoverOld{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph* _unitigGraph;
        u_int64_t _maxLength;
        unordered_set<UnitigGraph::Node*> _isUnitigBubble;


        

        //std::set<UnitigGraph::Node*, SuperbubbleComparator> _queue;
        //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , SuperbubbleComparator> _queue;
        //u_int64_t _nextReadIndexWriter;
        //u_int64_t _nbRemoved;
        //vector<UnitigGraph::Node*> _queue;

        SuperbubbleRemoverOld(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph = _progressiveAbundanceFilter->_unitigGraph;
            _maxLength = maxLength;
        }

        bool execute(){

            

            
            

            _progressiveAbundanceFilter->_debugOuputFile = ofstream("/home/gats/workspace/tmp/lala1.txt");
            //bool isModification = removeSuperbubble();

            //cout << removeSuperbubble() << endl;
            //getchar();
            
            bool isModification = false;

            while(true){

                /*
                for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                    for(u_int32_t nodeIndex : node->_nodes){
                        if(nodeIndex == 11231202){

                            cout << "haha " << node->_unitigIndex << endl;
                            for(u_int32_t nodeIndex : node->_nodes){
                                if(nodeIndex == 11231202){
                                    cout << nodeIndex << " *" << endl;
                                }
                                else{
                                    cout << nodeIndex << endl;
                                }
                            }

                            for(UnitigGraph::Node* nn : node->_successors){
                                cout << "S: " << nn->startNode() << endl;
                            }
                            for(UnitigGraph::Node* nn : node->_predecessors){
                                cout << "P: " << nn->startNode() << endl;
                            }
                        }
                    }
                }
                */

                /*
                for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){
                    if(_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()){
                        if(_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)->_unitigIndex) == _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()){
                            cout << "outch" << endl;
                            getchar();
                        }    
                    }
                }
                */


                auto start = high_resolution_clock::now();
                
                //_nbRemoved = 0;
                bool isMod = removeSuperbubble();


                cout << "Nb superbubble removed: " << _progressiveAbundanceFilter->_nbBubbleRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
                /*

                int nbUnitigs = 0;
                for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                    if(node->_unitigIndex == -1) continue;
                    nbUnitigs += 1;
                }

                cout << "Nb unitigs: " << nbUnitigs << endl;
                cout << "Checksum ab: " << _unitigGraph->computeChecksum_abundance() << endl;
                */

                /*
                cout << "Checksum ab: " << _unitigGraph->computeChecksum_abundance() << endl;
                _progressiveAbundanceFilter->_debugOuputFile = ofstream("/home/gats/workspace/tmp/lala1.txt");

                vector<UnitigGraph::Node*> nodes;
                for(UnitigGraph::Node* node : _unitigGraph->_nodes){
                    if(node->_unitigIndex == -1) continue;
                    if(node->isCircular()) continue;
                    nodes.push_back(node);
                }

                std::sort(nodes.begin(), nodes.end(), UnitigGraph::NodeComparator);
                for(UnitigGraph::Node* node : nodes){

                    _progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << node->_length << " " << node->_abundance << endl;
                    

                }

                _progressiveAbundanceFilter->_debugOuputFile.close();
                

                for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                    for(u_int32_t nodeIndex : node->_nodes){
                        if(nodeIndex == 11822266){

                            cout << node->isCircular() << endl;
                            cout << "haha " << node->_unitigIndex << endl;
                            for(u_int32_t nodeIndex : node->_nodes){
                                if(nodeIndex == 11822266){
                                    cout << nodeIndex << " *" << endl;
                                }
                                else{
                                    cout << nodeIndex << endl;
                                }
                            }
                        }
                    }
                }
                

                _progressiveAbundanceFilter->_debugOuputFile.close();
                cout << "up" << endl;
                getchar();
                */

                //for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                //    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                //    if(node->startNode() == 1845619) cout << "yes" << endl;
                //}

                if(isMod) isModification = true;
                if(!isMod) break;
            }
            




            

                    


            

            
            //exit(1);

            /*
            _progressiveAbundanceFilter->_debugOuputFile = ofstream("/home/gats/workspace/tmp/lala1.txt");

            vector<UnitigGraph::Node*> nodes = _unitigGraph->_nodes;
            std::sort(nodes.begin(), nodes.end(), UnitigGraph::NodeComparator);
            for(size_t i=0; i<nodes.size(); i++){

                UnitigGraph::Node* node = nodes[i];
                if(node->_unitigIndex == -1){
                    _progressiveAbundanceFilter->_debugOuputFile << i << " " << node->startNode() << " " << node->endNode() << " " << " removed" << endl;
                }
                else{
                    _progressiveAbundanceFilter->_debugOuputFile << i << " " << node->startNode() << " " << node->endNode() << " " << node->_length << " " << node->_abundance << endl;
                }

            }

            _progressiveAbundanceFilter->_debugOuputFile.close();
            cout << "up" << endl;
            //exit(1);
            */

            return isModification;
        }

        bool removeSuperbubble(){
           
            vector<UnitigGraph::Node*> queue;

            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){

                if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1) continue;

                if(_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()){
                    continue;
                }

                queue.push_back(node);
            }
            

            std::sort(queue.begin(), queue.end(), SuperbubbleComparator);

            _isUnitigBubble.clear();
            _progressiveAbundanceFilter->_nextReadIndexWriter = 0;
            _progressiveAbundanceFilter->_nbBubbleRemoved = 0;

            BubbleFunctor functor(*this);
            size_t i = 0;

            #pragma omp parallel num_threads(_progressiveAbundanceFilter->_nbCores)
            {

                BubbleFunctor functorSub(functor);
                SuperbubbleWriter sw;

                while(true){

                    UnitigGraph::Node* unitig = nullptr;
                    bool isEof = false;

                    #pragma omp critical
                    {
                        
                        if(i >= queue.size()){
                            isEof = true;
                        }
                        else{
                            unitig = queue[i];
                            //_logFile << nodeIndex << endl;
                        }

                        //_logFile << unitigIndex << " " << _unitigs[unitigIndex]._startNode << endl;
                        //_logFile << _unitigs[unitigIndex]._length << " " << _unitigs[unitigIndex]._startNode << endl;
                        //getchar();
                        //_logFile << i << endl;
                        //sw = {i, unitig, {}};
                        sw._readIndex = i;
                        sw._unitigSource = unitig;
                        sw._unitigExit = nullptr;
                        i += 1;
                    }

                    if(isEof) break;
                    functorSub(sw);

                }
                

            }


            /*
            bool isModification = false;
            unordered_set<UnitigGraph::Node*> isUnitigBubble;

            bool lala = false;
            
            for(UnitigGraph::Node* node : queue){


                if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1) continue;

                if(isUnitigBubble.find(node) != isUnitigBubble.end()) continue;
                if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)) != isUnitigBubble.end()) continue;
                
                //if(node->startNode() == 11231202){
                //    cout << "2" << endl;
                //}

                //cout << node->_nodes.size() << endl;

                //cout << "is bubble? " << node->_unitigIndex << " " << (_lol.find(node->_unitigIndex) != _lol.end()) << " " << (_lol.find(_unitigGraph->unitigIndex_toReverseDirection(node->_unitigIndex)) != _lol.end()) << endl;
                vector<u_int32_t> seenNodes;
                UnitigGraph::Node* nodeExit = isSuperbubble(node, isUnitigBubble, seenNodes);



                if(nodeExit == nullptr){
                    indexNotSuperbubble(node, seenNodes);
                    continue;
                }
                
                if(nodeExit == _unitigGraph->unitigIndex_toReverseDirection(node)){
                    indexNotSuperbubble(node, seenNodes);
                    continue; //loop side of an inverse repeat
                }
                if(nodeExit == node){
                    if(node->_length < 400000){
                        indexNotSuperbubble(node, seenNodes);
                        continue;
                    }
                }
                if(isUnitigBubble.find(nodeExit) != isUnitigBubble.end()){
                    indexNotSuperbubble(node, seenNodes);
                    continue;
                }
                if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(nodeExit)) != isUnitigBubble.end()){
                    indexNotSuperbubble(node, seenNodes);
                    continue;
                }

                //if(node->startNode() == 150812){
                //    cout << "3    " << nodeExit->startNode() << " " << nodeExit->endNode() << endl;
                //}

                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << nodeExit->startNode() << " " << nodeExit->endNode() << endl;
                //cout << "\tRemove superbubble: " << node->_unitigIndex << " " << nodeExit->_unitigIndex << endl;

                isModification = true;

                //if(node->startNode() == 1845619){
                //    cout << "4" << endl;
                //}

                collapseSuperbubble(node, nodeExit, isUnitigBubble);
                _nbRemoved += 1;


            }
            */

            unordered_set<UnitigGraph::Node*> recompactNodes;

            //cout << "lala 1" << endl;
            for(UnitigGraph::Node* node : _isUnitigBubble){
                if(node->_unitigIndex == -1) continue;

                vector<UnitigGraph::Node*> succs = node->_successors;
                vector<UnitigGraph::Node*> preds = node->_predecessors;
                _unitigGraph->removeNode(node);

                for(UnitigGraph::Node* predecessor : preds){
                    if(predecessor->_unitigIndex == -1) continue;
                    recompactNodes.insert(predecessor);
                }

                for(UnitigGraph::Node* successor : succs){
                    if(successor->_unitigIndex == -1) continue;
                    recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(successor));
                }

            }
            
            //cout << "lala 2" << endl;
            for(UnitigGraph::Node* node : recompactNodes){
                if(node->_unitigIndex == -1) continue;
                _unitigGraph->recompact(node);
            }

            //cout << "lala 3" << endl;

            //exit(1);

            return recompactNodes.size() > 0;
        }



        class BubbleFunctor {

            public:

            SuperbubbleRemoverOld& _parent;

            BubbleFunctor(SuperbubbleRemoverOld& parent) : _parent(parent){
            }

            BubbleFunctor(const BubbleFunctor& copy) : _parent(copy._parent){
            }

            ~BubbleFunctor(){
            }

            void operator () (SuperbubbleWriter sw) {

                UnitigGraph::Node* node = sw._unitigSource;



                //if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1){
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return;
                }


                bool exist = false;
                #pragma omp critical(superbubble)
                {
                    if(_parent._isUnitigBubble.find(node) != _parent._isUnitigBubble.end() || _parent._isUnitigBubble.find(_parent._unitigGraph->unitigIndex_toReverseDirection(node)) != _parent._isUnitigBubble.end()) exist = true;
                }
                
                if(exist){
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return;
                }
                //if(node->startNode() == 11231202){
                //    cout << "2" << endl;
                //}

                //cout << node->_nodes.size() << endl;

                //cout << "is bubble? " << node->_unitigIndex << " " << (_lol.find(node->_unitigIndex) != _lol.end()) << " " << (_lol.find(_unitigGraph->unitigIndex_toReverseDirection(node->_unitigIndex)) != _lol.end()) << endl;
                vector<u_int32_t> seenNodes;
                UnitigGraph::Node* nodeExit = isSuperbubble(node, _parent._isUnitigBubble, seenNodes);

                //sw._unitigExit = nodeExit;
                sw._seenNodes = seenNodes;


                if(nodeExit == nullptr){
                    //_parent._progressiveAbundanceFilter->indexNotSuperbubble(node, seenNodes);
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return;
                }
                
                if(nodeExit == _parent._unitigGraph->unitigIndex_toReverseDirection(node)){
                    //_parent._progressiveAbundanceFilter->indexNotSuperbubble(node, seenNodes);
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return; //loop side of an inverse repeat
                }
                if(nodeExit == node){
                    if(node->_length < 400000){
                        //_parent._progressiveAbundanceFilter->indexNotSuperbubble(node, seenNodes);
                        _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                        return;
                    }
                }

                /*
                if(isUnitigBubble.find(nodeExit) != isUnitigBubble.end()){
                    _parent._progressiveAbundanceFilter->indexNotSuperbubble(node, seenNodes);
                    return;
                }
                if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(nodeExit)) != isUnitigBubble.end()){
                    _parent._progressiveAbundanceFilter->indexNotSuperbubble(node, seenNodes);
                    return;
                }
                */

                //if(node->startNode() == 150812){
                //    cout << "3    " << nodeExit->startNode() << " " << nodeExit->endNode() << endl;
                //}

                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << nodeExit->startNode() << " " << nodeExit->endNode() << endl;
                //cout << "\tRemove superbubble: " << node->_unitigIndex << " " << nodeExit->_unitigIndex << endl;

                //isModification = true;

                //if(node->startNode() == 1845619){
                //    cout << "4" << endl;
                //}

                collapseSuperbubble(node, nodeExit, _parent._isUnitigBubble, sw);
                //_nbRemoved += 1;

                /*
                //if(_parent._isUnitigBubble.find(node) != _parent._isUnitigBubble.end()) return;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)) != isUnitigBubble.end()) continue;

                UnitigGraph::Node* bubbleSide = isBubble(node);


                if(bubbleSide == nullptr){
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return;
                }
                //if(isUnitigBubble.find(bubbleSide) != isUnitigBubble.end()) continue;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(bubbleSide)) != isUnitigBubble.end()) continue;

                //if(_parent._isUnitigBubble.find(bubbleSide) != _parent._isUnitigBubble.end()) return;

                //isModification = true;
                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << bubbleSide->startNode() << " " << bubbleSide->endNode() << endl;

                //_unitigGraph->removeNode(bubbleSide);
                //_unitigGraph->recompact(node);
                //_parent._isUnitigBubble.insert(bubbleSide);
                //_parent._isUnitigBubble.insert(_parent._unitigGraph->unitigIndex_toReverseDirection(bubbleSide));
                //_nbRemoved += 1;
                sw._unitigRemoved.push_back(bubbleSide);
                sw._unitigRemoved.push_back(_parent._unitigGraph->unitigIndex_toReverseDirection(bubbleSide));
                _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                */
            }


            UnitigGraph::Node* isSuperbubble(UnitigGraph::Node* nodeSource, unordered_set<UnitigGraph::Node*>& isUnitigBubble, vector<u_int32_t>& seenNodes){

                

                seenNodes.clear();

                //vector<ProgressiveAbundanceFilter::SuperbubbleNode> nodes;
                unordered_set<u_int32_t> isVisited;
                unordered_set<u_int32_t> seen;
                unordered_map<u_int32_t, u_int64_t> pathLength;
                vector<UnitigGraph::Node*> queue;

                queue.push_back(nodeSource);
                pathLength[nodeSource->_unitigIndex] = 0;
                //nodes.push_back({nodeSource, (u_int16_t)nodeSource->_nodes.size(), (u_int64_t)nodeSource->successorSum(), (u_int64_t)nodeSource->predecessorSum()});
                seenNodes.push_back(nodeSource->_unitigIndex);

                while(queue.size() > 0){

                    UnitigGraph::Node* vNode = queue[queue.size()-1];
                    u_int32_t v = vNode->_unitigIndex;

                    
                    //if(nodeSource->startNode() == 150812){
                    //    cout << vNode->startNode() << " " << pathLength[v] << " " << vNode->_successors.size() << endl;
                    //}

                    //nodes.push_back({vNode, (u_int16_t)vNode->_nodes.size(), (u_int16_t)vNode->_successors.size(), (u_int16_t)vNode->_predecessors.size()});

                    queue.pop_back();

                    if(pathLength[v] > _parent._maxLength){
                        //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                        return nullptr;
                    }

                    isVisited.insert(v);
                    if(seen.find(v) != seen.end()){
                        seen.erase(v);
                    }

                    //vector<u_int32_t> successors;
                    //_graph->getSuccessors_unitig(v, 0, successors);

                    if(vNode->_successors.size() == 0){
                        //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                        return nullptr; //abort tip
                    }


                    //std::sort(successors.begin(),successors.end(), [this](const u_int32_t& unitigIndex1, const u_int32_t& unitigIndex2) {
                    //    return _graph->_unitigs[unitigIndex1]._startNode < _graph->_unitigs[unitigIndex2]._startNode;
                    //});

                    //for(const u_int32_t& u : successors){
                    for(UnitigGraph::Node* uNode : vNode->_successors){

                        seenNodes.push_back(uNode->_unitigIndex);
                        /*
                        bool exist = false;
                        for(const ProgressiveAbundanceFilter::SuperbubbleNode& sNode : nodes){
                            if(sNode._node == uNode){
                                exist = true;
                                break;
                            }
                        }

                        if(!exist){
                            nodes.push_back({uNode, (u_int16_t)uNode->_nodes.size(), (u_int64_t)uNode->successorSum(), (u_int64_t)uNode->predecessorSum()});
                        }
                        */

                        /*
                        if(isUnitigBubble.find(uNode) != isUnitigBubble.end()){
                            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                            return nullptr;
                        }
                        if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(uNode)) != isUnitigBubble.end()){
                            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                            return nullptr;
                        }
                        */


                        u_int32_t u = uNode->_unitigIndex;
                        

                        //if(isBubble[_unitigs[u]._startNode]) return -1;
                        //if(isBubble[nodeIndex_toReverseDirection(_unitigs[u]._startNode)]) return -1;
                        if(uNode == nodeSource){
                            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                            return nullptr; //cycle including s
                        }

                        if(isVisited.find(u) == isVisited.end()){
                            seen.insert(u);
                            //long length = uNode->_length - _unitigGraph->_kminmerOverlapMean;
                            //_logFile << _unitigs[u]._length << " " << node._overlap << endl;
                            //if(length < 0) length = _unitigs[u]._length;
                            //pathLength[u] = pathLength[v] + length;
                            //pathLength[u] = pathLength[v] + _unitigs[u]._nbNodes;


                            if(vNode == nodeSource){
                                pathLength[u] = uNode->_length;
                            }
                            else{
                                pathLength[u] = pathLength[v] + (uNode->_nodes.size() * _parent._unitigGraph->_kminmerLengthNonOverlap);
                            }
                            
                        }
                        else{

                            //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                            //_logFile << "exit 4" << endl;
                            return nullptr; //Cycle within superbubble
                        }

                    }


                    for(UnitigGraph::Node* uNode : vNode->_successors){

                        u_int32_t u = uNode->_unitigIndex;
                        //u_int32_t u = node._index;
                        //_logFile << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

                        //vector<u_int32_t> predecessors;
                        //_graph->getPredecessors_unitig(u, 0, predecessors);
                        bool allPredecessorsAreVisited = true;
                        //for(u_int32_t p : predecessors){
                        for(UnitigGraph::Node* pNode : uNode->_predecessors){
                            if(isVisited.find(pNode->_unitigIndex) == isVisited.end()){
                                allPredecessorsAreVisited = false;
                                break;
                            }
                        }

                        if(allPredecessorsAreVisited){
                            //_logFile << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                            queue.push_back(uNode);
                        }

                        //_logFile << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                        if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]->_unitigIndex) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                            
                            UnitigGraph::Node* tNode = queue[0];
                            u_int32_t t = tNode->_unitigIndex;
                            //vector<u_int32_t> successors_t;
                            //_graph->getSuccessors_unitig(t, 0, successors_t);
                            

                            if(std::find(tNode->_successors.begin(), tNode->_successors.end(), nodeSource) == tNode->_successors.end()){
                                return tNode;
                            }
                            else{
                                //_logFile << "exit 5" << endl;
                                //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                                return nullptr; // cycle including s
                            }

                        }


                    }

                }


                //_progressiveAbundanceFilter->_superBubbleIndex[nodeSource->_unitigIndex] = nodes;
                return nullptr;

            }
            
            void collapseSuperbubble(UnitigGraph::Node* sourceNode, UnitigGraph::Node* exitNode, unordered_set<UnitigGraph::Node*>& isUnitigBubble, SuperbubbleWriter sw){

                vector<UnitigGraph::Node*> superbubbleNodes;
                collectSuperbubbleNodes(sourceNode, exitNode, superbubbleNodes);

                unordered_set<u_int32_t> keepNodes;
                //vector<u_int32_t> nodes;

                //Choose one path in the superbubble
                //u_int32_t unitigIndex = source_unitigIndex;
                UnitigGraph::Node* nodeCurrent = sourceNode;
                //vector<UnitigGraph::Node*> path = {nodeCurrent};

                /*
                bool print = false;
                for(UnitigGraph::Node* node : superbubbleNodes){

                    for(u_int32_t nodeIndex : node->_nodes){
                        if(nodeIndex == 11822266 || nodeIndex == 11822267){
                            cout << "-----" << endl;
                            print = true;
                        }
                    }
                }
                */


                while(true){
                    //vector<u_int32_t> successors;
                    //_graph->getSuccessors_unitig(unitigIndex, 0, successors);

                    /*
                    if(print){
                        cout << "----------------------" << endl;
                        for(u_int32_t nodeIndex : nodeCurrent->_nodes){
                            cout << nodeIndex << endl;
                        }
                    }
                    */



                    vector<UnitigGraph::Node*> maxAbNodes;
                    float maxAbundance = 0;
                    u_int32_t maxV = -1;

                    for(UnitigGraph::Node* nn : nodeCurrent->_successors){
                        
                        if(nn->_abundance == maxAbundance){
                            maxAbNodes.push_back(nn);
                        }
                        else if(nn->_abundance > maxAbundance){
                            maxAbundance = nn->_abundance;
                            maxAbNodes.clear();
                            maxAbNodes.push_back(nn);
                        }
                    }

                    if(maxAbNodes.size() == 1){
                        keepNodes.insert(maxAbNodes[0]->_unitigIndex);
                        //path.push_back(maxAbNodes[0]);
                        nodeCurrent = maxAbNodes[0];
                    }
                    else{

                        
                        vector<BubbleSide> bubbleSides;
                        for(UnitigGraph::Node* nn : maxAbNodes){

                            
                            //u_int64_t sum = 0;
                            //for(u_int32_t nodeIndex : _unitigs[unitigIndex]._nodes){
                            //    sum += nodeIndex;
                            //}
                            
                            bubbleSides.push_back({nn, nn->_nodes.size(), nn->startNode()}); //_unitigs[unitigIndex]._quality
                        }


                        std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);
                        //path.push_back(bubbleSides[0]._unitig);
                        keepNodes.insert(bubbleSides[0]._unitig->_unitigIndex);
                        //_logFile << "Keep: " << bubbleSides[0]._unitigIndex << endl;

                        nodeCurrent = bubbleSides[0]._unitig;

                        //_logFile << "Keep node: " << _unitigs[unitigIndex]._startNode << endl;
                        //_logFile << _unitigs[unitigIndex]._nbNodes;
                        //getchar();
                    }
                    

                    if(nodeCurrent == exitNode) break;
                }


                for(UnitigGraph::Node* node : superbubbleNodes){

                    //for(u_int32_t nodeIndex : node->_nodes){
                    //    if(nodeIndex == 11822266){
                    //        cout << "lala" << endl;
                    //    }
                    //}

                    if(keepNodes.find(node->_unitigIndex) != keepNodes.end()) continue;
                    
                    sw._unitigRemoved.push_back(node);
                    //sw._unitigRemoved.push_back(_parent._unitigGraph->unitigIndex_toReverseDirection(node));
                    //isUnitigBubble.insert(node);
                    //isUnitigBubble.insert(_unitigGraph->unitigIndex_toReverseDirection(node));
                    //_unitigGraph->removeNode(node);
                }
                
                sw._unitigExit = exitNode;
                _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);

            }

            void collectSuperbubbleNodes(UnitigGraph::Node* sourceNode, UnitigGraph::Node* exitNode, vector<UnitigGraph::Node*>& nodes){

                nodes.clear();

                unordered_set<u_int32_t> isVisited;
                queue<UnitigGraph::Node*> queue;

                isVisited.insert(sourceNode->_unitigIndex);
                isVisited.insert(exitNode->_unitigIndex);
                queue.push(sourceNode);

                while (!queue.empty()){

                    UnitigGraph::Node* currentNode = queue.front();
                    u_int32_t unitigIndex_current = currentNode->_unitigIndex;
                    queue.pop();

                    for(UnitigGraph::Node* nn : currentNode->_successors){

                        if (isVisited.find(nn->_unitigIndex) != isVisited.end()) continue;

                        queue.push(nn);
                        isVisited.insert(nn->_unitigIndex);
                        nodes.push_back(nn);

                    }
                }

            }

        };
        

    };
    

    class BubbleRemover{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph* _unitigGraph;
        u_int64_t _maxLength;
        //u_int64_t _nbRemoved;


        BubbleRemover(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph = _progressiveAbundanceFilter->_unitigGraph;
            _maxLength = maxLength;
        }

        bool execute(){

            bool isModification = false;

            while(true){


                //_progressiveAbundanceFilter->_debugOuputFile = ofstream("/home/gats/workspace/tmp/lala1.txt");
                auto start = high_resolution_clock::now();
                
                //_nbRemoved = 0;
                bool isMod = removeBubbleOld();

                cout << "Nb bubble removed: " << _progressiveAbundanceFilter->_nbBubbleRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
                //_progressiveAbundanceFilter->_debugOuputFile.close();
                //cout << "up" << endl;
                //getchar();

                if(isMod) isModification = true;
                if(!isMod) break;
            }

            return isModification;
        }

        /*
        bool removeBubble(){


            //if(isUnitigBubble.find(node) != isUnitigBubble.end()) continue;
            //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)) != isUnitigBubble.end()) continue;

            bool isModification = false;
            unordered_set<UnitigGraph::Node*> isUnitigBubble;

            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){

                if(node == nullptr) continue;
                if(node->_unitigIndex == -1) continue;
                

                //cout << "is bubble? " << node->_unitigIndex << " " << (_lol.find(node->_unitigIndex) != _lol.end()) << " " << (_lol.find(_unitigGraph->unitigIndex_toReverseDirection(node->_unitigIndex)) != _lol.end()) << endl;
                UnitigGraph::Node* bubbleSide = isBubble(node);

                if(bubbleSide != nullptr){

                    //cout << "\tRemove bubble: " << node->_unitigIndex << " " << bubbleSide->_unitigIndex << endl;

                    isModification = true;

                    _unitigGraph->removeNode(bubbleSide);
                    _unitigGraph->recompact(node);

                    //getchar();
                    //for(UnitigGraph::Node* predecessor : node->_predecessors){
                    //    _unitigGraph->removePredecessor(node, predecessor);
                    //    _unitigGraph->recompact(predecessor);
                    //}

                    //return false;
                }
                else{
                    //cout << "\tno" << endl;
                }
            }

            return isModification;
        }
        */

        unordered_set<UnitigGraph::Node*> _isUnitigBubble;
        //bool _isModification;

        bool removeBubbleOld(){

            vector<UnitigGraph::Node*> queue;

            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){

                if(node == nullptr) continue;
                if(node->_unitigIndex == -1) continue;
                if(node->_successors.size() <= 1 || node->_successors.size() > 5) continue;
                
                queue.push_back(node);
            }

            std::sort(queue.begin(), queue.end(), SuperbubbleComparator);

            _isUnitigBubble.clear();
            _progressiveAbundanceFilter->_nextReadIndexWriter = 0;
            _progressiveAbundanceFilter->_nbBubbleRemoved = 0;
            

            /*

            for(UnitigGraph::Node* node : queue){
                
                if(isUnitigBubble.find(node) != isUnitigBubble.end()) continue;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)) != isUnitigBubble.end()) continue;

                UnitigGraph::Node* bubbleSide = isBubble(node, isUnitigBubble);


                if(bubbleSide == nullptr) continue;
                //if(isUnitigBubble.find(bubbleSide) != isUnitigBubble.end()) continue;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(bubbleSide)) != isUnitigBubble.end()) continue;

                if(isUnitigBubble.find(bubbleSide) != isUnitigBubble.end()) continue;

                isModification = true;
                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << bubbleSide->startNode() << " " << bubbleSide->endNode() << endl;

                //_unitigGraph->removeNode(bubbleSide);
                //_unitigGraph->recompact(node);
                isUnitigBubble.insert(bubbleSide);
                isUnitigBubble.insert(_unitigGraph->unitigIndex_toReverseDirection(bubbleSide));
                _nbRemoved += 1;
            }
            */

            BubbleFunctor functor(*this);
            size_t i = 0;

            #pragma omp parallel num_threads(_progressiveAbundanceFilter->_nbCores)
            {

                BubbleFunctor functorSub(functor);
                SuperbubbleWriter sw;

                while(true){

                    UnitigGraph::Node* unitig = nullptr;
                    bool isEof = false;

                    #pragma omp critical
                    {
                        
                        if(i >= queue.size()){
                            isEof = true;
                        }
                        else{
                            unitig = queue[i];
                            //_logFile << nodeIndex << endl;
                        }

                        //_logFile << unitigIndex << " " << _unitigs[unitigIndex]._startNode << endl;
                        //_logFile << _unitigs[unitigIndex]._length << " " << _unitigs[unitigIndex]._startNode << endl;
                        //getchar();
                        //_logFile << i << endl;
                        sw._readIndex = i;
                        sw._unitigSource = unitig;
                        sw._unitigExit = nullptr;
                        
                        i += 1;
                    }

                    if(isEof) break;
                    functorSub(sw);

                }
                

            }

            unordered_set<UnitigGraph::Node*> recompactNodes;

            //cout << "lala 1" << endl;
            for(UnitigGraph::Node* node : _isUnitigBubble){
                if(node->_unitigIndex == -1) continue;

                vector<UnitigGraph::Node*> succs = node->_successors;
                vector<UnitigGraph::Node*> preds = node->_predecessors;
                _unitigGraph->removeNode(node);

                for(UnitigGraph::Node* predecessor : preds){
                    if(predecessor->_unitigIndex == -1) continue;
                    recompactNodes.insert(predecessor);
                }

                for(UnitigGraph::Node* successor : succs){
                    if(successor->_unitigIndex == -1) continue;
                    recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(successor));
                }

            }
            
            //cout << "lala 2" << endl;
            for(UnitigGraph::Node* node : recompactNodes){
                if(node->_unitigIndex == -1) continue;
                _unitigGraph->recompact(node);
            }

            return _isUnitigBubble.size() > 0;
        }




        class BubbleFunctor {

            public:

            BubbleRemover& _parent;

            BubbleFunctor(BubbleRemover& parent) : _parent(parent){
            }

            BubbleFunctor(const BubbleFunctor& copy) : _parent(copy._parent){
            }

            ~BubbleFunctor(){
            }

            void operator () (SuperbubbleWriter sw) {

                UnitigGraph::Node* node = sw._unitigSource;

                //if(_parent._isUnitigBubble.find(node) != _parent._isUnitigBubble.end()) return;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(node)) != isUnitigBubble.end()) continue;

                UnitigGraph::Node* bubbleExit = nullptr;
                UnitigGraph::Node* bubbleSide = isBubble(node, bubbleExit);


                if(bubbleSide == nullptr){
                    _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
                    return;
                }
                //if(isUnitigBubble.find(bubbleSide) != isUnitigBubble.end()) continue;
                //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(bubbleSide)) != isUnitigBubble.end()) continue;

                //if(_parent._isUnitigBubble.find(bubbleSide) != _parent._isUnitigBubble.end()) return;

                //isModification = true;
                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << bubbleSide->startNode() << " " << bubbleSide->endNode() << endl;

                //_unitigGraph->removeNode(bubbleSide);
                //_unitigGraph->recompact(node);
                //_parent._isUnitigBubble.insert(bubbleSide);
                //_parent._isUnitigBubble.insert(_parent._unitigGraph->unitigIndex_toReverseDirection(bubbleSide));
                //_nbRemoved += 1;
                sw._unitigRemoved.push_back(bubbleSide);
                //sw._unitigRemoved.push_back(_parent._unitigGraph->unitigIndex_toReverseDirection(bubbleSide));
                sw._unitigExit = bubbleExit;

                _parent._progressiveAbundanceFilter->simplifySuperbubble(sw, _parent._isUnitigBubble);
            }

            UnitigGraph::Node* isBubble(UnitigGraph::Node* utg_1, UnitigGraph::Node* unitigExit){


                bool isBubble = false;
                unitigExit = nullptr;


                for(size_t i=0; i<utg_1->_successors.size(); i++) {
                    
                    UnitigGraph::Node* utg_2 = utg_1->_successors[i]; //_unitigGraph->_nodes[successors[i]];
                    //if(isUnitigBubble.find(utg_2) != isUnitigBubble.end()) continue;

                    for(size_t j=i+1; j<utg_1->_successors.size(); j++){

                        UnitigGraph::Node* utg_3 = utg_1->_successors[j]; //_unitigGraph->_nodes[successors[j]];
                        //if(isUnitigBubble.find(utg_3) != isUnitigBubble.end()) continue;
                        
                        if(utg_2->_predecessors.size() != 1) continue;
                        if(utg_2->_successors.size() != 1) continue;

                        if(utg_3->_predecessors.size() != 1) continue;
                        if(utg_3->_successors.size() != 1) continue;

                        if(utg_2->_length > _parent._maxLength || utg_3->_length > _parent._maxLength) continue;

                        if(utg_2->_successors[0] != utg_3->_successors[0]) continue;
                        UnitigGraph::Node* utg_4 = utg_2->_successors[0];
                        //if(isUnitigBubble.find(utg_4) != isUnitigBubble.end()) continue;

                        if(BiGraph::nodeIndex_to_nodeName(utg_1->endNode()) == BiGraph::nodeIndex_to_nodeName(utg_4->startNode())) continue; //Repeated unitig with a cycle on on side

                        unitigExit = utg_4;

                        if(utg_2->_abundance > utg_3->_abundance){
                            return utg_3;
                        }
                        else if(utg_2->_abundance < utg_3->_abundance){
                            return utg_2;
                        }
                        else{
                            
                            vector<BubbleSide> bubbleSides;
                            bubbleSides.push_back({utg_2, utg_2->_nodes.size(), utg_2->startNode()}); //_unitigs[utg_2._index]._quality
                            bubbleSides.push_back({utg_3, utg_3->_nodes.size(), utg_3->startNode()}); //_unitigs[utg_3._index]._quality

                            std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);

                            return bubbleSides[1]._unitig;
                        }
                        
                    }
                }

                return nullptr;
            }

        };


    };

    void simplifySuperbubble(SuperbubbleWriter swNext, unordered_set<UnitigGraph::Node*>& isUnitigBubble){


        #pragma omp critical(superbubble)
        {

            _readWriterQueue.push(swNext);

            while(!_readWriterQueue.empty()){

                const SuperbubbleWriter& sw = _readWriterQueue.top();

                if(sw._readIndex == _nextReadIndexWriter){



                    bool isValid = true;
                    
                    for(UnitigGraph::Node* node : sw._unitigRemoved){
                        if(isUnitigBubble.find(node) != isUnitigBubble.end()){
                            isValid = false;
                            break;
                        }
                    }
                    if(isUnitigBubble.find(sw._unitigSource) != isUnitigBubble.end()){
                        isValid = false;
                    }
                    //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(sw._unitigSource)) != isUnitigBubble.end()){
                    //    isValid = false;
                    //}
                    if(sw._unitigExit != nullptr){
                        if(isUnitigBubble.find(sw._unitigExit) != isUnitigBubble.end()){
                            isValid = false;
                        }
                    }
                    //if(isUnitigBubble.find(_unitigGraph->unitigIndex_toReverseDirection(sw._unitigExit)) != isUnitigBubble.end()){
                    //    isValid = false;
                    //}
                    //}

                    if(sw._unitigRemoved.size() == 0){
                        isValid = false;
                    }

                    if(!isValid && sw._seenNodes.size() > 0){
                        indexNotSuperbubble(sw._unitigSource, sw._seenNodes);
                    }

                    if(isValid){

                        _nbBubbleRemoved += 1;

                        for(UnitigGraph::Node* node : sw._unitigRemoved){
                            
                            isUnitigBubble.insert(node);
                            isUnitigBubble.insert(_unitigGraph->unitigIndex_toReverseDirection(node));
                        }




                    }
                    


                    _readWriterQueue.pop();
                    _nextReadIndexWriter += 1;
                }
                else{
                    break;
                }
            }

        }
    }

    void indexNotSuperbubble(UnitigGraph::Node* nodeSource, const vector<u_int32_t>& seenNodes){
        for(u_int32_t unitigIndex : seenNodes){
            _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[unitigIndex].push_back(nodeSource->_unitigIndex);
        }
        _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.insert(nodeSource->_unitigIndex);
    }

    
    class TipRemover{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph* _unitigGraph;
        u_int64_t _maxLength;

        //struct Tip{
        //    UnitigGraph::Node* _unitig;
        //};

        u_int64_t _nbTipRemoved;

        struct TipComparator {
            bool operator()(UnitigGraph::Node* p1, UnitigGraph::Node* p2) const{
                if(p1->_nodes.size() == p2->_nodes.size()){
                    if(p1->_abundance == p2->_abundance){
                        return p1->startNode() > p2->startNode();
                    }
                    return p1->_abundance > p2->_abundance;
                }
                return p1->_nodes.size() < p2->_nodes.size();
            }
        };


        //vector<UnitigGraph::Node*> _queue;
        //vector<UnitigGraph::Node*> _nextQueue;
        std::set<UnitigGraph::Node*, TipComparator> _queue;
        //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator> _queue;
        u_int64_t _nextReadIndexWriter;

        TipRemover(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph = _progressiveAbundanceFilter->_unitigGraph;
            _nbTipRemoved = 0;
            _maxLength = maxLength;
	        //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator>();
        }


        u_int64_t execute(){


            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){


                //if(node == nullptr) continue;
                //if(node->_length > _maxLength) continue;
                //if(node->_unitigIndex % 2 == 1) continue;

                if(isTipAny(node)){
                    _queue.insert(node);
                    //_nextQueue.push_back(node);
                }
            }  



            //std::sort(queue.begin(), queue.end(), SuperbubbleComparator);

            
            auto start = high_resolution_clock::now();
            bool isModification = removeTips(false, _maxLength);
            cout << "Nb Tip removed: " << _nbTipRemoved << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)" << endl;
                
            return isModification;
            
            

            
            //return nbremoved;
        }

        bool removeTips(bool useK, u_int64_t maxLength){

            //_queue = _nextQueue;
            //_nextQueue.clear();
            //u_int64_t nbRemoved = 0;




            //unordered_set<UnitigGraph::Node*> recompactNodes;
            bool isModification = false;
            bool lala = false;

            //for(UnitigGraph::Node* node : _queue){
            while(!_queue.empty()){

            

                auto top = _queue.begin();
                UnitigGraph::Node* node = *top;
                _queue.erase(top);

                if(!isTipAny(node)){
                    //_nextQueue.push_back(node);
                    continue;
                }

                //if(node->_unitigIndex == -1){
                //    continue;
                //}

                //_progressiveAbundanceFilter->_debugOuputFile << node->startNode() << " " << node->endNode() << " " << node->_nodes.size() << " " << node->_predecessors.size() << " " << node->_successors.size() << endl;
                
           
                //if(node->_unitigIndex == 7186 || node->_unitigIndex == 7187){
                //    cout << "omg" << endl;
                //    getchar();
                //}
                //cout << "Remove tip: " << node->_nodes.size() << endl;
                //cout << "Remove tip: " << node->_unitigIndex << " " << node->_abundance << " " << node->_length << endl;
                //cout << node->_length << endl;
                //nbRemoved += 1;
                isModification = true;
                _nbTipRemoved += 1;

                vector<UnitigGraph::Node*> preds = node->_predecessors;
                vector<UnitigGraph::Node*> nodeToAdd;




                for(UnitigGraph::Node* predecessor : preds){
                    
                    if(predecessor->_unitigIndex == -1) continue;
                
                    UnitigGraph::Node* predecessor_rc = _unitigGraph->unitigIndex_toReverseDirection(predecessor);
                    //_queue.erase(predecessor);
                    //_queue.erase(predecessor_rc);
                    //_progressiveAbundanceFilter->_debugOuputFile << "    R" << predecessor->startNode() << " " << predecessor->_nodes.size()  << endl;

                    //cout << predecessor->_unitigIndex << endl;


                

                    _unitigGraph->removePredecessor(node, predecessor);
                    //recompactNodes.insert(predecessor);
                    _unitigGraph->recompact(predecessor);


                    if(isTipAny(predecessor)){
                        _queue.insert(predecessor);
                    }
                    if(isTipAny(predecessor_rc)){
                        _queue.insert(predecessor_rc);
                    }
                    

                }
                

                //for(UnitigGraph::Node* node : nodeToAdd){
                    //_queue.push(node);
                //}

            }

         
            
            //return nbRemoved > 0;
            return isModification;
        }

        bool isTipAny(const UnitigGraph::Node* node){

            //if(node == nullptr) return false;
            if(node->_unitigIndex == -1) return false;
            
            if(node->_length > _maxLength) return false;
            
            //if(node->_unitigIndex % 2 == 1) return false;
            if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() != 1) return false;
            if(node->_predecessors.size() == 0) return false;

            if(node->_isPalindrome) return false;

            return true;
        }

        bool isTip(const UnitigGraph::Node* node, bool useK, u_int64_t maxLength){

            //if(node == nullptr) return false;
            if(node->_unitigIndex == -1) return false;

            
            if(useK){
                if(node->_nodes.size() > maxLength) return false;
            }
            else{
                if(node->_length > maxLength) return false;
            }
            
            //if(node->_length > 50000) return false;
            
            //if(node->_unitigIndex % 2 == 1) return false;
            if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() != 1) return false;
            if(node->_predecessors.size() == 0) return false;

            return true;
        }
    };
    


    /*
    class TipRemover{

        public:

        ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
        UnitigGraph* _unitigGraph;
        u_int64_t _maxLength;

        //struct Tip{
        //    UnitigGraph::Node* _unitig;
        //};

        u_int64_t _nbTipRemoved;

        struct TipComparator {
            bool operator()(UnitigGraph::Node* p1, UnitigGraph::Node* p2) const{
                if(p1->_nodes.size() == p2->_nodes.size()){
                    if(p1->_abundance == p2->_abundance){
                        return p1->startNode() > p2->startNode();
                    }
                    return p1->_abundance > p2->_abundance;
                }
                return p1->_nodes.size() < p2->_nodes.size();
            }
        };


        vector<UnitigGraph::Node*> _queue;
        vector<UnitigGraph::Node*> _nextQueue;
        //std::set<UnitigGraph::Node*, TipComparator> _queue;
        //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator> _queue;
        u_int64_t _nextReadIndexWriter;

        TipRemover(ProgressiveAbundanceFilter* progressiveAbundanceFilter, u_int64_t maxLength){
            _progressiveAbundanceFilter = progressiveAbundanceFilter;
            _unitigGraph = _progressiveAbundanceFilter->_unitigGraph;
            _nbTipRemoved = 0;
            _maxLength = maxLength;
	        //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , TipComparator>();
        }


        u_int64_t execute(){


            for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){


                //if(node == nullptr) continue;
                //if(node->_length > _maxLength) continue;
                //if(node->_unitigIndex % 2 == 1) continue;

                if(isTipAny(node)){
                    //_queue.insert(node);
                    _nextQueue.push_back(node);
                }
            }  



            //std::sort(_nextQueue.begin(), _nextQueue.end(), SuperbubbleComparator);


            
            
            
            bool isModification = false;
            
            while(true){
                
                bool isModSub = false;

                while(true){
                    _nbTipRemoved = 0;
                    bool isModif = removeTips(true, _progressiveAbundanceFilter->_kminmerSize);
                    if(isModif){
                        isModification = true;
                        isModSub = true;
                    }

                    if(!isModif) break;

                    cout << "Nb removed tips: " << _nbTipRemoved << endl;// << "    " << _unitigGraph->computeChecksum_abundance() << endl;


                }

                
                if(isModSub) continue;

                while(true){
                    _nbTipRemoved = 0;
                    bool isModif = removeTips(true, _progressiveAbundanceFilter->_kminmerSize*2);
                    if(isModif){
                        isModification = true;
                        isModSub = true;
                    }

                    if(!isModif) break;
                    cout << "Nb removed tips: " << _nbTipRemoved << endl;//<< "    " << _unitigGraph->computeChecksum_abundance() << endl;
                }

                if(isModSub) continue;

                while(true){
                    _nbTipRemoved = 0;
                    bool isModif = removeTips(false, _maxLength);
                    if(isModif){
                        isModification = true;
                        isModSub = true;
                    }

                    if(!isModif) break;
                    cout << "Nb removed tips: " << _nbTipRemoved << endl;//<< "    " << _unitigGraph->computeChecksum_abundance() << endl;
                }
                

                if(!isModSub) break;

            }


            //cout << "Nb tips removed: " << _nbTipRemoved << endl;
            return isModification;
            
        }

        bool removeTips(bool useK, u_int64_t maxLength){

            _queue = _nextQueue;
            _nextQueue.clear();
            //u_int64_t nbRemoved = 0;



            unordered_set<UnitigGraph::Node*> recompactNodes;
            bool isModification = false;
            bool lala = false;

            for(UnitigGraph::Node* node : _queue){


                if(!isTipAny(node)){
                    _nextQueue.push_back(node);
                    continue;
                }

                isModification = true;
                _nbTipRemoved += 1;

                //cout << "Tip: " << node->_unitigIndex << " " << node->_length << " " << node->_abundance << endl;

                
                //if(node->_abundance > 1000){
                //    _unitigGraph->save(_progressiveAbundanceFilter->_tmpDir + "/debugGraph.gfa", 0);
                //    getchar();
                //}
                
                vector<UnitigGraph::Node*> preds = node->_predecessors;

                for(UnitigGraph::Node* predecessor : preds){
                    
                    if(predecessor->_unitigIndex == -1) continue;
                
                    UnitigGraph::Node* predecessor_rc = _unitigGraph->unitigIndex_toReverseDirection(predecessor);
                    //_queue.erase(predecessor);
                    //_queue.erase(predecessor_rc);
                    //_progressiveAbundanceFilter->_debugOuputFile << "    R" << predecessor->startNode() << " " << predecessor->_nodes.size()  << endl;

                    //cout << predecessor->_unitigIndex << endl;




                    _unitigGraph->removePredecessor(node, predecessor);
                    recompactNodes.insert(predecessor);
                    //_unitigGraph->recompact(predecessor);


                    //if(isTipAny(predecessor)){
                    //    _queue.insert(predecessor);
                    //}
                    //if(isTipAny(predecessor_rc)){
                    //    _queue.insert(predecessor_rc);
                    //}
                    

                }
                

                //for(UnitigGraph::Node* node : nodeToAdd){
                    //_queue.push(node);
                //}

            }


            
            for(UnitigGraph::Node* node : recompactNodes){
                if(node->_unitigIndex == -1) continue;
                _unitigGraph->recompact(node);

                if(isTipAny(node)){
                    _nextQueue.push_back(node);
                }

                UnitigGraph::Node* node_rc = _unitigGraph->unitigIndex_toReverseDirection(node);

                if(isTipAny(node_rc)){
                    _nextQueue.push_back(node_rc);
                }

            }
            
            
            //return nbRemoved > 0;
            return isModification;
        }

        bool isTipAny(const UnitigGraph::Node* node){

            //if(node == nullptr) return false;
            if(node->_unitigIndex == -1) return false;
            
            if(node->_length > _maxLength) return false;
            
            //if(node->_unitigIndex % 2 == 1) return false;
            if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() != 1) return false;
            if(node->_predecessors.size() == 0) return false;

            if(node->_isPalindrome) return false;

            return true;
        }

        bool isTip(const UnitigGraph::Node* node, bool useK, u_int64_t maxLength){

            //if(node == nullptr) return false;
            if(node->_unitigIndex == -1) return false;

            
            if(useK){
                if(node->_nodes.size() > maxLength) return false;
            }
            else{
                if(node->_length > maxLength) return false;
            }
            
            //if(node->_length > 50000) return false;
            
            //if(node->_unitigIndex % 2 == 1) return false;
            if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() != 1) return false;
            if(node->_predecessors.size() == 0) return false;

            return true;
        }
    };
    */
    UnitigGraph* _unitigGraph;
    vector<SaveState> _cachedGraphStates;
    SaveState _currentSaveState;
    size_t _kminmerSize;

    /*
    struct SuperbubbleNode{
        UnitigGraph::Node* _node;
        u_int16_t _nbNodes;
        //u_int16_t _nbSuccessors;
        //u_int16_t _nbPredecessors; 
        u_int64_t _succSum; 
        u_int64_t _predSum;
    };
    */

    struct CutoffIndex{
        u_int32_t _cutoffIndex;
        float _cutoffValue;
    };


    
    //unordered_map<u_int32_t, vector<SuperbubbleNode>> _superBubbleIndex;
    u_int32_t _cutoffIndex;
    vector<CutoffIndex> _cutoffIndexes;
    string _tmpDir;
    ofstream _debugOuputFile;
    bool _removeBubble;
    int _nbCores;
	unordered_set<u_int32_t> _isNodeNamePalindrome;


    ProgressiveAbundanceFilter(UnitigGraph* unitigGraph, const string& tmpDir, size_t kminmerSize, int nbCores, bool removeBubble){
        _unitigGraph = unitigGraph;
        _tmpDir = tmpDir;
        _kminmerSize = kminmerSize;
        _removeBubble = removeBubble;
        _nbCores = nbCores;
    }

    vector<UnitigGraph::Node*> _validNodes;

	template<typename Functor = X>
	//void parse(const Functor& functor){
    void execute(Functor& functor){


        //_unitigGraph->save(_tmpDir + "/lala.gfa", 0);

        _checksumAbundanceTotal = 0;
        _checksumNodeTotal = 0;
        _cutoffIndex = 0;
        _debugOuputFile = ofstream("/home/gats/workspace/tmp/lala_" + to_string(_kminmerSize) + ".txt");

        u_int32_t maxAbundance = 0;
        for(UnitigGraph::Node* node : _unitigGraph->_nodes){
            if(node->_abundance > maxAbundance){
                maxAbundance = node->_abundance;
            }
        }
        

        cout << "Max abundance: " << maxAbundance << endl;

        simplifyProgressive(functor);

        cout << "Checksum: " << _checksumNodeTotal << endl;
        cout << "Checksum abundnace: " << _checksumAbundanceTotal << endl;
        //getchar();

        _debugOuputFile.close();

        //if(_kminmerSize > 4){
        //string command = "diff /home/gats/workspace/tmp/bug/lala_" + to_string(_kminmerSize) + ".txt" + " /home/gats/workspace/tmp/lala_" + to_string(_kminmerSize-1) + ".txt";
        //system(command.c_str());
        //getchar();
        //}
    }

    unordered_set<float> isCutoffProcessed;

	template<typename Functor = X>
    void simplifyProgressive(Functor& functor){

        float maxAbundance = 2000;
        float currentCutoff = 0;
        _currentSaveState = {0, {}, {}};

        while(true){

            bool isModification = false;

            bool isMod = simplify();
            if(isMod){
                isModification = true;
            }

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

            //cout << "saving state: " << currentCutoff << endl;

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
                _unitigGraph->save(_tmpDir + "/minimizer_graph_u_cleaned.gfa", 0);
                //exit(1);

                
                //exit(1);
            }
            //compact(true, unitigDatas);
            //if(doesSaveUnitigGraph && currentCutoff == 0) saveUnitigGraph(_outputDir + "/minimizer_graph_u_cleaned.gfa", mdbg, minimizerSize, nbCores, true);
            

            //for(u_int32_t nodeName : currentSaveState._nodeNameRemoved_tmp){
            //    currentSaveState._nodeNameRemoved.push_back(nodeName);
            //}
            //currentSaveState._nodeNameRemoved_tmp.clear();

            int nbUnitigs = 0;
            for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                if(node->_unitigIndex == -1) continue;
                nbUnitigs += 1;
            }

            dumpUnitigs(currentCutoff);

            _currentSaveState._abundanceCutoff_min = currentCutoff;
            _cachedGraphStates.push_back(_currentSaveState);
            _currentSaveState = {0, {}, {}};


        
        }

    }
    

    bool simplify(){

        u_int64_t tipPass = 0;
        bool isModification = false;
        static u_int64_t maxTipLength = 50000;
        static u_int64_t maxBubbleLength = 50000;


        while(true){

            bool isModificationSub = false;
            
            _validNodes.clear();
            for(UnitigGraph::Node* node : _unitigGraph->_nodes){

                if(node->_unitigIndex == -1) continue;

                _validNodes.push_back(node);
            }

            /*
            cout << "||||||||||||||||||||||||||||||" << _unitigGraph->_nodes[236530]->_successors.size() << " " << _unitigGraph->_nodes[236530]->_predecessors.size() << endl;
            for(UnitigGraph::Node* nn : _unitigGraph->_nodes[236530]->_successors){
                cout << "\tS:" << nn->_unitigIndex << " " << nn->_abundance << " " << nn->_length << endl;
            }
            for(UnitigGraph::Node* nn : _unitigGraph->_nodes[236530]->_predecessors){
                cout << "\tP:" << nn->_unitigIndex << " " << nn->_abundance << " " << nn->_length << endl;
            }
            */
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

                    
                    cout << "Superbubble" << endl;
                    SuperbubbleRemoverOld superbubbleRemover(this, maxBubbleLength);
                    bool isModificationSuperbubble = superbubbleRemover.execute();
                    if(isModificationSuperbubble){
                        isModification = true;
                        isModBubble = true;
                        isModificationSub = true;
                    }
                    
                    
                    
                    
                    //cout << "Nb superbubble removed: " << superbubbleRemover._nbRemovedTotal << endl;

                    //cout << _unitigGraph->computeChecksum_nodeIndex() << endl;
                    //cout << _unitigGraph->computeChecksum_abundance() << endl;
                    //cout << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s" << endl;
                    //getchar();
                    
                    
                    cout << "Bubble" << endl;
                    auto start = high_resolution_clock::now();
                    
                    BubbleRemover bubbleRemover(this, maxBubbleLength);
                    bool isModificationBubble = bubbleRemover.execute();
                    if(isModificationBubble){
                        isModification = true;
                        isModBubble = true;
                        isModificationSub = true;
                    }
                    
                    //cout << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s" << endl;
                    

                    if(!isModBubble) break;
                }
            }

            
            
            cout << "Tip " << tipPass << endl;
            //auto start = high_resolution_clock::now();
            
            TipRemover tipRemover(this, maxTipLength);
            u_int64_t nbRemoved = tipRemover.execute();
            if(nbRemoved){
                isModification = true;
                isModificationSub = true;
            }

            tipPass += 1;
            
            //_unitigGraph->save(_tmpDir + "/debugGraph.gfa", 0);
            //getchar();

            //cout << "Nb unitigs: " << _unitigGraph->nbUnitigs() << endl;
            

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

        return isModification;
    }  

    struct AbundanceComparator {
        bool operator()(UnitigGraph::Node* p1, UnitigGraph::Node* p2) const{
            if(p1->_nodes.size() == p2->_nodes.size()){
                return p1->startNode() > p2->startNode();
            }
            return p1->_nodes.size() < p2->_nodes.size();
        }
    };


	template<typename Functor = X>
    u_int64_t removeAbundanceNoQueue(float abundanceCutoff_min, float& currentCutoff, Functor& functor){

        cout << "Remove Abundance" << endl;
        _debugOuputFile = ofstream("/home/gats/workspace/tmp/lala1.txt");
        cout << "Checksum ab: " << _unitigGraph->computeChecksum_abundance() << endl;
        /*
        std::set<UnitigGraph::Node*, TipComparatorLala> _queueLala;
        for(UnitigGraph::Node* node : _validNodes){
            _queueLala.insert(node);
        }  

        while(!_queueLala.empty()){

            auto top = _queueLala.begin();
            UnitigGraph::Node* node = *top;
            _queueLala.erase(top);

            
            if(node->startNode() == 11838786){
                int abSum = 0;
                for(float ab : node->_abundances){
                    abSum += ab;
                }
                cout << abSum << endl;
                cout << node->startNode() << " " << node->_abundance << " " << node->_nodes.size() << " " << node->_abundances.size() << endl;
                for(u_int32_t nodeIndex : node->_nodes){
                    cout << "\t" << nodeIndex << " " << _unitigGraph->_nodeDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]._abundance << endl;
                }
                getchar();
            }
            
            
            //cout << node->startNode() << " " << node->_abundance << endl;
            //for(u_int32_t nodeIndex : node->_nodes){
            //    cout << "\t" << nodeIndex << " " << _unitigGraph->_nodeDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]._abundance << endl;
            //}
            //getchar();
            if(node->isCircular()) continue;
            //_debugOuputFile << node->startNode() << " " << node->_nodes.size() << " " << node->_abundance << endl;
        }
        */
        
        _debugOuputFile.close();
        //exit(1);

        u_int64_t nbErrorsRemoved = 0;
        float localCutoffRate = 0.5;
        float cutoffGlobal = 1;
        float aplha = 0.1;
        float t=1.1;
        currentCutoff = min(t, abundanceCutoff_min);



        while(t < abundanceCutoff_min){ 
            
            //if(saveAllState) checkSaveState(currentCutoff, unitigDatas, detectRoundabout, maxBubbleLength, saveState, insertBubble, doesSaveUnitigGraph, mdbg, minimizerSize, nbCores);

            currentCutoff = t;
            //unordered_set<u_int32_t> removedUnitigs;

            //cout << t << " " << abundanceCutoff_min << endl;

            /*
            for(UnitigGraph::Node* node : nodes){

                if(node->_abundance >= t) continue;

                vector<UnitigGraph::Node*> preds = node->_predecessors;
                _unitigGraph->removeNode(node);
                nbErrorsRemoved += 1;


                for(UnitigGraph::Node* predecessor : preds){
                    
                    if(predecessor->_unitigIndex == -1) continue;
                    recompactNodes.insert(predecessor);
                    recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(predecessor));
                }
            }
            */


            //compact(true, unitigDatas);

            
            _queue = std::set<UnitigGraph::Node*, AbundanceComparator>();
            //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , AbundanceComparator>();
            //for(UnitigGraph::Node* node : _unitigGraph->_nodes){
            //    if(node->_unitigIndex == -1) continue;
            //    if(node->)
            //    _queue.push(node);
            //}
            /*
            for(size_t i=0; i<_unitigGraph->_nodes.size(); i+=2){

                UnitigGraph::Node* node1 = _unitigGraph->_nodes[i];
                UnitigGraph::Node* node2 = _unitigGraph->_nodes[i+1];
                UnitigGraph::Node* node = nullptr;

                if(node1->_unitigIndex == -1) continue;
                if(node2->_unitigIndex == -1) continue;

                if(node1->startNode() < node2->startNode()){
                    node = node1;
                }
                else{
                    node = node2;
                }

                if(node->_abundance >= t) continue;
            */
            unordered_set<UnitigGraph::Node*> recompactNodes;


            for(size_t i=0; i<_validNodes.size(); i+=2){

                UnitigGraph::Node* node1 = _validNodes[i];
                UnitigGraph::Node* node2 = _validNodes[i+1];
                UnitigGraph::Node* node = nullptr;

                if(node1->_unitigIndex == -1) continue;
                if(node2->_unitigIndex == -1) continue;

                if(node1->startNode() < node2->startNode()){
                    node = node1;
                }
                else{
                    node = node2;
                }

                if(node->_abundance >= t) continue;

                vector<UnitigGraph::Node*> succs = node->_successors;
                vector<UnitigGraph::Node*> preds = node->_predecessors;
                _unitigGraph->removeNode(node);

                nbErrorsRemoved += 1;


                for(UnitigGraph::Node* predecessor : preds){
                    
                    if(predecessor->_unitigIndex == -1) continue;
                    //if(predecessor->_abundance < t) continue;

                    recompactNodes.insert(predecessor);
                    //recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(predecessor));

                }


                for(UnitigGraph::Node* successor : succs){
                    
                    if(successor->_unitigIndex == -1) continue;
                    //if(predecessor->_abundance < t) continue;

                    //recompactNodes.insert(successor);
                    recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(successor));

                }
            }


            

            for(UnitigGraph::Node* node : recompactNodes){

                if(node->_unitigIndex == -1) continue;
        
                _unitigGraph->recompact(node);
            }
            

            t = t * (1+aplha);

            if(nbErrorsRemoved > 0){
                cout << "Nb errors removed: " << nbErrorsRemoved << endl;

                int nbUnitigs = 0;
                for(size_t i=0; i<_unitigGraph->_nodes.size(); i++){

                    UnitigGraph::Node* node = _unitigGraph->_nodes[i];
                    if(node->_unitigIndex == -1) continue;
                    nbUnitigs += 1;
                }

                cout << "Nb unitigs: " << nbUnitigs << endl;
                //getchar();
            }

            //if(!_removeBubble){
                if(isCutoffProcessed.find(t) == isCutoffProcessed.end()){
                    isCutoffProcessed.insert(t);
                    functor(t);
                }
            //}

            if(nbErrorsRemoved > 0) break;
        }
        

        return nbErrorsRemoved;

    }

    //priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , AbundanceComparator> _queue;
    std::set<UnitigGraph::Node*, AbundanceComparator> _queue;

    struct TipComparatorLala {
        bool operator()(UnitigGraph::Node* p1, UnitigGraph::Node* p2) const{
            return p1->startNode() > p2->startNode();
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

            /*
            for(UnitigGraph::Node* node : nodes){

                if(node->_abundance >= t) continue;

                vector<UnitigGraph::Node*> preds = node->_predecessors;
                _unitigGraph->removeNode(node);
                nbErrorsRemoved += 1;


                for(UnitigGraph::Node* predecessor : preds){
                    
                    if(predecessor->_unitigIndex == -1) continue;
                    recompactNodes.insert(predecessor);
                    recompactNodes.insert(_unitigGraph->unitigIndex_toReverseDirection(predecessor));
                }
            }
            */


            //compact(true, unitigDatas);

            
            _queue = std::set<UnitigGraph::Node*, AbundanceComparator>();
            //_queue = priority_queue<UnitigGraph::Node*, vector<UnitigGraph::Node*> , AbundanceComparator>();
            //for(UnitigGraph::Node* node : _unitigGraph->_nodes){
            //    if(node->_unitigIndex == -1) continue;
            //    if(node->)
            //    _queue.push(node);
            //}
            /*
            for(size_t i=0; i<_unitigGraph->_nodes.size(); i+=2){

                UnitigGraph::Node* node1 = _unitigGraph->_nodes[i];
                UnitigGraph::Node* node2 = _unitigGraph->_nodes[i+1];
                UnitigGraph::Node* node = nullptr;

                if(node1->_unitigIndex == -1) continue;
                if(node2->_unitigIndex == -1) continue;

                if(node1->startNode() < node2->startNode()){
                    node = node1;
                }
                else{
                    node = node2;
                }

                if(node->_abundance >= t) continue;
            */
            for(size_t i=0; i<_validNodes.size(); i+=2){

                UnitigGraph::Node* node1 = _validNodes[i];
                UnitigGraph::Node* node2 = _validNodes[i+1];
                UnitigGraph::Node* node = nullptr;

                if(node1->_unitigIndex == -1) continue;
                if(node2->_unitigIndex == -1) continue;

                if(node1->startNode() < node2->startNode()){
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

                            if(predecessor->startNode() < predecessor_rc->startNode()){
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
        /*
        for(UnitigGraph::Node* node : _unitigGraph->_nodes){
            if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;
            _unitigGraph->recompact(node);
        }*/
        
        

        return nbErrorsRemoved;

    }

    u_int64_t _checksumAbundanceTotal;
    u_int64_t _checksumNodeTotal;

    void dumpUnitigs(float currentCutoff){

        u_int64_t checksumTotal = 0;
        u_int64_t checksumAbundanceTotal = 0;
        _cutoffIndexes.push_back({_cutoffIndex, currentCutoff});

        //cout << "dump unitigs: " << currentCutoff << endl;
        ofstream outputFile(_tmpDir + "/filter/unitigs_" + to_string(_cutoffIndex) + ".bin");

        //for(UnitigGraph::Node* node : _unitigGraph->_nodes){
        for(size_t i=0; i<_unitigGraph->_nodes.size(); i+=2){

            UnitigGraph::Node* node1 = _unitigGraph->_nodes[i];
            UnitigGraph::Node* node2 = _unitigGraph->_nodes[i+1];
            UnitigGraph::Node* node = nullptr;

            if(node1->startNode() < node2->startNode()){
                node = node1;
            }
            else{
                node = node2;
            }

            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;
            //if(node->_unitigIndex % 2 == 1) continue;
            if(node->_successors.size() == 0 && node->_predecessors.size() == 0 && node->_abundance == 1) continue;

            u_int32_t size = node->_nodes.size();
            outputFile.write((const char*)&size, sizeof(size));
            

            bool isCircular = node->isCircular();
            //if(isCircular) continue;
            //isCircular = false;
            outputFile.write((const char*)&isCircular, sizeof(isCircular));
            
            bool isRepeatSide_ = isRepeatSide(node);
            outputFile.write((const char*)&isRepeatSide_, sizeof(isRepeatSide_));

            //if(node->_nodes.size() > 1000){
            //    cout << isCircular << endl;
            //    getchar();
            //}

            float abundance = node->_abundance;
            outputFile.write((const char*)&abundance, sizeof(abundance));

            outputFile.write((const char*)&node->_nodes[0], size * sizeof(u_int32_t));
            
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
        }


        outputFile.close();

        _cutoffIndex += 1;
    }

	bool isRepeatSide(UnitigGraph::Node* node){

		//if(node->_length > 50000) return false;
        if(node->_nodes.size() > _kminmerSize*2) return false;
		if(node->_successors.size() == 0) return false;
		if(node->_predecessors.size() == 0) return false;

		for(UnitigGraph::Node* nodeS : node->_successors){
			if(nodeS == node) continue;

			for(UnitigGraph::Node* nodeP : node->_predecessors){
				if(nodeP == node) continue;

				if(nodeS == nodeP){
					return true;
				}
			}
		}

		return false;
	}

};



#endif
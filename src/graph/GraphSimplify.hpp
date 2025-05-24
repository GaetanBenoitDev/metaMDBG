


/*
#ifndef MDBG_METAG_GRAPHSIMPLIFY
#define MDBG_METAG_GRAPHSIMPLIFY


#include "Commons.hpp"



class GraphSimplify{

public:


    BiGraph* _graphSuccessors;
    //BiGraph* _graphPredecessors;
    string _outputDir;
   
    vector<bool> _isNodeIndexIndexed;
    u_int32_t _nbNodes;
    int _nbCores;
    
    unordered_map<u_int32_t, u_int32_t> _endNode_to_unitigIndex;
    unordered_map<u_int32_t, u_int32_t> _startNode_to_unitigIndex;

    GraphSimplify(BiGraph* graphSuccessors, const string& outputDir, int nbCores){

        _graphSuccessors = graphSuccessors;
        _outputDir = outputDir;
        _nbCores = nbCores;
        
        _nbNodes = _graphSuccessors->_nbNodes;

    }



    ~GraphSimplify(){
    }


    u_int32_t _nextUnitigIndex;


    void addUnitigGraphEdges(){

        Logger::get().debug() << "Add unitig graph edges";

        //cout << "a " << _startNode_to_unitigIndex.size() << endl;
        //_logFile << "Unitig graph nb nodes: " << _unitigGraph->_nodes.size() << endl;
       

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


    void compact(){


        Logger::get().debug() << "\tCompacting";
        

        _nextUnitigIndex = 0;


        
        UnitigFunctor functor(this);
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

        Logger::get().debug() << "\tdone";
        
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

	class UnitigFunctor {

		public:

		GraphSimplify* _graph;

		UnitigFunctor(GraphSimplify* graph) : _graph(graph){
		}

		UnitigFunctor(const UnitigFunctor& copy) : _graph(copy._graph){
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
    

	ofstream _unitigGraphFile_nodes;
	ofstream _unitigGraphFile_edges_successors;
	ofstream _unitigGraphFile_edges_predecessors;

    void debug_writeGfaErrorfree(){

        _unitigGraphFile_nodes = ofstream(_outputDir + "/unitigGraph.nodes.bin");
        _unitigGraphFile_edges_successors = ofstream(_outputDir + "/unitigGraph.edges.successors.bin");
        _unitigGraphFile_edges_predecessors = ofstream(_outputDir + "/unitigGraph.edges.predecessors.bin");

        _nextUnitigIndex = 0;
        _isNodeIndexIndexed = vector<bool>(_graphSuccessors->_nbNodes, false);

        compact();

        _unitigGraphFile_nodes.close();
        _unitigGraphFile_edges_successors.close();
        _unitigGraphFile_edges_predecessors.close();

    }




    

};


#endif
*/
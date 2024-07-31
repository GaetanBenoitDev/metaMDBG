

#ifndef MDBG_METAG_UNITIGEXPANDER
#define MDBG_METAG_UNITIGEXPANDER

#include "Commons.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"

struct SuccessorData{
	u_int32_t _nodeIndex;
	u_int32_t _abundance;
	u_int32_t _sourceAbundance;
	int _prevRank;
	int _prevUnitigRank;
	bool _prevRankFinished;
	u_int32_t _nodeIndexSuccessor;
	u_int32_t _distanceFromSource;
	u_int32_t _prevUnitigIndex;
	int _nbSharedReads;
	vector<u_int32_t> _path;
	u_int32_t _backtrackPathLength;
	u_int32_t _backtrackPathLength_noFilter;
	//int _nbSharedReadsPrev;
	//vector<u_int32_t> _processedNodeIndex;
};

/*
class PathExplorer{

public:

	GraphSimplify* _graph;
	vector<u_int32_t> _prevNodes;
	vector<UnitigData>& _unitigDatas;
	bool _forward;
	unordered_set<u_int32_t>& _visitedNodes;

	PathExplorer(GraphSimplify* graph, vector<u_int32_t>& prevNodes, const vector<u_int32_t>& unitigDatas, bool forward, unordered_set<u_int32_t>& visitedNodes) : _prevNodes(prevNodes), _unitigDatas(unitigDatas): _visitedNodes(visitedNodes){
		_graph = graph;
		_forward = forward;
	}



	//u_int32_t _index;
	//unordered_set<DbgEdge, hash_pair> isEdgeVisited;
	vector<u_int32_t> _prevNodes;
	u_int32_t _source_abundance;
	u_int32_t _source_nodeIndex;
	u_int32_t _start_nodeIndex;
	float _currentAbundance;
	unordered_set<u_int32_t>& _visitedNodes;
	//unordered_map<u_int32_t, bool>& _isNodeImproved;
	vector<UnitigData>& _unitigDatas;

	vector<u_int32_t> _exploredNodes;
	//unordered_set<u_int32_t> _isPathAlreadyVisitedSourceNodes;
	u_int64_t _currentPathLength;
	//unordered_set<u_int32_t>& _solvedUnitigs;
	bool _canCheckScc;
	u_int32_t _cutoff_backtrackLength;
	AssemblyState& _assemblyState;
	GraphSimplify* _graph;

	PathExplorer(const vector<u_int32_t>& prevNodes, u_int32_t source_abundance, u_int32_t source_nodeIndex, u_int32_t start_nodeIndex, float currentAbundance, unordered_set<u_int32_t>& visitedNodes, vector<UnitigData>& unitigDatas, u_int64_t currentPathLength, bool canCheckScc, AssemblyState& assemblyState) : _visitedNodes(visitedNodes), _unitigDatas(unitigDatas), _assemblyState(assemblyState){
		_prevNodes = prevNodes;
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_start_nodeIndex = start_nodeIndex;
		//_abundanceCutoff_min = abundanceCutoff_min;
		_currentAbundance = currentAbundance;
		//_isPathAlreadyVisitedSourceNodes = isPathAlreadyVisitedSourceNodes;
		_currentPathLength = currentPathLength;
		_canCheckScc = canCheckScc;
	}

	bool containsSuccessor(const SuccessorData& successor, const vector<SuccessorData>& successorData){
		for(const SuccessorData& s : successorData){
			if(s._nodeIndex == successor._nodeIndex) return true;
		}
		return false;
	}

	u_int32_t getNextNode(u_int32_t current_nodeIndex, GraphSimplify* graph, bool forward, u_int32_t currentDepth, u_int32_t& resultType, vector<SuccessorData>& nextNodes, bool usePathSuccessors, bool canCheckTips){

		_graph = graph;

		nextNodes.clear();
		resultType = 0;

		bool orient_dummy = false;


		u_int32_t current_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
		u_int32_t current_abundance = graph->_nodeAbundances[current_nodeName]; //_unitigDatas[current_unitigIndex]._meanAbundance;


		vector<SuccessorData> data_successors;

		
		if(usePathSuccessors){
			if(isInInfiniteCycle(current_nodeIndex, graph, _unitigDatas, forward)){
				return -1; //Infinite simple path without exit
			}
		}


		vector<u_int32_t> successors;
		if(forward){
			graph->getSuccessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
		}
		else{
			graph->getPredecessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
		}

		//vector<u_int32_t> successors;
		//collectSolvedSuccessors(current_nodeIndex, directSuccessors, successors);
		//cout << directSuccessors.size() << " " << successors.size() << endl;

		if(successors.size() == 1 || !usePathSuccessors){
			


			for(u_int32_t utg_n : successors){
				//if(!includeVisitedSuccessors && (_visitedNodes.find(utg_n) != _visitedNodes.end())) continue;

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {utg_n, successor_abundance, 0, 0, 0, false, utg_n, 0, 0, -1, {}, 0};
				//successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				successor._sourceAbundance = abs(successor._abundance - _currentAbundance);
				data_successors.push_back(successor);

			}
			




		}
		else{



			unordered_map<u_int32_t, DataSuccessorPath> successorPaths;

			if(successors.size() == 0 && _assemblyState._currentPathIndex != -1){
				_assemblyState._currentPathIndex = -1;
			}

			//cout << "hisdfsdf: " << _assemblyState._currentPathIndex << endl;
			

			collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, true, -1, _visitedNodes);



			for(auto& it : successorPaths){

				u_int32_t nodeIndex = it.first;
				DataSuccessorPath& pathData = it.second;

				int distance = 0;


				u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {nodeIndex, successor_abundance, 0, 0, 0, false, nodeIndex, distance, 0, -1, pathData._path, 0};
				//successor._sourceAbundance = abs((long)successor._abundance - (long)_source_abundance);
				successor._sourceAbundance = abs(successor._abundance - _currentAbundance);
				data_successors.push_back(successor);

			}
			

			if(currentDepth == 0){
				if(successors.size() > 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "\t----------- " << endl;
					cout << "\tNb prev nodes: " << _prevNodes.size() << endl;
					cout << "\tNb possible successors: " << successorPaths.size() << endl;
					cout << "\tLocal abundance: " << graph->getNodeUnitigAbundance(current_nodeIndex) << endl;
 					//for(u_int32_t utg_n : successors){
					for(SuccessorData& successor : data_successors){
						u_int32_t utg_n = successor._nodeIndex;
						//for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "\t" << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << graph->getNodeUnitigAbundance(utg_n) << " " << Utils::computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
					
					}
				}
			}


		}

	
		






		if(data_successors.size() == 0){
			if(currentDepth == 0){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "No successors" << endl;
			}
			return -1;
		}
		else if(data_successors.size() == 1){


			if(data_successors[0]._path.size() == 0){ //Si useSuccessorPath, le path est deja complet
				data_successors[0]._path.push_back(current_nodeIndex);
				current_nodeIndex = data_successors[0]._nodeIndexSuccessor;
				data_successors[0]._path.push_back(current_nodeIndex);
			}

			updateNodeChosen(current_nodeIndex, currentDepth);
			nextNodes.push_back(data_successors[0]);
			return current_nodeIndex;
		}
		else{


			resultType = 2;

			

			if(usePathSuccessors){
				//cout << "allo " << data_successors.size() << endl;
				vector<SuccessorData> dummy;
				computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, 0, dummy, 0, true, true, forward, usePathSuccessors); //ATTENTION NE PAS ENLEVER CA, on s'en sert pour calculer la _backtrackLength_noFilter (juste desactiver print_debug si besoin)
			}

			//removeUnsupportingReads(graph, _unitigDatas, data_successors, forward);
			
			SuccessorData bestSuccessorIndex;
			vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
			if(usePathSuccessors){
				computeBestSuccessors_byUnitigRank_all2(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
				
				if(_assemblyState._cutoff_backtrackLength == 0){	
					std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byNbSharedReads);
					bestSuccessorIndex = data_successors[0];

				}


			}
			else{
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, false, usePathSuccessors);
				//cout << successors_bestPrevUnitigRank_noCutoff.size() << endl;
			}

			if(usePathSuccessors){
			

			}

			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

				//if(currentDepth == 0){
					//cout << "\tNode chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				//}

				//updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
				if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
				
				return current_nodeIndex;
			}




			if(usePathSuccessors){
				//exit(1);
				//cout << "lalala " << successors_bestPrevUnitigRank_noCutoff.size() << endl;
				if(successors_bestPrevUnitigRank_noCutoff.size() >= 2){
					bool isBubble = true;
					for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){
						if(!graph->_isBubble[successor._nodeIndexSuccessor]){
							isBubble = false;
						}
					}
					cout << "\tIs bubble ?: " <<  isBubble << endl;
					if(isBubble){
						
						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byPathLength_noFilter);

						//if(successors_bestPrevUnitigRank_noCutoff[0]._backtrackPathLength == successors_bestPrevUnitigRank_noCutoff[1]._backtrackPathLength){
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//}

						//if(currentDepth == 0){
						cout << "\t\tTake bubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
						//}
						//exit(1);
						//updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
						//if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
						return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;
					}

					
					//cout << currentDepth << " " << "detect" << endl;
					u_int32_t superbubbleSink = detectSuperbubble(current_nodeIndex, maxLength, graph, _unitigDatas, forward);
					if(superbubbleSink != -1){

						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byPathLength_noFilter);

						//if(successors_bestPrevUnitigRank_noCutoff[0]._backtrackPathLength == successors_bestPrevUnitigRank_noCutoff[1]._backtrackPathLength){
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//}

						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//if(currentDepth == 0){
							cout << "\t\tTake superbubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
						//getchar();
						//}
						//exit(1);
						updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
						return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

					}
					//cout << "lala: " << BiGraph::nodeIndex_to_nodeName(superbubbleSink) << endl;
				}

				//}

				if(successors_bestPrevUnitigRank_noCutoff.size() < 5){

				
				
					//if(currentDepth == 0){

					if(currentDepth == 0){
						for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "\tCheck simple cycle" << endl;
					}

					for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){

						if(isSmallCycle(successor._nodeIndexSuccessor, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
							updateNodeChosen(successor._nodeIndexSuccessor, currentDepth);
							nextNodes.push_back(successor);
							if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
							return successor._nodeIndexSuccessor;
						}
						
					}



					for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
						u_int32_t to_nodeIndex = successors_bestPrevUnitigRank_noCutoff[i]._nodeIndexSuccessor;
						for(size_t j=0; j<successors_bestPrevUnitigRank_noCutoff.size(); j++){
							if(i == j) continue;
							
							u_int32_t from_nodeIndex = successors_bestPrevUnitigRank_noCutoff[j]._nodeIndexSuccessor;

							if(currentDepth == 0){
								for(size_t i=0; i<currentDepth; i++) cout << "  ";
								cout << "Check simple cycle from: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << "    to:    " << graph->_graphSuccessors->nodeToString(to_nodeIndex) << endl;
							}

							if(isSmallCycle(from_nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
								//exit(1);
								updateNodeChosen(from_nodeIndex, currentDepth);
								nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[j]);
								if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
								return from_nodeIndex;
							}

						}	
					}
				}
			}
			
			for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[i]);
			}
			//if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);


			//if(!usePathSuccessors && currentDepth == 0){
			//	return getNextNode(current_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, true);
			//}

			return -1;
				
			//}

		}


	}

	
	bool isInInfiniteCycle(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		//cout << "Check infinite cycle" << endl;
		Unitig& unitig = graph->nodeIndex_to_unitig(current_nodeIndex);
		if(unitig._startNode != current_nodeIndex) return false;

		unordered_map<u_int32_t, DataSuccessorPath> successorPaths;
		collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, true, -1, _visitedNodes);

		//cout << "Check infinite cycle: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " <<  successorPaths.size() << endl;
		if(successorPaths.size() == 0){
			cout << "Infinite cycle" << endl;
			//getchar();
		} 

		return successorPaths.size() == 0;
	}

	void removeUnsupportingReads(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors, bool forward){
		
		graph->_removedReadIndex.clear();
		
		for(SuccessorData& successorData : data_successors){
			
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex);

			for(size_t i=0; i<_prevNodes.size()-1; i++){

				u_int32_t prevNodeIndex = _prevNodes[i];

				vector<u_int32_t> predecessors;
				if(forward){
					graph->getPredecessors(prevNodeIndex, 0, predecessors);
				}
				else{
					graph->getSuccessors(prevNodeIndex, 0, predecessors);
				}

				for(u_int32_t predecessor : predecessors){

					vector<u_int64_t> sharedReads;
					Utils::collectSharedReads(_unitigDatas[nodeName], _unitigDatas[BiGraph::nodeIndex_to_nodeName(predecessor)], sharedReads);	
					//if(sharedReads.size() == 0) continue;

					//cout << distance << " " << sharedReads.size() << endl;
					if(std::find(_prevNodes.begin(), _prevNodes.end(), predecessor) == _prevNodes.end()){
						for(u_int64_t readIndex : sharedReads){
							//if(graph->_removedReadIndex.find(readIndex) == graph->_removedReadIndex.end()){
								//cout << BiGraph::nodeIndex_to_nodeName(predecessor._index) << endl;
								graph->_removedReadIndex.insert(readIndex);
							//}
							//cout << readIndex << endl;
						}
					}
				}

			}
		}
		
	}

    u_int32_t detectSuperbubble(u_int32_t nodeIndex_source, u_int64_t maxLength, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		cout << "\tStart detection: " << BiGraph::nodeIndex_to_nodeName(nodeIndex_source) << endl;



        unordered_set<u_int32_t> isVisited;
        unordered_set<u_int32_t> seen;
        unordered_map<u_int32_t, u_int64_t> pathLength;
        vector<u_int32_t> queue;

        queue.push_back(nodeIndex_source);
        pathLength[nodeIndex_source] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

			//cout << "Visit: " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << " " << pathLength[v] << endl;
            //cout << "\tVisited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._endNode) << endl;
            queue.pop_back();

			//cout << "Visit 1 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

            if(pathLength[v] > 10000) continue;

			//cout << "Visit 2 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

			if(isVisited.find(v) != isVisited.end()){
				return -1; // ?? quick patch
			}

            isVisited.insert(v);
            if(seen.find(v) != seen.end()){
                seen.erase(v);
            }

            
			
			//u_int32_t current_nodeIndex = graph->_unitigs[v]._endNode;
			//PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, 0, _solvedUnitigs);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);
			
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			getNextNode(v, graph, forward, 0, resultType, nextNodes, false, true);

			//cout << "\tNb succ: " << nextNodes.size() << endl;

            //vector<u_int32_t> successors;
			//getSuccessors_unitig(v, successors);
            if(nextNodes.size() == 0) return -1; //abort tip

			for(SuccessorData& successorData : nextNodes){
				u_int32_t u = successorData._nodeIndex;
				
				//cout << "\tSucc : " << BiGraph::nodeIndex_to_nodeName(u) << " " << queue.size() << endl;

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(v, u);
				//bool isReachable = pathExplorer.isReachable2(successorData._nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)]-overlapLength), allowReverseDirection);
				//if(isReachable) return true;
				//}

            	//for(u_int32_t u : successors){
                if(u == nodeIndex_source) return -1; //cycle including s

                if(isVisited.find(u) == isVisited.end()){

                    	seen.insert(u);
                    	pathLength[u] = pathLength[v] + (graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(u)]-overlapLength);
					



                }

            }

            //for(u_int32_t u : successors){
			for(SuccessorData& successorData : nextNodes){
				u_int32_t u = successorData._nodeIndex;
                //cout << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

				u_int32_t resultType;
				vector<SuccessorData> prevNodes;
				getNextNode(u, graph, !forward, 0, resultType, prevNodes, false, true);
				//getPredecessors_unitig(u, predecessors);
                bool allPredecessorsAreVisited = true;
				for(SuccessorData& predecessorData : prevNodes){
					u_int32_t p = predecessorData._nodeIndex;
					//cout << "\tPred:" << BiGraph::nodeIndex_to_nodeName(p) << endl;
                	//for(u_int32_t p : predecessors){
                    if(isVisited.find(p) == isVisited.end()){
                        allPredecessorsAreVisited = false;
                        break;
                    }
                }

                if(allPredecessorsAreVisited){
					//cout << "\tPush: " << BiGraph::nodeIndex_to_nodeName(u) << endl;
                    //cout << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                    queue.push_back(u);
                }

                //cout << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                    u_int32_t t = queue[0];
					
					vector<SuccessorData> successors_t;
					getNextNode(t, graph, forward, 0, resultType, successors_t, false, true);

					bool foundSource = false;
					for(SuccessorData& successorData_t : successors_t){
						if(successorData_t._nodeIndex == nodeIndex_source){
							foundSource = true;
							break;
						}
					}

                    //vector<u_int32_t> successors_t;
                    //getSuccessors_unitig(t, successors_t);

                    if(!foundSource){
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

	bool isTip(u_int32_t source_nodeIndex, u_int64_t maxTipLength, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		//cout << "-----" << endl;
		u_int32_t targetNodeIndex = -1;
		bool print_debug = false;
		unordered_set<u_int32_t> visitedNodes = {};


		u_int32_t minNbJokers = -1;
		if(print_debug) cout << "\tCollect possible successors: " << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		u_int32_t currentDepth = 0;

		priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator> queue;
		queue.push({source_nodeIndex, UINT32_MAX, {}, 0, 0, _prevNodes, {}, {}, {}, _currentAbundance});


		while(true){
			
			if(queue.size() == 0) break;

			//cout << currentDepth << " " << queue.size() << endl;
			if(currentDepth > 100){
				//successorPaths.clear();
				return false;
			}

			DataSuccessorPath dataSuccessorPath = queue.top();
			float currentAbundance = dataSuccessorPath._currentAbundance;
        	queue.pop();

			

			//if(dataSuccessorPath._nbJokers > minNbJokers) continue;

			u_int32_t current_nodeIndex = dataSuccessorPath._currentNodeIndex; //dataSuccessorPath._path[dataSuccessorPath._path.size()-1];
			//cout << queue.size() << " " <<  BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._pathLength << endl;
			if(visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) continue;
			
			//cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			//if(print_debug) cout << "\tStart: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << minNbJokers << " " << currentDepth << " " << (visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) << endl;
			//getchar();
			//u_int32_t currentLength = dataSuccessorPath._pathLength;

			PathExplorer pathExplorer(dataSuccessorPath._prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, currentAbundance, visitedNodes, _unitigDatas, dataSuccessorPath._pathLength, false, _assemblyState);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);

			bool continueVisiting = isTip_visitSuccessor(current_nodeIndex, dataSuccessorPath._prevNodeIndex, dataSuccessorPath, pathExplorer, graph, minNbJokers, print_debug, targetNodeIndex, currentAbundance, maxTipLength);
			if(!continueVisiting) return false;
			visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));

			
			while(true){
				
				//cout << "miam" << endl;
				if(currentDepth > 100){
					return false;
				}

				vector<u_int32_t> allSuccessors;

				if(forward){
					graph->getSuccessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				else{
					graph->getPredecessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				
				//cout << currentDepth.size() << " " << allSuccessors.size() << endl;

				//cout << allSuccessors.size() << endl;
				if(allSuccessors.size() == 0){
					break;
				}
				else if(allSuccessors.size() == 1){
					
					u_int32_t prevNodeIndex = current_nodeIndex;
					current_nodeIndex = allSuccessors[0];
					if(visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) break;
					bool continueVisiting = isTip_visitSuccessor(current_nodeIndex, prevNodeIndex, dataSuccessorPath, pathExplorer, graph, minNbJokers, print_debug, targetNodeIndex, currentAbundance, maxTipLength);
					if(!continueVisiting) return false;
					visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));

					//cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._pathLength << endl;
				}
				else{

					u_int32_t resultType;
					vector<SuccessorData> nextNodes;
					pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, false);

					//cout << nextNodes.size() << endl;
					//cout << "\t\tNb valid successors: " << nextNodes.size() << endl;

					for(SuccessorData& successorData : nextNodes){
						//cout << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorData._nodeIndex, current_nodeIndex, queue, dataSuccessorPath, false, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					for(u_int32_t successorNodeIndex : allSuccessors){
						bool exist = false;
						for(SuccessorData& s : nextNodes){
							if(s._nodeIndex == successorNodeIndex){
								exist = true;
								break;
							}
						}
						if(exist) continue;

						//cout << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorNodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd invalid successor: " << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorNodeIndex, current_nodeIndex, queue, dataSuccessorPath, true, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					break;

				}


			}

		}

		return true;
	}
	
	bool isTip_visitSuccessor(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, DataSuccessorPath& dataSuccessorPath, PathExplorer& pathExplorer, GraphSimplify* graph, u_int32_t& minNbJokers, bool print_debug, u_int32_t targetNodeIndex, float& currentAbundance, u_int32_t maxTipLength){
		
		//cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;

		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			

		dataSuccessorPath._path.push_back(nodeIndex);


		//cout << dataSuccessorPath._pathLength << " " << maxTipLength << endl;
		if(dataSuccessorPath._pathLength > maxTipLength){//|| currentDepth > 50){
			//cout << "\t\t\tStop max length" << endl; 
			return false;
		}

		if(nodeIndexPrev != -1){ //First iteration
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
			dataSuccessorPath._pathLength += (graph->_nodeLengths[nodeName] - overlapLength);
		}
		pathExplorer.nodeExplored(nodeIndex, graph);
		//currentAbundance = PathExplorer::updateCurrentAbundance(nodeIndex, currentAbundance, graph, _assemblyState, _kminmerSize);

		return true;
	}

	unordered_set<u_int32_t> __beatenNodeIndex;

	void computeBestSuccessors_byUnitigRank_all2(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors2, int currentDepth, bool forward, vector<SuccessorData>& validSuccessors, bool includeVisitedSuccessors, bool usePathSuccessors){
	
		//cout << "1" << endl;
		__beatenNodeIndex.clear();



		for(size_t i=0; i<data_successors2.size(); i++){
			for(size_t j=i+1; j<data_successors2.size(); j++){

				if(__beatenNodeIndex.find(data_successors2[i]._nodeIndex) != __beatenNodeIndex.end()) continue;
				if(__beatenNodeIndex.find(data_successors2[j]._nodeIndex) != __beatenNodeIndex.end()) continue;


				//cout << i << " " << j << " " << _assemblyState._cutoff_backtrackLength << " " << BiGraph::nodeIndex_to_nodeName(data_successors2[i]._nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(data_successors2[j]._nodeIndex) << endl;

				//cout << i << " vs " << j << endl;
				vector<SuccessorData> data_successors;
				data_successors.push_back(data_successors2[i]);
				data_successors.push_back(data_successors2[j]);


				//cout << data_successors[0]._path.size() << " " << data_successors[1]._path.size() << endl;
				if(_assemblyState._cutoff_backtrackLength > 0){

					unordered_set<u_int32_t> ind;
					for(u_int32_t nodeIndex : data_successors[0]._path){
						ind.insert(nodeIndex);
					}

					u_int32_t sharedNodexIndex = -1;
					//u_int32_t lastSharedRank = 0;
					//size_t index=0;
					for(u_int32_t nodeIndex : data_successors[1]._path){
						if(ind.find(nodeIndex) != ind.end() || find(_prevNodes.begin(), _prevNodes.end(), nodeIndex) != _prevNodes.end()){
							//lastSharedRank = index;
							//cout << "shared rank: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
							sharedNodexIndex = nodeIndex;
						}
						//index += 1;
					}

					//Find LAST index of sharedNodexIndex (this node can be repeated in the path)
					u_int32_t index1 = 0;
					for(size_t i=0; i<data_successors[0]._path.size(); i++){
						if(data_successors[0]._path[i] == sharedNodexIndex){
							index1 = i;
						}
					}
					u_int32_t index2 = 0;
					for(size_t i=0; i<data_successors[1]._path.size(); i++){
						if(data_successors[1]._path[i] == sharedNodexIndex){
							index2 = i;
						}
					}

					//cout << "Shared node: " << sharedNodexIndex << endl;

					u_int32_t rank1 = index1 + 1;
					u_int32_t rank2 = index2 + 1;


					data_successors[0]._nodeIndexSuccessor = data_successors[0]._path[rank1];
					data_successors[1]._nodeIndexSuccessor = data_successors[1]._path[rank2];

					//cout << "NodeIndexSuccessor " << BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndex) << ": " << BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndexSuccessor) << endl;
					//cout << "NodeIndexSuccessor " << BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndex) << ": " << BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndexSuccessor) << endl;
				}



				u_int32_t prev_nodeIndex = _prevNodes[_prevNodes.size()-1];
				u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);

				u_int32_t nodeName1 = BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndexSuccessor);
				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndexSuccessor);

				vector<SuccessorData> bestSuccessors;
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, bestSuccessors, includeVisitedSuccessors, usePathSuccessors);
				

					
					if(bestSuccessors.size() == 1){
						cout << "\tWinner: " << bestSuccessors.size() << " " << BiGraph::nodeIndex_to_nodeName(bestSuccessors[0]._nodeIndex) << endl;
						if(bestSuccessors[0]._nodeIndex == data_successors[0]._nodeIndex){
							__beatenNodeIndex.insert(data_successors[1]._nodeIndex);
						}
						else{
							__beatenNodeIndex.insert(data_successors[0]._nodeIndex);
						}
					}
				//}
				

			}
		}

 		validSuccessors.clear();
		for(SuccessorData& s : data_successors2){
			if(__beatenNodeIndex.find(s._nodeIndex) != __beatenNodeIndex.end()) continue;
			validSuccessors.push_back(s);
		}


		//cout << validSuccessors.size() << " " << __beatenNodeIndex.size() << endl;
		if(validSuccessors.size() == 0) return;
		if(validSuccessors.size() > 10) return;
		//if(currentDepth == 0 && abundanceCutoff_min != 0){

		//cout << "3" << endl;
		unordered_set<u_int32_t> reachableSuccessors;
		
		if(_assemblyState._cutoffType == CutoffType::STRAIN_HIGH){
			
			//unordered_set<u_int32_t> reachableSuccessors;
			for(size_t i=0; i<validSuccessors.size(); i++){
				for(size_t j=0; j<validSuccessors.size(); j++){
					if(i == j) continue; 
					//if(reachableSuccessors.find(successors_bestPrevUnitigRank[j]._nodeIndex) != reachableSuccessors.end()) continue; //Already reached
					cout << "\tCheck reachable (" << currentDepth  << "): " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
					
					u_int32_t totalIter = 0;
					if(isReachable2(validSuccessors[i]._nodeIndex, validSuccessors[j]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, false, totalIter)){//} && isReachable2(validSuccessors[j]._nodeIndex, validSuccessors[i]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, false)){
						cout << "\t\tIs reachable: " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
						reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
					}
				}
			}
				
					
					
					
			for(size_t i=0; i<validSuccessors.size(); i++){
				for(size_t j=i+1; j<validSuccessors.size(); j++){
					//if(i == j) continue; 
					if(reachableSuccessors.find(validSuccessors[j]._nodeIndex) != reachableSuccessors.end()) continue; //Already reached
					
					u_int32_t totalIter = 0;
					if(isReachable2(validSuccessors[i]._nodeIndex, validSuccessors[j]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, true, totalIter) && isReachable2(validSuccessors[j]._nodeIndex, validSuccessors[i]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, true, totalIter)){
						cout << "\tIs reachable (RC): " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
						//cout << validSuccessors[i]._nbSharedReads << " " << validSuccessors[j]._nbSharedReads << endl;
						if(validSuccessors[i]._nbSharedReads > validSuccessors[j]._nbSharedReads){
							reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
						}
						else if(validSuccessors[j]._nbSharedReads > validSuccessors[i]._nbSharedReads){
							reachableSuccessors.insert(validSuccessors[i]._nodeIndex);
						}
						else{
							cout << "TODO egalite nb shared reads, on peut prendre en compte la distance au node source apr exemple" << endl;
								
							if(validSuccessors[i]._backtrackPathLength > validSuccessors[j]._backtrackPathLength){
								reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
							}
							else{
								reachableSuccessors.insert(validSuccessors[i]._nodeIndex);
							}

							//exit(1);
						}
					}
				}
			}
		}
		

		vector<SuccessorData> chosenSuccessors;
		for(SuccessorData& successorData : validSuccessors){
			if(reachableSuccessors.find(successorData._nodeIndex) != reachableSuccessors.end()) continue;
			chosenSuccessors.push_back(successorData);
			//cout << "lalala: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
		}

		if(chosenSuccessors.size() == 0 && reachableSuccessors.size() > 0){

			vector<SuccessorData> chosenSuccessors_tmp;
			for(SuccessorData& successorData : validSuccessors){
				chosenSuccessors_tmp.push_back(successorData);
			}

			sort(chosenSuccessors_tmp.begin(), chosenSuccessors_tmp.end(), SuccessorComparator_byDistanceFromSource);
			//sort(chosenSuccessors_tmp.begin(), chosenSuccessors_tmp.end(), SuccessorComparator_byNbSharedReads);
			chosenSuccessors.push_back(chosenSuccessors_tmp[0]);
		}

		validSuccessors.clear();
		for(SuccessorData& successorData : chosenSuccessors){
			validSuccessors.push_back(successorData);
		}



		//}

	}

	void computeBestSuccessors_byUnitigRank_all(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors2, int currentDepth, bool forward, vector<SuccessorData>& validSuccessors, bool includeVisitedSuccessors, bool usePathSuccessors){

		//cout << "Include visited successors: " << includeVisitedSuccessors << endl;
		//cout << currentDepth << " " << usePathSuccessors << endl;
		float cutoff_min = PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType) / 2.0; //max(1.0, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0)/4.0);
		cutoff_min = (int)(cutoff_min + 1); //max(cutoff_min, 1);
		//(!usePathSuccessors){
		//	cutoff_min = 1;
		//}

		//cout << "Cutoff min: " << cutoff_min << endl;

		
		vector<SuccessorData> data_successors;
		for(SuccessorData& s : data_successors2){
			//cout << "\t" << BiGraph::nodeIndex_to_nodeName(s._nodeIndex) << " " << (_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(s._nodeIndex)) != _visitedNodes.end()) << endl;
			if(!includeVisitedSuccessors && _visitedNodes.find(BiGraph::nodeIndex_to_nodeName(s._nodeIndex)) != _visitedNodes.end()) continue;
			data_successors.push_back(s);
			

			//if(usePathSuccessors){
			//	cout << "Possible succ: " << BiGraph::nodeIndex_to_nodeName(s._nodeIndex) << "   " << BiGraph::nodeIndex_to_nodeName(s._path[1]) << endl;
			//}
		}
		

		if(data_successors.size() == 0) return;
		//Print debug prev ranks


		vector<SuccessorData> successors_bestPrevUnitigRank_cutoff;
		//if(usePathSuccessors){
			//cout << cutoff_min << endl;
		computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_cutoff, cutoff_min, true, false, forward, usePathSuccessors);
		if(successors_bestPrevUnitigRank_cutoff.size() == 1){
			u_int32_t current_nodeIndex = successors_bestPrevUnitigRank_cutoff[0]._nodeIndexSuccessor;

			//if(currentDepth == 1){
			//cout << "\tNode chosen 1: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;// << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
			//}

			//updateNodeChosen(current_nodeIndex, currentDepth);
			//nextNodes.push_back(current_nodeIndex);
			//return current_nodeIndex;

			for(SuccessorData& s : validSuccessors){
				if(s._nodeIndex == current_nodeIndex) return;
			}

			//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << currentDepth << " " << (_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != _visitedNodes.end()) << endl; 
			validSuccessors.push_back(successors_bestPrevUnitigRank_cutoff[0]);
			return;
		}

		if(usePathSuccessors){
			for(SuccessorData& successorData : successors_bestPrevUnitigRank_cutoff){
				
				bool exist = false;
				for(SuccessorData& s : validSuccessors){
					if(s._nodeIndex == successorData._nodeIndex){
						exist = true;
						break;
					}
				}
				if(exist) continue;

				//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl; 
				validSuccessors.push_back(successorData);
			}
		}
		else{
			for(SuccessorData& successorData : successors_bestPrevUnitigRank_cutoff){
				
				bool exist = false;
				for(SuccessorData& s : validSuccessors){
					if(s._nodeIndex == successorData._nodeIndex){
						exist = true;
						break;
					}
				}
				if(exist) continue;

				//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl; 
				validSuccessors.push_back(successorData);
			}
		}
	}

	void computeBestSuccessors_byUnitigRank(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors, int currentDepth, vector<SuccessorData>& successors_bestPrevUnitigRank, float abundanceCutoff_min, bool useUnitigRank, bool print_debug, bool forward, bool usePathSuccessors){



		#ifdef PRINT_PREV_RANK
			//if(print_debug && usePathSuccessors){
			if(usePathSuccessors){
				cout << endl << "\tCutoff: " << abundanceCutoff_min << endl;
			}
		#endif
		#ifdef PRINT_PREV_RANK_ALL
		//if(print_debug && currentDepth == 0){
			cout << endl << "\tCutoff: " << abundanceCutoff_min << endl;
		//}
		#endif

		for(SuccessorData& successor : data_successors){
			successor._prevRankFinished = false;
			successor._prevUnitigRank = -1;
			successor._prevRank = -1;
			//successor._nbSharedReadsPrev = -1;
			successor._backtrackPathLength = 0;
			successor._backtrackPathLength_noFilter = 0;
		}

		successors_bestPrevUnitigRank.clear();

		u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

		u_int32_t prevRank = 0;
		u_int32_t prevRank_unitig = 0;
		size_t i = 0;

		u_int32_t lastNode = -1;
		u_int32_t pathLength = 0;

		while(true){
			

			bool isFinished = true;
			for(size_t i=0; i<data_successors.size(); i++){
				if(!data_successors[i]._prevRankFinished){
					isFinished = false;
				}
			}
			//cout << "        " << isFinished << endl;
			if(isFinished) break;

			int prevIndex = _prevNodes.size() - prevRank - 1;
			if(prevIndex < 0 ) break;

			u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
			u_int32_t prev_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(prev_nodeIndex);

			if(i == 0){
				pathLength = graph->_nodeLengths[prev_nodeName];
			}
			else{
				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(prev_nodeIndex, lastNode);
				//cout << "LALZALsdf    " << overlapLength << endl;
				pathLength += (graph->_nodeLengths[prev_nodeName] - overlapLength);
			}

			//if(!print_debug){
			if(pathLength > 9000 && abundanceCutoff_min != 0){
				break;
			}
			//}
			//cout << prevIndex << " " << prev_nodeIndex << endl;
			//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);

			bool currentUnitigChanged = false;
			if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
				prevRank_unitig += 1;
				currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
				currentUnitigChanged = true;
			}

			//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
			string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

			size_t nbUnfinishedSuccessors = 0;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRankFinished){
					str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
					continue;
				}

				//if(std::find(successor._processedNodeIndex.begin(), successor._processedNodeIndex.end(), prev_nodeIndex) != successor._processedNodeIndex.end()){
				//	successor._prevRankFinished = true;
				//	str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
				//	continue;
				//} 

				u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor._nodeIndexSuccessor);
				
				bool isContigNode = _unitigDatas[prev_nodeName]._readIndexes.size() == 0;
				//u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//vector<u_int64_t> sharedReads;
				u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//cout << prev_nodeName << " " << successor_nodeName << " " << nbSharedReads << endl;
				//if(nbSharedReads > _abundanceCutoff_min/2){
				//if(nbSharedReads > successor._abundance/5){
				//if(nbSharedReads > 0){
				if(prevRank < 200 && (isContigNode || (nbSharedReads > 0 && nbSharedReads > abundanceCutoff_min))){

					successor._backtrackPathLength_noFilter = pathLength;


						if(successor._nbSharedReads == -1){
							successor._nbSharedReads = nbSharedReads;
						}
						//if(currentUnitigChanged || successor._nbSharedReads == -1){
							//successor._nbSharedReads = nbSharedReads;
						//}
						successor._prevUnitigIndex = currentUnitigIndex;
						successor._prevUnitigRank = prevRank_unitig; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						nbUnfinishedSuccessors += 1;
						successor._backtrackPathLength = pathLength;




				}
				else{
					successor._prevRankFinished = true;
				}



				str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
				//cout << "    " << prevRank << ": " <<  graph->nodeToString(current_nodeIndex) << "     " << graph->nodeToString(successor._nodeIndex)  << ": " << nbSharedReads;
				
				//if(foundPath) break;

				lastNode = prev_nodeIndex;
				i += 1;
			}

			str_debug += "    " + to_string(pathLength);


			#ifdef PRINT_PREV_RANK
				//if(print_debug && usePathSuccessors){
				if(usePathSuccessors){
					cout << "\t" << str_debug << endl;
				}
			#endif
			#ifdef PRINT_PREV_RANK_ALL
			//if(print_debug && currentDepth == 0){
				cout << str_debug << endl;
			//}
			#endif

			prevRank += 1;

			//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


		}

		if(print_debug) return;


		std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPathLength);
		int maxPathLength = data_successors[0]._backtrackPathLength;
		//cout << "MAX 2: " << maxPrevRank << endl;
		//vector<SuccessorData> successors_bestPrevRank;
		for(SuccessorData& successor : data_successors){
			if(successor._backtrackPathLength >= maxPathLength-_assemblyState._cutoff_backtrackLength){
				successors_bestPrevUnitigRank.push_back(successor);
			}
		}
		



	}

	bool isReachable(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, bool forward, GraphSimplify* graph, bool allowReverseDirection){

		//unordered_set<u_int32_t> isVisited;
        unordered_map<u_int32_t, u_int64_t> distance;
        vector<u_int32_t> queue;

        queue.push_back(from_nodeIndex);
        distance[from_nodeIndex] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];
			
			//cout << queue.size() << endl;
            queue.pop_back();

            if(distance[v] > maxLength) continue;

            //isVisited.insert(v);
            //if(seen.find(v) != seen.end()){
            //    seen.erase(v);
            //}

			vector<u_int32_t> successors;
			if(forward){
				graph->getSuccessors(v, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
			}
			else{
				graph->getPredecessors(v, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
			}

            for(u_int32_t u : successors){
				if(u == to_nodeIndex) return true;

				if(allowReverseDirection){
					if(BiGraph::nodeIndex_to_nodeName(u) == BiGraph::nodeIndex_to_nodeName(to_nodeIndex)) return true;
				}

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(v, u);
				distance[u] = distance[v] + (graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(u)] - overlapLength);
				queue.push_back(u);
			}

            //vector<u_int32_t> successors;
            //getSuccessors_unitig(v, successors);
            //if(successors.size() == 0) return -1; //abort tip
		}

		return false;

	}

	void updateNodeChosen(u_int32_t nodeIndex, u_int64_t currentDepth){

	}



	u_int64_t maxLength = 100000;
	u_int64_t maxLength_reachable = 20000;



	unordered_map<u_int32_t, u_int16_t> _nbVisitedTimesLala;

	void collectPossibleSuccessors(u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, bool print_debug, u_int32_t targetNodeIndex, unordered_set<u_int32_t> visitedNodes){
		
		//float currentAbundance = _currentAbundance;
		_nbVisitedTimesLala.clear();
	
		//unordered_set<u_int32_t> visitedNodes;
		//if(extractSubgraph){
		//	cout << "------"  << endl;
		//	file_test = ofstream("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		//	file_test << "Name,Colour" << endl;
		//}

		//vector<u_int32_t> succ;
    	//graph->extractSubGraph(source_nodeIndex, succ, 100000, forward);
		//exit(1);

		u_int32_t minNbJokers = -1;
		if(print_debug) cout << "\tCollect possible successors: " << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		u_int32_t currentDepth = 0;

		priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator> queue;
		queue.push({source_nodeIndex, UINT32_MAX, {}, 0, 0, _prevNodes, {}, {}, {}, _currentAbundance});



		while(true){
			
			if(queue.size() == 0){
				//cout << "exit 0" << endl;
				break;
			}

			//cout << currentDepth << endl;
			if(currentDepth > 100){
				//successorPaths.clear();
				return;
			}

			DataSuccessorPath dataSuccessorPath = queue.top();
			float currentAbundance = dataSuccessorPath._currentAbundance;
        	queue.pop();

			if(dataSuccessorPath._nbJokers > minNbJokers) continue;

			u_int32_t current_nodeIndex = dataSuccessorPath._currentNodeIndex; //dataSuccessorPath._path[dataSuccessorPath._path.size()-1];
			if(print_debug) cout << "\t\tStart: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << minNbJokers << " " << currentDepth << " " << (visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) << endl;
			//getchar();
			//u_int32_t currentLength = dataSuccessorPath._pathLength;

			PathExplorer pathExplorer(dataSuccessorPath._prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, currentAbundance, visitedNodes, _unitigDatas, dataSuccessorPath._pathLength, false, _assemblyState);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);

			bool continueVisiting = visitSuccessor(current_nodeIndex, dataSuccessorPath._prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
			if(!continueVisiting) continue;

			while(true){


				
				if(currentDepth > 100){
					//successorPaths.clear();
					return;
				}

				vector<u_int32_t> allSuccessors;
				if(forward){
					graph->getSuccessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				else{
					graph->getPredecessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				


				//cout << "\t\tTotal successors: " << allSuccessors.size() << endl;
				if(allSuccessors.size() == 0){
					break;
				}
				else if(allSuccessors.size() == 1){
					
					u_int32_t prevNodeIndex = current_nodeIndex;
					current_nodeIndex = allSuccessors[0];
					bool continueVisiting = visitSuccessor(current_nodeIndex, prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
					if(!continueVisiting) break;

				}
				else{

					u_int32_t resultType;
					vector<SuccessorData> nextNodes;
					pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, true);

					//cout << nextNodes.size() << endl;
					//cout << "\t\tNb valid successors: " << nextNodes.size() << endl;

					for(SuccessorData& successorData : nextNodes){
						//cout << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorData._nodeIndex, current_nodeIndex, queue, dataSuccessorPath, false, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					for(u_int32_t successorNodeIndex : allSuccessors){
						bool exist = false;
						for(SuccessorData& s : nextNodes){
							if(s._nodeIndex == successorNodeIndex){
								exist = true;
								break;
							}
						}
						if(exist) continue;

						//cout << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorNodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd invalid successor: " << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorNodeIndex, current_nodeIndex, queue, dataSuccessorPath, true, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					break;

				}


			}

		}

		
	}

	void addSuccessorPath(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator>& queue, const DataSuccessorPath& dataSuccessorPath, bool addJoker, GraphSimplify* graph, const PathExplorer& pathExplorer, u_int32_t& currentDepth, bool print_debug, float currentAbundance){

		//currentDepth += 1;
		if(currentDepth > 1000) return;
		//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

		u_int32_t nbJokers = dataSuccessorPath._nbJokers;
		if(addJoker) nbJokers += 1;

		//u_int32_t pathLength = dataSuccessorPath._pathLength;
		//u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
		//pathLength += (graph->_nodeLengths[nodeName] - overlapLength);

		DataSuccessorPath dataSuccessorPath_next = {nodeIndex, nodeIndexPrev, dataSuccessorPath._path, dataSuccessorPath._pathLength, nbJokers, pathExplorer._prevNodes, {}, {}, {}, currentAbundance};
		//dataSuccessorPath_next._path.push_back(nodeIndex);

		queue.push(dataSuccessorPath_next);
	}

	
	bool visitSuccessor(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, DataSuccessorPath& dataSuccessorPath, PathExplorer& pathExplorer, GraphSimplify* graph, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, u_int32_t& minNbJokers, bool print_debug, u_int32_t targetNodeIndex, float& currentAbundance){
		
		//cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;

		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		
		
		if(_nbVisitedTimesLala.find(nodeName) == _nbVisitedTimesLala.end()){
			_nbVisitedTimesLala[nodeName] = 0;
			//cout << "reset: " << nodeName << endl;
			for(auto it : _nbVisitedTimesLala){
				_nbVisitedTimesLala[it.first] = 0;
			}
		}

		//cout << "\t\tVisit successor: " << BiGraph::nodeIndex_to_nodeName(nodeIndex)  << " " << (_visitedNodes.find(nodeName) != _visitedNodes.end()) << " " << (dataSuccessorPath._nbVisitedTimes[nodeName]) << endl;
		
		_nbVisitedTimesLala[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
		if(_nbVisitedTimesLala[nodeName] > maxVisitable){
			//cout << "\t\t\tStop max visitables" << endl; 
			return false;
		}
		

		dataSuccessorPath._path.push_back(nodeIndex);

		//if(!isNodeVisited(nodeIndex, nodeIndexPrev)){
		//if(!extractSubgraph){
		if(pathExplorer._visitedNodes.find(nodeName) == pathExplorer._visitedNodes.end()){
			if(successorPaths.find(nodeIndex) == successorPaths.end()){
				if(targetNodeIndex == -1){
					successorPaths[nodeIndex] = dataSuccessorPath;
				}
				if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
			}
			else{
				if(dataSuccessorPath._nbJokers <= successorPaths[nodeIndex]._nbJokers){
					if(dataSuccessorPath._pathLength < successorPaths[nodeIndex]._pathLength){
						successorPaths[nodeIndex] = dataSuccessorPath;
						if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
					}
				}
			}
			if(print_debug) cout << "\t\t\t" << "Not visited: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << dataSuccessorPath._pathLength << endl;
			//if(BiGraph::nodeIndex_to_nodeName(nodeIndex) == 1847) getchar();
			//cout << "lala? " << targetNodeIndex << endl;
			if(targetNodeIndex == -1) return false;
		}
		//}

		if(dataSuccessorPath._pathLength > maxLength){//|| currentDepth > 50){
			//cout << "\t\t\tStop max length" << endl; 
			return false;
		}

		if(nodeIndexPrev != -1){ //First iteration
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
			dataSuccessorPath._pathLength += (graph->_nodeLengths[nodeName] - overlapLength);
		}
		pathExplorer.nodeExplored(nodeIndex, graph);
		//currentAbundance = PathExplorer::updateCurrentAbundance(nodeIndex, currentAbundance, graph, _assemblyState, _kminmerSize);

		return true;
	}


	bool isReachable2(u_int32_t current_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, bool allowReverseDirection, u_int32_t& totalIter){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		if(allowReverseDirection){
			if(BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == BiGraph::nodeIndex_to_nodeName(to_nodeIndex)) return true;
		}
		else{
			if(current_nodeIndex == to_nodeIndex) return true;
		}




		totalIter += 1;
		if(totalIter > 1000) return false;

		if(currentLength > maxLength_reachable || currentDepth > 50){
			return false;
		}


		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _currentAbundance, visitedNodes,_unitigDatas, currentLength, false, _assemblyState);
		pathExplorer.nodeExplored(current_nodeIndex, graph);
		//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(current_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);
		
		u_int32_t resultType;
		vector<SuccessorData> nextNodes;
		pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, false);

		for(SuccessorData& successorData : nextNodes){
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(current_nodeIndex, successorData._nodeIndex);
			bool isReachable = pathExplorer.isReachable2(successorData._nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)]-overlapLength), allowReverseDirection, totalIter);
			if(isReachable) return true;
		}

		return false;

		
	}



	bool isSmallCycle(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		//bool lala = true;


		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is simple cycle: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << " ?    ";
		}

		unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, from_nodeIndex, from_nodeIndex, _currentAbundance, visitedNodes, _unitigDatas, 0, false, _assemblyState);


		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(from_nodeIndex, graph);
		pathExplorer._visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(from_nodeIndex));
		//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(from_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			from_nodeIndex = pathExplorer.getNextNode(from_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, false, true);
			//cout << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;
			pathExplorer._visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(from_nodeIndex));

			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(from_nodeIndex == to_nodeIndex){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " Yes" << endl;
				}
				return true;
			}
			if(from_nodeIndex == -1){ //dead end or multiple braching path
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}



			pathExplorer.nodeExplored(from_nodeIndex, graph);
			//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(from_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);

			if(iter > maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " No" << endl;
		}

		//if(!lala) return false;
		return false;
	}


	void nodeExplored(u_int32_t nodeIndex, GraphSimplify* graph){
		
		bool orient_dummy;
		_prevNodes.push_back(nodeIndex);
		//_exploredNodes.push_back(nodeIndex);

		PathExplorer::clampPrevNodes(_prevNodes, _unitigDatas);
		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		_currentPathLength += graph->_nodeLengths[nodeName];

		_visitedNodes.insert(nodeName);

	}



	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

	static bool SuccessorComparator_byPrevUnitigRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevUnitigRank > b._prevUnitigRank;
	}

	static bool SuccessorComparator_byDistance(const SuccessorData &a, const SuccessorData &b){
		return a._distanceFromSource < b._distanceFromSource;
	}

	static bool SuccessorComparator_byAbundanceToSource(const SuccessorData &a, const SuccessorData &b){
		return a._sourceAbundance < b._sourceAbundance;
	}

	static bool SuccessorComparator_byNbSharedReads(const SuccessorData &a, const SuccessorData &b){
		return a._nbSharedReads > b._nbSharedReads;
	}

	static bool SuccessorComparator_byPathLength(const SuccessorData &a, const SuccessorData &b){
		return a._backtrackPathLength > b._backtrackPathLength;
	}

	static bool SuccessorComparator_byPathLength_noFilter(const SuccessorData &a, const SuccessorData &b){
		return a._backtrackPathLength_noFilter > b._backtrackPathLength_noFilter;
	}

	static bool SuccessorComparator_byDistanceFromSource(const SuccessorData &a, const SuccessorData &b){
		return a._path.size() < b._path.size();
	}

	static void clampPrevNodes(vector<u_int32_t>& prevNodes, vector<UnitigData>& _unitigDatas){

		//cout << "clamp" << endl;
		if(prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1]);
			UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(prevNodes.size() == 0) break;
				u_int32_t nodeIndex = prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(prevNodes[0]);

				const UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				//cout << nodeName_current << " " << nodeName_first << " " << Utils::shareAnyRead(unitigData_current, unitigData_first) << " " << unitigData_current._readIndexes.size() << " " << unitigData_first._readIndexes.size() << endl;
				if(Utils::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					prevNodes.erase(prevNodes.begin());
				}

			}
		}

	}



};
*/












































class PathExplorer{

	GraphSimplify* _graph;
	vector<u_int32_t> _prevNodes;
	vector<UnitigData>& _unitigDatas;
	bool _forward;
	unordered_set<u_int32_t>& _visitedNodes;

	PathExplorer(GraphSimplify* graph, vector<u_int32_t>& prevNodes, const vector<u_int32_t>& unitigDatas, bool forward, unordered_set<u_int32_t>& visitedNodes) : _prevNodes(prevNodes), _unitigDatas(unitigDatas): _visitedNodes(visitedNodes){
		_graph = graph;
		_forward = forward;
	}

	bool containsSuccessor(const SuccessorData& successor, const vector<SuccessorData>& successorData){
		for(const SuccessorData& s : successorData){
			if(s._nodeIndex == successor._nodeIndex) return true;
		}
		return false;
	}



	
};




























class PathExplorer_SuccessorCollector : PathExplorer{


public:


	u_int64_t _currentPathLength;

	PathExplorer_SuccessorCollector(GraphSimplify* graph, vector<u_int32_t>& prevNodes, const vector<u_int32_t>& unitigDatas, bool forward, unordered_set<u_int32_t>& visitedNodes) : PathExplorer(graph, prevNodes, unitigDatas, forward, visitedNodes){

	}

	void execute(u_int32_t source_nodeIndex){

	}





	u_int64_t maxLength = 100000;
	u_int64_t maxLength_reachable = 20000;



	unordered_map<u_int32_t, u_int16_t> _nbVisitedTimesLala;

	void collectPossibleSuccessors(u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, bool print_debug, u_int32_t targetNodeIndex, unordered_set<u_int32_t> visitedNodes){
		
		//float currentAbundance = _currentAbundance;
		_nbVisitedTimesLala.clear();
	
		//unordered_set<u_int32_t> visitedNodes;
		//if(extractSubgraph){
		//	cout << "------"  << endl;
		//	file_test = ofstream("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		//	file_test << "Name,Colour" << endl;
		//}

		//vector<u_int32_t> succ;
    	//graph->extractSubGraph(source_nodeIndex, succ, 100000, forward);
		//exit(1);

		u_int32_t minNbJokers = -1;
		if(print_debug) cout << "\tCollect possible successors: " << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		u_int32_t currentDepth = 0;

		priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator> queue;
		queue.push({source_nodeIndex, UINT32_MAX, {}, 0, 0, _prevNodes, {}, {}, {}, _currentAbundance});



		while(true){
			
			if(queue.size() == 0){
				//cout << "exit 0" << endl;
				break;
			}

			//cout << currentDepth << endl;
			if(currentDepth > 100){
				//successorPaths.clear();
				return;
			}

			DataSuccessorPath dataSuccessorPath = queue.top();
			float currentAbundance = dataSuccessorPath._currentAbundance;
        	queue.pop();

			if(dataSuccessorPath._nbJokers > minNbJokers) continue;

			u_int32_t current_nodeIndex = dataSuccessorPath._currentNodeIndex; //dataSuccessorPath._path[dataSuccessorPath._path.size()-1];
			if(print_debug) cout << "\t\tStart: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << minNbJokers << " " << currentDepth << " " << (visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) << endl;
			//getchar();
			//u_int32_t currentLength = dataSuccessorPath._pathLength;

			PathExplorer pathExplorer(dataSuccessorPath._prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, currentAbundance, visitedNodes, _unitigDatas, dataSuccessorPath._pathLength, false, _assemblyState);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);

			bool continueVisiting = visitSuccessor(current_nodeIndex, dataSuccessorPath._prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
			if(!continueVisiting) continue;

			while(true){

				//cout << currentDepth << endl;
				/*
            	if(extractSubgraph){
					cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
					file_test << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << "," << "red" << endl;
				}
				*/
				
				if(currentDepth > 100){
					//successorPaths.clear();
					return;
				}

				vector<u_int32_t> allSuccessors;
				if(forward){
					graph->getSuccessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				else{
					graph->getPredecessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				


				//cout << "\t\tTotal successors: " << allSuccessors.size() << endl;
				if(allSuccessors.size() == 0){
					break;
				}
				else if(allSuccessors.size() == 1){
					
					u_int32_t prevNodeIndex = current_nodeIndex;
					current_nodeIndex = allSuccessors[0];
					bool continueVisiting = visitSuccessor(current_nodeIndex, prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
					if(!continueVisiting) break;

				}
				else{

					u_int32_t resultType;
					vector<SuccessorData> nextNodes;
					pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, true);

					//cout << nextNodes.size() << endl;
					//cout << "\t\tNb valid successors: " << nextNodes.size() << endl;

					for(SuccessorData& successorData : nextNodes){
						//cout << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorData._nodeIndex, current_nodeIndex, queue, dataSuccessorPath, false, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					for(u_int32_t successorNodeIndex : allSuccessors){
						bool exist = false;
						for(SuccessorData& s : nextNodes){
							if(s._nodeIndex == successorNodeIndex){
								exist = true;
								break;
							}
						}
						if(exist) continue;

						//cout << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorNodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd invalid successor: " << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorNodeIndex, current_nodeIndex, queue, dataSuccessorPath, true, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					break;

				}


			}

		}

		
	}

	void addSuccessorPath(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator>& queue, const DataSuccessorPath& dataSuccessorPath, bool addJoker, GraphSimplify* graph, const PathExplorer& pathExplorer, u_int32_t& currentDepth, bool print_debug, float currentAbundance){

		//currentDepth += 1;
		if(currentDepth > 1000) return;
		//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

		u_int32_t nbJokers = dataSuccessorPath._nbJokers;
		if(addJoker) nbJokers += 1;

		//u_int32_t pathLength = dataSuccessorPath._pathLength;
		//u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
		//pathLength += (graph->_nodeLengths[nodeName] - overlapLength);

		DataSuccessorPath dataSuccessorPath_next = {nodeIndex, nodeIndexPrev, dataSuccessorPath._path, dataSuccessorPath._pathLength, nbJokers, pathExplorer._prevNodes, {}, {}, {}, currentAbundance};
		//dataSuccessorPath_next._path.push_back(nodeIndex);

		queue.push(dataSuccessorPath_next);
	}

	
	bool visitSuccessor(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, DataSuccessorPath& dataSuccessorPath, PathExplorer& pathExplorer, GraphSimplify* graph, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, u_int32_t& minNbJokers, bool print_debug, u_int32_t targetNodeIndex, float& currentAbundance){
		
		//cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;

		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		
		
		if(_nbVisitedTimesLala.find(nodeName) == _nbVisitedTimesLala.end()){
			_nbVisitedTimesLala[nodeName] = 0;
			//cout << "reset: " << nodeName << endl;
			for(auto it : _nbVisitedTimesLala){
				_nbVisitedTimesLala[it.first] = 0;
			}
		}

		//cout << "\t\tVisit successor: " << BiGraph::nodeIndex_to_nodeName(nodeIndex)  << " " << (_visitedNodes.find(nodeName) != _visitedNodes.end()) << " " << (dataSuccessorPath._nbVisitedTimes[nodeName]) << endl;
		
		_nbVisitedTimesLala[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
		if(_nbVisitedTimesLala[nodeName] > maxVisitable){
			//cout << "\t\t\tStop max visitables" << endl; 
			return false;
		}
		

		dataSuccessorPath._path.push_back(nodeIndex);

		//if(!isNodeVisited(nodeIndex, nodeIndexPrev)){
		//if(!extractSubgraph){
		if(pathExplorer._visitedNodes.find(nodeName) == pathExplorer._visitedNodes.end()){
			if(successorPaths.find(nodeIndex) == successorPaths.end()){
				if(targetNodeIndex == -1){
					successorPaths[nodeIndex] = dataSuccessorPath;
				}
				if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
			}
			else{
				if(dataSuccessorPath._nbJokers <= successorPaths[nodeIndex]._nbJokers){
					if(dataSuccessorPath._pathLength < successorPaths[nodeIndex]._pathLength){
						successorPaths[nodeIndex] = dataSuccessorPath;
						if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
					}
				}
			}
			if(print_debug) cout << "\t\t\t" << "Not visited: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << dataSuccessorPath._pathLength << endl;
			//if(BiGraph::nodeIndex_to_nodeName(nodeIndex) == 1847) getchar();
			//cout << "lala? " << targetNodeIndex << endl;
			if(targetNodeIndex == -1) return false;
		}
		//}

		if(dataSuccessorPath._pathLength > maxLength){//|| currentDepth > 50){
			//cout << "\t\t\tStop max length" << endl; 
			return false;
		}

		if(nodeIndexPrev != -1){ //First iteration
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
			dataSuccessorPath._pathLength += (graph->_nodeLengths[nodeName] - overlapLength);
		}
		pathExplorer.nodeExplored(nodeIndex, graph);
		//currentAbundance = PathExplorer::updateCurrentAbundance(nodeIndex, currentAbundance, graph, _assemblyState, _kminmerSize);

		return true;
	}




	u_int32_t getNextNode(u_int32_t current_nodeIndex, GraphSimplify* graph, bool forward, u_int32_t currentDepth, u_int32_t& resultType, vector<SuccessorData>& nextNodes, bool usePathSuccessors, bool canCheckTips){



		nextNodes.clear();
		resultType = 0;



		u_int32_t current_nodeName = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);


		vector<SuccessorData> data_successors;

		
		//if(usePathSuccessors){
		//	if(isInInfiniteCycle(current_nodeIndex, graph, _unitigDatas, forward)){
		//		return -1; //Infinite simple path without exit
		//	}
		//}


		vector<u_int32_t> successors;
		if(forward){
			graph->getSuccessors(current_nodeIndex, 0, successors);
		}
		else{
			graph->getPredecessors(current_nodeIndex, 0, successors);
		}

		//vector<u_int32_t> successors;
		//collectSolvedSuccessors(current_nodeIndex, directSuccessors, successors);
		//cout << directSuccessors.size() << " " << successors.size() << endl;

		if(successors.size() == 1){
			


			for(u_int32_t utg_n : successors){
				//if(!includeVisitedSuccessors && (_visitedNodes.find(utg_n) != _visitedNodes.end())) continue;

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {utg_n, successor_abundance, 0, 0, 0, false, utg_n, 0, 0, -1, {}, 0};
				//successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				successor._sourceAbundance = abs(successor._abundance - _currentAbundance);
				data_successors.push_back(successor);

			}
			


			/*
			vector<u_int32_t> path;
			u_int32_t sink = checkScc(current_nodeIndex, data_successors, graph, forward, false, path);
			//if(sink == -1) sink = checkScc(current_nodeIndex, data_successors, graph, forward, true);

			if(sink != -1){
				u_int32_t utg_n = sink;
				data_successors.clear();
				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {utg_n, successor_abundance, 0, 0, 0, false, utg_n, 0, 0, -1, path, 0};
				successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				data_successors.push_back(successor);

				//for(u_int32_t nodeIndex : successor._path){
				//	cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				//}
				//getchar();
				//exit(1);
			}*/

		}
		else{

			/*
			//cout << "--------------" << endl;
			//cout << "Current node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			vector<u_int32_t> pathDummy;
			unordered_map<u_int32_t, DataSuccessorPath> successorPaths;
			
			//for(u_int32_t utg_n : successors){
			unordered_map<u_int32_t, u_int16_t> nbVisitedTimes;
			u_int64_t lalalala = 0;
			collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, 0, 0, pathDummy, successorPaths, nbVisitedTimes, lalalala);
			//collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, 0, 0, possibleSuccessors);
			*/

			unordered_map<u_int32_t, DataSuccessorPath> successorPaths;

			if(successors.size() == 0 && _assemblyState._currentPathIndex != -1){
				_assemblyState._currentPathIndex = -1;
			}

			//cout << "hisdfsdf: " << _assemblyState._currentPathIndex << endl;
			

			collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, true, -1, _visitedNodes);



			for(auto& it : successorPaths){

				u_int32_t nodeIndex = it.first;
				DataSuccessorPath& pathData = it.second;

				int distance = 0;


				u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {nodeIndex, successor_abundance, 0, 0, 0, false, nodeIndex, distance, 0, -1, pathData._path, 0};
				//successor._sourceAbundance = abs((long)successor._abundance - (long)_source_abundance);
				successor._sourceAbundance = abs(successor._abundance - _currentAbundance);
				data_successors.push_back(successor);

			}
			

			if(currentDepth == 0){
				if(successors.size() > 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "\t----------- " << endl;
					cout << "\tNb prev nodes: " << _prevNodes.size() << endl;
					cout << "\tNb possible successors: " << successorPaths.size() << endl;
					cout << "\tLocal abundance: " << graph->getNodeUnitigAbundance(current_nodeIndex) << endl;
 					//for(u_int32_t utg_n : successors){
					for(SuccessorData& successor : data_successors){
						u_int32_t utg_n = successor._nodeIndex;
						//for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "\t" << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << graph->getNodeUnitigAbundance(utg_n) << " " << Utils::computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
					
					}
				}
			}


		}

	
		






		if(data_successors.size() == 0){
			if(currentDepth == 0){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "No successors" << endl;
			}
			return -1;
		}
		else if(data_successors.size() == 1){


			if(data_successors[0]._path.size() == 0){ //Si useSuccessorPath, le path est deja complet
				data_successors[0]._path.push_back(current_nodeIndex);
				current_nodeIndex = data_successors[0]._nodeIndexSuccessor;
				data_successors[0]._path.push_back(current_nodeIndex);
			}

			updateNodeChosen(current_nodeIndex, currentDepth);
			nextNodes.push_back(data_successors[0]);
			return current_nodeIndex;
		}
		else{


			resultType = 2;

			

			if(usePathSuccessors){
				//cout << "allo " << data_successors.size() << endl;
				vector<SuccessorData> dummy;
				computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, 0, dummy, 0, true, true, forward, usePathSuccessors); //ATTENTION NE PAS ENLEVER CA, on s'en sert pour calculer la _backtrackLength_noFilter (juste desactiver print_debug si besoin)
			}

			//removeUnsupportingReads(graph, _unitigDatas, data_successors, forward);
			
			SuccessorData bestSuccessorIndex;
			vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
			if(usePathSuccessors){
				computeBestSuccessors_byUnitigRank_all2(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
				
				if(_assemblyState._cutoff_backtrackLength == 0){	
					std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byNbSharedReads);
					bestSuccessorIndex = data_successors[0];

				}


			}
			else{
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, false, usePathSuccessors);
				//cout << successors_bestPrevUnitigRank_noCutoff.size() << endl;
			}

			if(usePathSuccessors){
			

			}

			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

				//if(currentDepth == 0){
					//cout << "\tNode chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				//}

				//updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
				if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
				
				return current_nodeIndex;
			}




			if(usePathSuccessors){
				//exit(1);
				//cout << "lalala " << successors_bestPrevUnitigRank_noCutoff.size() << endl;
				if(successors_bestPrevUnitigRank_noCutoff.size() >= 2){
					bool isBubble = true;
					for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){
						if(!graph->_isBubble[successor._nodeIndexSuccessor]){
							isBubble = false;
						}
					}
					cout << "\tIs bubble ?: " <<  isBubble << endl;
					if(isBubble){
						
						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byPathLength_noFilter);

						//if(successors_bestPrevUnitigRank_noCutoff[0]._backtrackPathLength == successors_bestPrevUnitigRank_noCutoff[1]._backtrackPathLength){
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//}

						//if(currentDepth == 0){
						cout << "\t\tTake bubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
						//}
						//exit(1);
						//updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
						//if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
						return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;
					}

					
					//cout << currentDepth << " " << "detect" << endl;
					u_int32_t superbubbleSink = detectSuperbubble(current_nodeIndex, maxLength, graph, _unitigDatas, forward);
					if(superbubbleSink != -1){

						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byPathLength_noFilter);

						//if(successors_bestPrevUnitigRank_noCutoff[0]._backtrackPathLength == successors_bestPrevUnitigRank_noCutoff[1]._backtrackPathLength){
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//}

						//sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//if(currentDepth == 0){
							cout << "\t\tTake superbubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
						//getchar();
						//}
						//exit(1);
						updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
						return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

					}
					//cout << "lala: " << BiGraph::nodeIndex_to_nodeName(superbubbleSink) << endl;
				}

				//}

				if(successors_bestPrevUnitigRank_noCutoff.size() < 5){

				
				
					//if(currentDepth == 0){

					if(currentDepth == 0){
						for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "\tCheck simple cycle" << endl;
					}

					for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){

						if(isSmallCycle(successor._nodeIndexSuccessor, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
							updateNodeChosen(successor._nodeIndexSuccessor, currentDepth);
							nextNodes.push_back(successor);
							if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
							return successor._nodeIndexSuccessor;
						}
						
					}



					for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
						u_int32_t to_nodeIndex = successors_bestPrevUnitigRank_noCutoff[i]._nodeIndexSuccessor;
						for(size_t j=0; j<successors_bestPrevUnitigRank_noCutoff.size(); j++){
							if(i == j) continue;
							
							u_int32_t from_nodeIndex = successors_bestPrevUnitigRank_noCutoff[j]._nodeIndexSuccessor;

							if(currentDepth == 0){
								for(size_t i=0; i<currentDepth; i++) cout << "  ";
								cout << "Check simple cycle from: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << "    to:    " << graph->_graphSuccessors->nodeToString(to_nodeIndex) << endl;
							}

							if(isSmallCycle(from_nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
								//exit(1);
								updateNodeChosen(from_nodeIndex, currentDepth);
								nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[j]);
								if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);
								return from_nodeIndex;
							}

						}	
					}
				}
			}
			
			for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[i]);
			}
			//if(_assemblyState._solvingComplexArea) if(!containsSuccessor(bestSuccessorIndex, nextNodes)) nextNodes.push_back(bestSuccessorIndex);


			//if(!usePathSuccessors && currentDepth == 0){
			//	return getNextNode(current_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, true);
			//}

			return -1;
				
			//}

		}


	}























	
	bool isInInfiniteCycle(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		//cout << "Check infinite cycle" << endl;
		Unitig& unitig = graph->nodeIndex_to_unitig(current_nodeIndex);
		if(unitig._startNode != current_nodeIndex) return false;

		unordered_map<u_int32_t, DataSuccessorPath> successorPaths;
		collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, true, -1, _visitedNodes);

		//cout << "Check infinite cycle: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " <<  successorPaths.size() << endl;
		if(successorPaths.size() == 0){
			cout << "Infinite cycle" << endl;
			//getchar();
		} 

		return successorPaths.size() == 0;
	}

    u_int32_t detectSuperbubble(u_int32_t nodeIndex_source, u_int64_t maxLength, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		cout << "\tStart detection: " << BiGraph::nodeIndex_to_nodeName(nodeIndex_source) << endl;



        unordered_set<u_int32_t> isVisited;
        unordered_set<u_int32_t> seen;
        unordered_map<u_int32_t, u_int64_t> pathLength;
        vector<u_int32_t> queue;

        queue.push_back(nodeIndex_source);
        pathLength[nodeIndex_source] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

			//cout << "Visit: " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << " " << pathLength[v] << endl;
            //cout << "\tVisited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._endNode) << endl;
            queue.pop_back();

			//cout << "Visit 1 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

            if(pathLength[v] > 10000) continue;

			//cout << "Visit 2 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

			if(isVisited.find(v) != isVisited.end()){
				return -1; // ?? quick patch
			}

            isVisited.insert(v);
            if(seen.find(v) != seen.end()){
                seen.erase(v);
            }

            
			
			//u_int32_t current_nodeIndex = graph->_unitigs[v]._endNode;
			//PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, 0, _solvedUnitigs);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);
			
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			getNextNode(v, graph, forward, 0, resultType, nextNodes, false, true);

			//cout << "\tNb succ: " << nextNodes.size() << endl;

            //vector<u_int32_t> successors;
			//getSuccessors_unitig(v, successors);
            if(nextNodes.size() == 0) return -1; //abort tip

			for(SuccessorData& successorData : nextNodes){
				u_int32_t u = successorData._nodeIndex;
				
				//cout << "\tSucc : " << BiGraph::nodeIndex_to_nodeName(u) << " " << queue.size() << endl;

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(v, u);
				//bool isReachable = pathExplorer.isReachable2(successorData._nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)]-overlapLength), allowReverseDirection);
				//if(isReachable) return true;
				//}

            	//for(u_int32_t u : successors){
                if(u == nodeIndex_source) return -1; //cycle including s

                if(isVisited.find(u) == isVisited.end()){

                    	seen.insert(u);
                    	pathLength[u] = pathLength[v] + (graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(u)]-overlapLength);
					



                }

            }

            //for(u_int32_t u : successors){
			for(SuccessorData& successorData : nextNodes){
				u_int32_t u = successorData._nodeIndex;
                //cout << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

				u_int32_t resultType;
				vector<SuccessorData> prevNodes;
				getNextNode(u, graph, !forward, 0, resultType, prevNodes, false, true);
				//getPredecessors_unitig(u, predecessors);
                bool allPredecessorsAreVisited = true;
				for(SuccessorData& predecessorData : prevNodes){
					u_int32_t p = predecessorData._nodeIndex;
					//cout << "\tPred:" << BiGraph::nodeIndex_to_nodeName(p) << endl;
                	//for(u_int32_t p : predecessors){
                    if(isVisited.find(p) == isVisited.end()){
                        allPredecessorsAreVisited = false;
                        break;
                    }
                }

                if(allPredecessorsAreVisited){
					//cout << "\tPush: " << BiGraph::nodeIndex_to_nodeName(u) << endl;
                    //cout << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                    queue.push_back(u);
                }

                //cout << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                    u_int32_t t = queue[0];
					
					vector<SuccessorData> successors_t;
					getNextNode(t, graph, forward, 0, resultType, successors_t, false, true);

					bool foundSource = false;
					for(SuccessorData& successorData_t : successors_t){
						if(successorData_t._nodeIndex == nodeIndex_source){
							foundSource = true;
							break;
						}
					}

                    //vector<u_int32_t> successors_t;
                    //getSuccessors_unitig(t, successors_t);

                    if(!foundSource){
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


	unordered_set<u_int32_t> __beatenNodeIndex;

	void computeBestSuccessors_byUnitigRank_all2(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors2, int currentDepth, bool forward, vector<SuccessorData>& validSuccessors, bool includeVisitedSuccessors, bool usePathSuccessors){
	
		//cout << "1" << endl;
		__beatenNodeIndex.clear();



		for(size_t i=0; i<data_successors2.size(); i++){
			for(size_t j=i+1; j<data_successors2.size(); j++){

				if(__beatenNodeIndex.find(data_successors2[i]._nodeIndex) != __beatenNodeIndex.end()) continue;
				if(__beatenNodeIndex.find(data_successors2[j]._nodeIndex) != __beatenNodeIndex.end()) continue;


				//cout << i << " " << j << " " << _assemblyState._cutoff_backtrackLength << " " << BiGraph::nodeIndex_to_nodeName(data_successors2[i]._nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(data_successors2[j]._nodeIndex) << endl;

				//cout << i << " vs " << j << endl;
				vector<SuccessorData> data_successors;
				data_successors.push_back(data_successors2[i]);
				data_successors.push_back(data_successors2[j]);


				//cout << data_successors[0]._path.size() << " " << data_successors[1]._path.size() << endl;
				if(_assemblyState._cutoff_backtrackLength > 0){

					unordered_set<u_int32_t> ind;
					for(u_int32_t nodeIndex : data_successors[0]._path){
						ind.insert(nodeIndex);
					}

					u_int32_t sharedNodexIndex = -1;
					//u_int32_t lastSharedRank = 0;
					//size_t index=0;
					for(u_int32_t nodeIndex : data_successors[1]._path){
						if(ind.find(nodeIndex) != ind.end() || find(_prevNodes.begin(), _prevNodes.end(), nodeIndex) != _prevNodes.end()){
							//lastSharedRank = index;
							//cout << "shared rank: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
							sharedNodexIndex = nodeIndex;
						}
						//index += 1;
					}

					//Find LAST index of sharedNodexIndex (this node can be repeated in the path)
					u_int32_t index1 = 0;
					for(size_t i=0; i<data_successors[0]._path.size(); i++){
						if(data_successors[0]._path[i] == sharedNodexIndex){
							index1 = i;
						}
					}
					u_int32_t index2 = 0;
					for(size_t i=0; i<data_successors[1]._path.size(); i++){
						if(data_successors[1]._path[i] == sharedNodexIndex){
							index2 = i;
						}
					}

					//cout << "Shared node: " << sharedNodexIndex << endl;

					u_int32_t rank1 = index1 + 1;
					u_int32_t rank2 = index2 + 1;


					data_successors[0]._nodeIndexSuccessor = data_successors[0]._path[rank1];
					data_successors[1]._nodeIndexSuccessor = data_successors[1]._path[rank2];

					//cout << "NodeIndexSuccessor " << BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndex) << ": " << BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndexSuccessor) << endl;
					//cout << "NodeIndexSuccessor " << BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndex) << ": " << BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndexSuccessor) << endl;
				}



				u_int32_t prev_nodeIndex = _prevNodes[_prevNodes.size()-1];
				u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);

				u_int32_t nodeName1 = BiGraph::nodeIndex_to_nodeName(data_successors[0]._nodeIndexSuccessor);
				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(data_successors[1]._nodeIndexSuccessor);

				vector<SuccessorData> bestSuccessors;
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, bestSuccessors, includeVisitedSuccessors, usePathSuccessors);
				

					
					if(bestSuccessors.size() == 1){
						cout << "\tWinner: " << bestSuccessors.size() << " " << BiGraph::nodeIndex_to_nodeName(bestSuccessors[0]._nodeIndex) << endl;
						if(bestSuccessors[0]._nodeIndex == data_successors[0]._nodeIndex){
							__beatenNodeIndex.insert(data_successors[1]._nodeIndex);
						}
						else{
							__beatenNodeIndex.insert(data_successors[0]._nodeIndex);
						}
					}
				//}
				

			}
		}

 		validSuccessors.clear();
		for(SuccessorData& s : data_successors2){
			if(__beatenNodeIndex.find(s._nodeIndex) != __beatenNodeIndex.end()) continue;
			validSuccessors.push_back(s);
		}


		//cout << validSuccessors.size() << " " << __beatenNodeIndex.size() << endl;
		if(validSuccessors.size() == 0) return;
		if(validSuccessors.size() > 10) return;
		//if(currentDepth == 0 && abundanceCutoff_min != 0){

		//cout << "3" << endl;
		unordered_set<u_int32_t> reachableSuccessors;
		
		if(_assemblyState._cutoffType == CutoffType::STRAIN_HIGH){
			
			//unordered_set<u_int32_t> reachableSuccessors;
			for(size_t i=0; i<validSuccessors.size(); i++){
				for(size_t j=0; j<validSuccessors.size(); j++){
					if(i == j) continue; 
					//if(reachableSuccessors.find(successors_bestPrevUnitigRank[j]._nodeIndex) != reachableSuccessors.end()) continue; //Already reached
					cout << "\tCheck reachable (" << currentDepth  << "): " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
					
					u_int32_t totalIter = 0;
					if(isReachable2(validSuccessors[i]._nodeIndex, validSuccessors[j]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, false, totalIter)){//} && isReachable2(validSuccessors[j]._nodeIndex, validSuccessors[i]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, false)){
						cout << "\t\tIs reachable: " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
						reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
					}
				}
			}
				
					
					
					
			for(size_t i=0; i<validSuccessors.size(); i++){
				for(size_t j=i+1; j<validSuccessors.size(); j++){
					//if(i == j) continue; 
					if(reachableSuccessors.find(validSuccessors[j]._nodeIndex) != reachableSuccessors.end()) continue; //Already reached
					
					u_int32_t totalIter = 0;
					if(isReachable2(validSuccessors[i]._nodeIndex, validSuccessors[j]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, true, totalIter) && isReachable2(validSuccessors[j]._nodeIndex, validSuccessors[i]._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0, true, totalIter)){
						cout << "\tIs reachable (RC): " << BiGraph::nodeIndex_to_nodeName(validSuccessors[j]._nodeIndex) << " (From " << BiGraph::nodeIndex_to_nodeName(validSuccessors[i]._nodeIndex) << ")" << endl;
						//cout << validSuccessors[i]._nbSharedReads << " " << validSuccessors[j]._nbSharedReads << endl;
						if(validSuccessors[i]._nbSharedReads > validSuccessors[j]._nbSharedReads){
							reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
						}
						else if(validSuccessors[j]._nbSharedReads > validSuccessors[i]._nbSharedReads){
							reachableSuccessors.insert(validSuccessors[i]._nodeIndex);
						}
						else{
							cout << "TODO egalite nb shared reads, on peut prendre en compte la distance au node source apr exemple" << endl;
								
							if(validSuccessors[i]._backtrackPathLength > validSuccessors[j]._backtrackPathLength){
								reachableSuccessors.insert(validSuccessors[j]._nodeIndex);
							}
							else{
								reachableSuccessors.insert(validSuccessors[i]._nodeIndex);
							}

							//exit(1);
						}
					}
				}
			}
		}
		

		vector<SuccessorData> chosenSuccessors;
		for(SuccessorData& successorData : validSuccessors){
			if(reachableSuccessors.find(successorData._nodeIndex) != reachableSuccessors.end()) continue;
			chosenSuccessors.push_back(successorData);
			//cout << "lalala: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
		}

		if(chosenSuccessors.size() == 0 && reachableSuccessors.size() > 0){

			vector<SuccessorData> chosenSuccessors_tmp;
			for(SuccessorData& successorData : validSuccessors){
				chosenSuccessors_tmp.push_back(successorData);
			}

			sort(chosenSuccessors_tmp.begin(), chosenSuccessors_tmp.end(), SuccessorComparator_byDistanceFromSource);
			//sort(chosenSuccessors_tmp.begin(), chosenSuccessors_tmp.end(), SuccessorComparator_byNbSharedReads);
			chosenSuccessors.push_back(chosenSuccessors_tmp[0]);
		}

		validSuccessors.clear();
		for(SuccessorData& successorData : chosenSuccessors){
			validSuccessors.push_back(successorData);
		}



		//}

	}

	void computeBestSuccessors_byUnitigRank_all(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors2, int currentDepth, bool forward, vector<SuccessorData>& validSuccessors, bool includeVisitedSuccessors, bool usePathSuccessors){

		//cout << "Include visited successors: " << includeVisitedSuccessors << endl;
		//cout << currentDepth << " " << usePathSuccessors << endl;
		float cutoff_min = PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType) / 2.0; //max(1.0, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0)/4.0);
		cutoff_min = (int)(cutoff_min + 1); //max(cutoff_min, 1);
		//(!usePathSuccessors){
		//	cutoff_min = 1;
		//}

		//cout << "Cutoff min: " << cutoff_min << endl;

		
		vector<SuccessorData> data_successors;
		for(SuccessorData& s : data_successors2){
			//cout << "\t" << BiGraph::nodeIndex_to_nodeName(s._nodeIndex) << " " << (_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(s._nodeIndex)) != _visitedNodes.end()) << endl;
			if(!includeVisitedSuccessors && _visitedNodes.find(BiGraph::nodeIndex_to_nodeName(s._nodeIndex)) != _visitedNodes.end()) continue;
			data_successors.push_back(s);
			

			//if(usePathSuccessors){
			//	cout << "Possible succ: " << BiGraph::nodeIndex_to_nodeName(s._nodeIndex) << "   " << BiGraph::nodeIndex_to_nodeName(s._path[1]) << endl;
			//}
		}
		

		if(data_successors.size() == 0) return;
		//Print debug prev ranks


		vector<SuccessorData> successors_bestPrevUnitigRank_cutoff;
		//if(usePathSuccessors){
			//cout << cutoff_min << endl;
		computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_cutoff, cutoff_min, true, false, forward, usePathSuccessors);
		if(successors_bestPrevUnitigRank_cutoff.size() == 1){
			u_int32_t current_nodeIndex = successors_bestPrevUnitigRank_cutoff[0]._nodeIndexSuccessor;

			//if(currentDepth == 1){
			//cout << "\tNode chosen 1: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;// << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
			//}

			//updateNodeChosen(current_nodeIndex, currentDepth);
			//nextNodes.push_back(current_nodeIndex);
			//return current_nodeIndex;

			for(SuccessorData& s : validSuccessors){
				if(s._nodeIndex == current_nodeIndex) return;
			}

			//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << currentDepth << " " << (_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != _visitedNodes.end()) << endl; 
			validSuccessors.push_back(successors_bestPrevUnitigRank_cutoff[0]);
			return;
		}

		if(usePathSuccessors){
			for(SuccessorData& successorData : successors_bestPrevUnitigRank_cutoff){
				
				bool exist = false;
				for(SuccessorData& s : validSuccessors){
					if(s._nodeIndex == successorData._nodeIndex){
						exist = true;
						break;
					}
				}
				if(exist) continue;

				//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl; 
				validSuccessors.push_back(successorData);
			}
		}
		else{
			for(SuccessorData& successorData : successors_bestPrevUnitigRank_cutoff){
				
				bool exist = false;
				for(SuccessorData& s : validSuccessors){
					if(s._nodeIndex == successorData._nodeIndex){
						exist = true;
						break;
					}
				}
				if(exist) continue;

				//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl; 
				validSuccessors.push_back(successorData);
			}
		}
	}

	void computeBestSuccessors_byUnitigRank(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors, int currentDepth, vector<SuccessorData>& successors_bestPrevUnitigRank, float abundanceCutoff_min, bool useUnitigRank, bool print_debug, bool forward, bool usePathSuccessors){



		#ifdef PRINT_PREV_RANK
			//if(print_debug && usePathSuccessors){
			if(usePathSuccessors){
				cout << endl << "\tCutoff: " << abundanceCutoff_min << endl;
			}
		#endif
		#ifdef PRINT_PREV_RANK_ALL
		//if(print_debug && currentDepth == 0){
			cout << endl << "\tCutoff: " << abundanceCutoff_min << endl;
		//}
		#endif

		for(SuccessorData& successor : data_successors){
			successor._prevRankFinished = false;
			successor._prevUnitigRank = -1;
			successor._prevRank = -1;
			//successor._nbSharedReadsPrev = -1;
			successor._backtrackPathLength = 0;
			successor._backtrackPathLength_noFilter = 0;
		}

		successors_bestPrevUnitigRank.clear();

		u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

		u_int32_t prevRank = 0;
		u_int32_t prevRank_unitig = 0;
		size_t i = 0;

		u_int32_t lastNode = -1;
		u_int32_t pathLength = 0;

		while(true){
			

			bool isFinished = true;
			for(size_t i=0; i<data_successors.size(); i++){
				if(!data_successors[i]._prevRankFinished){
					isFinished = false;
				}
			}
			//cout << "        " << isFinished << endl;
			if(isFinished) break;

			int prevIndex = _prevNodes.size() - prevRank - 1;
			if(prevIndex < 0 ) break;

			u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
			u_int32_t prev_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(prev_nodeIndex);

			if(i == 0){
				pathLength = graph->_nodeLengths[prev_nodeName];
			}
			else{
				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(prev_nodeIndex, lastNode);
				//cout << "LALZALsdf    " << overlapLength << endl;
				pathLength += (graph->_nodeLengths[prev_nodeName] - overlapLength);
			}

			//if(!print_debug){
			if(pathLength > 9000 && abundanceCutoff_min != 0){
				break;
			}
			//}
			//cout << prevIndex << " " << prev_nodeIndex << endl;
			//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);

			bool currentUnitigChanged = false;
			if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
				prevRank_unitig += 1;
				currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
				currentUnitigChanged = true;
			}

			//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
			string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

			size_t nbUnfinishedSuccessors = 0;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRankFinished){
					str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
					continue;
				}

				//if(std::find(successor._processedNodeIndex.begin(), successor._processedNodeIndex.end(), prev_nodeIndex) != successor._processedNodeIndex.end()){
				//	successor._prevRankFinished = true;
				//	str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
				//	continue;
				//} 

				u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor._nodeIndexSuccessor);
				
				bool isContigNode = _unitigDatas[prev_nodeName]._readIndexes.size() == 0;
				//u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//vector<u_int64_t> sharedReads;
				u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//cout << prev_nodeName << " " << successor_nodeName << " " << nbSharedReads << endl;
				//if(nbSharedReads > _abundanceCutoff_min/2){
				//if(nbSharedReads > successor._abundance/5){
				//if(nbSharedReads > 0){
				if(prevRank < 200 && (isContigNode || (nbSharedReads > 0 && nbSharedReads > abundanceCutoff_min))){

					successor._backtrackPathLength_noFilter = pathLength;


						if(successor._nbSharedReads == -1){
							successor._nbSharedReads = nbSharedReads;
						}
						//if(currentUnitigChanged || successor._nbSharedReads == -1){
							//successor._nbSharedReads = nbSharedReads;
						//}
						successor._prevUnitigIndex = currentUnitigIndex;
						successor._prevUnitigRank = prevRank_unitig; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						nbUnfinishedSuccessors += 1;
						successor._backtrackPathLength = pathLength;




				}
				else{
					successor._prevRankFinished = true;
				}



				str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
				//cout << "    " << prevRank << ": " <<  graph->nodeToString(current_nodeIndex) << "     " << graph->nodeToString(successor._nodeIndex)  << ": " << nbSharedReads;
				
				//if(foundPath) break;

				lastNode = prev_nodeIndex;
				i += 1;
			}

			str_debug += "    " + to_string(pathLength);


			#ifdef PRINT_PREV_RANK
				//if(print_debug && usePathSuccessors){
				if(usePathSuccessors){
					cout << "\t" << str_debug << endl;
				}
			#endif
			#ifdef PRINT_PREV_RANK_ALL
			//if(print_debug && currentDepth == 0){
				cout << str_debug << endl;
			//}
			#endif

			prevRank += 1;

			//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


		}

		if(print_debug) return;


		std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPathLength);
		int maxPathLength = data_successors[0]._backtrackPathLength;
		//cout << "MAX 2: " << maxPrevRank << endl;
		//vector<SuccessorData> successors_bestPrevRank;
		for(SuccessorData& successor : data_successors){
			if(successor._backtrackPathLength >= maxPathLength-_assemblyState._cutoff_backtrackLength){
				successors_bestPrevUnitigRank.push_back(successor);
			}
		}
		



	}

	bool isReachable(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, bool forward, GraphSimplify* graph, bool allowReverseDirection){

		//unordered_set<u_int32_t> isVisited;
        unordered_map<u_int32_t, u_int64_t> distance;
        vector<u_int32_t> queue;

        queue.push_back(from_nodeIndex);
        distance[from_nodeIndex] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];
			
			//cout << queue.size() << endl;
            queue.pop_back();

            if(distance[v] > maxLength) continue;

            //isVisited.insert(v);
            //if(seen.find(v) != seen.end()){
            //    seen.erase(v);
            //}

			vector<u_int32_t> successors;
			if(forward){
				graph->getSuccessors(v, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
			}
			else{
				graph->getPredecessors(v, PathExplorer::computeAbundanceCutoff(_currentAbundance, 0, _assemblyState._cutoffType), successors);
			}

            for(u_int32_t u : successors){
				if(u == to_nodeIndex) return true;

				if(allowReverseDirection){
					if(BiGraph::nodeIndex_to_nodeName(u) == BiGraph::nodeIndex_to_nodeName(to_nodeIndex)) return true;
				}

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(v, u);
				distance[u] = distance[v] + (graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(u)] - overlapLength);
				queue.push_back(u);
			}

            //vector<u_int32_t> successors;
            //getSuccessors_unitig(v, successors);
            //if(successors.size() == 0) return -1; //abort tip
		}

		return false;

	}

	void updateNodeChosen(u_int32_t nodeIndex, u_int64_t currentDepth){

	}



	u_int64_t maxLength = 100000;
	u_int64_t maxLength_reachable = 20000;



	unordered_map<u_int32_t, u_int16_t> _nbVisitedTimesLala;

	void collectPossibleSuccessors(u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, bool print_debug, u_int32_t targetNodeIndex, unordered_set<u_int32_t> visitedNodes){
		
		//float currentAbundance = _currentAbundance;
		_nbVisitedTimesLala.clear();
	
		//unordered_set<u_int32_t> visitedNodes;
		//if(extractSubgraph){
		//	cout << "------"  << endl;
		//	file_test = ofstream("/home/gats/workspace/run/overlap_test_AD/subgraph.csv");
		//	file_test << "Name,Colour" << endl;
		//}

		//vector<u_int32_t> succ;
    	//graph->extractSubGraph(source_nodeIndex, succ, 100000, forward);
		//exit(1);

		u_int32_t minNbJokers = -1;
		if(print_debug) cout << "\tCollect possible successors: " << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		u_int32_t currentDepth = 0;

		priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator> queue;
		queue.push({source_nodeIndex, UINT32_MAX, {}, 0, 0, _prevNodes, {}, {}, {}, _currentAbundance});



		while(true){
			
			if(queue.size() == 0){
				//cout << "exit 0" << endl;
				break;
			}

			//cout << currentDepth << endl;
			if(currentDepth > 100){
				//successorPaths.clear();
				return;
			}

			DataSuccessorPath dataSuccessorPath = queue.top();
			float currentAbundance = dataSuccessorPath._currentAbundance;
        	queue.pop();

			if(dataSuccessorPath._nbJokers > minNbJokers) continue;

			u_int32_t current_nodeIndex = dataSuccessorPath._currentNodeIndex; //dataSuccessorPath._path[dataSuccessorPath._path.size()-1];
			if(print_debug) cout << "\t\tStart: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << minNbJokers << " " << currentDepth << " " << (visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != visitedNodes.end()) << endl;
			//getchar();
			//u_int32_t currentLength = dataSuccessorPath._pathLength;

			PathExplorer pathExplorer(dataSuccessorPath._prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, currentAbundance, visitedNodes, _unitigDatas, dataSuccessorPath._pathLength, false, _assemblyState);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);

			bool continueVisiting = visitSuccessor(current_nodeIndex, dataSuccessorPath._prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
			if(!continueVisiting) continue;

			while(true){

				//cout << currentDepth << endl;
				/*
            	if(extractSubgraph){
					cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
					file_test << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << "," << "red" << endl;
				}
				*/
				
				if(currentDepth > 100){
					//successorPaths.clear();
					return;
				}

				vector<u_int32_t> allSuccessors;
				if(forward){
					graph->getSuccessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				else{
					graph->getPredecessors(current_nodeIndex, PathExplorer::computeAbundanceCutoff(currentAbundance, 0, _assemblyState._cutoffType), allSuccessors);
				}
				


				//cout << "\t\tTotal successors: " << allSuccessors.size() << endl;
				if(allSuccessors.size() == 0){
					break;
				}
				else if(allSuccessors.size() == 1){
					
					u_int32_t prevNodeIndex = current_nodeIndex;
					current_nodeIndex = allSuccessors[0];
					bool continueVisiting = visitSuccessor(current_nodeIndex, prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug, targetNodeIndex, currentAbundance);
					if(!continueVisiting) break;

				}
				else{

					u_int32_t resultType;
					vector<SuccessorData> nextNodes;
					pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, true);

					//cout << nextNodes.size() << endl;
					//cout << "\t\tNb valid successors: " << nextNodes.size() << endl;

					for(SuccessorData& successorData : nextNodes){
						//cout << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorData._nodeIndex, current_nodeIndex, queue, dataSuccessorPath, false, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					for(u_int32_t successorNodeIndex : allSuccessors){
						bool exist = false;
						for(SuccessorData& s : nextNodes){
							if(s._nodeIndex == successorNodeIndex){
								exist = true;
								break;
							}
						}
						if(exist) continue;

						//cout << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//if(BiGraph::nodeIndex_to_nodeName(successorNodeIndex) == 840) getchar();
						//cout << "\t\t\tAdd invalid successor: " << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorNodeIndex, current_nodeIndex, queue, dataSuccessorPath, true, graph, pathExplorer, currentDepth, print_debug, currentAbundance);
					}

					break;

				}


			}

		}

		
	}

	void addSuccessorPath(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator>& queue, const DataSuccessorPath& dataSuccessorPath, bool addJoker, GraphSimplify* graph, const PathExplorer& pathExplorer, u_int32_t& currentDepth, bool print_debug, float currentAbundance){

		//currentDepth += 1;
		if(currentDepth > 1000) return;
		//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

		u_int32_t nbJokers = dataSuccessorPath._nbJokers;
		if(addJoker) nbJokers += 1;

		//u_int32_t pathLength = dataSuccessorPath._pathLength;
		//u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
		//pathLength += (graph->_nodeLengths[nodeName] - overlapLength);

		DataSuccessorPath dataSuccessorPath_next = {nodeIndex, nodeIndexPrev, dataSuccessorPath._path, dataSuccessorPath._pathLength, nbJokers, pathExplorer._prevNodes, {}, {}, {}, currentAbundance};
		//dataSuccessorPath_next._path.push_back(nodeIndex);

		queue.push(dataSuccessorPath_next);
	}

	
	bool visitSuccessor(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, DataSuccessorPath& dataSuccessorPath, PathExplorer& pathExplorer, GraphSimplify* graph, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, u_int32_t& minNbJokers, bool print_debug, u_int32_t targetNodeIndex, float& currentAbundance){
		
		//cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;

		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		
		
		if(_nbVisitedTimesLala.find(nodeName) == _nbVisitedTimesLala.end()){
			_nbVisitedTimesLala[nodeName] = 0;
			//cout << "reset: " << nodeName << endl;
			for(auto it : _nbVisitedTimesLala){
				_nbVisitedTimesLala[it.first] = 0;
			}
		}

		//cout << "\t\tVisit successor: " << BiGraph::nodeIndex_to_nodeName(nodeIndex)  << " " << (_visitedNodes.find(nodeName) != _visitedNodes.end()) << " " << (dataSuccessorPath._nbVisitedTimes[nodeName]) << endl;
		
		_nbVisitedTimesLala[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
		if(_nbVisitedTimesLala[nodeName] > maxVisitable){
			//cout << "\t\t\tStop max visitables" << endl; 
			return false;
		}
		

		dataSuccessorPath._path.push_back(nodeIndex);

		//if(!isNodeVisited(nodeIndex, nodeIndexPrev)){
		//if(!extractSubgraph){
		if(pathExplorer._visitedNodes.find(nodeName) == pathExplorer._visitedNodes.end()){
			if(successorPaths.find(nodeIndex) == successorPaths.end()){
				if(targetNodeIndex == -1){
					successorPaths[nodeIndex] = dataSuccessorPath;
				}
				if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
			}
			else{
				if(dataSuccessorPath._nbJokers <= successorPaths[nodeIndex]._nbJokers){
					if(dataSuccessorPath._pathLength < successorPaths[nodeIndex]._pathLength){
						successorPaths[nodeIndex] = dataSuccessorPath;
						if(dataSuccessorPath._nbJokers < minNbJokers) minNbJokers = dataSuccessorPath._nbJokers;
					}
				}
			}
			if(print_debug) cout << "\t\t\t" << "Not visited: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << dataSuccessorPath._pathLength << endl;
			//if(BiGraph::nodeIndex_to_nodeName(nodeIndex) == 1847) getchar();
			//cout << "lala? " << targetNodeIndex << endl;
			if(targetNodeIndex == -1) return false;
		}
		//}

		if(dataSuccessorPath._pathLength > maxLength){//|| currentDepth > 50){
			//cout << "\t\t\tStop max length" << endl; 
			return false;
		}

		if(nodeIndexPrev != -1){ //First iteration
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
			dataSuccessorPath._pathLength += (graph->_nodeLengths[nodeName] - overlapLength);
		}
		pathExplorer.nodeExplored(nodeIndex, graph);
		//currentAbundance = PathExplorer::updateCurrentAbundance(nodeIndex, currentAbundance, graph, _assemblyState, _kminmerSize);

		return true;
	}


	bool isReachable2(u_int32_t current_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, bool allowReverseDirection, u_int32_t& totalIter){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		if(allowReverseDirection){
			if(BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == BiGraph::nodeIndex_to_nodeName(to_nodeIndex)) return true;
		}
		else{
			if(current_nodeIndex == to_nodeIndex) return true;
		}




		totalIter += 1;
		if(totalIter > 1000) return false;

		if(currentLength > maxLength_reachable || currentDepth > 50){
			return false;
		}


		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _currentAbundance, visitedNodes,_unitigDatas, currentLength, false, _assemblyState);
		pathExplorer.nodeExplored(current_nodeIndex, graph);
		//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(current_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);
		
		u_int32_t resultType;
		vector<SuccessorData> nextNodes;
		pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false, false);

		for(SuccessorData& successorData : nextNodes){
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(current_nodeIndex, successorData._nodeIndex);
			bool isReachable = pathExplorer.isReachable2(successorData._nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)]-overlapLength), allowReverseDirection, totalIter);
			if(isReachable) return true;
		}

		return false;

		
	}



	bool isSmallCycle(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		//bool lala = true;


		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is simple cycle: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << " ?    ";
		}

		unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, from_nodeIndex, from_nodeIndex, _currentAbundance, visitedNodes, _unitigDatas, 0, false, _assemblyState);


		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(from_nodeIndex, graph);
		pathExplorer._visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(from_nodeIndex));
		//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(from_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			from_nodeIndex = pathExplorer.getNextNode(from_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, false, true);
			//cout << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;
			pathExplorer._visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(from_nodeIndex));

			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(from_nodeIndex == to_nodeIndex){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " Yes" << endl;
				}
				return true;
			}
			if(from_nodeIndex == -1){ //dead end or multiple braching path
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}



			pathExplorer.nodeExplored(from_nodeIndex, graph);
			//pathExplorer._currentAbundance = PathExplorer::updateCurrentAbundance(from_nodeIndex, pathExplorer._currentAbundance, graph, _assemblyState, _kminmerSize);

			if(iter > maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " No" << endl;
		}

		//if(!lala) return false;
		return false;
	}


	void nodeExplored(u_int32_t nodeIndex, GraphSimplify* graph){
		
		bool orient_dummy;
		_prevNodes.push_back(nodeIndex);
		//_exploredNodes.push_back(nodeIndex);

		PathExplorer::clampPrevNodes(_prevNodes, _unitigDatas);
		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		_currentPathLength += graph->_nodeLengths[nodeName];

		_visitedNodes.insert(nodeName);

	}



	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

	static bool SuccessorComparator_byPrevUnitigRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevUnitigRank > b._prevUnitigRank;
	}

	static bool SuccessorComparator_byDistance(const SuccessorData &a, const SuccessorData &b){
		return a._distanceFromSource < b._distanceFromSource;
	}

	static bool SuccessorComparator_byAbundanceToSource(const SuccessorData &a, const SuccessorData &b){
		return a._sourceAbundance < b._sourceAbundance;
	}

	static bool SuccessorComparator_byNbSharedReads(const SuccessorData &a, const SuccessorData &b){
		return a._nbSharedReads > b._nbSharedReads;
	}

	static bool SuccessorComparator_byPathLength(const SuccessorData &a, const SuccessorData &b){
		return a._backtrackPathLength > b._backtrackPathLength;
	}

	static bool SuccessorComparator_byPathLength_noFilter(const SuccessorData &a, const SuccessorData &b){
		return a._backtrackPathLength_noFilter > b._backtrackPathLength_noFilter;
	}

	static bool SuccessorComparator_byDistanceFromSource(const SuccessorData &a, const SuccessorData &b){
		return a._path.size() < b._path.size();
	}

	static void clampPrevNodes(vector<u_int32_t>& prevNodes, vector<UnitigData>& _unitigDatas){

		//cout << "clamp" << endl;
		if(prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1]);
			UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(prevNodes.size() == 0) break;
				u_int32_t nodeIndex = prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(prevNodes[0]);

				const UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				//cout << nodeName_current << " " << nodeName_first << " " << Utils::shareAnyRead(unitigData_current, unitigData_first) << " " << unitigData_current._readIndexes.size() << " " << unitigData_first._readIndexes.size() << endl;
				if(Utils::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					prevNodes.erase(prevNodes.begin());
				}

			}
		}

	}



};




















class UnitigExpanderPath{

	GraphSimplify* _graph;
	vector<u_int32_t> _prevNodes;
	vector<UnitigData>& _unitigDatas;
	bool _forward;

	vector<u_int32_t> _nodePath;
	unordered_set<u_int32_t> _visitedNodes;

	void UnitigExpanderPath(GraphSimplify* graph, const vector<UnitigData>& unitigDatas, bool forward) : _unitigDatas(unitigDatas){

		_graph = graph;
		_forward = forward;
	}

	bool execute(u_int32_t source_nodeIndex){

		
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(source_nodeIndex);

		cout << "Start expand: " << nodeName << endl;


		u_int32_t current_nodeIndex = source_nodeIndex;
		


		binNode(current_nodeIndex);
		_visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));
		
		u_int32_t lastNodeIndex = current_nodeIndex;


		
		while(true){

			PathExplorer pathExplorer(graph, _prevNodes, _unitigDatas, _forward, _visitedNodes);
			
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, resultType, nextNodes, true);
			
			if(resultType == 2){ //Banch solved
				cout << "Path size: " << pathData.nodePath.size() << endl;
			}


			if(current_nodeIndex == -1) return false;
			
			for(size_t i=1; i<nextNodes[0]._path.size(); i++){

				
				current_nodeIndex = nextNodes[0]._path[i];


				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(lastNodeIndex, current_nodeIndex);

				if(current_nodeIndex == source_nodeIndex){ //Path complete
					_nodePath.pop_back(); //if the path is solved, the source node exist as first and last element,thus we remove the last one

					//cout << "Path complete! (" << _nodePath.size() << ")" << endl;
					return true; 
				}

				
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
				//cout << "Add node: " << BiGraph::nodeToString(current_nodeIndex) << " " << (visitedNodes.find(nodeName) != visitedNodes.end()) << " " << pathData.nodePath.size() << " [Ab: " << currentAbundance << "]" << " " << "[Cutoff: " << PathExplorer::computeAbundanceCutoff(currentAbundance, 0, assemblyState._cutoffType) << "]";
				cout << "Add node: " << nodeName << endl;
				
				//if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
				//	cout << "    " << _evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName];
				//}
				//else{
				//	cout << "unknown pos" << endl;
				//	//getchar();
				//}

				
				
				//cout << endl;

				binNode(current_nodeIndex);
				

				lastNodeIndex = current_nodeIndex;
			}

		}


	}


	void binNode(u_int32_t nodeIndex){


		_prevNodes.push_back(nodeIndex);
		clampPrevNodes();

		_nodePath.push_back(nodeIndex);

		_visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		
		//u_int32_t current_unitigIndex = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		//file_groundTruth << current_unitigIndex << "," << pathIndex << endl;		
		
		//_binnedNodes.insert(current_unitigIndex);
		
	}

	void clampPrevNodes(){

		//cout << "clamp" << endl;
		if(_prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(_prevNodes[_prevNodes.size()-1]);
			const UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(_prevNodes.size() == 0) break;
				u_int32_t nodeIndex = _prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(_prevNodes[0]);

				const UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				//cout << nodeName_current << " " << nodeName_first << " " << Utils::shareAnyRead(unitigData_current, unitigData_first) << " " << unitigData_current._readIndexes.size() << " " << unitigData_first._readIndexes.size() << endl;
				if(Utils::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					_prevNodes.erase(_prevNodes.begin());
				}

			}
		}

	}

};

class UnitigExpander{

	GraphSimplify* _graph;
	vector<u_int32_t> _prevNodes;
	vector<UnitigData>& _unitigDatas;


	void UnitigExpander(GraphSimplify* graph, const vector<UnitigData>& unitigDatas) : _unitigDatas(unitigDatas){

		_graph = graph;
	}

	bool execute(u_int32_t source_nodeIndex){
		//Expand right
		//Expand left
	}

};

#endif 






#ifndef MDBG_METAG_ASSEMBLY
#define MDBG_METAG_ASSEMBLY

#include "Commons.hpp"

#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>
//#include "graph/Graph.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"


struct SuccessorData{
	u_int32_t _nodeIndex;
	u_int32_t _abundance;
	u_int32_t _sourceAbundance;
	u_int32_t _prevRank;
	u_int32_t _prevUnitigRank;
	bool _prevRankFinished;
	//vector<u_int32_t> _processedNodeIndex;
};






class PathExplorer{

public: 

	//u_int32_t _index;
	//unordered_set<DbgEdge, hash_pair> isEdgeVisited;
	vector<u_int32_t> _prevNodes;
	u_int32_t _source_abundance;
	u_int32_t _source_nodeIndex;
	u_int32_t _start_nodeIndex;
	float _abundanceCutoff_min;
	unordered_set<u_int32_t>& _visitedNodes;
	unordered_map<u_int32_t, bool>& _isNodeImproved;

	vector<u_int32_t> _exploredNodes;

	PathExplorer(const vector<u_int32_t>& prevNodes, u_int32_t source_abundance, u_int32_t source_nodeIndex, u_int32_t start_nodeIndex, float abundanceCutoff_min, unordered_set<u_int32_t>& visitedNodes, unordered_map<u_int32_t, bool>& isNodeImproved) : _visitedNodes(visitedNodes), _isNodeImproved(isNodeImproved){
		_prevNodes = prevNodes;
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_start_nodeIndex = start_nodeIndex;
		_abundanceCutoff_min = abundanceCutoff_min;
	}

	u_int32_t getNextNode(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int32_t& resultType){

		resultType = 0;
		if(currentDepth > 5) return -1;

		//u_int64_t iter = 0;
		bool orient_dummy = false;
		vector<u_int32_t> successors;

		//u_int32_t current_nodeIndex = _start_nodeIndex;

		u_int32_t current_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
		u_int32_t current_abundance = graph->_nodeAbundances[current_nodeName]; //_unitigDatas[current_unitigIndex]._meanAbundance;

		//if(_iter > 10000) return;

		//cout << "----------- " << iter << endl;
		//adjNode* node = graph->_nodes[utg_nodeIndex];
		vector<SuccessorData> data_successors;

		//if(!canExplorePath){
		//	cout << "HAAAAA " << current_nodeName << " " << graph->_graphSuccessors->nodeIndex_to_nodeName(_source_nodeIndex, orient_dummy) << endl;
		//}

		//cout << _prevNodes.size() << endl;
		//if(_prevNodes.size() > 1 && current_nodeIndex == _source_nodeIndex){
		//	cout << "Path complete! " << endl;
		//	//_pathDatas.push_back(pathData);
		//	return -2;
		//}
		//else if(iter > maxIter){
		//	return -1;
		//}


		if(forward){
			graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}
		else{
			graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}

		//bool isBranchingNode = successors.size() > 1;
		


		if(currentDepth == 0){
			if(successors.size() > 1){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "----------- " << endl;
				cout << "Nb prev nodes: " << _prevNodes.size() << endl;
				for(u_int32_t utg_n : successors){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
				
				}
			}
		}

		for(u_int32_t utg_n : successors){

			//u_int64_t utg_n = node->val;


			u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
			u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

			if(successors.size() > 1){
				//if(currentDepth == 0){
				if(isPathAlreadyExplored(utg_n, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 20)){
					//cout << "Already explored: " << current_nodeName << " " << successor_nodeName << endl;
					continue;
				}
				//}
			}


			SuccessorData successor = {utg_n, successor_abundance, 0, 0, false};
			successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
			data_successors.push_back(successor);

		}
		
		if(data_successors.size() == 0){
			if(currentDepth == 0){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "No successors" << endl;
			}
			return -1;
		}
		else if(data_successors.size() == 1){
			
			current_nodeIndex = data_successors[0]._nodeIndex;

			/*
			//Anti-bug temp
			if(_isNodeImproved.find(current_nodeIndex) != _isNodeImproved.end()){
				if(!_isNodeImproved[current_nodeIndex]){
					return -1;
				}
			}

			//updateNodeChosen(current_nodeIndex, currentDepth);
			*/
			return current_nodeIndex;
		}
		else{


			/*
			//-------------------------------------------------------------------------------
			//Solve multiple simple cycle
			vector<SuccessorData> successors_nonVisited;
			for(SuccessorData& successor : data_successors){
				if(isPathAlreadyExplored(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0)) continue;	
				successors_nonVisited.push_back(successor);
			}

			if(successors_nonVisited.size() == 1){
				for(SuccessorData& successor : successors_nonVisited){
					if(isSmallCycle(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
						return successor._nodeIndex;
					}
				}	
			}
	
			//-------------------------------------------------------------------------------
			*/

			resultType = 2;



			//Print debug prev ranks
			vector<SuccessorData> dummy;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, dummy, 0, true, true);



			vector<SuccessorData> successors_bestPrevUnitigRank_cutoff;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_cutoff, _abundanceCutoff_min/2, true, false);
			if(successors_bestPrevUnitigRank_cutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_cutoff[0]._nodeIndex;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 1: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				return current_nodeIndex;
			}

			vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_noCutoff, 0, true, false);
			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndex;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 2: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				return current_nodeIndex;
			}

			vector<SuccessorData> successors_bestPrevRank;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevRank, _abundanceCutoff_min/2, false, false);
			if(successors_bestPrevRank.size() == 1){
				current_nodeIndex = successors_bestPrevRank[0]._nodeIndex;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 3: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				return current_nodeIndex;
			}

			/*
			u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

			u_int32_t prevRank = 0;
			u_int32_t prevRank_unitig = 0;

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
				u_int32_t prev_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(prev_nodeIndex, orient_dummy);

				//cout << prevIndex << " " << prev_nodeIndex << endl;
				//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
				if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
					prevRank_unitig += 1;
					currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
				}

				//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
				//string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

				for(SuccessorData& successor : data_successors){
					if(successor._prevRankFinished){
						//str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
						continue;
					}

					//if(std::find(successor._processedNodeIndex.begin(), successor._processedNodeIndex.end(), prev_nodeIndex) != successor._processedNodeIndex.end()){
					//	successor._prevRankFinished = true;
					//	str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
					//	continue;
					//} 

					u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(successor._nodeIndex, orient_dummy);
					
					u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					//if(nbSharedReads > _abundanceCutoff_min/2){
					//if(nbSharedReads > successor._abundance/5){
					//if(nbSharedReads > 0){
					if(nbSharedReads > _abundanceCutoff_min/2){
						successor._prevUnitigRank = prevRank_unitig; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						//successor._processedNodeIndex.push_back(prev_nodeIndex);
					}
					else{
						successor._prevRankFinished = true;
					}


					//str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);

				}
				//cout << canExplorePath << endl;
				//cout << current_nodeIndex << " "  <<  _source_nodeIndex << endl;

				//if(currentDepth == 0){
				//	for(size_t i=0; i<currentDepth; i++) cout << "  ";
				//	cout << str_debug << endl;
				//}

				prevRank += 1;

				//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


			}


			
			//Check unitig diff
			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevUnitigRank);
			u_int32_t maxPrevUnitigRank = data_successors[0]._prevUnitigRank;
			vector<SuccessorData> successors_bestPrevUnitigRank;
			for(SuccessorData& successor : data_successors){
				if(successor._prevUnitigRank == maxPrevUnitigRank){
					successors_bestPrevUnitigRank.push_back(successor);
				}
			}
			if(successors_bestPrevUnitigRank.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank[0]._nodeIndex;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				return current_nodeIndex;
			}
			

			//Check node diff
			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevRank);
			u_int32_t maxPrevRank = data_successors[0]._prevRank;
			vector<SuccessorData> successors_bestPrevRank;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRank > maxPrevRank-5){
					successors_bestPrevRank.push_back(successor);
				}
			}


			if(successors_bestPrevRank.size() == 1){


				//DbgEdge edge = {current_nodeIndex, successors_bestPrevRank[0]._nodeIndex};
				//edge = edge.normalize();
				//isEdgeVisited.insert(edge);


				current_nodeIndex = successors_bestPrevRank[0]._nodeIndex;
				//nodeExplored(current_nodeIndex, graph);

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				return current_nodeIndex;
			}
			else{
			*/
				
			//cout << successors_bestPrevRank.size() << endl;
			//cout << "check bubble" << endl;
			if(successors_bestPrevUnitigRank_noCutoff.size() >= 2){
				bool isBubble = true;
				for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){
					if(!graph->_isBubble[successor._nodeIndex]){
						isBubble = false;
					}
				}
				if(isBubble){
					//if(currentDepth == 0){
						cout << "Take bubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndex) << endl;
					//}
					//exit(1);
					updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndex, currentDepth);
					return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndex;
				}

			}


			
			if(currentDepth == 0){

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Check simple cycle" << endl;
				}

				for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){

					if(isSmallCycle(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
						updateNodeChosen(successor._nodeIndex, currentDepth);
						return successor._nodeIndex;
					}
					
				}



				for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
					u_int32_t to_nodeIndex = successors_bestPrevUnitigRank_noCutoff[i]._nodeIndex;
					for(size_t j=0; j<successors_bestPrevUnitigRank_noCutoff.size(); j++){
						if(i == j) continue;
						
						u_int32_t from_nodeIndex = successors_bestPrevUnitigRank_noCutoff[j]._nodeIndex;

						if(currentDepth == 0){
							for(size_t i=0; i<currentDepth; i++) cout << "  ";
							cout << "Check simple cycle from: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << "    to:    " << graph->_graphSuccessors->nodeToString(to_nodeIndex) << endl;
						}

						if(isSmallCycle(from_nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
							//exit(1);
							updateNodeChosen(from_nodeIndex, currentDepth);
							return from_nodeIndex;
						}

					}	
				}
			}
			

			return -1;
				
			//}

		}


	}

	void computeBestSuccessors_byUnitigRank(GraphSimplify* graph, const vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors, int currentDepth, vector<SuccessorData>& successors_bestPrevUnitigRank, float abundanceCutoff_min, bool useUnitigRank, bool print_debug){


		for(SuccessorData& successor : data_successors){
			successor._prevRankFinished = false;
			successor._prevUnitigRank = 0;
			successor._prevRank = 0;
		}

		successors_bestPrevUnitigRank.clear();

		u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

		u_int32_t prevRank = 0;
		u_int32_t prevRank_unitig = 0;

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

			//cout << prevIndex << " " << prev_nodeIndex << endl;
			//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
			if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
				prevRank_unitig += 1;
				currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
			}

			//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
			string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

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

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(successor._nodeIndex);
				
				u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//if(nbSharedReads > _abundanceCutoff_min/2){
				//if(nbSharedReads > successor._abundance/5){
				//if(nbSharedReads > 0){
				if(nbSharedReads > abundanceCutoff_min){
					successor._prevUnitigRank = prevRank_unitig; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
					successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
					//successor._processedNodeIndex.push_back(prev_nodeIndex);
				}
				else{
					successor._prevRankFinished = true;
				}

				//if(nbSharedReads == 0) continue;
				//if(nbSharedReads > 0 && nbSharedReads >=  data_successors[j]/5){
				//	successor_foundPath[j] = true;
				//	foundPath = true;
				//}

				str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
				//cout << "    " << prevRank << ": " <<  graph->nodeToString(current_nodeIndex) << "     " << graph->nodeToString(successor._nodeIndex)  << ": " << nbSharedReads;
				
				//if(foundPath) break;
			}
			//cout << canExplorePath << endl;
			//cout << current_nodeIndex << " "  <<  _source_nodeIndex << endl;

			if(print_debug){
				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << str_debug << endl;
				}
			}


			prevRank += 1;

			//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


		}

		if(useUnitigRank){
			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevUnitigRank);
			u_int32_t maxPrevUnitigRank = data_successors[0]._prevUnitigRank;
			
			for(SuccessorData& successor : data_successors){
				if(successor._prevUnitigRank == maxPrevUnitigRank){
					successors_bestPrevUnitigRank.push_back(successor);
				}
			}
		}
		else{

			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevRank);
			u_int32_t maxPrevRank = data_successors[0]._prevRank;
			//vector<SuccessorData> successors_bestPrevRank;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRank > maxPrevRank-5){
					successors_bestPrevUnitigRank.push_back(successor);
				}
			}

		}


	}

	void updateNodeChosen(u_int32_t nodeIndex, u_int64_t currentDepth){
		if(currentDepth != 0) return;

		if(_visitedNodes.find(nodeIndex) == _visitedNodes.end()){
			for(auto& it : _isNodeImproved){
				_isNodeImproved[it.first] = true;
			}
		}

		_isNodeImproved[nodeIndex] = false;
	}

	bool isPathAlreadyExplored(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t maxIter){
		
		//bool lala = true;

		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		}

		//auto it = _isNodeImproved.find(current_nodeIndex);
		if(_isNodeImproved.find(current_nodeIndex) != _isNodeImproved.end()){
			if(!_isNodeImproved[current_nodeIndex]){
				if(currentDepth == 1){
					cout << " Yes (Is not improved)" << endl;
				}
				return true;
			}
		}

		if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
			if(currentDepth == 1){
				cout << " No" << endl;
			}
			//lala = false;
			return false;
		}

		if(maxIter == 0){
			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				if(currentDepth == 1){
					cout << " No" << endl;
				}
				return false;
			}
			else{
				if(currentDepth == 1){
					cout << " Yes" << endl;
				}
				return true;
			}
		}

		//u_int32_t source_nodeIndex = current_nodeIndex;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved);

		//u_int64_t maxIter = 10;


		u_int64_t iter = 0;
		//u_int32_t current_nodeIndex = _start_nodeIndex;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, currentDepth, resultType);
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(current_nodeIndex == source_nodeIndex) break;
			if(current_nodeIndex == -1){ //dead end or multiple braching path
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}

			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}

			pathExplorer.nodeExplored(current_nodeIndex, graph);

			if(iter >= maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " Yes" << endl;
		}

		//if(!lala) return false;
		return true;
	}

	/*
	bool isSimpleCycle(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){
		
		//bool lala = true;

		cout << "\tIs simple cycle: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";

		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes);

		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, false);
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(current_nodeIndex == source_nodeIndex){
				cout << " Yes" << endl;
				return true;
			}
			if(current_nodeIndex == -1){ //dead end or multiple braching path
				cout << " No" << endl;
				//lala = false;
				return false;
			}



			pathExplorer.nodeExplored(current_nodeIndex, graph);

			if(iter > maxIter) break;

			iter += 1;
		}
		
		cout << " No" << endl;

		//if(!lala) return false;
		return false;
	}
	*/

	bool isSmallCycle(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		//bool lala = true;


		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is simple cycle: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << " ?    ";
		}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, from_nodeIndex, from_nodeIndex, _abundanceCutoff_min, visitedNodes, _isNodeImproved);


		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(from_nodeIndex, graph);
		pathExplorer._visitedNodes.insert(from_nodeIndex);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			from_nodeIndex = pathExplorer.getNextNode(from_nodeIndex, graph, _unitigDatas, forward, currentDepth, resultType);
			//cout << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;
			pathExplorer._visitedNodes.insert(from_nodeIndex);

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

	/*
	bool isPathAlreadyExplored(u_int32_t start_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, u_int64_t maxIter, bool forward){

		//bool dummy = false;
		//if(graph->_graphSuccessors->nodeIndex_to_nodeName(start_nodeIndex, dummy) == 6307){
		//	cout << "hey" << endl;
		//	exit(1);
		//}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, start_nodeIndex, start_nodeIndex, _abundanceCutoff_min, _visitedNodes);
		pathExplorer.extend(graph, _unitigDatas, maxIter, forward, false);
		
		//cout << "explored" << endl;
		for(u_int32_t nodeIndex : pathExplorer._exploredNodes){
			//cout << "\t\t" << (_visitedNodes.find(nodeIndex) == _visitedNodes.end()) << endl;

			if(_visitedNodes.find(nodeIndex) == _visitedNodes.end()){
				cout << "No" << endl;
				return false;
			}
			//binNode(nodeIndex, pathData.prevNodes, graph, pathData._index);
			//visitedNodes.insert(nodeIndex);
		}

		cout << "Yes" << endl;
		//cout << "Is explored" << endl;
		return true;
	}
	*/

	void nodeExplored(u_int32_t nodeIndex, GraphSimplify* graph){
		
		bool orient_dummy;
		_prevNodes.push_back(nodeIndex);
		_exploredNodes.push_back(nodeIndex);

		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);


		//_binnedNodes.insert(current_unitigIndex);
		//cout << "Node explored: " << graph->_graphSuccessors->nodeToString(nodeIndex) << " " << graph->_nodeAbundances[nodeName]  << endl;

		//_nbVisitedTimes[current_unitigIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}

	static u_int64_t computeSharedReads(const UnitigData& utg1, const UnitigData& utg2){

		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){
			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t collectSharedReads(const UnitigData& utg1, const UnitigData& utg2, vector<u_int64_t>& sharedReads){

		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}

		sharedReads.clear();

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){
			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				sharedReads.push_back(utg1._readIndexes[i]);
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}


	static bool shareAnyRead(const UnitigData& utg1, const UnitigData& utg2){

		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}

		size_t i=0;
		size_t j=0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){

			//cout << i << " " << j << endl;
			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				return true;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return false;
	}

	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

	static bool SuccessorComparator_byPrevUnitigRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevUnitigRank > b._prevUnitigRank;
	}
};




class Assembly : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	vector<UnitigData> _unitigDatas;

	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	gzFile _outputContigFile;

	Assembly(): Tool ("asm"){


		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename", true));
		//getParser()->push_back (new OptionOneParam (STR_MINIM_SIZE, "minimizer length", false, "16"));
		//getParser()->push_back (new OptionOneParam (STR_KMINMER_SIZE, "k-min-mer length", false, "3"));
		//getParser()->push_back (new OptionOneParam (STR_DENSITY, "density of minimizers", false, "0.005"));
		//getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", false, ""));
		getParser()->push_back (new OptionOneParam (STR_INPUT_TRUTH, "", false, ""));

	}

	~Assembly(){

	}

	void execute (){
		parseArgs();


		
		
		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");

		execute_assembly();
		
		gzclose(_outputContigFile);
	}

	void parseArgs(){
		_inputDir = getInput()->getStr(STR_INPUT_DIR);
		_truthInputFilename = getInput()->get(STR_INPUT_TRUTH) ? getInput()->getStr(STR_INPUT_TRUTH) : "";

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_filename_readMinimizers = _inputDir + "/read_data.gz";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";

	}


	void execute_assembly(){

		//size_t globalAbundanceCutoff_min = 3;

		string gfa_filename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";



		/*
		cout << "Simplifying graph" << endl;
		GraphSimplify* graphSimplify2 = new GraphSimplify(gfa_filename, _inputDir);
		graphSimplify2->execute(5);
		//graphSimplify2->clear(1);
		//graphSimplify2->compact();
		graphSimplify2->extractComponent(2246262);
		*/




		
		cout << gfa_filename << endl;
		MDBG* mdbg = new MDBG(_kminmerSize);
		mdbg->load(mdbg_filename);

		cout << "Nb nodes: " <<  mdbg->_dbg_nodes.size() << endl;

		if(_truthInputFilename != ""){
			extract_truth_kminmers(mdbg);
		}




		//vector<u_int32_t> unitigLengths;
		//vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
		//GfaParser gfaParser;
		//BiGraph* graph = gfaParser.createBiGraph_lol(gfa_filename);
		//AdjGraph* graph = gfaParser.createGraph_lol(gfa_filename);
		
		//cout << "Nb nodes: " << graph->_nbNodes << endl;
		//cout << "Nb edges: " << graph->_nbEdges << endl;








		//unordered_set<u_int64_t> filteredMinimizers;



		//load_read_compositions();



		//vector<UnitigData> unitigDatas;
		_unitigDatas.resize(mdbg->_dbg_nodes.size());

		//for(auto it : mdbg->_dbg_nodes){

			//const KmerVec& vec = it.first;

			//u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			//u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			//if(unitigIndex == -1) continue;

			//unitigDatas[unitigIndex]._nbKminmers += 1;
			//unitigDatas[unitigIndex]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;
			//_unitigDatas[kminmer_index]._meanAbundance = mdbg->_dbg_nodes[vec]._abundance;
		//}
		/*
		for(size_t i=0; i<graph->_nbNodes; i++){
			//cout << i << endl;
			//unitigDatas.push_back({i, {}});
			//unitigDatas.push_back({i, {}, {}, 0});
			unitigDatas[i]._compositionMean.resize(_compositionManager->_compositionVectorSize);
		}*/
		removeUnsupportedEdges(mdbg, gfa_filename, gfa_filename_noUnsupportedEdges);
		cout << "done" << endl;
		//delete graph; //NEED THIS FOR DEBUG TO EXTRACT HIFIASM SUB GRAPH
		//BiGraph* graphBi_successors = gfaParser.createBiGraph_lol(gfa_filename_noUnsupportedEdges, true);
		//BiGraph* graphBi_predecessors = gfaParser.createBiGraph_lol(gfa_filename_noUnsupportedEdges, false);
		//cout << "Nb nodes: " << (graphBi_successors->_nbNodes/2) << endl;
		//cout << "Nb edges: " << graphBi_successors->_nbEdges << endl;



		/*

		string command = "python3 ~/workspace/scripts/assembly/simplify_gfa.py " + gfa_filename_noUnsupportedEdges + " -out " + gfa_filename_unitigs;
		cout << command << endl;
		int ret = system(command.c_str());
		if(ret != 0){
			cout << "ERROR IN GFA TOOLS" << endl;
			exit(ret);
		}

		
		//vector<u_int32_t> unitigLengths;
		_node_to_unitig.resize(mdbg->_dbg_nodes.size(), -1);
		GfaParser::getNodeToUnitig(gfa_filename_unitigs, _node_to_unitig);
		
		*/

		//vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		//UnitigGraph* graph_unitig = gfaParser.createGraph(gfa_filename_unitigs, _node_to_unitig, unitigLengths);
		//delete graph_unitig;



		/*
		cout << "---------------------" << endl;
		for(u_int32_t r1 : unitigDatas[857640]._readIndexes){
			for(u_int32_t r2 : unitigDatas[857640]._readIndexes){
				//float tnf_dist = computeDistanceTNF(_readDatas[r1], _readDatas[r2]);
				//cout << tnf_dist << " " << readIndex_1 << " " << readIndex_2 << " " << utg << " " << utg_n << endl;
				cout << euclidianDistance(_readDatas[r1]._composition, _readDatas[r2]._composition) << " " << r1 << " " << r2 << " " << endl;
				//cout << euclidianDistance(_readDatas[readIndex_1]._composition, _readDatas[readIndex_2]._composition) << " " << utg << " " << utg_n << endl;
			}
		}
		cout << "---------------------" << endl;
		*/


		/*
		for(UnitigData& data : unitigDatas){
			for(size_t i=0; i<data._compositionMean.size(); i++){
				data._compositionMean[i] /= data._compositionNb;
			}

		}
		*/


		/*
		//Simulation
		ofstream file_groundTruth_hifiasm_position(_inputDir + "/groundtruth_hifiasm_position.csv");
		ofstream file_groundTruth_hifiasm(_inputDir + "/groundtruth_hifiasm.csv");
		file_groundTruth_hifiasm << "Name,Colour" << endl;
		file_groundTruth_hifiasm_position << "Name,Order" << endl;

		//unordered_set<u_int32_t> groundTruth_kminmers;
		int founded = 0;
		for(auto it : mdbg->_dbg_nodes){
			const KmerVec& vec = it.first;

			//vec = vec.normalize();
			if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;

			founded += 1;
			//groundTruth_kminmers.insert(it.second._index);

			file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
			file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

			//vector<u_int32_t> neighbors;
			//graph->collectNeighbors(it.second._index, 100, neighbors, 100, visitedNodes);
			//for(u_int32_t nn : neighbors){
				//cout << nn << endl;
				//groundTruth_kminmers.insert(nn);
			//}
			//cout << "n: " << neighbors.size() << endl;
			//cout << n << endl;


			//cout << groundTruth_kminmers.size() << endl;
		}
		//cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
		cout << "Found: " << founded << endl;
		//gfaParser.rewriteGfa_withNodes(gfa_filename, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_kminmers);
		file_groundTruth_hifiasm.close();
		file_groundTruth_hifiasm_position.close();
		*/

		
		/*
		cout << "Collecting truth kminmers" << endl;
		//ofstream file_groundTruth_hifiasm_position(_outputDir + "/groundtruth_hifiasm_position.csv");
		ofstream file_groundTruth_hifiasm(_outputDir + "/groundtruth_hifiasm.csv");
		file_groundTruth_hifiasm << "Name,Colour" << endl;
		//file_groundTruth_hifiasm_position << "Name,Colour" << endl;

		unordered_set<u_int32_t> visitedNodes;
		for(auto it : mdbg->_dbg_nodes){
			if(_evaluation_hifiasmGroundTruth.find(it.first) == _evaluation_hifiasmGroundTruth.end()) continue;
			visitedNodes.insert(it.second._index);
		}

		unordered_set<u_int32_t> groundTruth_kminmers;
		int founded = 0;
		for(auto it : mdbg->_dbg_nodes){

			const KmerVec& vec = it.first;

			//vec = vec.normalize();
			if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;

			founded += 1;
			groundTruth_kminmers.insert(it.second._index);

			file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
			//file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

			vector<u_int32_t> neighbors;
			graph->collectNeighbors(it.second._index, 100, neighbors, 100, visitedNodes);
			for(u_int32_t nn : neighbors){
				//cout << nn << endl;
				groundTruth_kminmers.insert(nn);
			}
			//cout << "n: " << neighbors.size() << endl;
			//cout << n << endl;


			cout << groundTruth_kminmers.size() << endl;
		}
		cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
		cout << "Founded: " << founded << endl;
		gfaParser.rewriteGfa_withNodes(gfa_filename_noUnsupportedEdges, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_kminmers);
		file_groundTruth_hifiasm.close();
		//file_groundTruth_hifiasm_position.close();
		*/
		








		/*
		computeNodeAbundance(mdbg, gfa_filename_noUnsupportedEdges);
		getUnitigLengths(gfa_filename_unitigs);
		*/


		cout << "Simplifying graph" << endl;
		GraphSimplify* graphSimplify = new GraphSimplify(gfa_filename_noUnsupportedEdges, _inputDir);
		cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		graphSimplify->execute(5);

		
		//cout << "remettre cleaning" << endl;
		//graphSimplify->compact();
		
		
		//graphSimplify->extractComponent(1224833);





		//exit(1);
		
		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm.csv");
		file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		// mini genome test
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(16, true), graphSimplify->nodeIndex_to_unitig(4750)._abundance, graphSimplify, 0);

		// 201
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(4750, true), graphSimplify->nodeIndex_to_unitig(4750)._abundance, graphSimplify, 0);

  
		//562 (ecoli)
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(1224833, false), 35, graphSimplify, 0);
		//exit(1);


		//return;



		//file_groundTruth.close();

		auto rng = std::default_random_engine {};
		unordered_set<u_int32_t> binNodes;
		
		
		//unordered_set<u_int32_t> visitedNodes;
		vector<u_int32_t> startingNodesIndex;

		vector<u_int32_t> unitigLength_cutoffs = {100000, 50000, 30000, 10000};
		size_t binIndex=0;

		vector<vector<UnitigLength>> startingUnitigs;
		startingUnitigs.resize(unitigLength_cutoffs.size());

		unordered_set<u_int32_t> usedNodeNames;
		for(size_t i=0; i<unitigLength_cutoffs.size(); i++){
				
			u_int32_t unitigLength_cutoff_min = unitigLength_cutoffs[i];
			u_int32_t unitigLength_cutoff_max = -1;
			if(i > 0){
				unitigLength_cutoff_max = unitigLength_cutoffs[i-1];
			}

			vector<UnitigLength>& unitigLengths = startingUnitigs[i];

			for(const Unitig& unitig : graphSimplify->_unitigs){
				//if(unitig._index % 2 == 1) continue;
				if(unitig._length < unitigLength_cutoff_min) continue;
				if(unitig._length > unitigLength_cutoff_max) continue;
				if(unitig._abundance < 24) continue; //200


				vector<u_int32_t> nodes;
				//Unitig& unitig = graphSimplify->_unitigs[unitig._index];

				//unitigLength._abundance = unitig._abundance;

				//float abundanceCutoff_min = computeAbundanceCutoff_min(unitigLength._abundance);
				//if(unitigLength._abundance < 100) continue;


				
				graphSimplify->getUnitigNodes(unitig, nodes);

				/*
				cout << endl << "----" << endl;
				for(u_int32_t nIndex : nodes){
					u_int32_t nName = BiGraph::nodeIndex_to_nodeName(nIndex);
					cout << nName << " " << graphSimplify->_nodeAbundances[nName] << endl;
				}*/




				//std::shuffle(nodes.begin(), nodes.end(), rng);
				
				bool isValid = true;
				u_int32_t maxAbundance = 0;
				u_int32_t startNodeIndex = -1;
				for(u_int32_t nodeIndex : nodes){

					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

					if(usedNodeNames.find(nodeName) != usedNodeNames.end()){
						isValid = false;
						break;
					}


					usedNodeNames.insert(nodeName);
					//if(graphSimplify->_isBubble[nodeIndex]) continue;

					u_int32_t abundance = graphSimplify->_nodeAbundances[nodeName];
					if(abundance > maxAbundance){
						maxAbundance = abundance;
						startNodeIndex = nodeIndex;
					}
					//startNodeIndex = nodeIndex;
					//unitigLength._startNodeIndex = nodeIndex;
					//break;
				}

				if(!isValid) continue;

				/*
				if(BiGraph::nodeIndex_to_nodeName(startNodeIndex) == 2246258){
					cout << "-" << endl;
					cout << "OMG" << endl;
					cout << unitig._abundance << endl;

					for(u_int32_t nIndex : nodes){
						u_int32_t nName = BiGraph::nodeIndex_to_nodeName(nIndex);
						cout << nName << " " << graphSimplify->_nodeAbundances[nName] << endl;
					}
					exit(1);
				}
				*/

				unitigLengths.push_back({unitig._length, unitig._abundance, startNodeIndex});
			}

			std::sort(unitigLengths.begin(), unitigLengths.end(), UnitigComparator_ByAbundance);


			cout << "Length cutoff: " << unitigLength_cutoff_min << " " << unitigLength_cutoff_max << endl;
			cout << "Nb starting unitigs: " << unitigLengths.size() << endl;
		}
		usedNodeNames.clear();


		for(const vector<UnitigLength>& startingUnitig : startingUnitigs){	

			cout << "---- New cutoff ---- " << endl;

			for(const UnitigLength& unitigLength : startingUnitig){

				cout << unitigLength._length << " " << unitigLength._abundance << " " << BiGraph::nodeIndex_to_nodeName(unitigLength._startNodeIndex) << endl;
				//cout << "Try asm: " << unitigLength._startNodeIndex << endl;
				//if(unitigLength._index % 2 == 1) continue;

				//if(unitigLength._length < 30000) continue;
				//loat abundanceCutoff_min = computeAbundanceCutoff_min(unitigLength._abundance);
				//if(unitigLength._abundance < 100) continue;
				//if(visitedNodes.find(unitigLength._index) != visitedNodes.end()) continue;


				//cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;

				u_int32_t nodeIndex = unitigLength._startNodeIndex;
				
				u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
				if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;

				//visitedNodes.insert(unitigLength._index);

				

				//vector<u_int32_t> nodes;
				//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
				//graphSimplify->getUnitigNodes(unitig, nodes);
				//for(u_int32_t node : nodes){
				//	visitedNodes.insert(node);
				//}
				/*
				cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;
				
				if(unitigLength._length < 10000) continue;

				//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
				//if(graphSimplify->_nodeToUnitig[unitig._startNode] == -1) continue;

				vector<u_int32_t> nodes;
				Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
				graphSimplify->getUnitigNodes(unitig, nodes);

				u_int32_t nodeIndex = nodes[rand() % nodes.size()];
				u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
				if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;

				float abundanceCutoff_min = computeAbundanceCutoff_min(nodeIndex, graphSimplify);
				if(abundanceCutoff_min < 30) continue;
				//u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(unitigLength.startNode, orient_dummy);
				//if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;
					
				*/

				//for(size_t n=0; n<graphSimplify->_nodeToUnitig.size(); n++){
				//	if(n%2 == 1) continue;


					//if(graphSimplify->_nodeToUnitig[n] == unitigLength._index){

						
						//file_groundTruth = ofstream(_outputDir + "/binning_results.csv");
			
						bool isContigValid = solveBin(nodeIndex, unitigLength._abundance, graphSimplify, binIndex);
						if(isContigValid){
							binIndex += 1;
						} 
						else{
							bool isContigValid = solveBin(graphSimplify->nodeIndex_toReverseDirection(nodeIndex), unitigLength._abundance, graphSimplify, binIndex);
							if(isContigValid) binIndex += 1;
						}

						//file_groundTruth.close();

						//exit(1);
						
						//break;
					//}
				//}
				
			}
		}

		cout << "Nb bins: " << binIndex << endl;


		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
	}

		
	
				/*
			bool orient_dummy;
			//size_t binIndex=0;
			for(UnitigLength& unitigLength : unitigLengths){
				//cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;
				
				//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
				//if(graphSimplify->_nodeToUnitig[unitig._startNode] == -1) continue;



				//u_int32_t nodeIndex = nodes[rand() % nodes.size()];
				//u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
				//if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;



				//unitigLength._startNodeIndex = nodeIndex;
				//startingNodesIndex.push_back(nodeIndex);
			}
			*/


			//for(size_t n=0; n<graphSimplify->_nodeToUnitig.size(); n++){
			//	u_int32_t unitigIndex = graphSimplify->_nodeToUnitig[n];
			//	graphSimplify->_nodeAbundances[n] = graphSimplify->_unitigs[unitigIndex]._abundance;
			//}

			//genome3
			//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(1869, true), graphSimplify, 0);
			
			//genome_201_50x
			//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(115862, false), graphSimplify, 0);
			
			//exit(1);

			//for(u_int32_t nodeIndex : startingNodesIndex){


	struct PathData{
		u_int32_t _index;
		unordered_set<DbgEdge, hash_pair> isEdgeVisited;
		vector<u_int32_t> nodePath;
		vector<u_int32_t> prevNodes;
		u_int32_t source_abundance;
		u_int32_t source_nodeIndex;
		u_int32_t source_nodeIndex_path;
		float _abundanceCutoff_min;
		unordered_map<u_int32_t, bool> _isNodeImproved;
	};

	void removeUnsupportedEdges(MDBG* mdbg, const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges){
		cout << "Removing unsupported edges" << endl;

		//GraphSimplify* graphSimplify = new GraphSimplify(gfaFilename, _inputDir);
		//graphSimplify->compact();


		gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
		ReadIndexType readIndex = 0;

		while(true){
			
			//cout << readIndex << endl;
			//"reprise: essayer avec gfatools unitigs"
			u_int16_t size;
			vector<u_int64_t> minimizers;
			gzread(file_readData, (char*)&size, sizeof(size));

			if(gzeof(file_readData)) break;
			
			minimizers.resize(size);
			gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> minimizersPos; 
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0);

			/*
			if(readIndex == 75286){
				cout << "-------------------------- " << size << endl;
				for(u_int64_t m : minimizers){
					cout << m << endl;
				}
			}
			else if(readIndex == 75888){
				cout << "-------------------------- " << size << endl;
				for(u_int64_t m : minimizers){
					cout << m << endl;
				}
			}*/

			vector<ReadIndexType> unitigIndexex;
		
			
			for(KmerVec& vec : kminmers){
				//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

				//cout << mdbg->_dbg_nodes[vec]._index << endl;
				if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

				




				u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
				UnitigData& unitigData = _unitigDatas[nodeName];
				unitigData._readIndexes.push_back(readIndex);
				//u_int32_t nodeIndex1 = BiGraph::nodeName_to_nodeIndex(nodeName, true);
				//u_int32_t nodeIndex2 = BiGraph::nodeName_to_nodeIndex(nodeName, false);
				//cout << nodeName << " " << nodeIndex << endl;
				//u_int32_t unitigIndex1 = graphSimplify->nodeIndex_to_unitigIndex(nodeIndex1);   //kminmer_index;//node_to_unitig[kminmer_index];
				//u_int32_t unitigIndex2 = graphSimplify->nodeIndex_to_unitigIndex(nodeIndex2);   //kminmer_index;//node_to_unitig[kminmer_index];
				//cout << unitigIndex << endl;
				//if(unitigIndex == -1) continue;

				//UnitigData& unitigData = _unitigDatas[nodeName];
				//unitigData._readIndexes.push_back(readIndex);

				//if(kminmer_index == 55479) cout << "AAAAA" << endl;
				//cout << kminmer_index << " " << unitigIndex << endl;
				//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex1) != unitigIndexex.end()) continue;
				//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex2) != unitigIndexex.end()) continue;

				//unitigIndexex.push_back(unitigIndex1);
				//if(unitigIndexex.size() >= 2) break;
			}

			/*
			if(unitigIndexex.size() >= 2){
				for(KmerVec& vec : kminmers){


					if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

					u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
					//u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
					//u_int32_t unitigIndex = graphSimplify->nodeIndex_to_unitigIndex(nodeIndex);   //kminmer_index;//node_to_unitig[kminmer_index];
					//if(unitigIndex == -1) continue;

					UnitigData& unitigData = _unitigDatas[nodeName];
					
					//if(unitigData._readIndexes_exists.find(readIndex) != unitigData._readIndexes_exists.end()) continue;
					//cout << unitigData._readIndexes.size() << endl;
					//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
						
					//unitigData._readIndexes_exists.insert(readIndex);
					unitigData._readIndexes.push_back(readIndex);
					
					//cout << unitigData._readIndexes.size() << " " << unitigLengths[unitigIndex] << endl;
					//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue
					//unitigIndexex.push_back(unitigIndex);

				}
			}*/


			
			readIndex += 1;
		}

		/*
		unordered_set<DbgEdge, hash_pair> unsupportedEdges;
        vector<u_int32_t> neighbors;

        for(size_t nodeIndex=0; nodeIndex< graphSimplify->_graphSuccessors->_nbNodes; nodeIndex++){
           
		    u_int32_t nodeName1 = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex);
			
			graphSimplify->getPredecessors(nodeIndex, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(!PathExplorer::shareAnyRead(_unitigDatas[nodeName1], _unitigDatas[nodeName2])){
					unsupportedEdges.insert({nodeName1, nodeName2});
					unsupportedEdges.insert({nodeName2, nodeName1});
				}
			}

			graphSimplify->getSuccessors(nodeIndex, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(!PathExplorer::shareAnyRead(_unitigDatas[nodeName1], _unitigDatas[nodeName2])){
					unsupportedEdges.insert({nodeName1, nodeName2});
					unsupportedEdges.insert({nodeName2, nodeName1});
				}
			}

        }
		*/

		/*
		unordered_set<DbgEdge, hash_pair> unsupportedEdges;
		
        vector<u_int32_t> neighbors;
		
        for(Unitig& unitig : graphSimplify->_unitigs){

			u_int32_t startNodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);
			graphSimplify->getPredecessors(unitig._startNode, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

				if(unsupportedEdges.find({startNodeName, nodeName}) != unsupportedEdges.end() || unsupportedEdges.find({nodeName, startNodeName}) != unsupportedEdges.end()) {	
					//node = node->next;
					continue;
				}

				if(PathExplorer::shareAnyRead(_unitigDatas[startNodeName], _unitigDatas[nodeName])){
					continue;
				}

				unsupportedEdges.insert({startNodeName, nodeName});
				unsupportedEdges.insert({nodeName, startNodeName});
				
			}

			graphSimplify->getSuccessors(unitig._startNode, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

				if(unsupportedEdges.find({startNodeName, nodeName}) != unsupportedEdges.end() || unsupportedEdges.find({nodeName, startNodeName}) != unsupportedEdges.end()) {	
					//node = node->next;
					continue;
				}

				if(PathExplorer::shareAnyRead(_unitigDatas[startNodeName], _unitigDatas[nodeName])){
					continue;
				}

				unsupportedEdges.insert({startNodeName, nodeName});
				unsupportedEdges.insert({nodeName, startNodeName});
				
			}

			u_int32_t endNodeName = BiGraph::nodeIndex_to_nodeName(unitig._endNode);
            graphSimplify->getSuccessors(unitig._endNode, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(!PathExplorer::shareAnyRead(_unitigDatas[endNodeName], _unitigDatas[nodeName])){
					unsupportedEdges.insert({endNodeName, nodeName});
					unsupportedEdges.insert({nodeName, endNodeName});
				}
			}
			
			
            graphSimplify->getPredecessors(unitig._endNode, 0, neighbors);
			for(u_int32_t nodeIndex : neighbors){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(!PathExplorer::shareAnyRead(_unitigDatas[endNodeName], _unitigDatas[nodeName])){
					unsupportedEdges.insert({endNodeName, nodeName});
					unsupportedEdges.insert({nodeName, endNodeName});
				}
			}

		}

		delete graphSimplify;
		*/
		

		unordered_set<DbgEdge, hash_pair> unsupportedEdges;
		
		AdjGraph* graph = GfaParser::createGraph_lol(gfaFilename);


		for(size_t utg=0; utg<graph->_nbNodes; utg++){

			//u_int32_t unitigIndex = graphSimplify->nodeIndex_to_unitigIndex(nodeIndex);
			//cout << utg << " " << graph->_nbNodes << endl;

			//int nbNeighbors = 0;
			//adjNode* node = graph->_nodes[utg];
			//while (node != nullptr) {
			//	nbNeighbors += 1;
			//}

			adjNode* node = graph->_nodes[utg];
			while (node != nullptr) {
				
				ReadIndexType utg_n = node->val;

				if(unsupportedEdges.find({utg, utg_n}) != unsupportedEdges.end() || unsupportedEdges.find({utg_n, utg}) != unsupportedEdges.end()) {	
					node = node->next;
					continue;
				}

				if(PathExplorer::shareAnyRead(_unitigDatas[utg], _unitigDatas[utg_n])){
					node = node->next;
					continue;
				}

				//if(shareAnyRead(unitigDatas[utg], unitigDatas[utg_n])){
				unsupportedEdges.insert({utg, utg_n});
				unsupportedEdges.insert({utg_n, utg});

				//if(nbNeighbors == 2){
				//	lala += 1;
				//}
				//}

				node = node->next;
			}

		}
		
		
		delete graph;
		

		cout << "Nb unsupported edges: " << unsupportedEdges.size() << endl;
		GfaParser::rewriteGfa_withoutEdges(gfaFilename, gfa_filename_noUnsupportedEdges, unsupportedEdges);
		//unsupportedEdges.clear();
		//exit(1);
		//exit(1);

	}
	//vector<PathData> _pathDatas;
	//vector<int32_t> _node_to_unitig;

	//unordered_map<u_int32_t, u_int32_t> _nodeLabel;
	//unordered_set<u_int32_t> _globalVisitedNodes;

	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

	static bool SuccessorComparator_byAbundance(const SuccessorData &a, const SuccessorData &b){
		return a._sourceAbundance < b._sourceAbundance;
	}

	static bool UnitigComparator_ByLength(const UnitigLength &a, const UnitigLength &b){
		return a._length > b._length;
	}

	static bool UnitigComparator_ByAbundance(const UnitigLength &a, const UnitigLength &b){
		return a._abundance > b._abundance;
	}

	size_t _iter;
	ofstream file_groundTruth;
	ofstream file_groundTruth_hifiasmContigs;

	/*
	void getSuccessors(u_int32_t nodeIndex, const PathData& pathData, BiGraph* graph, vector<u_int32_t>& successors){

		bool orient_dummy = false;
		successors.clear();

		adjNode* node = graph->_nodes[nodeIndex];

		while(node != nullptr){

			u_int64_t utg_n = node->val;
			u_int32_t nodeName = graph->nodeIndex_to_nodeName(utg_n, orient_dummy);
			//cout << unitigIndex << " " << _node_to_unitig[unitigIndex] << endl;
			
			u_int32_t unitigIndex = _node_to_unitig[nodeName];

			//cout << _nodeAbundances[nodeName] << endl;
			if(unitigIndex == -1){ //Cleaned
				node = node->next;
				continue;
			}

			if(_nodeAbundances[nodeName] < pathData._abundanceCutoff_min){ //Abundance min cutoff
				node = node->next;
				continue;
			}

			//cout << _nodeAbundances[nodeName] << " " <<  pathData._abundanceCutoff_min << endl;
			successors.push_back(node->val);

			node = node->next;
		}

	}
	*/

	unordered_set<u_int32_t> _binnedNodes;
	unordered_map<u_int32_t, u_int32_t> _nbVisitedTimes;

	void binNode(u_int32_t nodeIndex, vector<u_int32_t>& prevNodes, vector<u_int32_t>& nodePath, GraphSimplify* graph, u_int32_t pathIndex){

		bool orient_dummy;
		//u_int32_t utg_nodeIndex = nodeIndex; //successors[0]._nodeIndex;
		prevNodes.push_back(nodeIndex);

		if(prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1]);
			UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(prevNodes.size() == 0) break;
				u_int32_t nodeIndex = prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(prevNodes[0]);

				UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				if(PathExplorer::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					prevNodes.erase(prevNodes.begin());
				}

			}
		}


		//cout << "Prev nodes size: " << prevNodes.size() << endl;
		//cout << "Complete path size: " << nodePath.size() << endl;

		nodePath.push_back(nodeIndex);

		//_nodeLabel[nodeIndex] = pathIndex;
		u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		file_groundTruth << current_unitigIndex << "," << pathIndex << endl;		

		_binnedNodes.insert(current_unitigIndex);
		//cout << "Add node: " << graph->_graphSuccessors->nodeToString(nodeIndex) << " " << graph->_nodeAbundances[current_unitigIndex]  << endl;

		_nbVisitedTimes[nodeIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}

	float computeAbundanceCutoff_min(u_int32_t abundance){
		return abundance / 4.0;
	}

	bool solveBin(u_int32_t source_nodeIndex, u_int32_t abundance, GraphSimplify* graph, int pathIndex){


		bool orient_dummy = false;
		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(source_nodeIndex, orient_dummy);

		cout << endl << endl << endl << endl << endl << endl << endl << endl << "----- Start solve bin -----------------------------------------------------------------------------------------------" << endl;
		cout << "Source: " << nodeName << " " << graph->_graphSuccessors->nodeToString(source_nodeIndex) << endl;
		cout << "Source abundance: " << abundance << endl;

		//cout << graph->_nodeToUnitig[source_nodeIndex] << endl;
		//cout << graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance << endl;
		//bool orient_dummy;
		//u_int32_t utg_nodeIndex = source_nodeIndex; //graph->nodeName_to_nodeIndex(utg, orient_dummy);
		
		//u_int32_t utg_nodeIndex = graph->nodeName_to_nodeIndex(utg, orient_dummy);
		

		//float abundanceCutoff_min = graph->_nodeAbundances[nodeName] / 5.0;
		float abundanceCutoff_min = computeAbundanceCutoff_min(abundance);
		cout << "Abundance cutoff min: " << abundanceCutoff_min << endl;
		//if(abundanceCutoff_min < 30) return;

		cout << "Simplifying graph local" << endl;
		//graph->clear();
		graph->debug_writeGfaErrorfree(abundance, abundanceCutoff_min, source_nodeIndex);
		//graph->clear();
		//graph->execute(abundanceCutoff_min);
		//graph = graph->clone();
		//cout << "clone done" << endl;
		//graph->execute(abundanceCutoff_min);


		//u_int32_t source_unitigIndex = graph->nodeIndex_to_nodeName(utg, orient_dummy);
		if(graph->_nodeToUnitig[source_nodeIndex] == -1){ //actuellement ça peut arriver si l'assemblage demarre sur une tip
			cout << "Source node removed :(" << endl;
			return false; //????
		}
		u_int32_t source_abundance = graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance;
		cout << graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance << endl;
		//cout << "PAS BON CA: utiliser successeur puis predeesseur" << endl;
			cout << endl << endl << endl << endl << endl << "----- Forward -------------------------------------------------------------------------------------------------------------------------------------" << endl;
		_iter = 0;
		//u_int32_t source_nodeIndex = graph->_graphPredecessors->nodeName_to_nodeIndex(utg, false);
		PathData pathData = {pathIndex, {}, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		bool pathSolved = solveBin_path(pathData, graph, true);
		if(pathSolved){
			cout << "Path is solve forward" << endl;
		}
		
		vector<u_int64_t> supportingReads_forward;
		vector<u_int32_t> nodePath_forward = pathData.nodePath;
		getSupportingReads(nodePath_forward, supportingReads_forward);


		vector<u_int32_t> nodePath_backward;
		vector<u_int64_t> supportingReads_backward;

		if(!pathSolved){
			cout << endl << endl << endl << endl << endl << "----- Backward -------------------------------------------------------------------------------------------------------------------------------------" << endl;
			_iter = 0;
			pathData = {pathIndex, {}, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
			pathSolved = solveBin_path(pathData, graph, false);
			if(pathSolved){ //Assez bizarre, l'algo peut echouer en forward mais résoudre en backward et inversement

				supportingReads_forward.clear();
				nodePath_forward.clear();
				cout << "Path is solve backward" << endl;
			}

			nodePath_backward = pathData.nodePath;
			getSupportingReads(nodePath_backward, supportingReads_backward);
		}
		
        
		vector<u_int32_t> nodePath;
		vector<u_int64_t> nodePath_supportingReads;

		if(nodePath_backward.size() > 1){
			std::reverse(nodePath_backward.begin(), nodePath_backward.end());
			std::reverse(supportingReads_backward.begin(), supportingReads_backward.end());
			nodePath_backward.pop_back(); //Remove source node
			supportingReads_backward.pop_back(); //Remove source node
			nodePath = nodePath_backward;
			nodePath_supportingReads = supportingReads_backward;
		}

		nodePath.insert(nodePath.end(), nodePath_forward.begin(), nodePath_forward.end());
		nodePath_supportingReads.insert(nodePath_supportingReads.end(), supportingReads_forward.begin(), supportingReads_forward.end());

		//if(!pathSolved && nodePath.size() < 3000) return false;
		if(!pathSolved) return false;
		
		if(nodePath.size() > 0){
			u_int64_t size = nodePath.size();
			gzwrite(_outputContigFile, (const char*)&size, sizeof(size));
			gzwrite(_outputContigFile, (const char*)&nodePath[0], size * sizeof(u_int32_t));
			gzwrite(_outputContigFile, (const char*)&nodePath_supportingReads[0], size * sizeof(u_int64_t));
		}

		if(_truthInputFilename != ""){
			for(u_int32_t nodeIndex : nodePath){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
					for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
						file_groundTruth_hifiasmContigs << unitigName << "," << pathIndex << endl;
					}
					//cout << _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName] << " " << pathIndex << endl;
				}
			}
		}




		//for(u_int64_t lala : nodePath_supportingReads){
		//	cout << lala << endl;
		//}
		cout << nodePath.size() << endl;
		cout << nodePath_supportingReads.size() << endl;

		/*
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[7836]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[1519]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[1520]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[1521]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[6641]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[8327]) << endl;

		
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[11316]) << endl;
		cout << PathExplorer::computeSharedReads(_unitigDatas[8326], _unitigDatas[1523]) << endl;
		*/
		//exit(1);
		/*
		//cout << nodePath_forward.size() << endl;
		//cout << nodePath_backward.size() << endl;

		if(nodePath_backward.size() > 1){
			
			std::reverse(nodePath_backward.begin(), nodePath_backward.end());
			nodePath_backward.pop_back(); //Remove source node

			u_int64_t size = nodePath_backward.size();
			gzwrite(_outputContigFile, (const char*)&size, sizeof(size));
			gzwrite(_outputContigFile, (const char*)&nodePath_backward[0], size * sizeof(u_int32_t));
		}
		else{
			u_int64_t size = 0;
			gzwrite(_outputContigFile, (const char*)&size, sizeof(size));
		}
		
		

		//u_int64_t pathLength_backward = nodePath_backward.size();
		//u_int64_t pathLength_forward = nodePath_forward.size();
		
		u_int64_t size = nodePath_forward.size();
		gzwrite(_outputContigFile, (const char*)&size, sizeof(size));
		gzwrite(_outputContigFile, (const char*)&nodePath_forward[0], size * sizeof(u_int32_t));
		*/



		/*
		//Skip starting node (nodePath_backward.rend()-1)), added from forward prevNodes
		for (vector<u_int32_t>::reverse_iterator i = nodePath_backward.rbegin(); i != (nodePath_backward.rend()-1); ++i ) { 
			u_int32_t nodeIndex = *i;
			cout << graph->_graphSuccessors->nodeToString(nodeIndex) << endl;
			
			gzwrite(_outputContigFile, (const char*)&nodeIndex, sizeof(nodeIndex));
		} 

		cout << "---------------------------------------------------------------------------" << endl;
		//for(size_t i = nodePath_backward.size() - 1; i >= 0; i--){
		//	cout << i << endl;
			//cout << nodePath_backward[i] << endl;
		//}

		for(size_t i=0; i<nodePath_forward.size(); i++){
			u_int32_t nodeIndex = nodePath_forward[i];
			cout << graph->_graphSuccessors->nodeToString(nodeIndex) << endl;
			
			gzwrite(_outputContigFile, (const char*)&nodeIndex, sizeof(nodeIndex));
			//cout << i << endl;
			//cout << nodePath_forward[i] << endl;
		}
		
		string nextLine = "\n";
		gzwrite(_outputContigFile, (const char*)&nodeIndex, sizeof(nodeIndex));
		*/
		//*/

		/*
		cout << "----- Backward ------" << endl;
		_iter = 0;
		//u_int32_t source_nodeIndex = graph->_graphPredecessors->nodeName_to_nodeIndex(utg, false);
		PathData pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		solveBin_path(pathData, graph);
		*/

		/*
		cout << "----- Start extending left ------" << endl;
		_iter = 0;
		source_nodeIndex = graph->nodeName_to_nodeIndex(utg, true);
		pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		solveBin_path(pathData, graph);
		//solveBin_path(pathData, graph, false);
		*/
		//cout << _pathDatas.size() << endl;

		return true;
	}

	void getSupportingReads(const vector<u_int32_t>& pathNodes, vector<u_int64_t>& supportingReads){

		supportingReads.clear();
		vector<u_int32_t> prevNodes;

		for(u_int32_t nodeIndex : pathNodes){

			//cout << nodeIndex << " " << prevNodes.size() <<  endl;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			UnitigData& unitigData = _unitigDatas[nodeName];

			if(prevNodes.size() == 0){
				prevNodes.push_back(nodeIndex);
				supportingReads.push_back(unitigData._readIndexes[0]);
				continue;
			}

			u_int32_t prevRank = 0;
			u_int64_t supportingRead = -1;

			while(true){
				
				int prevIndex = prevNodes.size() - prevRank - 1;
				//cout << "lalalala     " << prevIndex << endl;
				if(prevIndex < 0 ) break;
				
				u_int32_t prev_nodeIndex = prevNodes[prevIndex];
				u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);
				UnitigData& prev_unitigData = _unitigDatas[prev_nodeName];
				
				vector<u_int64_t> sharedReads;
				PathExplorer::collectSharedReads(unitigData, prev_unitigData, sharedReads);
				//cout << "sdfsdfsd     " << sharedReads.size() << endl;
				if(sharedReads.size() > 0){
					prevRank += 1;
					supportingRead = sharedReads[0];

					//cout << "lala: " << supportingRead << endl;
				}
				else{
					break;
				}
			}


			prevNodes.push_back(nodeIndex);
			supportingReads.push_back(supportingRead);

		}
	}

	bool solveBin_path(PathData& pathData, GraphSimplify* graph, bool forward){

		unordered_set<u_int32_t> visitedNodes;

		u_int32_t current_nodeIndex = pathData.source_nodeIndex;
		binNode(current_nodeIndex, pathData.prevNodes, pathData.nodePath, graph, pathData._index);

		while(true){
			PathExplorer pathExplorer(pathData.prevNodes, pathData.source_abundance, pathData.source_nodeIndex, current_nodeIndex, pathData._abundanceCutoff_min, visitedNodes, pathData._isNodeImproved);
			//u_int32_t nextNodeIndex = pathExplorer.getNextNode( graph, _unitigDatas, 100000, forward, true);
			
			u_int32_t resultType;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, 0, resultType);
			
			if(resultType == 2){ //Banch solved
				cout << "Path size: " << pathData.nodePath.size() << endl;
			}
			//cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

			//cout << pathData.nodePath.size() << endl;
			//cout << current_nodeIndex << endl;
			//if(_binnedNodes.size() > 100) return false; //DEBUG assemble small fragment
			if(current_nodeIndex == -1) return false; //No more successors, or no branching solution
			
			if(current_nodeIndex == pathData.source_nodeIndex){ //Path complete
				pathData.nodePath.pop_back(); //if the path is solved, the source node exist as first and last element,thus we remove the last one

				cout << "Path complete!" << endl;
				return true; 
			}
			//if(current_nodeIndex == -2) return true; //Path complete

			binNode(current_nodeIndex, pathData.prevNodes, pathData.nodePath, graph, pathData._index);
			visitedNodes.insert(current_nodeIndex);

			//anti-bug, infinite loop
			if(_nbVisitedTimes.find(current_nodeIndex) != _nbVisitedTimes.end()){
				if(_nbVisitedTimes[current_nodeIndex] > 1000){
					return false;
				}
			}
			/*
			for(size_t i=1; i<pathExplorer._exploredNodes.size(); i++){
				u_int32_t nodeIndex = pathExplorer._exploredNodes[i];
				current_nodeIndex = nodeIndex;
			}

			if(!foundPath) break;
			*/
		}
	}



	void extract_truth_kminmers(MDBG* mdbg){

		
		u_int64_t maxHashValue = -1;

		u_int64_t _hash_otpt[2];
		int _seed = 42;
		setDispatcher (new SerialDispatcher());


		IBank* inbank = Bank::open(_truthInputFilename);

		
		Iterator<Sequence>* itSeq = createIterator<Sequence> (
															inbank->iterator(),
															inbank->estimateNbItems(),
															"Parsing reads"
															);

		LOCAL (itSeq);
			
		std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();
		u_int32_t readIndex = 0;
		u_int32_t datasetID = 0;

		

		ModelCanonical model (_minimizerSize);
		ModelCanonical::Iterator itKmer (model);

		hash<KmerVec> h;

		for (size_t i=0; i<itBanks.size(); i++)
		{
			itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

			u_int32_t position = 0;

			for (itSeq->first(); !itSeq->isDone(); itSeq->next()){


				Sequence& sequence = itSeq->item();



				string rleSequence;
				vector<u_int64_t> rlePositions;
				Encoder::encode_rle(sequence.getDataBuffer(), sequence.getDataSize(), rleSequence, rlePositions);


				Data buf((char*)rleSequence.c_str());
				itKmer.setData (buf);

				/*
				string sequence_str;

				char lastChar = '0';
				for(size_t i=0; i<sequence.getDataSize(); i++){
					if(readseq[i] == lastChar) continue;
					sequence_str += readseq[i];
					lastChar = readseq[i];
				}


				size_t nbMinimizersPerRead = 0;

				Data buf((char*)sequence_str.c_str());



				itKmer.setData (buf);
				*/

				//u_int64_t lastMinimizer = -1;
				vector<u_int64_t> minimizers;
				//vector<u_int64_t> minimizers_pos;
				//u_int64_t nbMinizersRead = 0;

				//vector<MinimizerPair> minimizerPairs;
				

				//u_int64_t pos = 0;
				//u_int32_t lastMinimizerPos = -1;
				for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

					kmer_type kmerMin = itKmer->value();
					u_int64_t kmerValue = kmerMin.getVal();
					u_int64_t minimizer;
					MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
					minimizer = _hash_otpt[0];



					//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
					double kmerHashed_norm = ((double) minimizer) / maxHashValue;
					if(kmerHashed_norm < _minimizerDensity){


						minimizers.push_back(minimizer);
						//minimizers_pos.push_back(pos);

						//cout << pos << endl;

						//minimizerCounts[minimizer] += 1;
						

					}

					//cout << kmerHashed << endl;
					//pos += 1;
				}

				
				int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
				for(int i=0; i<i_max; i++){

					//cout << minimizers[i] << endl;
					KmerVec vec;

					for(int j=i; j<i+_kminmerSize; j++){
						u_int64_t minimizer = minimizers[j];

						vec._kmers.push_back(minimizer);
					}

					vec = vec.normalize();

					if(mdbg->_dbg_nodes.find(vec) != mdbg->_dbg_nodes.end()){
						_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[mdbg->_dbg_nodes[vec]._index].push_back(sequence.getComment());
						//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
					}


					if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
						//cout << position << endl;
						continue;
					}
					_evaluation_hifiasmGroundTruth[vec] = datasetID;
					_evaluation_hifiasmGroundTruth_position[vec] = position;
					position += 1;



					//cout << position << endl;
					//if(position == 2930) break;
					//cout << position << endl;
					//gzwrite(file, (const char*)&vec._kmers[0], _kminmerSize * sizeof(u_int64_t));
					//gzwrite(file, (const char*)&datasetID, sizeof(datasetID));
					
				}


				//for(size_t i=0; i<minimizers.size(); i++){
				//	cout << minimizers[i] << endl;
				//}

				/*
				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				vector<u_int64_t> minimizersPos; 
				//vector<u_int64_t> rlePositions;
				MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0);
				

				for(size_t i=0; i<kminmers.size(); i++){
					KmerVec& vec = kminmers[i];

					if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
						//cout << endl << "HI: " << position << endl;
						continue;
					}
					_evaluation_hifiasmGroundTruth[vec] = datasetID;
					_evaluation_hifiasmGroundTruth_position[vec] = position;

				
					
					position += 1;
				}
				*/

				readIndex += 1;
			}



			datasetID += 1;
		}
		
		cout << "Nb minimizers groundtruth: " << _evaluation_hifiasmGroundTruth.size() << endl;
		//gzclose(file);
		
		//exit(1);
	}

	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
};


#endif
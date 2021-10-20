
//	"reprise: retester IsImproved pendant pathExploration, si on revient au noeud source sans amelioration on empeche ce path
//	hifiasm ground truth: utiliser contig plutot que unitig pour avoir un meilleur ordering
//	"

//PB: quand zone complexe tres fragmenté en unitig (succession d'unitig tres court), le bestPrevUnitigRank n'est pas très faible si ça se joue à un unitig prêt, ajouter une notion 
//collectPossibleSuccessors: ici on ne veut prune que les chemins dont on est sûr à 100% qu'il ne peuvent pas etre successeur, donc plutot une methode permissive, les loop font qu'on peut rater le bon successeur si on est trop stringent
//Idée: Utiliser strongly connected component -> si solved, on parcours le chemin résolu et on vérifie qu'il n'y a pas d'alternative dû a un génome incomplet (grande tips probablement), si chemin alternatif, on cancel le chemin et on l'assemble avec un connected component classique

/*
		CollectPossibleSuccessor: Mieux gérer la récursivité et le currentDepth
		"reprise: isReachable2: ajouter le maxVisitable?
		Opti: CollectPossibleSUccessors: reset du nbVisitedTimes quand node non visité: on peut faire ça uniquement si on a changé d'unitigs par rapport au node précédent 
		Detection de superbubble:
			- Ajouter un quick patch :( boucle infini)
			- ajouter le isReachable pour prendre en compte les cycle interne, peut etre tricky
		Detecter la source et sink des patterns complexes en les delimitant par les long unitigs, puis developper une methode pour trouer un chemin dedans
		Low abundant genome / gap: parvenir à lier les gaps, il y a deux type de gap:
			- gap qui n'existe pas dans les données (low abundant genome avec une partie non séquencé), Il faut bridge deux partie du graphe avec le nbSharedReads (potentiellement des tips)
			- gap dû a un cutoff min trop élevé: ici le successeur est juste masqué par le isNodeValid2
		"
		Graph::getOverlap(): pas d'overlap parfois ou dans le mauvais sens, a check
		- ToBasespace: mdbg_nodes: on en a besoin pour obtenir le nodeName des kminmer, mais dans ce MDBG, on peut enlever tous les nodes qui ne font pas parti d'un contig pour sauver de la mémoire
*/

#define PRINT_DEBUT_ASM

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
	int _prevRank;
	int _prevUnitigRank;
	bool _prevRankFinished;
	u_int32_t _nodeIndexSuccessor;
	u_int32_t _distanceFromSource;
	u_int32_t _prevUnitigIndex;
	int _nbSharedReads;
	vector<u_int32_t> _path;
	u_int32_t _backtrackPathLength;
	//vector<u_int32_t> _processedNodeIndex;
};

struct SuccessorDistance{
	u_int32_t _nodeIndex;
	u_int32_t _distance;
};

struct DataSuccessorPath{
	u_int32_t _currentNodeIndex;
	u_int32_t _prevNodeIndex;
	vector<u_int32_t> _path;
	u_int32_t _pathLength;
	u_int32_t _nbJokers;
	vector<u_int32_t> _prevNodes;
	unordered_map<u_int32_t, u_int16_t> _nbVisitedTimes;
};

struct DataSuccessorPath_Comparator {
    bool operator()(DataSuccessorPath const& p1, DataSuccessorPath const& p2)
    {
        return p1._nbJokers > p2._nbJokers;
    }
};

/*
class SourceSinkSolver{

	GraphSimplify* _graph;

	SourceSinkSolver(GraphSimplify* graph){
		_graph = graph;
	}

	void solve_unitig(u_int32_t source_unitigIndex, u_int32_t sink_unitigIndex, vector<u_int32_t>& path){
		shortestPath_unitig(source_unitigIndex, sink_unitigIndex, path);
	}


    u_int32_t shortestPath_unitig(u_int32_t source_unitigIndex, u_int32_t sink_unitigIndex, vector<u_int32_t>& path){

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
			_graph->getSuccessors_unitig(unitigIndex_current, successors);


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

            path.clear();
            u_int32_t n = sink_unitigIndex;
            while(n != source_unitigIndex){
                path.push_back(n);
                n = prev[n];
            }
            path.push_back(source_unitigIndex);

            return distance[sink_unitigIndex];
        }

        return -1;
    }

};
*/

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
	vector<UnitigData>& _unitigDatas;

	vector<u_int32_t> _exploredNodes;
	unordered_set<u_int32_t> _isPathAlreadyVisitedSourceNodes;
	u_int64_t _currentPathLength;
	unordered_set<u_int32_t>& _solvedUnitigs;

	PathExplorer(const vector<u_int32_t>& prevNodes, u_int32_t source_abundance, u_int32_t source_nodeIndex, u_int32_t start_nodeIndex, float abundanceCutoff_min, unordered_set<u_int32_t>& visitedNodes, unordered_map<u_int32_t, bool>& isNodeImproved, unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes, vector<UnitigData>& unitigDatas, u_int64_t currentPathLength, unordered_set<u_int32_t>& solvedUnitigs) : _visitedNodes(visitedNodes), _isNodeImproved(isNodeImproved), _unitigDatas(unitigDatas), _solvedUnitigs(solvedUnitigs){
		_prevNodes = prevNodes;
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_start_nodeIndex = start_nodeIndex;
		_abundanceCutoff_min = abundanceCutoff_min;
		_isPathAlreadyVisitedSourceNodes = isPathAlreadyVisitedSourceNodes;
		_currentPathLength = currentPathLength;
	}

	u_int32_t getNextNode(u_int32_t current_nodeIndex, GraphSimplify* graph, bool forward, u_int32_t currentDepth, u_int32_t& resultType, vector<SuccessorData>& nextNodes, bool usePathSuccessors){

		nextNodes.clear();
		resultType = 0;
		//if(currentDepth > 20) return -1;

		//u_int64_t iter = 0;
		bool orient_dummy = false;

		//u_int32_t current_nodeIndex = _start_nodeIndex;

		u_int32_t current_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
		u_int32_t current_abundance = graph->_nodeAbundances[current_nodeName]; //_unitigDatas[current_unitigIndex]._meanAbundance;

		cout << "Current node: " << current_nodeName << endl;
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

		/*
			vector<u_int32_t> successors;
			if(forward){
				graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, successors);
			}
			else{
				graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, successors);
			}

			


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
				//cout << utg_n << " " << _currentPathLength <<  endl;

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				
				if(successors.size() > 1){

					//u_int64_t depth = 0;
					//u_int64_t pathLength = 0;
					//if(currentDepth > 0){
					//	depth = currentDepth + 1;
					//	pathLength = currentPathLength + graph->nodeLengths[successor_nodeName];
					//}

					if(currentDepth == 0){
						if(isPathAlreadyExplored2(utg_n, current_nodeIndex, graph, _unitigDatas, forward, 0, 0)){
							//cout << "Already explored: " << current_nodeName << " " << successor_nodeName << endl;
							continue;
						}
					}

					//}
				}


				SuccessorData successor = {utg_n, successor_abundance, 0, 0, false};
				successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				data_successors.push_back(successor);

			}
			*/
		

		vector<u_int32_t> successors;
		if(forward){
			graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}
		else{
			graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}

		if(successors.size() <= 1 || !usePathSuccessors){
			
			for(u_int32_t utg_n : successors){
				//if(!includeVisitedSuccessors && (_visitedNodes.find(utg_n) != _visitedNodes.end())) continue;

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {utg_n, successor_abundance, 0, 0, 0, false, utg_n, 0, 0, -1, {}, 0};
				successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				data_successors.push_back(successor);

			}

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
			collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, true);

			for(auto& it : successorPaths){

				u_int32_t nodeIndex = it.first;
				DataSuccessorPath& pathData = it.second;
				/*
				nodeIndex = path[1]; //Next node after source_nodeIndex = node

				bool exist = false;
				for(SuccessorData& s : data_successors){
					if(s._nodeIndex == nodeIndex){
						exist= true;
						break;
					}
				}
				if(exist) continue;
				*/
			
				//cout << "Successor: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				//for(u_int32_t node : path){
				//	cout << "\t" << BiGraph::nodeIndex_to_nodeName(node) << endl;
				//}
				int distance = 0;
				/*
				bool exist = false;
				for(SuccessorData& suc: data_successors){
					if(suc._nodeIndex == nodeIndex){
						exist = true;
						break;
					}
				}
				if(exist) continue;
				*/
				//u_int64_t utg_n = node->val;
				//cout << utg_n << " " << _currentPathLength <<  endl;

				u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				SuccessorData successor = {nodeIndex, successor_abundance, 0, 0, 0, false, nodeIndex, distance, 0, -1, pathData._path, 0};
				successor._sourceAbundance = abs((long)successor._abundance - (long)_source_abundance);
				data_successors.push_back(successor);

			}
			//exit(1);
			//}

			if(currentDepth == 0){
				if(successors.size() > 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "\t----------- " << endl;
					cout << "\tNb prev nodes: " << _prevNodes.size() << endl;
					cout << "\tNb possible successors: " << successorPaths.size() << endl;
 					//for(u_int32_t utg_n : successors){
					for(SuccessorData& successor : data_successors){
						u_int32_t utg_n = successor._nodeIndex;
						//for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "\t" << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
					
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

			if(usePathSuccessors){
				if(isInInfiniteCycle(current_nodeIndex, graph, _unitigDatas, forward)) return -1; //Infinite simple path without exit
			}
			
			//cout << "ONE SUCC: " << data_successors[0]._path.size() << endl;
			//for(u_int32_t nodeIndex: data_successors[0]._path){
			//	cout << "\t" << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
			//}
			if(data_successors[0]._path.size() == 0){ //Si useSuccessorPath, le path est deja complet
				data_successors[0]._path.push_back(current_nodeIndex);
				current_nodeIndex = data_successors[0]._nodeIndexSuccessor;
				data_successors[0]._path.push_back(current_nodeIndex);
			}
			//_solvedUnitigs.insert(graph->nodeIndex_to_unitigIndex(current_nodeIndex));
			/*
			//Anti-bug temp
			if(_isNodeImproved.find(current_nodeIndex) != _isNodeImproved.end()){
				if(!_isNodeImproved[current_nodeIndex]){
					return -1;
				}
			}

			*/
			updateNodeChosen(current_nodeIndex, currentDepth);
			nextNodes.push_back(data_successors[0]);
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

			/*
			if(usePathSuccessors){
			
				//cout << "---------------" << endl;
				//cout << successors_bestPrevUnitigRank_noCutoff.size() << endl;
				cout << "\tPossible succesors path:" << endl;
				for(SuccessorData& s : data_successors){
					cout << "\t\t" << BiGraph::nodeToString(s._nodeIndex) << " " << (_visitedNodes.find(s._nodeIndex) != _visitedNodes.end()) << endl;
				}
			}*/
			
			
			if(usePathSuccessors){
				//cout << "allo " << data_successors.size() << endl;
				vector<SuccessorData> dummy;
				computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, 0, dummy, 0, true, true, forward);
			}

			vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
			if(usePathSuccessors){
				computeBestSuccessors_byUnitigRank_all2(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
			}
			else{
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, true, usePathSuccessors);
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, successors_bestPrevUnitigRank_noCutoff, false, usePathSuccessors);
			}

			if(usePathSuccessors){
			
				//cout << "---------------" << endl;
				//cout << successors_bestPrevUnitigRank_noCutoff.size() << endl;
				//cout << "\tPossible succesors path:" << endl;
				//for(SuccessorData& s : successors_bestPrevUnitigRank_noCutoff){
				//	cout << "\t\t" << BiGraph::nodeToString(s._nodeIndex) << " " << (_visitedNodes.find(s._nodeIndex) != _visitedNodes.end()) << endl;
				//}


				/*
				vector<SuccessorData> dummy;
				computeBestSuccessors_byUnitigRank(graph, _unitigDatas, successors_bestPrevUnitigRank_noCutoff, 0, dummy, 0, true, true, forward);

				vector<SuccessorData> pathSuccessors;
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, successors_bestPrevUnitigRank_noCutoff, currentDepth, forward, pathSuccessors, true);
			
				successors_bestPrevUnitigRank_noCutoff = pathSuccessors;
				*/

				//if(current_nodeIndex == BiGraph::nodeName_to_nodeIndex(539588, false)){
				//	exit(1);
				//}
			}

			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

				if(currentDepth == 0){
					cout << "\tNode chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				//updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
				return current_nodeIndex;
			}

			/*
			//Print debug prev ranks
			vector<SuccessorData> dummy;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, dummy, 0, true, true, forward);



			vector<SuccessorData> successors_bestPrevUnitigRank_cutoff;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_cutoff, _abundanceCutoff_min/2, true, false, forward);
			if(successors_bestPrevUnitigRank_cutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_cutoff[0]._nodeIndexSuccessor;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 1: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(current_nodeIndex);
				return current_nodeIndex;
			}

			vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_noCutoff, 0, true, false, forward);
			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 2: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(current_nodeIndex);
				return current_nodeIndex;
			}

			vector<SuccessorData> successors_bestPrevRank;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevRank, _abundanceCutoff_min/2, false, false, forward);
			if(successors_bestPrevRank.size() == 1){
				current_nodeIndex = successors_bestPrevRank[0]._nodeIndexSuccessor;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen 3: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				updateNodeChosen(current_nodeIndex, currentDepth);
				nextNodes.push_back(current_nodeIndex);
				return current_nodeIndex;
			}
			*/

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




			if(usePathSuccessors){
				//cout << "lalala " << successors_bestPrevUnitigRank_noCutoff.size() << endl;
				if(successors_bestPrevUnitigRank_noCutoff.size() >= 2){
					bool isBubble = true;
					for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){
						if(!graph->_isBubble[successor._nodeIndexSuccessor]){
							isBubble = false;
						}
					}
					//cout << isBubble << endl;
					if(isBubble){
						
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);

						//if(currentDepth == 0){
							cout << "\tTake bubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
						//}
						//exit(1);
						updateNodeChosen(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
						return successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;
					}

					//cout << currentDepth << " " << "detect" << endl;
					u_int32_t superbubbleSink = detectSuperbubble(current_nodeIndex, maxLength, graph, _unitigDatas, forward);
					if(superbubbleSink != 1){
						sort(successors_bestPrevUnitigRank_noCutoff.begin(), successors_bestPrevUnitigRank_noCutoff.end(), SuccessorComparator_byAbundanceToSource);
						//if(currentDepth == 0){
							cout << "\tTake superbubble: " << graph->_graphSuccessors->nodeToString(successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor) << endl;
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


			
			//if(currentDepth == 0){

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "\tCheck simple cycle" << endl;
				}

				for(SuccessorData& successor : successors_bestPrevUnitigRank_noCutoff){

					if(isSmallCycle(successor._nodeIndexSuccessor, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
						updateNodeChosen(successor._nodeIndexSuccessor, currentDepth);
						nextNodes.push_back(successor);
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
							return from_nodeIndex;
						}

					}	
				}
			}
			
			
			for(size_t i=0; i<successors_bestPrevUnitigRank_noCutoff.size(); i++){
				nextNodes.push_back(successors_bestPrevUnitigRank_noCutoff[i]);
			}


			//if(!usePathSuccessors && currentDepth == 0){
			//	return getNextNode(current_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, true);
			//}

			return -1;
				
			//}

		}


	}

	
	bool isInInfiniteCycle(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){
		Unitig& unitig = graph->nodeIndex_to_unitig(current_nodeIndex);
		if(unitig._startNode != current_nodeIndex) return false;

		unordered_map<u_int32_t, DataSuccessorPath> successorPaths;
		collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, successorPaths, false);

		//cout << "Check infinite cycle: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " <<  successorPaths.size() << endl;
		if(successorPaths.size() == 0){
			cout << "Inifinite cycle" << endl;
			//getchar();
		} 

		return successorPaths.size() == 0;
	}

    u_int32_t detectSuperbubble(u_int32_t nodeIndex_source, u_int64_t maxLength, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){

		cout << "Start detection: " << BiGraph::nodeIndex_to_nodeName(nodeIndex_source) << endl;
		/*
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
		nbVisitedTimes[nodeName] += 1;

		u_int32_t maxVisitable = (graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
	
		//cout << nodeName << " " << maxVisitable << " " << graph->getNodeUnitigAbundance(current_nodeIndex) << " " << _source_abundance << endl;

		if(nbVisitedTimes[nodeName] > maxVisitable){
			return;
		}*/


        unordered_set<u_int32_t> isVisited;
        unordered_set<u_int32_t> seen;
        unordered_map<u_int32_t, u_int64_t> pathLength;
        vector<u_int32_t> queue;

        queue.push_back(nodeIndex_source);
        pathLength[nodeIndex_source] = 0;

        while(queue.size() > 0){
            u_int32_t v = queue[queue.size()-1];

			cout << "Visit: " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;
            //cout << "\tVisited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[v]._endNode) << endl;
            queue.pop_back();

			cout << "Visit 1 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

            if(pathLength[v] > 10000) continue;

			cout << "Visit 2 : " << BiGraph::nodeIndex_to_nodeName(v) << " " << queue.size() << endl;

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
			getNextNode(v, graph, forward, 0, resultType, nextNodes, false);

			cout << "Nb succ: " << nextNodes.size() << endl;

            //vector<u_int32_t> successors;
			//getSuccessors_unitig(v, successors);
            if(nextNodes.size() == 0) return -1; //abort tip

			for(SuccessorData& successorData : nextNodes){
				u_int32_t u = successorData._nodeIndex;
				
				cout << "\tSucc : " << BiGraph::nodeIndex_to_nodeName(u) << " " << queue.size() << endl;

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
				getNextNode(u, graph, !forward, 0, resultType, prevNodes, false);
				//getPredecessors_unitig(u, predecessors);
                bool allPredecessorsAreVisited = true;
				for(SuccessorData& predecessorData : prevNodes){
					u_int32_t p = predecessorData._nodeIndex;
                	//for(u_int32_t p : predecessors){
                    if(isVisited.find(p) == isVisited.end()){
                        allPredecessorsAreVisited = false;
                        break;
                    }
                }

                if(allPredecessorsAreVisited){
					cout << "Push: " << BiGraph::nodeIndex_to_nodeName(u) << endl;
                    //cout << "\t\tAll predecessors visited: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << endl;
                    queue.push_back(u);
                }

                //cout << "\t\t\tQueue size: " << queue.size() << " " << seen.size() << endl;
                if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                    u_int32_t t = queue[0];
					
					vector<SuccessorData> successors_t;
					getNextNode(t, graph, forward, 0, resultType, successors_t, false);

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
	
		__beatenNodeIndex.clear();

		for(size_t i=0; i<data_successors2.size(); i++){
			for(size_t j=i+1; j<data_successors2.size(); j++){

				if(__beatenNodeIndex.find(data_successors2[i]._nodeIndex) != __beatenNodeIndex.end()) continue;
				if(__beatenNodeIndex.find(data_successors2[j]._nodeIndex) != __beatenNodeIndex.end()) continue;

				//cout << i << " vs " << j << endl;
				vector<SuccessorData> data_successors;
				data_successors.push_back(data_successors2[i]);
				data_successors.push_back(data_successors2[j]);

				size_t rank = 1;
				while(true){
					if(data_successors[0]._path[rank] != data_successors[1]._path[rank]) break;
					rank +=1;
				}

				data_successors[0]._nodeIndexSuccessor = data_successors[0]._path[rank];
				data_successors[1]._nodeIndexSuccessor = data_successors[1]._path[rank];
				vector<SuccessorData> bestSuccessors;
				computeBestSuccessors_byUnitigRank_all(graph, _unitigDatas, data_successors, currentDepth, forward, bestSuccessors, includeVisitedSuccessors, usePathSuccessors);
			
				//cout << "Winner: " << bestSuccessors.size() << " " << BiGraph::nodeIndex_to_nodeName(bestSuccessors[0]._nodeIndex) << endl;
				if(bestSuccessors.size() == 1){
					if(bestSuccessors[0]._nodeIndex == data_successors[0]._nodeIndex){
						__beatenNodeIndex.insert(data_successors[1]._nodeIndex);
					}
					else{
						__beatenNodeIndex.insert(data_successors[0]._nodeIndex);
					}
				}
			}
		}

		validSuccessors.clear();
		for(SuccessorData& s : data_successors2){
			if(__beatenNodeIndex.find(s._nodeIndex) != __beatenNodeIndex.end()) continue;
			validSuccessors.push_back(s);
		}

		if(validSuccessors.size() == 0) return;

		//if(currentDepth == 0 && abundanceCutoff_min != 0){

			unordered_set<u_int32_t> reachableSuccessors;
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
		float cutoff_min = _abundanceCutoff_min/4;
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
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_cutoff, cutoff_min, true, false, forward);
			if(successors_bestPrevUnitigRank_cutoff.size() == 1){
				u_int32_t current_nodeIndex = successors_bestPrevUnitigRank_cutoff[0]._nodeIndexSuccessor;

				if(currentDepth == 1){
					cout << "\tNode chosen 1: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;// << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

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
		//}

		
		vector<SuccessorData> successors_bestPrevUnitigRank_noCutoff;
		if(usePathSuccessors){
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevUnitigRank_noCutoff, 0, true, false, forward);
			if(successors_bestPrevUnitigRank_noCutoff.size() == 1){
				u_int32_t current_nodeIndex = successors_bestPrevUnitigRank_noCutoff[0]._nodeIndexSuccessor;

				if(currentDepth == 1){
					cout << "\tNode chosen 2: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;//<< " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				for(SuccessorData& s : validSuccessors){
					if(s._nodeIndex == current_nodeIndex) return;
				}

				//updateNodeChosen(current_nodeIndex, currentDepth);
				//nextNodes.push_back(current_nodeIndex);
				//return current_nodeIndex;
				//cout << "Add valid successor: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << currentDepth << (_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != _visitedNodes.end()) << endl; 
				validSuccessors.push_back(successors_bestPrevUnitigRank_noCutoff[0]);
				return;
			}
		}

		/*
		if(usePathSuccessors){
			vector<SuccessorData> successors_bestPrevRank;
			computeBestSuccessors_byUnitigRank(graph, _unitigDatas, data_successors, currentDepth, successors_bestPrevRank, cutoff_min, false, false, forward);
			if(successors_bestPrevRank.size() == 1){
				u_int32_t current_nodeIndex = successors_bestPrevRank[0]._nodeIndexSuccessor;

				if(currentDepth == 1){
					cout << "\tNode chosen 3: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;//<< " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				//updateNodeChosen(current_nodeIndex, currentDepth);
				//nextNodes.push_back(current_nodeIndex);
				//return current_nodeIndex;
				validSuccessors.push_back(successors_bestPrevRank[0]);
				return;
			}
		}
		*/

		//cout << "ALLO: " << usePathSuccessors << " " << includeVisitedSuccessors << endl;

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

	void computeBestSuccessors_byUnitigRank(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, vector<SuccessorData>& data_successors, int currentDepth, vector<SuccessorData>& successors_bestPrevUnitigRank, float abundanceCutoff_min, bool useUnitigRank, bool print_debug, bool forward){

		//if(print_debug && currentDepth == 0){
			//cout << endl << "Cutoff: " << abundanceCutoff_min << endl;
		//}

		for(SuccessorData& successor : data_successors){
			successor._prevRankFinished = false;
			successor._prevUnitigRank = -1;
			successor._prevRank = -1;
			successor._backtrackPathLength = 0;
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

				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(successor._nodeIndexSuccessor);
				
				u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
				//if(nbSharedReads > _abundanceCutoff_min/2){
				//if(nbSharedReads > successor._abundance/5){
				//if(nbSharedReads > 0){
				if(nbSharedReads > 0 && nbSharedReads > abundanceCutoff_min){

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

				lastNode = prev_nodeIndex;
				i += 1;
			}

			str_debug += "    " + to_string(pathLength);
			/*
			if(!print_debug && currentDepth == 0){
				if(nbUnfinishedSuccessors == 1){ //We record solved unitigs
					//if(_solvedUnitigs.find(currentUnitigIndex) == _solvedUnitigs.end()) cout << "Unitig solved: " << currentUnitigIndex << endl;
					_solvedUnitigs.insert(currentUnitigIndex);
				}
				
				bool onlySolvedNodeRemains = true;
				for(SuccessorData& successor : data_successors){
					if(std::find(chosenSuccessors.begin(), chosenSuccessors.end(), successor._nodeIndex) == chosenSuccessors.end()){
						onlySolvedNodeRemains = false;
						break;
					}
				}

				if(onlySolvedNodeRemains){
					_solvedUnitigs.insert(currentUnitigIndex);
				}
			}*/


			//cout << canExplorePath << endl;
			//cout << current_nodeIndex << " "  <<  _source_nodeIndex << endl;

			//if(print_debug && currentDepth == 0){
				//cout << str_debug << endl;
			//}

			prevRank += 1;

			//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


		}

		if(print_debug) return;

		std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPathLength);
		int maxPathLength = data_successors[0]._backtrackPathLength;
		//cout << "MAX 2: " << maxPrevRank << endl;
		//vector<SuccessorData> successors_bestPrevRank;
		for(SuccessorData& successor : data_successors){
			if(successor._backtrackPathLength > maxPathLength-2000){
				successors_bestPrevUnitigRank.push_back(successor);
			}
		}


		/*
		if(useUnitigRank){
			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevUnitigRank);
			int maxPrevUnitigRank = data_successors[0]._prevUnitigRank;
			//cout << "MAX 1: " << maxPrevUnitigRank << endl;
			
			for(SuccessorData& successor : data_successors){
				//cout << BiGraph::nodeIndex_to_nodeName(successor._nodeIndex) << " " <<  successor._prevUnitigRank << " " << maxPrevUnitigRank << endl;
				if(successor._prevUnitigRank == maxPrevUnitigRank){
					successors_bestPrevUnitigRank.push_back(successor);
				}
			}
		}
		else{

			
			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevRank);
			int maxPrevRank = data_successors[0]._prevRank;
			//cout << "MAX 2: " << maxPrevRank << endl;
			//vector<SuccessorData> successors_bestPrevRank;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRank > maxPrevRank-5){
					successors_bestPrevUnitigRank.push_back(successor);
				}
			}

		}
		*/

		

		

		/*
		vector<SuccessorData> chosenSuccessors;
		std::sort(successors_bestPrevUnitigRank.begin(), successors_bestPrevUnitigRank.end(), SuccessorComparator_byNbSharedReads);
		for(SuccessorData& successorData : successors_bestPrevUnitigRank){
			if(successorData._distanceFromSource == successors_bestPrevUnitigRank[0]._distanceFromSource){
				chosenSuccessors.push_back(successorData);
			}
		}


		*/

		/*
		if(!print_debug && useUnitigRank && successors_bestPrevUnitigRank.size() > 1){

			//cout << (_solvedUnitigs.find(79) != _solvedUnitigs.end()) << endl;
			bool containsOnlySolvedUnitigs = true;
			for(SuccessorData& successor : successors_bestPrevUnitigRank){
				//cout << "lala: " << BiGraph::nodeIndex_to_nodeName(successor._nodeIndex) << " " << successor._prevUnitigIndex << endl;
				if(_solvedUnitigs.find(successor._prevUnitigIndex) == _solvedUnitigs.end()){
					containsOnlySolvedUnitigs = false;
				}
			}

			//cout << "loulou: " << containsOnlySolvedUnitigs << endl;
			if(containsOnlySolvedUnitigs){

				vector<SuccessorData> chosenSuccessors;

				std::sort(successors_bestPrevUnitigRank.begin(), successors_bestPrevUnitigRank.end(), SuccessorComparator_byDistance);
				for(SuccessorData& successorData : successors_bestPrevUnitigRank){
					if(successorData._distanceFromSource == successors_bestPrevUnitigRank[0]._distanceFromSource){
						chosenSuccessors.push_back(successorData);
					}
				}

				std::sort(chosenSuccessors.begin(), chosenSuccessors.end(), SuccessorComparator_byNbSharedReads);
				//SuccessorData& bestSuccessorData = successors_bestPrevUnitigRank[0];

				successors_bestPrevUnitigRank.clear();
				successors_bestPrevUnitigRank.push_back(chosenSuccessors[0]);

				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! successors: " << successors_bestPrevUnitigRank.size() << endl;
				//exit(1);
			}



		}
		*/



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
				graph->getSuccessors(v, 0, successors);
			}
			else{
				graph->getPredecessors(v, 0, successors);
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

			/*
            for(u_int32_t u : successors){
                if(u == unitigIndex_source) return -1; //cycle including s

                if(isVisited.find(u) == isVisited.end()){
                    seen.insert(u);
                    pathLength[u] = pathLength[v] + _unitigs[u]._length;
                }

            }

            for(u_int32_t u : successors){

                //cout << "\t\tVisiting: " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._startNode) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[u]._endNode) << endl;

                vector<u_int32_t> predecessors;
                getPredecessors_unitig(u, predecessors);
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
                    getSuccessors_unitig(t, successors_t);
                    if(std::find(successors_t.begin(), successors_t.end(), unitigIndex_source) == successors_t.end()){
                        return t;
                    }
                    else{
                        return -1; // cycle including s
                    }

                }


            }
			*/


	}

	void updateNodeChosen(u_int32_t nodeIndex, u_int64_t currentDepth){
		/*
		if(currentDepth != 0) return;

		if(_visitedNodes.find(nodeIndex) == _visitedNodes.end()){
			for(auto& it : _isNodeImproved){
				_isNodeImproved[it.first] = true;
			}
		}

		_isNodeImproved[nodeIndex] = false;
		*/
	}

	/*
	bool isPathAlreadyExplored(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t maxIter){
		
		//if(_isPathAlreadyVisitedSourceNodes.find(source_nodeIndex) != _isPathAlreadyVisitedSourceNodes.end()){
		//	return true;
		//}
		

		//bool lala = true;

		if(currentDepth == 1){
			//unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes;
			_isPathAlreadyVisitedSourceNodes.clear();
			_isPathAlreadyVisitedSourceNodes.insert(source_nodeIndex);
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		}
		else{
			_isPathAlreadyVisitedSourceNodes.insert(source_nodeIndex);
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
		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas);

		//u_int64_t maxIter = 10;


		u_int64_t iter = 0;
		//u_int32_t current_nodeIndex = _start_nodeIndex;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			vector<u_int32_t> nextNodes;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth, resultType, nextNodes);
			//cout << currentDepth << " " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
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

			//for(u_int32_t nodeIndex : nextNodes){
			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}
			//}


			//for(u_int32_t nodeIndex : nextNodes){
			pathExplorer.nodeExplored(current_nodeIndex, graph);
			//}

			if(iter >= maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " Yes" << endl;
		}

		//if(!lala) return false;
		return true;
	}*/


	u_int64_t maxLength = 100000;
	u_int64_t maxLength_reachable = 20000;

	/*
	bool isPathAlreadyExplored2(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength){
		

		if(currentDepth == 0){	
			cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		}

		if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
			if(currentDepth == 0) cout << " No" << endl;
			return false;
		}


		if(currentLength > maxLength){
			if(currentDepth == 0) cout << " Yes (1)" << endl;
			return true;
		}

		//u_int32_t source_nodeIndex = current_nodeIndex;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, currentLength, _solvedUnitigs);

		//u_int64_t maxIter = 10;


		//u_int64_t length = 0;
		//u_int32_t current_nodeIndex = _start_nodeIndex;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			vector<u_int32_t> nextNodes;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);

			cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << pathExplorer._currentPathLength << " " << nextNodes.size() << endl;
			//cout << currentDepth << " " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			//if(current_nodeIndex == source_nodeIndex) break;

			if(nextNodes.size() == 0){ //dead end or multiple braching path
				if(currentDepth == 0){
					cout << " No" << endl;
				}
				return false;
			}
			else if(nextNodes.size() == 1){ //dead end or multiple braching path
				pathExplorer.nodeExplored(current_nodeIndex, graph);
				current_nodeIndex = nextNodes[0];
			}
			else{
				for(u_int32_t nodeIndex  : nextNodes){
					cout << "HHAAAAA " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
					//pathExplorer.nodeExplored(nodeIndex, graph);
					if(!isPathAlreadyExplored2(nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, pathExplorer._currentPathLength)){
						return false;
					}
				}
			}


			//for(u_int32_t nodeIndex : nextNodes){
			//}

			//if(currentLength >= maxLength) break;

			//currentLength += 1;
			
			if(pathExplorer._currentPathLength > maxLength){
				cout << "max length exceed: " << pathExplorer._currentPathLength << " " << maxLength << endl;
				break;
			}
		}
		
		if(currentDepth == 0){
			cout << " Yes (2)" << endl;
		}

		return true;
	}
	*/

	void collectPossibleSuccessors(u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, bool print_debug){
		
		u_int32_t minNbJokers = -1;
		if(print_debug) cout << "Collect possible successors" << endl;
		u_int32_t currentDepth = 0;

		priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator> queue;
		queue.push({source_nodeIndex, UINT32_MAX, {}, 0, 0, _prevNodes, {}});



		while(queue.size() > 0){

			DataSuccessorPath dataSuccessorPath = queue.top();
        	queue.pop();

			if(dataSuccessorPath._nbJokers > minNbJokers) continue;

			u_int32_t current_nodeIndex = dataSuccessorPath._currentNodeIndex; //dataSuccessorPath._path[dataSuccessorPath._path.size()-1];
			if(print_debug) cout << "\tStart: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << dataSuccessorPath._nbJokers << " " << minNbJokers << " " << currentDepth << endl;
			//getchar();
			//u_int32_t currentLength = dataSuccessorPath._pathLength;

			PathExplorer pathExplorer(dataSuccessorPath._prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, dataSuccessorPath._pathLength, _solvedUnitigs);
			//pathExplorer.nodeExplored(current_nodeIndex, graph);

			bool continueVisiting = visitSuccessor(current_nodeIndex, dataSuccessorPath._prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug);
			if(!continueVisiting) continue;

			while(true){

				vector<u_int32_t> allSuccessors;
				if(forward){
					graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, allSuccessors);
				}
				else{
					graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, allSuccessors);
				}

				//cout << "\t\tTotal successors: " << allSuccessors.size() << endl;
				if(allSuccessors.size() == 0){
					break;
				}
				else if(allSuccessors.size() == 1){
					
					u_int32_t prevNodeIndex = current_nodeIndex;
					current_nodeIndex = allSuccessors[0];
					bool continueVisiting = visitSuccessor(current_nodeIndex, prevNodeIndex, dataSuccessorPath, pathExplorer, graph, successorPaths, minNbJokers, print_debug);
					if(!continueVisiting) break;

				}
				else{

					u_int32_t resultType;
					vector<SuccessorData> nextNodes;
					pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);

					//cout << "\t\tNb valid successors: " << nextNodes.size() << endl;

					for(SuccessorData& successorData : nextNodes){
						//cout << "\t\t\tAdd valid successor: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorData._nodeIndex, current_nodeIndex, queue, dataSuccessorPath, false, graph, pathExplorer, currentDepth, print_debug);
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

						//cout << "\t\t\tAdd invalid successor: " << BiGraph::nodeIndex_to_nodeName(successorNodeIndex) << endl;
						//getchar();
						currentDepth += 1;
						addSuccessorPath(successorNodeIndex, current_nodeIndex, queue, dataSuccessorPath, true, graph, pathExplorer, currentDepth, print_debug);
					}

					break;

				}


			}

		}

		/*
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
		
		if(nbVisitedTimes.find(nodeName) == nbVisitedTimes.end()){
			for(auto it : nbVisitedTimes){
				nbVisitedTimes[it.first] = 0;
			}
		}

		nbVisitedTimes[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
	
		cout << lalalala << "    " << nodeName << " " << maxVisitable << " " << nbVisitedTimes[nodeName]  << " " << currentDepth << endl;//<< "        " << graph->getNodeUnitigAbundance(current_nodeIndex) << " " << _source_abundance << endl;

		lalalala += 1;
		if(nbVisitedTimes[nodeName] > maxVisitable){
			return;
		}

		if(lalalala > 1000) return;

		//cout << currentLength << ":     " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
		path.push_back(current_nodeIndex);
		//if(currentDepth == 0){	
		//	cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		//}

		if(_visitedNodes.find(nodeName) == _visitedNodes.end()){
			cout << "Not visited: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			if(successorPaths.find(current_nodeIndex) == successorPaths.end()){
				successorPaths[current_nodeIndex] = {path, currentLength};
			}
			else{
				if(currentLength < successorPaths[current_nodeIndex]._pathLength){
					successorPaths[current_nodeIndex] = {path, currentLength};
				}
			}
			return;
		}

		//cout << currentDepth << " " << currentLength << endl;
		//cout << nbBranchs << endl;
		if(currentLength > maxLength){//|| currentDepth > 50){
			return;
		}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, currentLength, _solvedUnitigs);
		pathExplorer.nodeExplored(current_nodeIndex, graph);
		
		while(true){
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);

			if(nextNodes.size() == 0){
				return;
			}
			else if(nextNodes.size() == 1){
				u_int32_t prev_nodeIndex = current_nodeIndex;
				current_nodeIndex = nextNodes[0]._nodeIndex;
				cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;


				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
				
				if(nbVisitedTimes.find(nodeName2) == nbVisitedTimes.end()){
					for(auto it : nbVisitedTimes){
						nbVisitedTimes[it.first] = 0;
					}
				}

				nbVisitedTimes[nodeName2] += 1;

				if(nbVisitedTimes[nodeName2] > maxVisitable){
					return;
				}

				path.push_back(current_nodeIndex);

				if(_visitedNodes.find(nodeName2) == _visitedNodes.end()){
					if(successorPaths.find(current_nodeIndex) == successorPaths.end()){
						successorPaths[current_nodeIndex] = {path, currentLength};
					}
					else{
						if(currentLength < successorPaths[current_nodeIndex]._pathLength){
							successorPaths[current_nodeIndex] = {path, currentLength};
						}
					}
					return;
				}

				if(currentLength > maxLength){//|| currentDepth > 50){
					return;
				}

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(prev_nodeIndex, current_nodeIndex);
				currentLength += (graph->_nodeLengths[nodeName2] - overlapLength);
				pathExplorer.nodeExplored(current_nodeIndex, graph);
			}
			else{
				for(SuccessorData& successorData : nextNodes){
					cout << "\t\tStart new iter: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
					u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(current_nodeIndex, successorData._nodeIndex);
					pathExplorer.collectPossibleSuccessors(successorData._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)] - overlapLength), path, successorPaths, nbVisitedTimes, lalalala);
				}
				return;
			}
		}
		*/

		
	}

	void addSuccessorPath(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, priority_queue<DataSuccessorPath, vector<DataSuccessorPath>, DataSuccessorPath_Comparator>& queue, const DataSuccessorPath& dataSuccessorPath, bool addJoker, GraphSimplify* graph, const PathExplorer& pathExplorer, u_int32_t& currentDepth, bool print_debug){

		//currentDepth += 1;
		if(currentDepth > 1000) return;
		//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

		u_int32_t nbJokers = dataSuccessorPath._nbJokers;
		if(addJoker) nbJokers += 1;

		//u_int32_t pathLength = dataSuccessorPath._pathLength;
		//u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
		//pathLength += (graph->_nodeLengths[nodeName] - overlapLength);

		DataSuccessorPath dataSuccessorPath_next = {nodeIndex, nodeIndexPrev, dataSuccessorPath._path, dataSuccessorPath._pathLength, nbJokers, pathExplorer._prevNodes};
		//dataSuccessorPath_next._path.push_back(nodeIndex);

		queue.push(dataSuccessorPath_next);
	}

	
	bool visitSuccessor(u_int32_t nodeIndex, u_int32_t nodeIndexPrev, DataSuccessorPath& dataSuccessorPath, PathExplorer& pathExplorer, GraphSimplify* graph, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, u_int32_t& minNbJokers, bool print_debug){
		
		//cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;

		//cout << "\t\tVisit successor: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		
		
		if(dataSuccessorPath._nbVisitedTimes.find(nodeName) == dataSuccessorPath._nbVisitedTimes.end()){
			for(auto it : dataSuccessorPath._nbVisitedTimes){
				dataSuccessorPath._nbVisitedTimes[it.first] = 0;
			}
		}

		dataSuccessorPath._nbVisitedTimes[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
		if(dataSuccessorPath._nbVisitedTimes[nodeName] > maxVisitable){
			//cout << "\t\t\tStop max visitables" << endl; 
			return false;
		}
		

		dataSuccessorPath._path.push_back(nodeIndex);

		
		if(_visitedNodes.find(nodeName) == _visitedNodes.end()){
			if(successorPaths.find(nodeIndex) == successorPaths.end()){
				successorPaths[nodeIndex] = dataSuccessorPath;
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
			if(print_debug) cout << "\t\t\t" << "Not visited: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << dataSuccessorPath._nbJokers << endl;
			return false;
		}

		if(dataSuccessorPath._pathLength > maxLength){//|| currentDepth > 50){
			//cout << "\t\t\tStop max length" << endl; 
			return false;
		}

		if(nodeIndexPrev != -1){ //First iteration
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(nodeIndexPrev, nodeIndex);
			dataSuccessorPath._pathLength += (graph->_nodeLengths[nodeName] - overlapLength);
		}
		pathExplorer.nodeExplored(nodeIndex, graph);

		return true;
	}

	/*
	void collectPossibleSuccessors(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, vector<u_int32_t> path, unordered_map<u_int32_t, DataSuccessorPath>& successorPaths, unordered_map<u_int32_t, u_int16_t> nbVisitedTimes, u_int64_t& lalalala){
		
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
		
		if(nbVisitedTimes.find(nodeName) == nbVisitedTimes.end()){
			for(auto it : nbVisitedTimes){
				nbVisitedTimes[it.first] = 0;
			}
		}

		nbVisitedTimes[nodeName] += 1;

		u_int32_t maxVisitable = 2; //(graph->getNodeUnitigAbundance(current_nodeIndex) / (float) _source_abundance) * 20;
		
	
		cout << lalalala << "    " << nodeName << " " << maxVisitable << " " << nbVisitedTimes[nodeName]  << " " << currentDepth << endl;//<< "        " << graph->getNodeUnitigAbundance(current_nodeIndex) << " " << _source_abundance << endl;

		lalalala += 1;
		if(nbVisitedTimes[nodeName] > maxVisitable){
			return;
		}

		if(lalalala > 1000) return;

		//cout << currentLength << ":     " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
		path.push_back(current_nodeIndex);
		//if(currentDepth == 0){	
		//	cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		//}

		if(_visitedNodes.find(nodeName) == _visitedNodes.end()){
			cout << "Not visited: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			if(successorPaths.find(current_nodeIndex) == successorPaths.end()){
				successorPaths[current_nodeIndex] = {path, currentLength};
			}
			else{
				if(currentLength < successorPaths[current_nodeIndex]._pathLength){
					successorPaths[current_nodeIndex] = {path, currentLength};
				}
			}
			return;
		}

		//cout << currentDepth << " " << currentLength << endl;
		//cout << nbBranchs << endl;
		if(currentLength > maxLength){//|| currentDepth > 50){
			return;
		}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, currentLength, _solvedUnitigs);
		pathExplorer.nodeExplored(current_nodeIndex, graph);
		
		while(true){
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);

			if(nextNodes.size() == 0){
				return;
			}
			else if(nextNodes.size() == 1){
				u_int32_t prev_nodeIndex = current_nodeIndex;
				current_nodeIndex = nextNodes[0]._nodeIndex;
				cout << "\t\tSimple node: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;


				u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
				
				if(nbVisitedTimes.find(nodeName2) == nbVisitedTimes.end()){
					for(auto it : nbVisitedTimes){
						nbVisitedTimes[it.first] = 0;
					}
				}

				nbVisitedTimes[nodeName2] += 1;

				if(nbVisitedTimes[nodeName2] > maxVisitable){
					return;
				}

				path.push_back(current_nodeIndex);

				if(_visitedNodes.find(nodeName2) == _visitedNodes.end()){
					if(successorPaths.find(current_nodeIndex) == successorPaths.end()){
						successorPaths[current_nodeIndex] = {path, currentLength};
					}
					else{
						if(currentLength < successorPaths[current_nodeIndex]._pathLength){
							successorPaths[current_nodeIndex] = {path, currentLength};
						}
					}
					return;
				}

				if(currentLength > maxLength){//|| currentDepth > 50){
					return;
				}

				u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(prev_nodeIndex, current_nodeIndex);
				currentLength += (graph->_nodeLengths[nodeName2] - overlapLength);
				pathExplorer.nodeExplored(current_nodeIndex, graph);
			}
			else{
				for(SuccessorData& successorData : nextNodes){
					cout << "\t\tStart new iter: " << BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex) << endl;
					u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(current_nodeIndex, successorData._nodeIndex);
					pathExplorer.collectPossibleSuccessors(successorData._nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)] - overlapLength), path, successorPaths, nbVisitedTimes, lalalala);
				}
				return;
			}
		}






	
		
	}
	*/

	bool isReachable2(u_int32_t current_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, bool allowReverseDirection, u_int32_t& totalIter){
		

		if(allowReverseDirection){
			if(BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == BiGraph::nodeIndex_to_nodeName(to_nodeIndex)) return true;
		}
		else{
			if(current_nodeIndex == to_nodeIndex) return true;
		}

		//cout << currentLength << ":     " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
		//path.push_back(current_nodeIndex);
		//if(currentDepth == 0){	
		//	cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		//}

		//if(_visitedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) == _visitedNodes.end()){
		//	successorPaths[current_nodeIndex] = path;
		//	return;
		//}


		totalIter += 1;
		if(totalIter > 1000) return false;

		if(currentLength > maxLength_reachable || currentDepth > 50){
			return false;
		}


		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, currentLength, _solvedUnitigs);
		pathExplorer.nodeExplored(current_nodeIndex, graph);
		
		u_int32_t resultType;
		vector<SuccessorData> nextNodes;
		pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);

		for(SuccessorData& successorData : nextNodes){
			u_int16_t overlapLength = graph->_graphSuccessors->getOverlap(current_nodeIndex, successorData._nodeIndex);
			bool isReachable = pathExplorer.isReachable2(successorData._nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, currentLength+(graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex)]-overlapLength), allowReverseDirection, totalIter);
			if(isReachable) return true;
		}

		return false;
		/*
		while(true){
		

			cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << pathExplorer._currentPathLength << " " << nextNodes.size() << endl;
			//cout << currentDepth << " " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			//if(current_nodeIndex == source_nodeIndex) break;

			if(nextNodes.size() == 0){ //dead end or multiple braching path
				if(currentDepth == 0){
					cout << " No" << endl;
				}
				return false;
			}
			else if(nextNodes.size() == 1){ //dead end or multiple braching path
				pathExplorer.nodeExplored(current_nodeIndex, graph);
				current_nodeIndex = nextNodes[0];
			}
			else{
				for(u_int32_t nodeIndex  : nextNodes){
					cout << "HHAAAAA " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
					//pathExplorer.nodeExplored(nodeIndex, graph);
					if(!isPathAlreadyExplored2(nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, pathExplorer._currentPathLength)){
						return false;
					}
				}
			}
			
		}*/
		
	}

	/*
	bool isChimericPath(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, vector<u_int32_t> path, unordered_map<u_int32_t, vector<u_int32_t>>& successorPaths){
		
		u_int32_t unitigIndex_source = graph->nodeIndex_to_unitigIndex(current_nodeIndex);
		u_int32_t nodeIndex_start = current_nodeIndex;
		
		while(true){
			vector<u_int32_t> successors;
			if(forward){
				graph->getSuccessors(current_nodeIndex, 0, successors);
			}
			else{
				graph->getPredecessors(current_nodeIndex, 0, successors);
			}

			if(successors.size() == 0) return false;
			
			current_nodeIndex = successors[0];
			if(graph->nodeIndex_to_unitigIndex(current_nodeIndex) == unitigIndex_source){
				nodeIndex_start = current_nodeIndex;
			}
		}

		current_nodeIndex = nodeIndex_start;

		vector<u_int32_t> pathDummy;
		unordered_map<u_int32_t, vector<u_int32_t>> successorPaths;
		collectPossibleSuccessors(current_nodeIndex, graph, _unitigDatas, forward, 0, 0, pathDummy, successorPaths);

		return successorPaths.size() == 0;
		//PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes, _isNodeImproved, _isPathAlreadyVisitedSourceNodes, _unitigDatas, currentLength, _solvedUnitigs);
		
		//u_int32_t resultType;
		//vector<u_int32_t> nextNodes;
		//pathExplorer.getNextNode(current_nodeIndex, graph, forward, currentDepth+1, resultType, nextNodes, false);
		
	}*/


	/*
    void collectPossibleSuccessors(u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t currentLength, unordered_map<u_int32_t, u_int32_t>& possibleSuccessors){

		//cout << "Source: " << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << endl;
		possibleSuccessors.clear();


		unordered_set<u_int32_t> isVisited;
		unordered_map<u_int32_t, u_int32_t> distance;
        queue<u_int32_t> queue;


        //distanceNode[source_nodeIndex] = 0;
        distance[source_nodeIndex] = 0;
        isVisited.insert(source_nodeIndex);

        queue.push(source_nodeIndex);

		while (!queue.empty()){

            u_int32_t current_nodeIndex = queue.front();
            queue.pop();

			//cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				//cout << "\tPossible: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
				possibleSuccessors[current_nodeIndex] = distance[current_nodeIndex];
				continue;
			}

			
			u_int32_t currentLength = distance[current_nodeIndex];
			if(currentLength > maxLength) continue;

			vector<u_int32_t> successors;
			if(forward){
				graph->getSuccessors(current_nodeIndex, 0, successors);
			}
			else{
				graph->getPredecessors(current_nodeIndex, 0, successors);
			}

			//getSuccessors_unitig(unitigIndex_current, successors);


			for(u_int32_t nodeIndex_successor : successors){

                if (isVisited.find(nodeIndex_successor) != isVisited.end()) continue;

                distance[nodeIndex_successor] = distance[current_nodeIndex] + graph->_nodeLengths[BiGraph::nodeIndex_to_nodeName(nodeIndex_successor)];
                queue.push(nodeIndex_successor);
                isVisited.insert(nodeIndex_successor);

				//cout << "Visited: " << BiGraph::nodeIndex_to_nodeName(nodeIndex_successor) << endl;

            }
			



        }

    }
	*/


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

		unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, from_nodeIndex, from_nodeIndex, _abundanceCutoff_min, visitedNodes, _isNodeImproved, isPathAlreadyVisitedSourceNodes, _unitigDatas, 0, _solvedUnitigs);


		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(from_nodeIndex, graph);
		pathExplorer._visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(from_nodeIndex));

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;

		while(true){
		
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			from_nodeIndex = pathExplorer.getNextNode(from_nodeIndex, graph, forward, currentDepth, resultType, nextNodes, false);
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
		//_exploredNodes.push_back(nodeIndex);

		PathExplorer::clampPrevNodes(_prevNodes, _unitigDatas);
		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		_currentPathLength += graph->_nodeLengths[nodeName];


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

	static bool SuccessorComparator_byDistanceFromSource(const SuccessorData &a, const SuccessorData &b){
		return a._path.size() < b._path.size();
	}

	static void clampPrevNodes(vector<u_int32_t>& prevNodes, vector<UnitigData>& _unitigDatas){

		if(prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1]);
			UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(prevNodes.size() == 0) break;
				u_int32_t nodeIndex = prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(prevNodes[0]);

				const UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				if(PathExplorer::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					prevNodes.erase(prevNodes.begin());
				}

			}
		}

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
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;

	Assembly(): Tool ("asm"){


		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename", true));
		//getParser()->push_back (new OptionOneParam (STR_MINIM_SIZE, "minimizer length", false, "16"));
		//getParser()->push_back (new OptionOneParam (STR_KMINMER_SIZE, "k-min-mer length", false, "3"));
		//getParser()->push_back (new OptionOneParam (STR_DENSITY, "density of minimizers", false, "0.005"));
		//getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", false, ""));
		getParser()->push_back (new OptionOneParam (STR_INPUT_TRUTH, "", false, ""));
		getParser()->push_back (new OptionNoParam (STR_HIFIASM_DEBUG, "", false, false));

	}

	~Assembly(){

	}

	void execute (){

		/*
        BiGraph* graph = new BiGraph(8);
		graph->addEdge_debug(0, 1); //a, b
		graph->addEdge_debug(1, 2); //b, c
		graph->addEdge_debug(1, 4); //b, e
		graph->addEdge_debug(1, 5); //b, f
		graph->addEdge_debug(4, 0); //e, a
		graph->addEdge_debug(4, 5); //e, f
		graph->addEdge_debug(5, 6); //f, g
		graph->addEdge_debug(6, 5); //g, f
		graph->addEdge_debug(2, 6); //c, g
		graph->addEdge_debug(2, 3); //c, d
		graph->addEdge_debug(3, 2); //d, c
		graph->addEdge_debug(3, 7); //d, h
		graph->addEdge_debug(7, 3); //h, d
		graph->addEdge_debug(7, 6); //h, g

		GraphSimplify* g = new GraphSimplify(graph);
		vector<vector<u_int32_t>> sccs;
		g->getStronglyConnectedComponents(sccs);
		cout << sccs.size() << endl;

		for(vector<u_int32_t>& scc : sccs){
			cout << "----" << endl;
			for(u_int32_t unitigIndex : scc){
				vector<u_int32_t> nodes;
				g->getUnitigNodes(g->_unitigs[unitigIndex], nodes);
				for(u_int32_t nodeIndex : nodes){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				}
			}
		}

		exit(1);
		*/

		parseArgs();

		if(getInput()->get(STR_HIFIASM_DEBUG)){
			asmDebug();
			exit(0);
		}
		
		
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
        BiGraph* lol = GfaParser::createBiGraph_lol(_inputDir + "/minimizer_graph_noUnsupportedEdges.gfa_errorFree.gfa", true);
		cout << "lala" << endl;
		GraphSimplify* g = new GraphSimplify(lol);
		cout << "lala" << endl;
		g->debug_writeGfaErrorfree(0, 0, g->_graphSuccessors->nodeName_to_nodeIndex(582450, true));
		exit(1);
		*/

		/*
		GraphSimplify* g = new GraphSimplify(gfa_filename, _inputDir);
		for(size_t i=0; i<1000000000; i++){
			cout << "lala" << endl;
		}
		exit(1);
		*/
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

		
		GraphSimplify* graphSimplify = new GraphSimplify(gfa_filename_noUnsupportedEdges, _inputDir, 0);
		
		
		/*
		graphSimplify->clear(0);
		
		cout << "Collecting truth kminmers" << endl;
		ofstream file_groundTruth_hifiasm_position(_inputDir + "/groundtruth_hifiasm_position.csv");
		ofstream file_groundTruth_hifiasm(_inputDir + "/groundtruth_hifiasm.csv");
		file_groundTruth_hifiasm << "Name,Colour" << endl;
		file_groundTruth_hifiasm_position << "Name,Position" << endl;

		//unordered_set<u_int32_t> visitedNodes;
		//for(auto it : mdbg->_dbg_nodes){
		//	if(_evaluation_hifiasmGroundTruth.find(it.first) == _evaluation_hifiasmGroundTruth.end()) continue;
		//	visitedNodes.insert(it.second._index);
		//}

		unordered_set<u_int32_t> groundTruth_kminmers;
		int founded = 0;
		for(auto it : mdbg->_dbg_nodes){

			const KmerVec& vec = it.first;

			//vec = vec.normalize();
			if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;

			founded += 1;
			groundTruth_kminmers.insert(it.second._index);

			file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
			file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

			unordered_set<u_int32_t> neighbors;
    		graphSimplify->collectNeighbors(BiGraph::nodeName_to_nodeIndex(it.second._index, true), 100, 100, neighbors);
			//graph->collectNeighbors(it.second._index, 100, neighbors, 100, visitedNodes);
			for(u_int32_t nn : neighbors){
				//cout << nn << endl;
				groundTruth_kminmers.insert(BiGraph::nodeIndex_to_nodeName(nn));
			}
			//cout << "n: " << neighbors.size() << endl;
			//cout << n << endl;

    		graphSimplify->collectNeighbors(BiGraph::nodeName_to_nodeIndex(it.second._index, false), 100, 100, neighbors);
			cout << neighbors.size() << endl;
			for(u_int32_t nn : neighbors){
				groundTruth_kminmers.insert(BiGraph::nodeIndex_to_nodeName(nn));
			}

			cout << groundTruth_kminmers.size() << endl;
		}
		cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
		cout << "Found: " << founded << endl;
		GfaParser::rewriteGfa_withNodes(gfa_filename_noUnsupportedEdges, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_kminmers);
		file_groundTruth_hifiasm.close();
		file_groundTruth_hifiasm_position.close();
		
		gzFile file_groundTruth_hifiasm_data = gzopen((_inputDir + "/groundtruth_hifiasm_data.gz").c_str(),"wb");
		
		u_int32_t nbNodes = _unitigDatas.size();
		gzwrite(file_groundTruth_hifiasm_data, (const char*)&nbNodes, sizeof(nbNodes));

		for(u_int32_t nodeName : groundTruth_kminmers){
			gzwrite(file_groundTruth_hifiasm_data, (const char*)&nodeName, sizeof(nodeName));
			u_int32_t nbReads = _unitigDatas[nodeName]._readIndexes.size();
			gzwrite(file_groundTruth_hifiasm_data, (const char*)&nbReads, sizeof(nbReads));
			gzwrite(file_groundTruth_hifiasm_data, (const char*)&_unitigDatas[nodeName]._readIndexes[0], nbReads * sizeof(u_int32_t));
			
		}
		
		gzclose(file_groundTruth_hifiasm_data);

		u_int32_t startingNodeName = -1;
		for(auto& it: _evaluation_hifiasmGroundTruth_nodeNamePosition){
			if(it.second == 0){
				startingNodeName = it.first;
				break;
			}
		}
		cout << "Position 0: " << startingNodeName << endl;
		u_int32_t nodeIndex = graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(startingNodeName, true);
		graphSimplify->debug_writeGfaErrorfree(5, 5, nodeIndex);
		u_int32_t abundance = graphSimplify->getNodeUnitigAbundance(nodeIndex);
		solveBin(nodeIndex, abundance, graphSimplify, 0, true);
		gzclose(_outputContigFile);
		exit(1);
		*/
		









		/*
		computeNodeAbundance(mdbg, gfa_filename_noUnsupportedEdges);
		getUnitigLengths(gfa_filename_unitigs);
		*/


		cout << "Simplifying graph" << endl;
		//graphSimplify->debug_writeGfaErrorfree(0, 0, BiGraph::nodeName_to_nodeIndex(5643, true));
		//exit(1);


		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		graphSimplify->execute(5);

		
		//cout << "remettre initial cleaning" << endl;
		//graphSimplify->compact();
		
		
		//graphSimplify->extractComponent(1224833);





		//exit(1);
		
		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm.csv");
		file_groundTruth_hifiasmContigs << "Name,Colour" << endl;


		//562 (ecoli)
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(558906, false), 87, graphSimplify, 0, true);
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
			
						bool isContigValid = solveBin(nodeIndex, unitigLength._abundance, graphSimplify, binIndex, true);
						if(isContigValid){
							binIndex += 1;
						} 
						else{
							bool isContigValid = solveBin(graphSimplify->nodeIndex_toReverseDirection(nodeIndex), unitigLength._abundance, graphSimplify, binIndex, true);
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
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, readIndex);

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

				/*
				if(nodeName == 1839895){
					cout << nodeName << " " << readIndex << endl;
					cout << vec._kmers[0] << endl;
					cout << vec._kmers[1] << endl;
					cout << vec._kmers[2] << endl;
					exit(1);
				}*/
				
				//if(nodeName == 460064){
				//	cout << nodeName << " " << readIndex << endl;
				//}
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
		
		//cout << "------------------" << endl;
		//cout << PathExplorer::computeSharedReads(_unitigDatas[460067], _unitigDatas[460064]) << endl;
		//cout << "------------------" << endl;


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
		PathExplorer::clampPrevNodes(prevNodes, _unitigDatas);



		//cout << "Prev nodes size: " << prevNodes.size() << endl;
		//cout << "Complete path size: " << nodePath.size() << endl;

		nodePath.push_back(nodeIndex);

		//_nodeLabel[nodeIndex] = pathIndex;
		u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		file_groundTruth << current_unitigIndex << "," << pathIndex << endl;		

		_binnedNodes.insert(current_unitigIndex);
		
		_nbVisitedTimes[nodeIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}

	float computeAbundanceCutoff_min(u_int32_t abundance){
		return abundance / 4.0;
	}

	bool solveBin(u_int32_t source_nodeIndex, u_int32_t abundance, GraphSimplify* graph, int pathIndex, bool performCleaning){


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

		//if(performCleaning){
			cout << "Simplifying graph local" << endl;
			//graph->clear();
			graph->debug_writeGfaErrorfree(abundance, abundanceCutoff_min, source_nodeIndex);
		//}
		//else{
			//graph->clear(0);
			//graph->compact();
		//}
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
		//cout << source_abundance << endl;
		cout << graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance << endl;
		//cout << "PAS BON CA: utiliser successeur puis predeesseur" << endl;
			cout << endl << endl << endl << endl << endl << "----- Forward -------------------------------------------------------------------------------------------------------------------------------------" << endl;
		_iter = 0;
		//u_int32_t source_nodeIndex = graph->_graphPredecessors->nodeName_to_nodeIndex(utg, false);
		PathData pathData = {pathIndex, {}, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		bool pathSolved = solveBin_path(pathData, graph, true);
		if(pathSolved){
			cout << "Path is solve forward (" << pathData.nodePath.size() << ")" << endl;
		}
		
		vector<u_int64_t> supportingReads_forward;
		vector<u_int32_t> nodePath_forward = pathData.nodePath;
		getSupportingReads(nodePath_forward, supportingReads_forward);

		//cout << "to remvoe exit" << endl;
		//exit(1);
		
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
		visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));
		unordered_set<u_int32_t> solvedUnitigs;
		u_int32_t lastNodeIndex = current_nodeIndex;

		int truthPath_indexDirection = 1;
		long truth_pathIndex = -1;
		if(_evaluation_hifiasmGroundTruth_path.size() > 0){
			auto it = find(_evaluation_hifiasmGroundTruth_path.begin(), _evaluation_hifiasmGroundTruth_path.end(), BiGraph::nodeIndex_to_nodeName(lastNodeIndex));
			truth_pathIndex = it - _evaluation_hifiasmGroundTruth_path.begin();
		}
		
			//if(it != _evaluation_hifiasmGroundTruth_path.end()){
			//	index += 1;
			//}
		
		while(true){

			unordered_set<u_int32_t> isPathAlreadyVisitedSourceNodes;
			PathExplorer pathExplorer(pathData.prevNodes, pathData.source_abundance, pathData.source_nodeIndex, current_nodeIndex, pathData._abundanceCutoff_min, visitedNodes, pathData._isNodeImproved, isPathAlreadyVisitedSourceNodes, _unitigDatas, 0, solvedUnitigs);
			//u_int32_t nextNodeIndex = pathExplorer.getNextNode( graph, _unitigDatas, 100000, forward, true);
			
			u_int32_t resultType;
			vector<SuccessorData> nextNodes;
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, forward, 0, resultType, nextNodes, true);
			
			if(resultType == 2){ //Banch solved
				cout << "Path size: " << pathData.nodePath.size() << endl;
			}
			//cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

			//cout << pathData.nodePath.size() << endl;
			//cout << current_nodeIndex << endl;
			//if(_binnedNodes.size() > 100) return false; //DEBUG assemble small fragment
			if(current_nodeIndex == -1) return false; //No more successors, or no branching solution
			
			//if(current_nodeIndex == pathData.source_nodeIndex){ //Path complete
			//	pathData.nodePath.pop_back(); //if the path is solved, the source node exist as first and last element,thus we remove the last one

			//	cout << "Path complete!" << endl;
			//	return true; 
			//}
			//if(current_nodeIndex == -2) return true; //Path complete


			
			for(size_t i=1; i<nextNodes[0]._path.size(); i++){
				current_nodeIndex = nextNodes[0]._path[i];

				if(current_nodeIndex == pathData.source_nodeIndex){ //Path complete
					pathData.nodePath.pop_back(); //if the path is solved, the source node exist as first and last element,thus we remove the last one

					cout << "Path complete!" << endl;
					return true; 
				}

				//cout << "\t" << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
				
				u_int32_t nodeName =BiGraph::nodeIndex_to_nodeName(current_nodeIndex);
				cout << "Add node: " << BiGraph::nodeToString(current_nodeIndex) << " " << (visitedNodes.find(nodeName) != visitedNodes.end());
				
				if(BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == 403217){
					cout << BiGraph::nodeToString(current_nodeIndex) << endl;
					getchar();
				}
				/*
				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					cout << "    " << _evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName];
				}


				if(_evaluation_hifiasmGroundTruth_path.size() > 0){
					long truth_pathIndex_tmp = truth_pathIndex + truthPath_indexDirection;
					if(truth_pathIndex_tmp >= _evaluation_hifiasmGroundTruth_path.size()){
						truth_pathIndex_tmp = 0;
					}
					else if(truth_pathIndex_tmp < 0){
						truth_pathIndex_tmp = _evaluation_hifiasmGroundTruth_path.size()-1;
					}

					if(_evaluation_hifiasmGroundTruth_path[truth_pathIndex_tmp] == nodeName){
						truth_pathIndex = truth_pathIndex_tmp;
						cout << "    OK";
					}
					else{
						if(std::find(_evaluation_hifiasmGroundTruth_path.begin(), _evaluation_hifiasmGroundTruth_path.end(), nodeName) != _evaluation_hifiasmGroundTruth_path.end()){
							while(true){
								//cout << truth_pathIndex << endl;
								truth_pathIndex += truthPath_indexDirection;
								//cout << truth_pathIndex << endl;
								if(truth_pathIndex >= (long)_evaluation_hifiasmGroundTruth_path.size()){
									//cout << "lala" << endl;
									truth_pathIndex = 0;
								}
								else if(truth_pathIndex < 0){
									truth_pathIndex = _evaluation_hifiasmGroundTruth_path.size()-1;
								}

								u_int32_t nodeNameTruth = _evaluation_hifiasmGroundTruth_path[truth_pathIndex];
								if(graph->_isNodeValid2.find(BiGraph::nodeName_to_nodeIndex(nodeNameTruth, true)) != graph->_isNodeValid2.end()){
									cout << "    SKIPPED " << nodeNameTruth << endl;
									getchar();
								}

								//cout << truth_pathIndex << endl;
								if(_evaluation_hifiasmGroundTruth_path[truth_pathIndex] == nodeName){
									cout << "    OK";
									break;
								}


							}
						}
						else{
							cout << "    NOK";
							//getchar();
						}
					}
				}
				*/
				
				cout << endl;
				//cout << "\t\t" << _evaluation_hifiasmGroundTruth_path[truth_pathIndex] << endl;

				binNode(current_nodeIndex, pathData.prevNodes, pathData.nodePath, graph, pathData._index);
				visitedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));

				
				//anti-bug, infinite loop
				//if(_nbVisitedTimes.find(current_nodeIndex) != _nbVisitedTimes.end()){
				//	if(_nbVisitedTimes[current_nodeIndex] > 1000){
				//		return false;
				//	}
				//}

				lastNodeIndex = current_nodeIndex;
			}
			//for(u_int32_t nodeIndex : nextNodes[0]._path){
			//}
			//if(current_nodeIndex == BiGraph::nodeName_to_nodeIndex(539588, true) || BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == 539597){
				//char str [80];
				//scanf ("%79s",str);  
				//cout << "haaaa " << BiGraph::nodeToString(current_nodeIndex) << endl;
				//exit(1);
			//}
			


			//if(BiGraph::nodeIndex_to_nodeName(current_nodeIndex) == 10026){
			//	exit(1);
			//}

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
					
					//cout << (mdbg->_dbg_nodes.find(vec) != mdbg->_dbg_nodes.end()) << endl;
					//cout << (mdbg->_dbg_nodes.find(vec) != mdbg->_dbg_nodes.end()) << endl;
					//if(mdbg->_dbg_nodes[vec]._index == 0){
					//	cout << "lala " << position << endl;
					//}


					if(mdbg->_dbg_nodes.find(vec) != mdbg->_dbg_nodes.end()){
						_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[mdbg->_dbg_nodes[vec]._index].push_back(sequence.getComment());
						_evaluation_hifiasmGroundTruth_path.push_back(mdbg->_dbg_nodes[vec]._index);

						if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(mdbg->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
							_evaluation_hifiasmGroundTruth_nodeNamePosition[mdbg->_dbg_nodes[vec]._index] = position;
						}
						//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
					}
					else{
						//cout << "Not found position: " << position << endl;
						_evaluation_hifiasmGroundTruth_path.push_back(-1);
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
		//exit(1);
		//gzclose(file);
		
		//exit(1);
	}

	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;

	void asmDebug(){
		


		if(_truthInputFilename != ""){
			string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
			MDBG* mdbg = new MDBG(_kminmerSize);
			mdbg->load(mdbg_filename);

			extract_truth_kminmers(mdbg);

			delete mdbg;
		}


		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		gfa_filename_noUnsupportedEdges += "_errorFree.gfa"; //Comment to redo cleaning process
		cout << gfa_filename_noUnsupportedEdges << endl;



		string gfa_filename = _inputDir + "/minimizer_graph.gfa_groundTruth_hifiasm.gfa";

		gzFile file_groundTruth_hifiasm_data = gzopen((_inputDir + "/groundtruth_hifiasm_data.gz").c_str(),"rb");
		
		u_int32_t nbNodes;
		gzread(file_groundTruth_hifiasm_data, (char*)&nbNodes, sizeof(nbNodes));
		_unitigDatas.resize(nbNodes);

		while(true){
			
			u_int32_t nodeName;
			gzread(file_groundTruth_hifiasm_data, (char*)&nodeName, sizeof(nodeName));

			if(gzeof(file_groundTruth_hifiasm_data)) break;

			u_int32_t nbReads;
			vector<u_int64_t> readIndexes;
			gzread(file_groundTruth_hifiasm_data, (char*)&nbReads, sizeof(nbReads));

			
			_unitigDatas[nodeName]._readIndexes.resize(nbReads);
			gzread(file_groundTruth_hifiasm_data, (char*)&_unitigDatas[nodeName]._readIndexes[0], nbReads * sizeof(u_int32_t));
			//_unitigDatas[nodeName] = readIndexes;

		}

		gzclose(file_groundTruth_hifiasm_data);


		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		GraphSimplify* graphSimplify = new GraphSimplify(gfa_filename_noUnsupportedEdges, _inputDir, nbNodes);

		//C6
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(954857, false), 40, graphSimplify, 0, false);
		
		//C2
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(28286, true), 40, graphSimplify, 0, false);

		//c3
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(7685, true), 83, graphSimplify, 0, false);

		//c4 = spaced high complexity area (need to identify long unitig + long contig, then check if BFS leeds to single sink) 
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(25697, true), 650, graphSimplify, 0, false);

		//c4 k10
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(19971, true), 460, graphSimplify, 0, false);

		//c5 (low abundant genome with gap)
		//29885 (10 abundance max)

		//c7 (complex region, 2 strains)
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(2038, true), 70, graphSimplify, 0, true);
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(302618, false), 70, graphSimplify, 0, true);
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(537443, true), 70, graphSimplify, 0, true);

		//c8
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(426182, true), 160, graphSimplify, 0, true);

		//c9 (low abundant genome ~10)
		//Start node name: 42693

		//c10
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(236892, true), 20, graphSimplify, 0, true);

		//c11
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(83861, true), 41, graphSimplify, 0, true);

		//c13
		//Source: 3080 3080+
		//Source abundance: 238

		//c14
		//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(163348, true), 500, graphSimplify, 0, true);

		//c15
		//Source: 91818 91818+
		//Source abundance: 7



		//c16
		solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(18437, true), 31, graphSimplify, 0, true);

		file_groundTruth.close();
		gzclose(_outputContigFile);
	}

};



#endif
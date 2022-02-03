
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
		- SUperbubble et COmpelx area detection: on veut enlever les branches qui n'appartiennent pas au contexte actuelle (actuellement on verifie juste que les noeuds de l'area partage des reads avec le dernier noeud visité mais ça fait qu'on est limité à une certaine taille d'area)
		- Take bubble: si on prend un chemin arbitraire dans une bubble, il faut marqué le/les autre chemins comme visité pour ne qu'ils soient détecté comme successeurs ensuite
		- Complex area: quand on skip une complex area, il faut ajouté ses nodes dans _binnedNodes pour ne pa que l'assemblage puissent partir de'un de ses unitigs
		- cutoff plus intelligent: utiliser l'abondance de tous les long unitigs pour filtrer les erreurs, plutot qu'un seul unitig, vu que l'abondance varie sur l'ensemble du genome
			- ou utiliser un filtre de base assez bas, puis filtrer on the fly, mais ça va reduire la taille des unitigs de base
		- getOverlap: on devrait bannir cette methode et toujour utiliser getSuccessors_overlap pour obtenir les infos des edge en plus des successeur
		- Complex area detection issue: quand on lance un ecoli avec 40 de coverage, il y a des pb d'assemblage a cause des erreurs,
Gros changement non testé:
	- attention on a ajouter les visitedNodes dans nodeExplored collectpossiblesuccessor (ça va etre a test) qu'on explore bien tout malgré ça"

	- cutoff onthe fly = assemblyAbundance / 2, remettre superbubble detection, isTip on he fly egalement ?"
	- dans detection superbubble classique, on pourrait jouter la gestion des cycle des la phase de cleaning
	- minimizer contig: pas besoin de recuperer/ecrire les supporting reads
	- ToBasespace: recuperer les sequenceModel et variant en une seule passes sur les data, peut etre plus encessaire avec l'approche multi k?
	- Kminmer / kmer parser: verifier qu'on a la taille de seuqnece suffisante (pour les kmer) et nb minimizers suffisant (pour kminmer) pour les générer (k+1)
	- Gfa format: créer notre propre format binaire, mais attention gfa est pratique pour debug
	- Palindrome: augmenter la densité au endroit des positions bannis, pour conserver la meme taille d'overlaps (idée seb)
	- essayer open syncmer a la place de minimzier ?
	- Superbubble qui n'a qu'une branche variante (l'autre coté va direct de la source vers sink), actuellement on supprime la branche variante mais on pourrait voir laquel garder en fonction du score checkm ou quast
	- Chimeric reads: a mettre dans readselection direct
	- graphSImplify: indexer des nodeNames plutot que des nodeIndex (_isNodeValid), car on ne conserve jamais qu'un seul des deux coté d'un nodename
	- metaflye: plein de bonne idée a check (identification superbubble, roundabout)
	- ToBasespace: correction: on peut appliquer la correction d'un kminmer dés qu'on a atteint les 20 variants, comme ça on peut free la mémoire de ces variants dés que possible 
	- Unsupported edges: on peut distinguer les fake edges (qui n'existe dans aucun reads), des bons (provenant des especes rares), en chackant l'abondance de l'edge. Les fakes edge devrait avoir une abondance extrement elevé

	REPRISE:
		- on veut detecter le plus de bubble possible
		- on veut verifier qu'un branchement est une bubble avant quoi que ce soit pour accelerer le processus d'assemblage
*/
//MetaBinner: https://www.biorxiv.org/content/10.1101/2021.07.25.453671v1

/*
Solve repeat: attention startNode et endNode d'un meme unitig pourrait etre un meme successor/predecessor si l'unitig est couvert par plusieurs long reads (il faut choisir le bon dans ce cas)
ToBasespaceOnTheFly: en théorie, on peut reconstruire la partie manquante des contigs basespace (les k-1 premiers nodes), en utilisant les predecesseurs de ce noeud de départ
*/

//:/mnt/chris-native/chris/Projects/AD_HiFi
//mnt/gpfs/Hackathon/meren_hifi

#ifndef MDBG_METAG_ASSEMBLY2
#define MDBG_METAG_ASSEMBLY2

//#define PRINT_ADDED_NODES
//#define PRINT_PREV_RANK

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
#include "eval/ContigStatistics.hpp"
#include "toBasespace/ToBasespaceOnTheFly.hpp"
#include "contigFeatures/ContigFeature.hpp"

const u_int32_t LONG_UNITIG_LENGTH = 10000;
const u_int32_t MIN_SUPPORTED_READS = 2;
const u_int32_t MIN_REPEAT_SIDE_LENGTH = 10000;









class Assembly2 : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	vector<UnitigData> _unitigDatas;

	string _filename_solidKminmers;
	string _filename_inputContigs;
	string _filename_outputContigs;
	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;
	gzFile _outputContigFile_complete;
	MDBG* _mdbg;
	GraphSimplify* _graph;
	ToBasespaceOnTheFly _toBasespace;
	ContigFeature _contigFeature;
	string _gfaFilename;

	Assembly2(): Tool (){

	}

	~Assembly2(){

	}

	void execute (){

		_contigFeature.loadAbundanceFile_nodename(_filename_abundance);
		//_contigFeature.computeFastaComposition("/home/gats/workspace/run/overlap_test_AD_k7/binByCutoff/component_16_unitigs.fasta");
		
		//exit(1);

		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		_outputContigFile_complete = gzopen(_outputFilename_complete.c_str(),"wb");

		loadGraph();
		execute_binning();
		//execute_detectSpecies_byCutoff();
		
		gzclose(_outputContigFile);
		gzclose(_outputContigFile_complete);
	}

	void parseArgs(int argc, char* argv[]){



		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_DEBUG, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_UNITIG_NT, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_UNITIG_CLUSTER, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""));



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			_inputFilename_unitigNt = result[ARG_INPUT_FILENAME_UNITIG_NT].as<string>();
			_inputFilename_unitigCluster = result[ARG_INPUT_FILENAME_UNITIG_CLUSTER].as<string>();
			_debug = result[ARG_DEBUG].as<bool>();
			_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}


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
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.gz";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;


		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}


		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0);
		_graph = graphSimplify;
		
		
		//Generate unitigs
		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		delete _mdbg;
		cout << "done" << endl;
	


		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		//graphSimplify->execute(5, _kminmerSize);
		//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);
		_graph->debug_writeGfaErrorfree(500, 500, -1, _kminmerSize, false, true, false, _unitigDatas);




		_toBasespace.create(_inputDir);


	}

	void execute_detectSpecies_byCutoff(){



		//size_t globalAbundanceCutoff_min = 3;
		
		

		//detectIsolatedSpecies(gfa_filename, toBasespace);
		//exit(1);
		u_int64_t nbUnitigs = 0;
		unordered_set<u_int32_t> writtenNodeNames;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_unitigName;
		

		string clusterDir = _inputDir + "/" + "binByCutoff";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;

		u_int32_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;


		for(Unitig& unitig : _graph->_startingUnitigstest){

			
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);

			/*
			bool isValid = false;
			for(u_int32_t nodeIndex : unitig._nodes){
				if(_hifiasm_startingNodenames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _hifiasm_startingNodenames.end()){
					isValid = true;
				}
			}
			if(!isValid) continue;
			*/

			float cutoff = (unitig._abundance*0.9)*0.2;


			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff || reloadState){
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
				//_graph->loadState2(0, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
				reloadState = false;
			}
			

			//float cutoff = unitig._abundance*0.1;
			//_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);

        	//unordered_set<u_int32_t> component;
			//for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			//	component.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
			//}

        	unordered_set<u_int32_t> component;
        	_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), component);

			if(component.size() <= 2) continue;

			vector<u_int32_t> cluster;
			//unordered_set<u_int32_t> validNodes;

			for(u_int32_t unitigIndex : component){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					cluster.push_back(nodeName);
					//validNodes.insert(nodeName);
				}
			}

			bool isNewCluster = true;
			std::sort(cluster.begin(), cluster.end());

			//cout << endl << "Unitig: " << unitig._index << endl;

			for(const vector<u_int32_t>& existingCluster : clusters){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << "Existing: " << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.33){
					isNewCluster = false;
					break;
				}
			}
			
			for(const vector<u_int32_t>& existingCluster : clustersValid){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << "Valid: " << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.6){
					isNewCluster = false;
					break;
				}
			}

			if(!isNewCluster) continue;
		

			//cout << clusterIndex << " " << component.size() << endl;
			
			/*
			if(validNodes.size() > 0){
				
				reloadState = true;

				cout << endl << "Start new cluster: " << endl;
				
				const string& filenameComponentNodes = clusterDir + "/component_" + to_string(clusterIndex) + "_nodes.csv";
				ofstream fileComponenetNodes(filenameComponentNodes);
				fileComponenetNodes << "Name,Colour" << endl;
				for(u_int32_t nodeName : visitedNodes){
					fileComponenetNodes << nodeName << "," << clusterIndex << endl;
				}
				fileComponenetNodes.close();
				
				//exit(1);
				clusters.push_back(cluster);
			*/

			clusters.push_back(cluster);
			reloadState = true;

			bool hasValidContamination = false;
			float componentCompletneness = 0;
			vector<u_int32_t> componentNodesValid;

			/*
			//for(float cutoffRate : {0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1}){			
			for(float cutoffRate : {0.5, 0.4, 0.3, 0.2, 0.1}){
					
				cout << "Cutoff rate: " << cutoffRate << endl;

				float cutoff = (unitig._abundance*0.9)*cutoffRate;
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);


				unordered_set<u_int32_t> component;
				_graph->getConnectedComponent(unitig._startNode, component);

				//unordered_set<u_int32_t> component;
				//for(u_int32_t nodeIndex : _graph->_isNodeValid2){
				//	component.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
				//}

				if(component.size() <= 1) continue;





				vector<u_int32_t> cluster;
				//unordered_set<u_int32_t> validNodes;

				for(u_int32_t unitigIndex : component){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						cluster.push_back(nodeName);
						//validNodes.insert(nodeName);
					}
				}

				std::sort(cluster.begin(), cluster.end());

				bool isNewCluster = true;
				for(const vector<u_int32_t>& existingCluster : clustersValid){
					float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
					cout << jaccardDistance << endl;
					//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
					//cout << sharedRate << endl;
					if(jaccardDistance < 0.6){
						isNewCluster = false;
						break;
					}
				}
				if(!isNewCluster) break;

				*/
				/*
				cout << "subgraph extract" << endl;
				//_graph->compact(false, _unitigDatas);
				_graph->extractReadpathSubgraph(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(963565, true)));
				exit(1);
				*/

				//cout << _graph->_unitigs.size() << endl;
				for(u_int32_t unitigIndex : component){
					const Unitig& u = _graph->_unitigs[unitigIndex];
					if(u._length < MIN_REPEAT_SIDE_LENGTH){
						for(u_int32_t nodeIndex : u._nodes){
							_graph->_repeatedNodenames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
						}
					}
				}
				//_graph->compact(false, _unitigDatas);
				//cout << _graph->_unitigs.size() << endl;
				//getchar();

				const string& filenameUnitigColor= clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
				ofstream fUnitigColor(filenameUnitigColor);
				fUnitigColor << "Name,Colour" << endl;

				cout << "Generate contigs" << endl;
				const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
				//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
				ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

				vector<u_int32_t> componentNodes;

				//u_int32_t contigIndex = 0;
				unordered_set<u_int32_t> writtenUnitigs;
				for(u_int32_t unitigIndex : component){
					
					const Unitig& u = _graph->_unitigs[unitigIndex];
					u_int32_t contigIndex = u._index;
					//if(u._length < 10000) continue;
					
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

					string unitigSequence;
					_toBasespace.createSequence(u._nodes, unitigSequence);

					string header = ">ctg" + to_string(contigIndex) + '\n';
					basespaceUnitigFile << header;
					//gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
					unitigSequence +=  '\n';
					basespaceUnitigFile << unitigSequence;
					//gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());
					
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						componentNodes.push_back(nodeName);
						fUnitigColor << nodeName << "," << contigIndex << endl;
					}

					//cout << unitigSequence.size() << endl;

					contigIndex += 1;
				}

				fUnitigColor.close();
				//gzclose(basespaceUnitigFile);
				basespaceUnitigFile.close();

				
				string command_annotation = "python3 /home/gats/workspace/tools/binner/src/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor;// + " --isSingleSpecies";
				cout << command_annotation << endl;
				int ret = system(command_annotation.c_str());
				if(ret != 0){
					cerr << "Command failed: " << ret << endl;
					exit(ret);
				}
				
				
				
				
				//char isSingleSpeciesResult;
				//ifstream fileIsSingleSpecies(basespaceUnitigFilename + "_isSingleSpecies.txt");
				//fileIsSingleSpecies.get(isSingleSpeciesResult);
				//fileIsSingleSpecies.close();

				//cout << isSingleSpeciesResult << endl;


				indexScgClusters(basespaceUnitigFilename);
				float contamination;
				float completeness;
				computeContamination(component, completeness, contamination);
				cout << "Completeness: " << completeness << endl;
				cout << "Contamination: " << contamination << endl;

				//computeScgCentrality(clusterDir + "/contaminated_" + to_string(contaminatedIndex) + ".gfa.centrality.csv", component);
				//exit(1);

				//if(contamination <= 0.01){

				//contamination = 1;

				if(contamination < 0.05){


					componentNodesValid = componentNodes;
					componentCompletneness = completeness;
					hasValidContamination = true;
				}
				else{

					

					unordered_set<u_int32_t> validNodes;


					indexScgClusters(basespaceUnitigFilename);
					writeScgDebug(basespaceUnitigFilename, component);

					//vector<UnitigCentrality> unitigCentralities;
					//computeScgCentrality(outputFilename, component, 0, unitigCentralities);
					
					u_int32_t removedUnitiIndex = -1;
					unordered_set<u_int32_t> solvedRepeats;

					size_t pass = 0;
					u_int32_t source_unitigIndex = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
					//keepUnitigs.insert(source_unitigIndex);
					//keepUnitigs.insert(_graph->unitigIndex_toReverseDirection(source_unitigIndex));

					_unsolvableRepeat.clear();

					while(true){
						
						source_unitigIndex = _graph->nodeIndex_to_unitigIndex(unitig._startNode);

						_graph->getConnectedComponent_unitig(source_unitigIndex, component);
						cout << "Rebloch: " << source_unitigIndex << endl;
						cout << "Component: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << component.size() << endl;
						//_graph->saveGraph2(_inputDir + "/" + "subgraph.gfa", component);
						//getchar();


						vector<u_int32_t> componentNodes;
						validNodes.clear();
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								validNodes.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
								componentNodes.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
							}
						}

						string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
						//string outputFilename = clusterDir + "/contaminated_" + to_string(contaminatedIndex) + "_pass_" + to_string(pass) + ".gfa";
						GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
						break;


						float contamination;
						float completeness;
						computeContamination(component, completeness, contamination);
						cout << "Completeness: " << completeness << endl;
						cout << "Contamination: " << contamination << endl;

						//contamination = 1;


						if(contamination < 0.05){
							
							componentNodesValid = componentNodes;
							componentCompletneness = completeness;
							hasValidContamination = true;
							break;
						}



						
						//if(unitigCentralities.size() == 0){
						//	cout << "No more unitig to process" << endl;
						//	break;
						//}

						//const UnitigCentrality& unitigCentrality = unitigCentralities[unitigCentralities.size()-1];
						//unitigCentralities.pop_back();

						//u_int32_t unitigIndexMostCentral = unitigCentrality._unitigIndex;
						

						u_int32_t unitigIndexMostCentral = -1;

						while(true){
							vector<UnitigCentrality> unitigCentralities;
							unitigIndexMostCentral = computeScgCentrality(outputFilename, component, pass, unitigCentralities);
							if(unitigIndexMostCentral == -1){
								cout << "No more unitig to process" << endl;
								break;
							}


							//if(component.find(unitigIndexMostCentral) == component.end()) continue;
							//if(keepUnitigs.find(unitigIndexMostCentral) != keepUnitigs.end()) continue;



							break;
						}

						if(unitigIndexMostCentral == -1){
							cout << "No more unitig to process" << endl;
							break;
						}

						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndexMostCentral]._nodes){
							_unsolvableRepeat.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
						}
						
						//cout << unitigCentrality._unitigIndex << " " << unitigCentrality._centrality << endl;

						u_int32_t unitigIndexMostCentralRev = _graph->unitigIndex_toReverseDirection(unitigIndexMostCentral);
						
						solveRepeat(unitigIndexMostCentral, component);


						pass += 1;



					}


					contaminatedIndex += 1;
					
				}

			//}

			//if(!hasValidContamination){
			//	clusters.push_back(cluster);
			//}

			if(hasValidContamination && componentCompletneness > 0.5){

				

				std::sort(componentNodesValid.begin(), componentNodesValid.end());
				clustersValid.push_back(componentNodesValid);



				cout << "Found valid component: " << componentCompletneness << endl;

				const string& filenameHifiasm = clusterDir + "/component_" + to_string(clusterIndex) + "_hifiasm.csv";
				unordered_set<string> hifiasmWrittenUnitigs;
				ofstream fileHifiasm(filenameHifiasm);
				fileHifiasm << "Name,Colour" << endl;

				unordered_set<u_int32_t> componentNodes;
				for(u_int32_t nodeName : componentNodesValid){
					componentNodes.insert(nodeName);

					fileComponentNodeAll << nodeName << "," << clusterIndex << endl;

					if(_truthInputFilename != ""){
						

						if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
							for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
								if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
								fileHifiasm << unitigName << "," << clusterIndex << endl;
								fileHifiasmAll << unitigName << "," << clusterIndex << endl;
								hifiasmWrittenUnitigs.insert(unitigName);
							}
						}

					}

				}

				fileHifiasm.close();

				string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
				GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, componentNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
		





				//clusterIndex += 1;
				//getchar();
			}

			
			clusterIndex += 1;
		}



		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();

	}


	unordered_set<u_int32_t> _unsolvableRepeat;

					/*
					vector<vector<u_int32_t>> components;
					_graph->getStronglyConnectedComponents(unitig._index, true, components);
					ofstream file_scc1(clusterDir + "/component_" + to_string(clusterIndex) + ".gfa" + ".scc_1.csv");
					file_scc1 << "Name,Colour" << endl;
					int sccCluster = 0;
					for(vector<u_int32_t>& component : components){
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								file_scc1 << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << sccCluster << endl;
							}
						}
						sccCluster += 1;
					}
					file_scc1.close();

					_graph->getStronglyConnectedComponents(_graph->unitigIndex_toReverseDirection(unitig._index), true, components);
					ofstream file_scc2(clusterDir + "/component_" + to_string(clusterIndex) + ".gfa" + ".scc_2.csv");
					file_scc2 << "Name,Colour" << endl;
					sccCluster = 0;
					for(vector<u_int32_t>& component : components){
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								file_scc2 << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << sccCluster << endl;
							}
						}
						sccCluster += 1;
					}
					file_scc2.close();
					*/


						/*
						size_t pass = 0;
						while(true){

							outputFilename = clusterDir + "/contaminated_" + to_string(contaminatedIndex) + "_pass_" + to_string(pass) + ".gfa";
							GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			

							cout << endl << "Pass: " << pass << endl;
							u_int32_t unitigIndexMostCentral = computeScgCentrality(outputFilename, component, pass);

							if(unitigIndexMostCentral == -1) break; //no more thing to do

							u_int32_t unitigIndexMostCentralRev = _graph->unitigIndex_toReverseDirection(unitigIndexMostCentral);
							
							cout << "Most central unitig: " << unitigIndexMostCentral << " " << unitigIndexMostCentralRev << endl;
							cout << (component.find(unitigIndexMostCentral) != component.end()) << " " << (component.find(unitigIndexMostCentralRev) != component.end()) << endl;
							cout << _graph->_unitigs[unitigIndexMostCentral]._startNode << endl;

							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndexMostCentral]._nodes){
								cout << nodeIndex << " " << (_graph->_isNodeValid2.find(nodeIndex) != _graph->_isNodeValid2.end()) << endl;

								if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
								_graph->_isNodeValid2.erase(nodeIndex);
								_graph->_nodeToUnitig.erase(nodeIndex);
								cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;

							}
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndexMostCentralRev]._nodes){
								cout << nodeIndex << " " << (_graph->_isNodeValid2.find(nodeIndex) != _graph->_isNodeValid2.end()) << endl;
								
								if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
								
								_graph->_isNodeValid2.erase(nodeIndex);
								_graph->_nodeToUnitig.erase(nodeIndex);
								cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
							}

							_graph->_unitigs[unitigIndexMostCentral]._startNode = -1;
							_graph->_unitigs[unitigIndexMostCentral]._endNode = -1;
							_graph->_unitigs[unitigIndexMostCentralRev]._startNode = -1;
							_graph->_unitigs[unitigIndexMostCentralRev]._endNode = -1;
							
							cout << "Disable unitig: " << unitigIndexMostCentral << " " << unitigIndexMostCentralRev << endl;
							
        					//unordered_set<u_int32_t> component;
							cout << "Get connected componenet from: " << unitig._startNode << endl;
        					//cout_graph->getConnectedComponent_unitig(unitig._index, component);
							cout << _graph->nodeIndex_to_unitigIndex(unitig._startNode) << " " << _graph->_unitigs.size() << endl;
        					_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), component);
							cout << "Component: " << component.size() << endl;
	
							//unordered_set<u_int32_t> validNodes;
							validNodes.clear();
							for(u_int32_t unitigIndex : component){
								for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
									validNodes.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
								}
							}



							float contamination;
							float completeness;
							computeContamination(component, completeness, contamination);
							cout << "Completeness: " << completeness << endl;
							cout << "Contamination: " << contamination << endl;

							//if(contamination < 0.05) break;

							pass += 1;

							//getchar();
						}
						*/


	struct UnitigScgCluster{
		u_int32_t _unitigIndex;
		u_int32_t _scgIndex;
		u_int32_t _scgCluster;
	};

	struct ScgCluster{
		u_int32_t _scgIndex;
		u_int32_t _scgCluster;
	};

	unordered_map<u_int32_t, vector<ScgCluster>> _scgClusters;

	void indexScgClusters(const string& contigFilename){

		_scgClusters.clear();

        vector<string>* fields = new vector<string>();

    	const string& scgClusterFilename = contigFilename + ".scgCluster.txt";

		string line;
		ifstream infile(scgClusterFilename);

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			u_int32_t unitigIndex = stoull((*fields)[0]);
			u_int32_t scgIndex = stoull((*fields)[1]);
			u_int32_t clusterIndex = stoull((*fields)[2]);

			_scgClusters[unitigIndex].push_back({scgIndex, clusterIndex});

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;
		}

		delete fields;
		infile.close();
	}

	void writeScgDebug(const string& contigFilename, unordered_set<u_int32_t>& component){

		unordered_set<u_int32_t> scgs;
		scgs.insert(0);
		/*
		for(const auto& it : _scgClusters){
			const vector<ScgCluster>& scgClusters = it.second;
			for(const ScgCluster& scgCluster : scgClusters){
				scgs.insert(scgCluster._scgIndex);
			}
		}
		*/
				
		ofstream file(contigFilename + ".scg.csv");
		string header = "Name";

		for(size_t i=0; i<scgs.size(); i++){
			header += ",Color";//",SCG_" + to_string(i);
		}
		file << header << endl;

		unordered_set<u_int32_t> writtenUnitigs;
				
		for(u_int32_t unitigIndex : component){
					
			const Unitig& u = _graph->_unitigs[unitigIndex];
			u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			if(_scgClusters.find(unitigIndex) == _scgClusters.end()){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					string line = "";
					line += to_string(nodeName);
					for(size_t i=0; i<scgs.size(); i++){
						line += ",-1";
					}
					file << line << endl;
				}

			}

			unordered_map<u_int32_t, u_int32_t> unitigScgs;
			for(const ScgCluster& scgCluster : _scgClusters[unitigIndex]){
				unitigScgs[scgCluster._scgIndex] = scgCluster._scgCluster;
			}

			for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				string line = "";
				line += to_string(nodeName);
				for(size_t i=0; i<scgs.size(); i++){
					if(unitigScgs.find(i) == unitigScgs.end()){
						line += ",-1";
					}
					else{
						line += "," + to_string(unitigScgs[i]);
					}
				}
				line += "\n";
				file << line << endl;
			}

		}

		file.close();

	}

	void computeContamination(unordered_set<u_int32_t>& unitigs, float& completeness, float& contamination){

		completeness = 0;
		contamination = 1;

		unordered_map<u_int32_t, vector<u_int32_t>> scgCounts;

		for(u_int32_t unitigIndex : unitigs){

			if(_scgClusters.find(unitigIndex) == _scgClusters.end()) continue;

			const vector<ScgCluster>& scgClusters = _scgClusters[unitigIndex];

			for(const ScgCluster& scgCluster : scgClusters){

				if(scgCounts.find(scgCluster._scgIndex) == scgCounts.end()){
					scgCounts[scgCluster._scgIndex].push_back(scgCluster._scgCluster);
				}
				else{
					vector<u_int32_t>& clusters = scgCounts[scgCluster._scgIndex];
					if(std::find(clusters.begin(), clusters.end(), scgCluster._scgCluster) == clusters.end()){
						clusters.push_back(scgCluster._scgCluster);
					}
				}
			}
		}

		u_int32_t nbContaminatedScg = 0;
		for(auto& it : scgCounts){

			u_int32_t scgIndex = it.first;
			const vector<u_int32_t>& scgCluters = it.second;

			if(scgCluters.size() > 1){
				nbContaminatedScg += 1;
			}

			cout << scgIndex << " " << scgCluters.size() << endl;
		}

		completeness = ((float) scgCounts.size()) / ((float) 107);
		contamination = ((float) nbContaminatedScg) / ((float) 107); //((float) scgCounts.size());
	}

	ofstream file_repeat;

	void solveRepeat(u_int32_t source_unitigIndex, unordered_set<u_int32_t>& component){

		bool isRepeatSolved = false;
		cout << "Solving repeat: " << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[source_unitigIndex]._startNode) << endl;

		_removedUnitigs.clear();

		file_repeat.open("/home/gats/workspace/run/repeat.csv");
		file_repeat << "Name,Colour" << endl;

		unordered_set<u_int32_t> possiblePredecessors;
		unordered_set<u_int32_t> possibleSuccessors;

		unordered_map<u_int32_t, vector<u_int64_t>> unitigDatas2;
		indexUnitiData2(source_unitigIndex, unitigDatas2);

        //cout << "-----" << endl;
        //const UnitigData& source_readIndexes = _unitigDatas2[source_unitigIndex];
		//ofstream file_scc("/home/gats/workspace/run/subgraph.csv");
		//file_scc << "Name,Colour" << endl;

		for(u_int32_t nodeIndex : _graph->_unitigs[source_unitigIndex]._nodes){
			file_repeat << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
		}

		solveRepeat_collectSides(source_unitigIndex, true, unitigDatas2, possibleSuccessors);
		solveRepeat_collectSides(source_unitigIndex, false, unitigDatas2, possiblePredecessors);

		file_repeat.close();

		_graph->saveGraph2(_inputDir + "/" + "subgraph.gfa", component);
		cout << "\tRepeat written to csv: " << possibleSuccessors.size() << " " << possiblePredecessors.size() << endl;
		getchar();

		//unordered_set<u_int32_t> processedNodeNames;
		const vector<u_int64_t>& source_readIndexes = unitigDatas2[source_unitigIndex];
		
		vector<RepeatData> repeatDatas;

		for(u_int32_t successor_nodeName : possibleSuccessors){


			for(u_int32_t predecessor_nodeName : possiblePredecessors){

				vector<u_int64_t> sharedReads_succ;
				vector<u_int64_t> sharedReads_pred;
				vector<u_int64_t> sharedReads;
				vector<u_int64_t> sharedReadsLala;
            	Utils::collectSharedReads(source_readIndexes, _unitigDatas[successor_nodeName]._readIndexes, sharedReads_succ);
            	Utils::collectSharedReads(source_readIndexes, _unitigDatas[predecessor_nodeName]._readIndexes, sharedReads_pred);
            	Utils::collectSharedReads(sharedReads_succ, sharedReads_pred, sharedReads);
            	Utils::collectSharedReads(_unitigDatas[successor_nodeName]._readIndexes, _unitigDatas[predecessor_nodeName]._readIndexes, sharedReadsLala);

				cout << "\t" << successor_nodeName << " " << predecessor_nodeName << " " << sharedReadsLala.size() << " " << sharedReads.size() << endl;
			
				if(successor_nodeName == predecessor_nodeName) continue;
				if(sharedReads.size() < MIN_SUPPORTED_READS) continue;
				//if(processedNodeNames.find(successor_nodeName) != processedNodeNames.end()) continue;
				//if(processedNodeNames.find(predecessor_nodeName) != processedNodeNames.end()) continue;

				repeatDatas.push_back({successor_nodeName, predecessor_nodeName, sharedReads.size()});

			}



		}

		if(repeatDatas.size() >= 1){

			sort(repeatDatas.begin(), repeatDatas.end(), RepeatComparator_bySupportedReads);

			for(const RepeatData& repeatData : repeatDatas){
				cout << "\t\t" << repeatData._nodeName_successor << " -> " << repeatData._nodeName_predecessor << " " << repeatData._nbSupportingReads << endl;
			}

			const RepeatData& repeatData = repeatDatas[0];
			u_int32_t nodeName_successor = repeatData._nodeName_successor;
			u_int32_t nodeName_predecessor = repeatData._nodeName_predecessor;

			
			cout << "\t\tSolving repeat: " << nodeName_successor << " -> " <<  BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[source_unitigIndex]._startNode) << " -> " << nodeName_predecessor << endl; //<< " " << repeatData._nbSupportingReads << endl;
			solveRepeat_disconnectNode(nodeName_successor, true);
			//solveRepeat_disconnectNode(nodeName_successor, false);
			solveRepeat_disconnectNode(nodeName_predecessor, true);
			//solveRepeat_disconnectNode(nodeName_predecessor, false);
			solveRepeat_connectNode(nodeName_successor, nodeName_predecessor);

			//const Unitig& unitigLala1 = _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(nodeName_successor, true));
			//processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(unitigLala1._startNode));
			//processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(unitigLala1._endNode));
			//const Unitig& unitigLala2 = _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(nodeName_predecessor, true));
			//processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(unitigLala2._startNode));
			//processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(unitigLala2._endNode));

			isRepeatSolved = true;
			//return;
		}

		
		unordered_map<u_int32_t, vector<ScgCluster>> scgIndexTmp;
		for(const auto& it : _scgClusters){

			u_int32_t unitigIndex = it.first;
			const vector<ScgCluster>& scgClusters = it.second;

			for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
				scgIndexTmp[BiGraph::nodeIndex_to_nodeName(nodeIndex)] = scgClusters;
			}
		}

		//cout << "Source unitig: " << source_unitigIndex << endl;
		for(u_int32_t unitigIndex : _removedUnitigs){
			cout << "Disabling unitig: " << unitigIndex << endl;
			_graph->disableUnitig(unitigIndex);
		}

		_graph->compact(true, _unitigDatas);


		cout << _scgClusters.size() << endl;
		
		_scgClusters.clear();
		unordered_set<u_int32_t> usedUnitigIndex;

		for(const auto& it : scgIndexTmp){

			u_int32_t nodeName = it.first;
			const vector<ScgCluster>& scgClusters = it.second;

			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(nodeName, true));
			u_int32_t unitigIndex_rev = _graph->unitigIndex_toReverseDirection(unitigIndex);

			if(usedUnitigIndex.find(unitigIndex) != usedUnitigIndex.end()){

			}
			else if(usedUnitigIndex.find(unitigIndex_rev) != usedUnitigIndex.end()){
				unitigIndex = unitigIndex_rev;
			}
			else{
				usedUnitigIndex.insert(unitigIndex);
			}

			vector<ScgCluster>& indexedScgClusters = _scgClusters[unitigIndex];

			for(const ScgCluster& scgCluster : scgClusters){

				bool hasScg = false;

				for(const ScgCluster& existingScgCluster : indexedScgClusters){
					if(scgCluster._scgIndex == existingScgCluster._scgIndex && scgCluster._scgCluster == existingScgCluster._scgCluster){
						hasScg = true;
						break;
					}
				}

				if(!hasScg){
					_scgClusters[unitigIndex].push_back(scgCluster);
				}
			}
		}

		cout << _scgClusters.size() << endl;

		getchar();

		if(isRepeatSolved){
			_unsolvableRepeat.clear();
		}
		else{
			
			for(u_int32_t successor_nodeName : possibleSuccessors){
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(successor_nodeName, false);
				u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
				_graph->_isUnitigInvalid_tmp.insert(unitigIndex);
				_graph->_isUnitigInvalid_tmp.insert(_graph->unitigIndex_toReverseDirection(unitigIndex));
			}

			for(u_int32_t predecessor_nodeName : possiblePredecessors){
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(predecessor_nodeName, false);
				u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
				_graph->_isUnitigInvalid_tmp.insert(unitigIndex);
				_graph->_isUnitigInvalid_tmp.insert(_graph->unitigIndex_toReverseDirection(unitigIndex));
			}


		}

	}

	void solveRepeat_disconnectNode(u_int32_t nodeName, bool forward){

		u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, forward);
		u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
		const Unitig& unitig = _graph->_unitigs[unitigIndex];

		vector<u_int32_t> successors;
		vector<u_int32_t> successors2;
		if(nodeIndex == unitig._startNode){
			_graph->getPredecessors(nodeIndex, 0, successors);

			//cout << "Nb succs: " << BiGraph::nodeToString(nodeIndex) << " " << successors.size() << endl;

			for(u_int32_t nodeIndex_succ : successors){
				_graph->_graphSuccessors->removeEdge(nodeIndex_succ, nodeIndex);
				_graph->_graphSuccessors->removeEdge(_graph->nodeIndex_toReverseDirection(nodeIndex), _graph->nodeIndex_toReverseDirection(nodeIndex_succ));
				//cout << "\tRemoving edge 1: " << BiGraph::nodeToString(nodeIndex_succ) << " " << BiGraph::nodeToString(nodeIndex) << endl;
				//cout << "\t" << _graph->_graphSuccessors->removeEdge(nodeIndex_succ, nodeIndex) << endl;
				//cout << "\t"<< _graph->_graphSuccessors->removeEdge(_graph->nodeIndex_toReverseDirection(nodeIndex), _graph->nodeIndex_toReverseDirection(nodeIndex_succ)) << endl;
				//cout << _graph->_graphSuccessors->removeEdge(nodeIndex_succ, nodeIndex) << endl;
			}

			//_graph->getSuccessors(nodeIndex, 0, successors2);
		}
		else{
			_graph->getSuccessors(nodeIndex, 0, successors);
			//_graph->getPredecessors(nodeIndex, 0, successors2);

			//cout << "Nb succs: " << BiGraph::nodeToString(nodeIndex) << " " << successors.size() << endl;
			for(u_int32_t nodeIndex_succ : successors){
				_graph->_graphSuccessors->removeEdge(nodeIndex, nodeIndex_succ);
				_graph->_graphSuccessors->removeEdge(_graph->nodeIndex_toReverseDirection(nodeIndex_succ), _graph->nodeIndex_toReverseDirection(nodeIndex));
				//cout << "\tRemoving edge 2: " << BiGraph::nodeToString(nodeIndex) << " " << BiGraph::nodeToString(nodeIndex_succ) << endl;
				//cout << "\t"<< _graph->_graphSuccessors->removeEdge(nodeIndex, nodeIndex_succ) << endl;
				//cout << "\t"<< _graph->_graphSuccessors->removeEdge(_graph->nodeIndex_toReverseDirection(nodeIndex_succ), _graph->nodeIndex_toReverseDirection(nodeIndex)) << endl;
				//cout << _graph->_graphSuccessors->removeEdge(nodeIndex_succ, nodeIndex) << endl;
			}

		}

		//cout << "\t" << successors.size() << " " << successors2.size() << endl;
		//for(u_int32_t nodeIndex_succ : successors){
		//	cout << "\t\tedge: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_succ) << endl;
		//}
		//for(u_int32_t nodeIndex_succ : successors2){
		//	cout << "\t\tedge: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_succ) << endl;
		//}

		//for(u_int32_t nodeIndex_succ : successors){
		//	cout << "\tRemoving edge: " << BiGraph::nodeToString(nodeIndex) << " " << BiGraph::nodeToString(nodeIndex_succ) << endl;
		//	cout << _graph->_graphSuccessors->removeEdge(nodeIndex, nodeIndex_succ) << endl;
			//cout << _graph->_graphSuccessors->removeEdge(nodeIndex_succ, nodeIndex) << endl;
		//}

	}

	
    unordered_set<u_int32_t> _removedUnitigs;

	void solveRepeat_connectNode(u_int32_t nodeName1, u_int32_t nodeName2){

		u_int32_t nodeIndex_1_noSucc;
		u_int32_t nodeIndex_1_noPred;
		u_int32_t nodeIndex_2_noSucc;
		u_int32_t nodeIndex_2_noPred;

		u_int32_t nodeIndex_1 = BiGraph::nodeName_to_nodeIndex(nodeName1, true);
		u_int32_t nodeIndex_1_rev = BiGraph::nodeName_to_nodeIndex(nodeName1, false);
		u_int32_t nodeIndex_2 = BiGraph::nodeName_to_nodeIndex(nodeName2, true);
		u_int32_t nodeIndex_2_rev = BiGraph::nodeName_to_nodeIndex(nodeName2, false);


		vector<u_int32_t> successors;
		_graph->getSuccessors(nodeIndex_1, 0, successors);
		if(successors.size() == 0){
			nodeIndex_1_noSucc = nodeIndex_1;
			nodeIndex_1_noPred = nodeIndex_1_rev;
		}
		else{
			nodeIndex_1_noSucc = nodeIndex_1_rev;
			nodeIndex_1_noPred = nodeIndex_1;
		}

		_graph->getSuccessors(nodeIndex_2, 0, successors);
		if(successors.size() == 0){
			nodeIndex_2_noSucc = nodeIndex_2;
			nodeIndex_2_noPred = nodeIndex_2_rev;
		}
		else{
			nodeIndex_2_noSucc = nodeIndex_2_rev;
			nodeIndex_2_noPred = nodeIndex_2;
		}

		_graph->_graphSuccessors->addEdge(nodeIndex_1_noSucc, nodeIndex_2_noPred);
		_graph->_graphSuccessors->addEdge(nodeIndex_2_noSucc, nodeIndex_1_noPred);

		//cout << "\t"<< "Add edge: " << BiGraph::nodeToString(nodeIndex_1_noSucc) << " -> " << BiGraph::nodeToString(nodeIndex_2_noPred) << endl;
		//cout << "\t"<< "Add edge: " << BiGraph::nodeToString(nodeIndex_2_noSucc) << " -> " << BiGraph::nodeToString(nodeIndex_1_noPred) << endl;
	
		u_int32_t unitigIndex1 = _graph->nodeIndex_to_unitigIndex(nodeIndex_1_noSucc);
		u_int32_t unitigIndex1_rev = _graph->nodeIndex_toReverseDirection(unitigIndex1);
		u_int32_t unitigIndex2 = _graph->nodeIndex_to_unitigIndex(nodeIndex_2_noSucc);
		u_int32_t unitigIndex2_rev = _graph->nodeIndex_toReverseDirection(unitigIndex2);

		//Unitig& unitig1 = _graph->_unitigs[unitigIndex1];

		
		_removedUnitigs.insert(unitigIndex1);
		_removedUnitigs.insert(unitigIndex1_rev);
		_removedUnitigs.insert(unitigIndex2);
		_removedUnitigs.insert(unitigIndex2_rev);



		//disableUnitig(unitigIndex1);
		//disableUnitig(unitigIndex1_rev);
		//disableUnitig(unitigIndex2);
		//disableUnitig(unitigIndex2_rev);
	
	}

	struct RepeatData{
		u_int32_t _nodeName_successor;
		u_int32_t _nodeName_predecessor;
		u_int32_t _nbSupportingReads;
	};

	static bool RepeatComparator_bySupportedReads(const RepeatData &a, const RepeatData &b){
		return a._nbSupportingReads > b._nbSupportingReads;
	}

	void solveRepeat_collectSides(u_int32_t source_unitigIndex, bool forward, unordered_map<u_int32_t, vector<u_int64_t>>& unitigDatas2, unordered_set<u_int32_t>& possibleSuccessors){

		const vector<u_int64_t>& source_readIndexes = unitigDatas2[source_unitigIndex];

      	unordered_set<u_int32_t> isVisited;
        queue<u_int32_t> queue;

        queue.push(source_unitigIndex);

        while (!queue.empty()){

            u_int64_t unitigIndex = queue.front();
            queue.pop();

            if (isVisited.find(unitigIndex) != isVisited.end()) continue;
            isVisited.insert(unitigIndex);

			
			if(unitigIndex != source_unitigIndex) indexUnitiData2(unitigIndex, unitigDatas2);

            const vector<u_int64_t>& readIndexes = unitigDatas2[unitigIndex];
            //cout << _unitigDatas2[unitigIndex]._readIndexes.size() << endl;

            u_int64_t nbSharedReads = Utils::computeSharedReads(source_readIndexes, readIndexes);
            //cout << nbSharedReads << endl;
            cout << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._endNode) << " " << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode) << " " << nbSharedReads << " " << _graph->_unitigs[unitigIndex]._length << endl;
            //cout << "Shared reads: " << nbSharedReads << " " << readIndexes._readIndexes.size() << " " << source_readIndexes._readIndexes.size() << endl;
            if(nbSharedReads < MIN_SUPPORTED_READS) continue;

			if(unitigIndex != source_unitigIndex){
				if(_graph->_unitigs[unitigIndex]._length >= MIN_REPEAT_SIDE_LENGTH){

					u_int32_t startNodeName = BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode);
					u_int32_t endNodeName = BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._endNode);

            		u_int64_t nbSharedReads_startNode = Utils::computeSharedReads(source_readIndexes, _unitigDatas[startNodeName]._readIndexes);
            		u_int64_t nbSharedReads_endNode = Utils::computeSharedReads(source_readIndexes, _unitigDatas[endNodeName]._readIndexes);

					if(nbSharedReads_startNode > nbSharedReads_endNode){
						possibleSuccessors.insert(startNodeName);
						if(forward){
							file_repeat << startNodeName << "," << "green" << endl;
						}
						else{
							file_repeat << startNodeName << "," << "blue" << endl;
						}
					}
					else if(nbSharedReads_endNode > 0){
						possibleSuccessors.insert(endNodeName);
						if(forward){
							file_repeat << endNodeName << "," << "green" << endl;
						}
						else{
							file_repeat << endNodeName << "," << "blue" << endl;
						}
					}

					cout << "Added repeat side: " << endNodeName << " " << startNodeName << " " << _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(startNodeName, true))._nodes.size() << " " << (_graph->_repeatedNodenames.find(startNodeName) != _graph->_repeatedNodenames.end()) << endl;

					/*
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						if(forward){
							file_repeat << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "green" << endl;
						}
						else{
							file_repeat << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "blue" << endl;
						}
					}
					*/

					continue;
				}
			}



            vector<u_int32_t> successors;
			if(forward){
            	_graph->getSuccessors_unitig(unitigIndex, 0, successors);
			}
			else{
            	_graph->getPredecessors_unitig(unitigIndex, 0, successors);
			}
            
            for(u_int32_t unitigIndex_nn : successors){
                queue.push(unitigIndex_nn);
            }



        }

        //file_scc.close();

	}

	void indexUnitiData2(u_int32_t unitigIndex, unordered_map<u_int32_t, vector<u_int64_t>>& unitigData2){
		
		//if(unitigData2.find(unitigIndex) != unitigData2.end()) return;

		//readIndexes.clear();

        unordered_set<u_int64_t> readIndexes_unique;
        for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
            const UnitigData& unitigData = _unitigDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)];
            //cout << unitigData._readIndexes.size() << endl;
            for(u_int64_t readIndex : unitigData._readIndexes){
                readIndexes_unique.insert(readIndex);
            }
        }

        vector<u_int64_t> readIndexes;
        for(u_int32_t readIndex : readIndexes_unique){
            readIndexes.push_back(readIndex);
        }
        std::sort(readIndexes.begin(), readIndexes.end());

		unitigData2[unitigIndex] = readIndexes;
	}


	void indexReads_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		//vector<ReadIndexType> unitigIndexex;

		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}



	}


	void removeUnsupportedEdges(const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges, GraphSimplify* graph){

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly2::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parse(fp);
		
		vector<DbgEdge> removedEdges;

		for(Unitig& unitig : graph->_unitigs){

			//bool nodeName_ori;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._endNode);

			vector<u_int32_t> successors;
			graph->getSuccessors(unitig._endNode, 0, successors);

			for(u_int32_t successor : successors){
				
				//bool nodeName_succ_ori;
				u_int32_t nodeName_succ = BiGraph::nodeIndex_to_nodeName(successor);

				if(Utils::shareAnyRead(_unitigDatas[nodeName], _unitigDatas[nodeName_succ])){
					continue;
				}

				Unitig& unitig_succ = graph->_unitigs[graph->nodeIndex_to_unitigIndex(successor)];

				//graph->_graphSuccessors->removeEdge(nodeName, nodeName_ori, nodeName_succ, nodeName_succ_ori);
				removedEdges.push_back({unitig._endNode, successor});
				
				
			}



		}


		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;

		for(const DbgEdge& edge: removedEdges){

			bool edgeExist = graph->_graphSuccessors->removeEdge(edge._from, edge._to);
			//if(!edgeExist) cout << "nop 1" << endl;
			edgeExist = graph->_graphSuccessors->removeEdge(GraphSimplify::nodeIndex_toReverseDirection(edge._to), GraphSimplify::nodeIndex_toReverseDirection(edge._from));
			//if(!edgeExist) cout << "nop 2" << endl;
			
			//graph->_isEdgeUnsupported.insert(edge);
		}
		
		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;
		
		cout << "Nb unsupported edges: " << graph->_isEdgeUnsupported.size() << endl;


	}



	size_t _iter;
	ofstream file_groundTruth;
	ofstream file_groundTruth_hifiasmContigs;

	ofstream file_kminmersContigs;
	

	float computeAbundanceCutoff_min(u_int32_t abundance){
		return abundance / 4.0;
	}

	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;

	void extract_truth_kminmers_read(kseq_t* read, u_int64_t readIndex){
		//ottalSize += strlen(read->seq.s);


		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//2262408
				if("utg009434l" == string(read->name.s)){
					cout << nodeName << endl;
					_hifiasm_startingNodenames.insert(nodeName);
				}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(string(read->name.s));
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);

				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(_mdbg->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					_evaluation_hifiasmGroundTruth_nodeNamePosition[_mdbg->_dbg_nodes[vec]._index] = _extract_truth_kminmers_read_position;
				}
				//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
				
				_file_groundTruth_hifiasm_position << nodeName << "," << _extract_truth_kminmers_read_position << endl;
			}
			else{
				//cout << "Not found position: " << position << endl;
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}


			if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
				//cout << position << endl;
				continue;
			}


			//_evaluation_hifiasmGroundTruth[vec] = datasetID;
			_evaluation_hifiasmGroundTruth_position[vec] = _extract_truth_kminmers_read_position;
			_extract_truth_kminmers_read_position += 1;

			
			//cout << _extract_truth_kminmers_read_position << endl;
		}


	}


	void extract_truth_kminmers(){

		_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&Assembly2::extract_truth_kminmers_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_truthInputFilename, true);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		delete _minimizerParser;
	}

	struct UnitigCentrality{
		u_int32_t _unitigIndex;
		float _centrality;
	};

	static bool UnitigComparator_byCentrality(const UnitigCentrality &a, const UnitigCentrality &b){
		return a._centrality < b._centrality;
	}



	u_int32_t computeScgCentrality(const string& outputFilenameBase, unordered_set<u_int32_t>& unitigs, u_int32_t pass, vector<UnitigCentrality>& unitigCentralities){


		unitigCentralities.clear();

		u_int32_t nbPaths = 0;
		unordered_map<u_int32_t, u_int32_t> nbPathsTrought;
		
    	unordered_set<DbgEdge, hash_pair> sourceDests;
    	unordered_set<DbgEdge, hash_pair> processedPaths;

		for(auto& it1 : _scgClusters){
			u_int32_t unitigIndex1 = it1.first;
			if(_graph->_unitigs[unitigIndex1]._startNode == -1) continue;
			if(_graph->_isNodeValid2.find(_graph->_unitigs[unitigIndex1]._startNode) == _graph->_isNodeValid2.end()) continue;
			if(unitigs.find(unitigIndex1) == unitigs.end()) continue;

			const vector<ScgCluster>& scgClusters1 = it1.second;
			for(const ScgCluster& scgCluster1: scgClusters1){
				
				
				for(auto& it2 : _scgClusters){
					u_int32_t unitigIndex2 = it2.first;
					if(_graph->_unitigs[unitigIndex2]._startNode == -1) continue;
					if(_graph->_isNodeValid2.find(_graph->_unitigs[unitigIndex2]._startNode) == _graph->_isNodeValid2.end()) continue;
					if(unitigs.find(unitigIndex2) == unitigs.end()) continue;
					const vector<ScgCluster>& scgClusters2 = it2.second;
					for(const ScgCluster& scgCluster2: scgClusters2){
						
						
						if(scgCluster1._scgIndex == scgCluster2._scgIndex && scgCluster1._scgCluster == scgCluster2._scgCluster) continue;
						if(scgCluster1._scgIndex != scgCluster2._scgIndex) continue;
						if(unitigIndex1 == unitigIndex2) continue;

						//cout << "Decomenter scgCluster1._scgIndex != scgCluster2._scgIndex" << endl;
						DbgEdge edge = {unitigIndex1, unitigIndex2};
						edge = edge.normalize();

						sourceDests.insert(edge);
					}

				}

				
			}

		}

		cout << sourceDests.size() << endl;

		unordered_set<u_int32_t> processedUnitigs;

		while(true){
		
			vector<DbgEdge> toProcess;
			u_int32_t sourceUnitigIndex = -1;

			for(const DbgEdge& sourceDest : sourceDests){

				u_int32_t unitigIndex1 = sourceDest._from;
				u_int32_t unitigIndex2 = sourceDest._to;


				if(sourceUnitigIndex == -1){

					if(processedUnitigs.find(unitigIndex1) == processedUnitigs.end()){
						sourceUnitigIndex = unitigIndex1;
					}
					else if(processedUnitigs.find(unitigIndex2) == processedUnitigs.end()){
						sourceUnitigIndex = unitigIndex2;
					}
				}

				if(sourceUnitigIndex == -1) continue;

				if(sourceUnitigIndex == unitigIndex1 || sourceUnitigIndex == unitigIndex2){
					toProcess.push_back(sourceDest);
				}

			}



			if(toProcess.size() == 0) break;
			processedUnitigs.insert(sourceUnitigIndex);

			//cout << "Source: " << toProcess.size() << " " << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[sourceUnitigIndex]._startNode) << endl;

			unordered_map<u_int32_t, vector<u_int32_t>> destPaths1;
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths2;
			vector<u_int32_t> unitigIndexDests;

			for(const DbgEdge& sourceDest : toProcess){

				if(sourceDest._from != sourceUnitigIndex){

					DbgEdge edge = {sourceDest._from, sourceUnitigIndex};
					edge = edge.normalize();
					if(processedPaths.find(edge)  != processedPaths.end()) continue;
					processedPaths.insert(edge);

					destPaths1[sourceDest._from] = {};
					destPaths2[sourceDest._from] = {};
					unitigIndexDests.push_back(sourceDest._from);

					//cout << "Path: " << sourceUnitigIndex << " " << sourceDest._from << endl;
				}
				else if(sourceDest._to != sourceUnitigIndex){

					DbgEdge edge = {sourceDest._to, sourceUnitigIndex};
					edge = edge.normalize();
					if(processedPaths.find(edge) != processedPaths.end()) continue;
					processedPaths.insert(edge);

					destPaths1[sourceDest._to] = {};
					destPaths2[sourceDest._to] = {};
					unitigIndexDests.push_back(sourceDest._to);
					//cout << "Path: " << sourceUnitigIndex << " " << sourceDest._to << endl;
				}
			}

			//ofstream file_shortest("/home/gats/workspace/run/shortestPath.csv");
			//file_shortest << "Name,Colour" << endl;

			//cout << destPaths1.size() << endl;
			//cout << destPaths2.size() << endl;
			_graph->shortestPath_unitig_weighted(sourceUnitigIndex, destPaths1, true, true, 0);
			//_graph->shortestPath_unitig_weighted(_graph->unitigIndex_toReverseDirection(sourceUnitigIndex), destPaths2, false, false, 1);

			for(u_int32_t unitigIndex : unitigIndexDests){

				for(u_int32_t unitigIndex : destPaths1[unitigIndex]){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
					}
				}
				nbPaths += 1;

				//cout << destPaths1[unitigIndex].size() << " " << destPaths2[unitigIndex].size() << endl;

			
				/*
				u_int32_t path1_nbNodes = 0;
				u_int32_t path2_nbNodes = 0;
				for(u_int32_t unitigIndex : destPaths1[unitigIndex]){
					path1_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;

					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						file_shortest << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "0" << endl;
					}
				}
				for(u_int32_t unitigIndex : destPaths2[unitigIndex]){
					path2_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;

					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						file_shortest << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "1" << endl;
					}
				}

				//file_shortest.close();

				cout << "\tDest: " << unitigIndex << " " <<  path1_nbNodes << " " << path2_nbNodes << "   " << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode)  << endl;
				//if(unitigIndex == 38) getchar();

				if(path1_nbNodes < path2_nbNodes){
					
					for(u_int32_t unitigIndex : destPaths1[unitigIndex]){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}
				else{

					for(u_int32_t unitigIndex : destPaths2[unitigIndex]){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}

				nbPaths += 1;

				*/
			}
		}

		//cout << nbPathsTrought.size() << endl;

		/*
		vector<UnitigScgCluster> allScgClusters;

		for(auto& it : _scgClusters){
			
			u_int32_t unitigIndex = it.first;
			const vector<ScgCluster>& scgClusters = it.second;

			for(const ScgCluster& scgCluster: scgClusters){
				allScgClusters.push_back({unitigIndex, scgCluster._scgIndex, scgCluster._scgCluster});
			}

		}

		u_int32_t nbPaths = 0;
		unordered_map<u_int32_t, u_int32_t> nbPathsTrought;


		for(size_t i=0; i<allScgClusters.size(); i++){

			cout << i << " " << allScgClusters.size() << endl;

			const UnitigScgCluster& scgCluster1 = allScgClusters[i];
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths1;
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths2;

			vector<u_int32_t> unitigIndexDests;

			for(size_t j=0; j<allScgClusters.size(); j++){
			
				const UnitigScgCluster& scgCluster2 = allScgClusters[j];

				if(scgCluster1._scgIndex == scgCluster2._scgIndex && scgCluster1._scgCluster == scgCluster2._scgCluster) continue;
				if(scgCluster1._unitigIndex == scgCluster2._unitigIndex) continue;

				destPaths1[scgCluster2._unitigIndex] = {};
				destPaths1[_graph->unitigIndex_toReverseDirection(scgCluster2._unitigIndex)] = {};
				destPaths2[scgCluster2._unitigIndex] = {};
				destPaths2[_graph->unitigIndex_toReverseDirection(scgCluster2._unitigIndex)] = {};

				unitigIndexDests.push_back(scgCluster2._unitigIndex);
			}

			_graph->shortestPath_unitig_weighted(scgCluster1._unitigIndex, destPaths1, false, false);
			_graph->shortestPath_unitig_weighted(_graph->unitigIndex_toReverseDirection(scgCluster1._unitigIndex), destPaths2, false, false);

			for(u_int32_t unitigIndex : unitigIndexDests){

				vector<u_int32_t> path1;
				const vector<u_int32_t>& path1_forward = destPaths1[unitigIndex];
				const vector<u_int32_t>& path1_reverse = destPaths1[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path1_forward.size() == 0){
					path1 = path1_reverse;
				}
				else{
					path1 = path1_forward;
				}

				vector<u_int32_t> path2;
				const vector<u_int32_t>& path2_forward = destPaths2[unitigIndex];
				const vector<u_int32_t>& path2_reverse = destPaths2[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path2_forward.size() == 0){
					path2 = path2_reverse;
				}
				else{
					path2 = path2_forward;
				}

				u_int32_t path1_nbNodes = 0;
				u_int32_t path2_nbNodes = 0;
				for(u_int32_t unitigIndex : path1){
					path1_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}
				for(u_int32_t unitigIndex : path2){
					path2_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}


				if(path1_nbNodes < path2_nbNodes){
					
					for(u_int32_t unitigIndex : path1){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}
				else{

					for(u_int32_t unitigIndex : path2){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}

				nbPaths += 1;

			}

		}

		*/


		
		float min_value = 2;
		float max_value = -1;

		for(auto& it : nbPathsTrought){
			//u_int32_t componentIndex = unitigComponentIndex[i];
			//u_int32_t nbPaths = componentNbPaths[componentIndex];
			
			//float min_value = 2;
			//float max_value = -1;

			float centrality = it.second / ((float)nbPaths);
			if(centrality < min_value) min_value = centrality;
			if(centrality > max_value) max_value = centrality;
		}

		//cout << min_value << " " << max_value << endl;
		double maxCentrality = 0;
		u_int32_t minLength = -1;
		u_int32_t maxUnitigIndex = -1;

		ofstream output_file(outputFilenameBase + ".centrality.csv");
		output_file << "Name,Colour" << endl;
		
		unordered_set<u_int32_t> writtenUnitigs;

		for(auto& it : nbPathsTrought){
			
			u_int32_t nodeName = it.first;
			u_int32_t nbPathNode = it.second;
			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(nodeName, true));
			
			const Unitig& u = _graph->_unitigs[unitigIndex];
			//u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));


			float centrality = nbPathNode / ((float)nbPaths);
			float centrality_normalized = (centrality-min_value) / (max_value-min_value);

			for(u_int32_t nodeIndex : u._nodes){
				//output_file << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << centrality_normalized << endl;
			}

			if(centrality_normalized >= maxCentrality){

				const Unitig& unitig = _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(nodeName, true));
				//if(keepUnitigs.find(unitig._index) != keepUnitigs.end()) continue;
				if(unitigs.find(unitig._index) == unitigs.end()) continue;

				bool isRepeat = false;
				for(u_int32_t nodeIndex : unitig._nodes){
					if(_unsolvableRepeat.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _unsolvableRepeat.end()){
						isRepeat = true;
						break;
					}
				}
				if(isRepeat) continue;

				if(centrality_normalized > maxCentrality){
					minLength = unitig._length;
					maxUnitigIndex = unitig._index;
				}
				else{
					if(unitig._length < minLength){
						minLength = unitig._length;
						maxUnitigIndex = unitig._index;
					}
				}

				maxCentrality = centrality_normalized;
				cout << "Centrality: " << centrality_normalized << " " << nodeName << " " << maxUnitigIndex << endl;

			}

			unitigCentralities.push_back({unitigIndex, centrality});
		}

		output_file.close();
		
		sort(unitigCentralities.begin(), unitigCentralities.end(), UnitigComparator_byCentrality);

		return maxUnitigIndex;
	}











	void execute_binning(){


		

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;

		u_int32_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;

		unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		for(Unitig& unitig : _graph->_startingUnitigstest){
			
			cout << processedUnitigs << " " << _graph->_startingUnitigstest.size() << endl;
			processedUnitigs += 1;

			bool isProcessed = false;
			for(u_int32_t nodeIndex : unitig._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
					isProcessed = true;
					break;
				}
			}
			if(isProcessed) continue;

			for(u_int32_t nodeIndex : unitig._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				processedNodeNames.insert(nodeName);
			}

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);

			//const Unitig& u = _graph->_unitigs[unitigIndex];
			//u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			//if(processedUnitigs.find(BiGraph::nodeIndex_to_nodeName(unitig._startNode)) != processedUnitigs.end()) continue;
			//if(processedUnitigs.find(BiGraph::nodeIndex_to_nodeName(unitig._endNode)) != processedUnitigs.end()) continue;

			//processedUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig._startNode));
			//processedUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig._endNode));

			cout << endl << "Processing unitig: " << unitig._index << endl;

			//cout << unitig._length << endl;
			//continue;


			float cutoff = unitig._abundance*0.1;


			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff || reloadState){
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
				//_graph->loadState2(0, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
				reloadState = false;
			}
			


        	unordered_set<u_int32_t> component;
        	_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), 20, component);

			if(component.size() <= 2) continue;

			reloadState = true;

			cout << "Component size: " << component.size() << endl;



			//const string& filenameUnitigColor= clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
			//ofstream fUnitigColor(filenameUnitigColor);
			//fUnitigColor << "Name,Colour" << endl;

			//cout << "Generate contigs" << endl;
			//const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
			//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
			//ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

			unordered_set<u_int32_t> writtenUnitigs;


			u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
			const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._endNode));

			string unitigSequence_model;
			_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model);

			vector<float> compositionModel;
			_contigFeature.sequenceToComposition(unitigSequence_model, compositionModel);

			vector<float> abundancesModel;
			_contigFeature.sequenceToAbundance(unitig_model._nodes, abundancesModel);

			ContigFeatures contigFeatureModel = {unitigIndex_model, compositionModel, abundancesModel};

			cout << "\t";
			for(float count : abundancesModel) cout << count << " ";
			cout << endl;


			if(_truthInputFilename != ""){
			
				for(u_int32_t nodeIndex : unitig._nodes){

					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

					if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
						for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
							fileHifiasmAll << unitigName << "," << clusterIndex << endl;
						}
					}
				}
			}


			//vector<u_int32_t> componentNodes;
			//vector<ContigFeatures> contigFeatures;

			//u_int32_t contigIndex = 0;
			for(u_int32_t unitigIndex : component){
				
				const Unitig& u = _graph->_unitigs[unitigIndex];


				bool isProcessed = false;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
						isProcessed = true;
						break;
					}
				}
				if(isProcessed) continue;


				u_int32_t contigIndex = u._index;
				if(u._length < 5000) continue;
				
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				string unitigSequence;
				_toBasespace.createSequence(u._nodes, unitigSequence);

				//string header = ">ctg" + to_string(contigIndex) + '\n';
				//basespaceUnitigFile << header;
				//gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
				//unitigSequence +=  '\n';
				//basespaceUnitigFile << unitigSequence;
				//gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());
				
				//for(u_int32_t nodeIndex : u._nodes){
				//	u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				//	//componentNodes.push_back(nodeName);
				//	fUnitigColor << nodeName << "," << contigIndex << endl;
				//}

				vector<float> composition;
				_contigFeature.sequenceToComposition(unitigSequence, composition);

				vector<float> abundances;
				_contigFeature.sequenceToAbundance(u._nodes, abundances);

				ContigFeatures contigFeature = {unitigIndex, composition, abundances};



				//contigFeatures.push_back({unitigIndex, composition, abundances});

				//cout << unitigSequence.size() << endl;

				if(_contigFeature.isIntra(contigFeatureModel, contigFeature)){

					cout << "\t";
					for(float count : abundances) cout << count << " ";
					cout << endl;

					cout << "\t>" << u._length << ": " << _contigFeature.computeProbability(contigFeatureModel, contigFeature) << endl;
				
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						processedNodeNames.insert(nodeName);
					}

		
					if(_truthInputFilename != ""){
					
						unordered_set<string> hifiasmWrittenUnitigs;

						for(u_int32_t nodeIndex : u._nodes){

							u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

							if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
								for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
									if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
									fileHifiasmAll << unitigName << "," << clusterIndex << endl;
									hifiasmWrittenUnitigs.insert(unitigName);
								}
							}
						}
					}



				}



			}

			//fUnitigColor.close();
			//gzclose(basespaceUnitigFile);
			//basespaceUnitigFile.close();

			cout << "Bin done" << endl;
			getchar();
			clusterIndex += 1;
		}



		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();

	}

};



#endif

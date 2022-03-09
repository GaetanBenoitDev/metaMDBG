
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
//#include "assembly/Assembly.hpp"

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


		//vector<float> lala1 = {0.0704131, 0, 0}; 
		//vector<float> lala2 = {4.79126, 6.7738, 2.99172};
	
		//ContigFeatures f1 = {0, {}, lala1, lala1};
		//ContigFeatures f2 = {0, {}, lala2, lala2};

		//cout << _contigFeature.isIntra(f1, f2) << endl;
		//exit(1);
		/*
		unordered_map<string, vector<u_int32_t>> contigToNodes;
		ifstream file_contigToNode(_inputDir + "/contigs.fasta.gz" + ".nodes.csv");
        vector<string>* fields = new vector<string>();
		string line;
		while (std::getline(file_contigToNode, line)){
			
			GfaParser::tokenize(line, fields, ';');

			u_int32_t nodeName = stoull((*fields)[0]);
			u_int32_t contigIndex = stoull((*fields)[1]);
			//cout << nodeName << " " << contigIndex << endl;
			contigToNodes["ctg" + to_string(contigIndex)].push_back(nodeName);
		}
		delete fields;
		file_contigToNode.close();

		_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodes);
		*/
	

		//_contigFeature.computeFastaComposition("/home/gats/workspace/run/overlap_test_AD_k7/binByCutoff/component_16_unitigs.fasta");
		
		//exit(1);

		//unordered_map<string, vector<u_int32_t>> contigToNodesLala;
		//_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodesLala);

		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		_outputContigFile_complete = gzopen(_outputFilename_complete.c_str(),"wb");

		loadGraph();

		//generateContigs2(_inputDir + "/contigs.min.gz", _inputDir + "/contigs.fasta.gz");
		//return;




		
		//_contigFeature.loadAbundanceFile_nodename(_filename_abundance);
		
		/*
		_graph->loadState2(0, -1, _unitigDatas);
		unordered_map<string, vector<u_int32_t>> contigToNodes;
		for(const Unitig& u : _graph->_unitigs){
			if(u._index%2 == 0){
				for(u_int32_t nodeIndex : u._nodes){
					contigToNodes["ctg" + to_string(u._index)].push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
			}
		}
		_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodes);
		*/



		execute_binning2();
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
		_filename_readMinimizers = _inputDir + "/read_data.txt";
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

		extractContigKminmers2(_inputDir + "/contigs.fasta.gz");

		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}

		//cout << "lala " << _evaluation_hifiasmGroundTruth_position.size() << endl;

		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		//file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		//file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize);
		_graph = graphSimplify;
		
		
		//Generate unitigs
		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);




		/*
		//ofstream file_groundTruth_hifiasm_position(_inputDir + "/groundtruth_hifiasm_position.csv");
		ofstream file_groundTruth_hifiasm(_inputDir + "/groundtruth_hifiasm.csv");
		file_groundTruth_hifiasm << "Name,Colour" << endl;
		//file_groundTruth_hifiasm_position << "Name,Order" << endl;

		unordered_set<u_int32_t> validNodes;

		//unordered_set<u_int32_t> groundTruth_kminmers;
		int founded = 0;
		for(auto it : _mdbg->_dbg_nodes){
			const KmerVec& vec = it.first;

			u_int32_t nodeName = it.second._index;
			
			//vec = vec.normalize();
			if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()) continue;

			founded += 1;
			//groundTruth_kminmers.insert(it.second._index);

			//file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName] << endl;
			//file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

			validNodes.insert(nodeName);

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
		//file_groundTruth_hifiasm_position.close();

		string outputFilename = _inputDir + "/minimizer_graph_truth.gfa";
		GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
		*/




		delete _mdbg;
		cout << "done" << endl;
	









			






		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		//graphSimplify->execute(5, _kminmerSize);
		//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);
		_graph->debug_writeGfaErrorfree(500, 500, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false);




		//_toBasespace.create(_inputDir);


		//_graph->loadState2(0, -1, _unitigDatas);
		//generateFasta(_inputDir + "/allContigs.fasta.gz");
		//generateFastaExpanded(_inputDir + "/allContigsExpanded.fasta.gz");

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
		//file_groundTruth_hifiasmContigs.close();
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


	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

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

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly2::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		if(_kminmerSize > 3) return;
		
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
	//ofstream file_groundTruth_hifiasmContigs;

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
		ReadParser readParser(_truthInputFilename, true, false);
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

		//Assembly assembly(_unitigDatas, _contigFeature);



		vector<vector<string>> bins;
		

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;

		u_int64_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;

		unordered_set<u_int32_t> processedNodeNames;
		unordered_set<u_int32_t> processedContigIndex;
		u_int64_t processedUnitigs = 0;


		ofstream file_bin_all(_inputDir + "/bin_all.csv");
		file_bin_all << "Name,Color" << endl;

		for(Unitig& unitig : _graph->_startingUnitigstest){
			
			cout << processedUnitigs << " " << _graph->_startingUnitigstest.size() << "     " << unitig._length << endl;
			processedUnitigs += 1;

			/*
			bool isProcessed = false;
			for(u_int32_t nodeIndex : unitig._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
					isProcessed = true;
					break;
				}
			}
			if(isProcessed) continue;
			*/
			if(isContigBinned(unitig._nodes, processedNodeNames)) continue;


            //if(unitig._length < 20000) continue; //A TEST

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);

			//const Unitig& u = _graph->_unitigs[unitigIndex];
			//u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			//if(processedUnitigs.find(BiGraph::nodeIndex_to_nodeName(unitig._startNode)) != processedUnitigs.end()) continue;
			//if(processedUnitigs.find(BiGraph::nodeIndex_to_nodeName(unitig._endNode)) != processedUnitigs.end()) continue;

			//processedUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig._startNode));
			//processedUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig._endNode));

			cout << endl << "Starting unitig: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._length << endl;

			//cout << unitig._length << endl;
			//continue;


			float cutoff = unitig._abundance*0.2;


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
			//reloadState = true;
			
			if(_graph->_isNodeValid2.find(unitig._startNode) == _graph->_isNodeValid2.end()) continue;

			u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
			const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];


			for(u_int32_t nodeIndex : unitig_model._nodes){
				processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}

			binByReadpath(unitig._startNode, processedNodeNames, processedContigIndex, clusterDir, filename_binStats, fileHifiasmAll, fileComponentNodeAll, clusterIndex);
			continue;

        	unordered_set<u_int32_t> component;
        	//_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), 50, component);
        	_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), component);


			/*
			unordered_set<u_int32_t> allowedNodeNames;
			for(u_int32_t unitigIndex : component){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					allowedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
			}

			component.clear();
			_graph->loadState2(0, -1, _unitigDatas);
			for(Unitig& unitig : _graph->_startingUnitigstest){
				if(allowedNodeNames.find(BiGraph::nodeIndex_to_nodeName(unitig._startNode)) == allowedNodeNames.end()) continue;
				component.insert(unitig._index);
				component.insert(_graph->unitigIndex_toReverseDirection(unitig._index));
			}
			*/
			//if(component.size() <= 2) continue;


			cout << "Component size: " << component.size() << endl;

			/*
			cout << "Compute reachable unitigs" << endl;
			unordered_set<u_int32_t> reachedUnitigs;
			PathExplorerManager pathExplorerManager(_graph, _unitigDatas, component, reachedUnitigs);
			pathExplorerManager.collectReachableNodes(unitig._startNode);

			component = reachedUnitigs;
			cout << "Component size (reachable): " << component.size() << endl; 
			getchar();
			*/
		
			//const string& filenameUnitigColor= clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
			//ofstream fUnitigColor(filenameUnitigColor);
			//fUnitigColor << "Name,Colour" << endl;

			//cout << "Generate contigs" << endl;
			//const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
			//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
			//ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

			unordered_set<u_int32_t> writtenUnitigs;

	
			//u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
			//const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._endNode));

			for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex_model]._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				processedNodeNames.insert(nodeName);
			}

			//vector<u_int32_t> nodePath_model;
			//assembly.solveBin2(unitig._startNode, 0, _graph, 0, 0, false, nodePath_model);
			//cout << nodePath_model.size() << endl;


			//string unitigSequence_model;
			//_toBasespace.createSequence(nodePath_model, unitigSequence_model);
			string unitigSequence_model_init;
			_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model_init);

			vector<float> compositionModel_init;
			_contigFeature.nodepathToComposition(unitig_model._nodes, compositionModel_init);
			//_contigFeature.sequenceToComposition(unitigSequence_model_init, compositionModel_init);
			//vector<float> compositionModel;
			//_contigFeature.sequenceToComposition(unitigSequence_model, compositionModel);

			vector<float> abundancesModel;
			vector<float> abundancesModel_var;
			//_contigFeature.sequenceToAbundance(nodePath_model_init, abundancesModel);
			bool isAbValid = _contigFeature.sequenceToAbundance(unitig_model._nodes, abundancesModel, abundancesModel_var);
			if(!isAbValid) continue;
			
			ContigFeatures contigFeatureModel = {unitigIndex_model, compositionModel_init, abundancesModel, abundancesModel_var};

			//cout << unitigSequence_model.size() << endl;

			vector<string> bin;
			bin.push_back(unitigSequence_model_init);






			/*
			ofstream file_coverages(_inputDir + "/coverages.csv");
        	file_coverages << "Name,Abundance" << endl;

			ofstream file_coverages_component(_inputDir + "/coverages_component.csv");
        	file_coverages_component << "Name,Color" << endl;

			_graph->_isNodeInvalid_tmp.clear();

			unordered_set<u_int32_t> unitigIndexCovLala;
			for(u_int32_t unitigIndex : component){
				
				const Unitig& u = _graph->_unitigs[unitigIndex];
				
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != unitigIndexCovLala.end()) continue;
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != unitigIndexCovLala.end()) continue;

				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));




				float abRatio = 1;
				vector<float> abundances;
				vector<float> abundances_var;
				bool isAbValid = _contigFeature.sequenceToAbundance(u._nodes, abundances, abundances_var);
				if(isAbValid){


					vector<float> means_f1 = abundancesModel;
					vector<float> means_f2 = abundances;

					float mean_ratio = 0;
					for(size_t i=0; i<means_f1.size(); i++){
						if(means_f2[i] > 0){
							float ratio = means_f1[i] / means_f2[i];
							mean_ratio += ratio;
							//cout << ratio << " ";
						}
						else{
							//cout << "0" << " ";
						}
					}
					//cout << endl;

					mean_ratio /= means_f1.size();

					for(size_t i=0; i<means_f1.size(); i++){
						means_f2[i] *= mean_ratio;
					}


					bool hasValidAb = true;

					for(size_t i=0; i<means_f2.size(); i++){
						if(abundancesModel[i] <= 2) continue;
						float abOffset = abundancesModel[i] * abRatio;
						if(means_f2[i] < abundancesModel[i]-abOffset || means_f2[i] > abundancesModel[i]+abOffset){
							hasValidAb = false;
							break;
						}
					}

					if(hasValidAb){
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",green" << endl;
						}
					}
					else{
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",red" << endl;
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(_graph->nodeIndex_toReverseDirection(nodeIndex)));
						}

					}

				}


				unordered_set<u_int32_t> contigCovs;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

					if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;
					
					u_int32_t contigIndex = _contigFeature._nodeName_to_contigIndex[nodeName];

					if(contigCovs.find(contigIndex) != contigCovs.end()) continue;

					const vector<float>& contigAbundances = _contigFeature._contigCoverages[contigIndex];

					contigCovs.insert(contigIndex);

					string str = "";
					for(u_int32_t ab : contigAbundances){
						str += to_string(ab) + ";";
					}

					file_coverages << nodeName << "," << str << endl;
				}

				//vector<float> abundances;
				//vector<float> abundancesVar;
				//_contigFeature.sequenceToAbundance(u._nodes, abundances, abundancesVar);



			}

			file_coverages.close();
			file_coverages_component.close();
			
			component.clear();
        	//_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), 50, component);
        	_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(unitig._startNode), component);
			cout << "Component size filtered: " << component.size() << endl;

			cout << "\t";
			for(u_int32_t count : abundancesModel) cout << count << " ";
			cout << endl;
			getchar();
			*/

			/*
            unordered_set<u_int32_t> validNodes;
            for (u_int32_t unitigIndex : component){
                for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
                    u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
                    //if(_debug_groundTruthNodeNames.find(nodeName) == _debug_groundTruthNodeNames.end()) continue;
                    validNodes.insert(nodeName);
                }
            }
            string outputFilename = _inputDir + "/minimizer_graph_sub_2.gfa";
            GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			*/

			/*
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
			*/

			vector<string> hifiasmUnitigNames;
			//vector<u_int32_t> componentNodes;
			//vector<ContigFeatures> contigFeatures;

			//u_int32_t contigIndex = 0;

			ofstream file_bin(_inputDir + "/bin.csv");
        	file_bin << "Name,Color" << endl;

			for(u_int32_t nodeIndex : unitig_model._nodes){
				file_bin << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "green" << endl;
				file_bin_all << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "green" << endl;
			}

			//unordered_set<u_int32_t> binnedUnitigs;
			//for(u_int32_t nodeIndex : unitig._nodes){
			//	binnedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			//}
			
			for(u_int32_t unitigIndex : component){
				
				const Unitig& u = _graph->_unitigs[unitigIndex];

				/*
				bool isProcessed = false;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
						isProcessed = true;
						break;
					}
				}
				if(isProcessed) continue;
				*/
				if(isContigBinned(u._nodes, processedNodeNames)) continue;
				

				u_int32_t contigIndex = u._index;
				if(u._length < 2500){
					for(u_int32_t nodeIndex : u._nodes){
						file_repeat << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "blue" << endl;
					}
					continue;
				}
				
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				vector<u_int32_t> nodePath;
				nodePath = u._nodes;
				//vector<u_int32_t> nodePath;
				vector<u_int64_t> nodePath_supportingReads;
				//assembly.solveBin2(u._startNode, 0, _graph, 0, 0, false, nodePath, nodePath_supportingReads, 100000);

				string unitigSequence;
				_toBasespace.createSequence(nodePath, unitigSequence);

				cout << endl << "\tUnitig: " << BiGraph::nodeIndex_to_nodeName(u._startNode) << " " << u._length << " " << nodePath.size() << " " << unitigSequence.size() << endl;
				
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
				bool hasComposition = _contigFeature.nodepathToComposition(nodePath, composition);
				//_contigFeature.sequenceToComposition(unitigSequence, composition);
				
				vector<float> abundances;
				vector<float> abundancesVar;
				bool hasAbundances = _contigFeature.sequenceToAbundance(nodePath, abundances, abundancesVar);

				ContigFeatures contigFeature = {unitigIndex, composition, abundances, abundancesVar};
				
				
				//string unitigSequence_init;
				//_toBasespace.createSequence(u._nodes, unitigSequence_init);
				//vector<float> composition_init;
				//_contigFeature.sequenceToComposition(unitigSequence_init, composition_init);
				//vector<float> abundances_init;
				//vector<float> abundancesVar_init;
				//_contigFeature.sequenceToAbundance(u._nodes, abundances_init, abundancesVar_init);

				//ContigFeatures contigFeature_init = {unitigIndex, composition_init, abundances_init, abundancesVar_init};

				if(_contigFeature.isIntra(contigFeatureModel, contigFeature, hasComposition, hasAbundances)){
					//cout << "\tComposition dist: " << -log10(_contigFeature.computeCompositionProbability(compositionModel_init, composition)) << " " << _contigFeature.cal_tnf_dist(compositionModel_init, composition) << endl; // << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances_init) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature_init)  << endl; 
					//cout << "\tComposition dist (extended): " << _contigFeature.computeEuclideanDistance(compositionModel_init, composition) << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature)  << endl; 
					int nnz = 0;
					cout << "\tMetabat Abudance: " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;// << " " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetabat Abudance new: " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl; // << " " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetacoag abundance: " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl; // << " " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl;
					cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl; // << " " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl;
					//cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances_init) << " " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances) << endl;
					cout << "\t";
					for(u_int32_t count : abundancesModel) cout << count << " ";
					cout << endl;
					cout << "\t";
					for(u_int32_t count : abundances) cout << count << " ";
					cout << endl;
					//cout << "\t";
					//for(float count : abundancesModel_var) cout << count << " ";
					//cout << endl;
					//cout << "\t";
					//for(float count : abundancesVar) cout << count << " ";
					//cout << endl;
				}

				//contigFeatures.push_back({unitigIndex, composition, abundances});

				//cout << unitigSequence.size() << endl;

				if(_contigFeature.isIntra(contigFeatureModel, contigFeature, hasComposition, hasAbundances)){

					cout << "\t>" << u._length << ": " << _contigFeature.computeProbability(contigFeatureModel, contigFeature) << endl;
				
					//binnedUnitigs.insert(u._index);


					string unitigSequence_init;
					_toBasespace.createSequence(u._nodes, unitigSequence_init);
					bin.push_back(unitigSequence_init);



					for(u_int32_t nodeIndex : u._nodes){
						processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
						file_bin << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
						file_bin_all << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;

						if(_truthInputFilename != ""){
						
							unordered_set<string> hifiasmWrittenUnitigs;
							u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

							if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
								for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
									if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
									//fileHifiasmAll << unitigName << "," << clusterIndex << endl;
									hifiasmWrittenUnitigs.insert(unitigName);
									hifiasmUnitigNames.push_back(unitigName);
								}
							}
						}
					}

					//for(u_int32_t nodeIndex : u._nodes){
						//file_repeat << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
						//binnedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					//}

					/*
					bin.push_back(unitigSequence_init);



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
									//fileHifiasmAll << unitigName << "," << clusterIndex << endl;
									hifiasmWrittenUnitigs.insert(unitigName);

									hifiasmUnitigNames.push_back(unitigName);
								}
							}
						}
					}
					*/



				}



			}

			/*
			writtenUnitigs.clear();

			//_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
			for(u_int32_t nodeNameBinned : binnedNodeNames){

				const Unitig& u = _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(nodeNameBinned, true));

				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				string unitigSequence_init;
				_toBasespace.createSequence(u._nodes, unitigSequence_init);
				bin.push_back(unitigSequence_init);


				for(u_int32_t nodeIndex : u._nodes){
					processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));

					if(_truthInputFilename != ""){
					
						unordered_set<string> hifiasmWrittenUnitigs;
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

						if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
							for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
								if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
								//fileHifiasmAll << unitigName << "," << clusterIndex << endl;
								hifiasmWrittenUnitigs.insert(unitigName);

								hifiasmUnitigNames.push_back(unitigName);
							}
						}
					}
				}

			}
			//fUnitigColor.close();
			//gzclose(basespaceUnitigFile);
			//basespaceUnitigFile.close();
			*/
			
			int ret = computeBinStats(clusterDir, bin, filename_binStats);

			if(ret == 0){

				string lastLine = "";
				string line = "";
				ifstream file_binStats(filename_binStats);
				while (std::getline(file_binStats, line)){
					if(line.empty()) continue;
					lastLine = line;
				}
				file_binStats.close();
				
				vector<string>* fields = new vector<string>();
				GfaParser::tokenize(lastLine, fields, ',');
				float completeness = stof((*fields)[0]);
				float contamination = stof((*fields)[1]);
				delete fields;
				cout << "Result: " << completeness << " " << contamination << endl;

				if(completeness > 0.4){
					for(const string& unitigName : hifiasmUnitigNames){
						fileHifiasmAll << unitigName << "," << clusterIndex << endl;
					}
				}

				if(contamination > 0.2){
					getchar();
				}

			}


			//bins.push_back(bin);
			cout << "Bin done" << endl;
			file_bin.close();
			//getchar();
			clusterIndex += 1;
		}



		file_groundTruth.close();
		//file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();
		file_binStats.close();
		file_bin_all.close();

	}

	int computeBinStats(const string& outputDir, const vector<string>& bin, const string& filename_binStats){

		int lala = 0;
		const string& basespaceUnitigFilename = outputDir + "/contigsTmp.fasta";
		ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

		for(const string& sequence : bin){



			string header = ">ctg" + to_string(lala);
			basespaceUnitigFile << header << endl;
			basespaceUnitigFile << sequence << endl;
			
			lala += 1;
		}

		basespaceUnitigFile.close();

		string command_annotation = "python3 /home/gats/workspace/scripts/annotation/annotation2/computeBinStats.py -t 4 " + basespaceUnitigFilename + " " + filename_binStats;
		cout << command_annotation << endl;
		int ret = system(command_annotation.c_str());
		if(ret != 0){
			cerr << "Command failed: " << ret << endl;
			//exit(ret);
		}

		return ret;
		
	}

	void binByReadpath(u_int32_t source_nodeIndex, unordered_set<u_int32_t>& processedNodeNames, unordered_set<u_int32_t>& processedContigIndex, const string& clusterDir, const string& filename_binStats, ofstream& fileHifiasmAll, ofstream& fileComponentNodeAll, u_int64_t& clusterIndex){


		u_int32_t source_unitigIndex = _graph->nodeIndex_to_unitigIndex(source_nodeIndex);
		const Unitig& unitig_model = _graph->_unitigs[source_unitigIndex];

		u_int32_t contigIndex_model;
		string unitigSequence_model;
		//_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model);
		_contigFeature.nodepathToContigSequence(unitig_model._nodes, unitigSequence_model, contigIndex_model);
		if(unitigSequence_model.empty()) return;				
		if(processedContigIndex.find(contigIndex_model) != processedContigIndex.end()) return;




		//unordered_map<u_int32_t, u_int32_t> lala;
		unordered_set<u_int32_t> componentContigIndex;

		vector<string> bin;

		unordered_set<u_int32_t> binnedContigIndex;

      	unordered_set<u_int32_t> isVisited;
        queue<u_int32_t> queue;

        queue.push(source_unitigIndex);





		bin.push_back(unitigSequence_model);
		processedContigIndex.insert(contigIndex_model);
		cout << "\tContig index model: " << contigIndex_model << endl;

		vector<float> compositionModel;
		_contigFeature.nodepathToComposition(unitig_model._nodes, compositionModel);

		vector<float> abundancesModel;
		vector<float> abundancesModel_var;
		bool isAbValid = _contigFeature.sequenceToAbundance(unitig_model._nodes, abundancesModel, abundancesModel_var);
		if(!isAbValid) return;
			
		componentContigIndex.insert(contigIndex_model);
		ContigFeatures contigFeatureModel = {contigIndex_model, compositionModel, abundancesModel, abundancesModel_var};

		unordered_set<u_int32_t> nonIntraUnitigs;

		unordered_set<u_int32_t> writtenUnitigs;
		//u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
		//const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];
		writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._startNode));
		writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._endNode));

		unordered_set<string> hifiasmUnitigNames;
		unordered_set<u_int32_t> allComponentNodenames;




		for(u_int32_t nodeIndex : unitig_model._nodes){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			allComponentNodenames.insert(nodeName);
		}

		if(_truthInputFilename != ""){
			for(u_int32_t nodeIndex : unitig_model._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

				if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
					for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
						hifiasmUnitigNames.insert(unitigName);
					}
				}
			}
		}



		while(!queue.empty()){

            u_int64_t unitigIndex = queue.front();
            queue.pop();

            if (isVisited.find(unitigIndex) != isVisited.end()) continue;

            isVisited.insert(unitigIndex);
            isVisited.insert(_graph->unitigIndex_toReverseDirection(unitigIndex));




			size_t componentNbNodes = 10;
			size_t componentLength = 20000;

			unordered_set<u_int32_t> component;
			unordered_set<u_int32_t> componentNbNodes2;
			//_graph->getConnectedComponent_readpath(unitigIndex, _unitigDatas, 1, component);
			_graph->getConnectedComponent_unitig_nt(unitigIndex, componentLength, component);
			_graph->getConnectedComponent_unitig(unitigIndex, componentNbNodes, componentNbNodes2);
			//cout << "Component size: " << componentNbNodes2.size() << " " << component.size() << endl;


			//component = componentNbNodes2;





			/*
			ofstream file_coverages(_inputDir + "/coverages.csv");
        	file_coverages << "Name,Abundance" << endl;

			ofstream file_coverages_component(_inputDir + "/coverages_component.csv");
        	file_coverages_component << "Name,Color" << endl;

			_graph->_isNodeInvalid_tmp.clear();

			
			unordered_set<u_int32_t> unitigIndexCovLala;
			for(u_int32_t unitigIndex : component){
				
				const Unitig& u = _graph->_unitigs[unitigIndex];
				
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != unitigIndexCovLala.end()) continue;
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != unitigIndexCovLala.end()) continue;

				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));




				float abRatio = 1;
				vector<float> abundances;
				vector<float> abundances_var;
				bool isAbValid = _contigFeature.sequenceToAbundance(u._nodes, abundances, abundances_var);
				if(isAbValid){


					vector<float> means_f1 = abundancesModel;
					vector<float> means_f2 = abundances;

					float mean_ratio = 0;
					for(size_t i=0; i<means_f1.size(); i++){
						if(means_f2[i] > 0){
							float ratio = means_f1[i] / means_f2[i];
							mean_ratio += ratio;
							//cout << ratio << " ";
						}
						else{
							//cout << "0" << " ";
						}
					}
					//cout << endl;

					mean_ratio /= means_f1.size();

					for(size_t i=0; i<means_f1.size(); i++){
						means_f2[i] *= mean_ratio;
					}


					bool hasValidAb = true;

					for(size_t i=0; i<means_f2.size(); i++){
						if(abundancesModel[i] <= 2) continue;
						float abOffset = abundancesModel[i] * abRatio;
						if(means_f2[i] < abundancesModel[i]-abOffset || means_f2[i] > abundancesModel[i]+abOffset){
							hasValidAb = false;
							break;
						}
					}

					if(hasValidAb){
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",green" << endl;
						}
					}
					else{
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",red" << endl;
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(_graph->nodeIndex_toReverseDirection(nodeIndex)));
						}

					}

				}


				unordered_set<u_int32_t> contigCovs;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

					if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;
					
					u_int32_t contigIndex = _contigFeature._nodeName_to_contigIndex[nodeName];

					if(contigCovs.find(contigIndex) != contigCovs.end()) continue;

					const vector<float>& contigAbundances = _contigFeature._contigCoverages[contigIndex];

					contigCovs.insert(contigIndex);

					string str = "";
					for(u_int32_t ab : contigAbundances){
						str += to_string(ab) + ";";
					}

					file_coverages << nodeName << "," << str << endl;
				}

				//vector<float> abundances;
				//vector<float> abundancesVar;
				//_contigFeature.sequenceToAbundance(u._nodes, abundances, abundancesVar);



			}

			file_coverages.close();
			file_coverages_component.close();
			
			component.clear();
        	_graph->getConnectedComponent_unitig_nt(unitigIndex, componentLength, component);
			*/
		
			//cout << "Component size filtered: " << component.size() << endl;

			
			//cout << "\t";
			//for(u_int32_t count : abundancesModel) cout << count << " ";
			//cout << endl;
			//getchar();
			

			//for(u_int32_t utgIndex : component) queue.push(utgIndex);




			//"Reprise: better shortest path from set of nodes, si on entre dans un noeud qui a une distance plus faible que celle actuelle, pas besoin de le visiter"
			//"Extraire les bins (contigs) et les evaluer avec checkm"
			//"evaluer metaflye sur Human data"
			//if(unitigIndex != source_unitigIndex && unitigIndex != _graph->unitigIndex_toReverseDirection(source_unitigIndex)){
			vector<u_int32_t> validUnitigs;
			
			for(u_int32_t unitigIndex : component){
				const Unitig& u = _graph->_unitigs[unitigIndex];
				
				if(nonIntraUnitigs.find(unitigIndex) != nonIntraUnitigs.end()) continue;

				/*
				bool isProcessed = false;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
						isProcessed = true;
						break;
					}
				}

				if(isProcessed) continue;
				*/

				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				//if(isContigBinned(u._nodes, processedNodeNames)) continue;
				if(u._length < 2500) continue;
				//if(u._abundance < _minUnitigAbundance) continue;
				//if(u._nbNodes <= _kminmerSize*2) continue;
				
				

				


				string unitigSequence;
				u_int32_t contigIndex;
				//_toBasespace.createSequence(u._nodes, unitigSequence);
				_contigFeature.nodepathToContigSequence(u._nodes, unitigSequence, contigIndex);
				if(unitigSequence.empty()) continue;

				if(componentContigIndex.find(contigIndex) != componentContigIndex.end()){
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						allComponentNodenames.insert(nodeName);

						if(!_truthInputFilename.empty()){
							if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
								for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
									hifiasmUnitigNames.insert(unitigName);
								}
							}
						}

					}
				}

				if(processedContigIndex.find(contigIndex) != processedContigIndex.end()){
					validUnitigs.push_back(unitigIndex);
					continue;
				}

				cout << endl << "\tUnitig: " << BiGraph::nodeIndex_to_nodeName(u._startNode) << " " << u._length << " " << u._nodes.size() << " " << unitigSequence.size() << endl;
				
				vector<float> composition;
				bool hasComposition = _contigFeature.nodepathToComposition(u._nodes, composition);
				//_contigFeature.sequenceToComposition(unitigSequence, composition);
				
				vector<float> abundances;
				vector<float> abundancesVar;
				bool hasAbundances = _contigFeature.sequenceToAbundance(u._nodes, abundances, abundancesVar);

				ContigFeatures contigFeature = {contigIndex, composition, abundances, abundancesVar};






				if(_contigFeature.isIntra(contigFeatureModel, contigFeature, hasComposition, hasAbundances)){

					cout << "\tComposition dist: " << -log10(_contigFeature.computeCompositionProbability(compositionModel, composition)) << " " << _contigFeature.cal_tnf_dist(compositionModel, composition, contigIndex_model, contigIndex) << endl; // << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances_init) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature_init)  << endl; 
					//cout << "\tComposition dist (extended): " << _contigFeature.computeEuclideanDistance(compositionModel_init, composition) << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature)  << endl; 
					int nnz = 0;
					cout << "\tMetabat Abudance: " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;// << " " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetabat Abudance new: " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl; // << " " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetacoag abundance: " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl; // << " " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl;
					cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl; // << " " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl;
					//cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances_init) << " " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances) << endl;
					cout << "\t";
					for(u_int32_t count : abundancesModel) cout << count << " ";
					cout << endl;
					cout << "\t";
					for(u_int32_t count : abundances) cout << count << " ";
					cout << endl;
					cout << "\t>" << u._length << ": " << contigIndex << endl;
				
					//lala[BiGraph::nodeIndex_to_nodeName(u._startNode)] += 1;
					//lala[BiGraph::nodeIndex_to_nodeName(u._endNode)] += 1;

					//if(lala[BiGraph::nodeIndex_to_nodeName(u._startNode)] > 1 || lala[BiGraph::nodeIndex_to_nodeName(u._endNode)] > 1){
					//	cout << "nani??" << endl;
					//	getchar();
					//}


					//binnedUnitigs.insert(u._index);

					processedContigIndex.insert(contigIndex);
					componentContigIndex.insert(contigIndex);

					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						processedNodeNames.insert(nodeName);
						allComponentNodenames.insert(nodeName);
					}

					if(_truthInputFilename != ""){
					
						for(u_int32_t nodeIndex : u._nodes){

							u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

							
							if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
								for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
									hifiasmUnitigNames.insert(unitigName);
								}
							}
						}
					}
					
					//string unitigSequence_init;
					//_toBasespace.createSequence(u._nodes, unitigSequence_init);
					bin.push_back(unitigSequence);

					validUnitigs.push_back(unitigIndex);
					//queue.push(unitigIndex);
				}
				else{
					nonIntraUnitigs.insert(u._index);
					nonIntraUnitigs.insert(_graph->unitigIndex_toReverseDirection(u._index));
				}
				
			}

			for(u_int32_t unitigIndex : validUnitigs){
				queue.push(unitigIndex);
			}


		}


		int ret = computeBinStats(clusterDir, bin, filename_binStats);

		if(ret == 0){

			string lastLine = "";
			string line = "";
			ifstream file_binStats(filename_binStats);
			while (std::getline(file_binStats, line)){
				if(line.empty()) continue;
				lastLine = line;
			}
			file_binStats.close();
			
			vector<string>* fields = new vector<string>();
			GfaParser::tokenize(lastLine, fields, ',');
			float completeness = stof((*fields)[0]);
			float contamination = stof((*fields)[1]);
			delete fields;
			cout << "Result: " << completeness << " " << contamination << endl;

			if(completeness > 0.3){

				for(u_int32_t nodeName : allComponentNodenames){
					fileComponentNodeAll << nodeName << "," << clusterIndex << endl;
				}
				fileComponentNodeAll.flush();

				for(const string& unitigName : hifiasmUnitigNames){
					fileHifiasmAll << unitigName << "," << clusterIndex << endl;
				}
				fileHifiasmAll.flush();

				/*
				unordered_set<u_int32_t> component;
            	_graph->getConnectedComponent(source_nodeIndex, component);
				unordered_set<u_int32_t> validNodes;
				for (u_int32_t unitigIndex : component){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						validNodes.insert(nodeName);
					}
				}
				string outputFilename = _inputDir + "/minimizer_graph_sub_4_" + to_string(clusterIndex) + ".gfa";
				GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
				*/
			


				//if(clusterIndex >= 4) getchar();
				clusterIndex += 1;
			}

			if(contamination > 0.2){
				getchar();
			}

		}

		cout << "bin done" << endl;
		//getchar();
	}

	void generateFasta(const string& outputFilename){

		gzFile basespaceUnitigFile = gzopen(outputFilename.c_str(), "wb");

		unordered_set<u_int32_t> writtenUnitigs;

		for(const Unitig& u : _graph->_unitigs){
			
			//u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			
			vector<u_int32_t> nodePath = u._nodes;
			u_int32_t nodeIndex = nodePath[0];
			for(int i=0; i<(_kminmerSize-1); i++){

				vector<u_int32_t> preds;
				_graph->getPredecessors(nodeIndex, 0, preds);

				if(preds.size() == 0) break;
				nodeIndex = preds[0];

				nodePath.insert(nodePath.begin(), nodeIndex);
			}
			

			string unitigSequence;
			_toBasespace.createSequence(nodePath, unitigSequence);

			string header = ">ctg" + to_string(u._index) + '\n';
			gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
			unitigSequence +=  '\n';
			gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());

		}

		gzclose(basespaceUnitigFile);
		
	}


	struct Contig{
		vector<u_int32_t> _nodePath;
		vector<u_int32_t> _nodePath_sorted;
	};

	void generateContigs(const string& outputFilename, const string& outputFilename_fasta){

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
		if(!fs::exists (path)){
			fs::create_directory(path);
		} 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(_inputDir + "/binning_results_hifiasm.csv");
		fileHifiasmAll << "Name,Colour" << endl;

		gzFile outputContigFile = gzopen(outputFilename.c_str(),"wb");
		gzFile outputContigFile_fasta = gzopen(outputFilename_fasta.c_str(),"wb");

		ofstream file_asmResult = ofstream(_inputDir + "/binning_results.csv");
		file_asmResult << "Name,Colour" << endl;

		//Assembly assembly(_unitigDatas, _contigFeature);

		float prevCutoff = -1;


		unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		u_int64_t contigIndex = 0;
		u_int64_t __loadState2_index = 0;

		vector<Contig> contigs;
		vector<vector<u_int32_t>> clusters;

		for(Unitig& u : _graph->_startingUnitigstest){
			
			//cout << processedNodeNames.size() << " " << _graph->_startingUnitigstes.size() << endl;

			//cout << processedUnitigs << " " << _graph->_startingUnitigstest.size() << endl;
			processedUnitigs += 1;

			bool isProcessed = false;
			for(u_int32_t nodeIndex : u._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
					isProcessed = true;
					break;
				}
			}
			if(isProcessed) continue;

			
			float cutoff = u._abundance*0.2;

			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			//string unitigSequenceStarting;
			//_toBasespace.createSequence(u._nodes, unitigSequenceStarting);
			cout << "------- " << BiGraph::nodeIndex_to_nodeName(u._startNode) << " " << u._length << endl; // << " " << unitigSequenceStarting.size() << endl;

			if(realCutoff != prevCutoff){
				_graph->loadState2(cutoff, u._startNode, _unitigDatas);
				prevCutoff = realCutoff;
			}
			//continue;

            unordered_set<u_int32_t> component;

			/*
			_graph->determineRepeats(_unitigDatas, cutoff);

			cout <<  "repeat" << endl;
			getchar();
			*/

			/*
            unordered_set<u_int32_t> component;
            _graph->getConnectedComponent(u._startNode, component);

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
				u_int64_t nbShared = Utils::computeSharedElements(existingCluster, cluster);
				//float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				//cout << "Existing: " << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				float similarity = nbShared / (float)cluster.size();
				cout << "Similarity: " << similarity << endl;
				if(similarity > 0.5){
					isNewCluster = false;
					break;
				}
			}

			if(!isNewCluster) continue;

			clusters.push_back(cluster);

            unordered_set<u_int64_t> readset;
            for (u_int32_t unitigIndex : component){
                for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
                    for(u_int64_t readIndex : _unitigDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]._readIndexes){
                        readset.insert(readIndex);
                    }
                }
            }

            const string& readFilename = "/home/gats/workspace/data/HiFi_AD/input.txt";
            //const string& readFilename = "/home/gats/workspace/data/overlap_test/genome_both/input.txt";
            ReadParser readParser(readFilename, false);
            //readParser.extractSubsample(_inputDir + "/readset_" + to_string(__loadState2_index) + ".fasta.gz", readset);
			__loadState2_index += 1;
			*/






			vector<u_int32_t> nodePath;

			
			//PathExplorerManager pathExplorerManager(_graph, _unitigDatas, component, _contigFeature, file_asmResult);
			//pathExplorerManager.execute(u._startNode, nodePath);
			//exit(1);

			vector<u_int64_t> nodePath_supportingReads;
			//assembly.solveBin2(u._startNode, u._abundance, _graph, 0, 0, false, nodePath, nodePath_supportingReads, 0);


			for(u_int32_t nodeIndex : nodePath){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				processedNodeNames.insert(nodeName);
			}


			string unitigSequenceExpanded;
			_toBasespace.createSequence(nodePath, unitigSequenceExpanded);
			cout << "\tExpanded contigs: " << nodePath.size() << " " << unitigSequenceExpanded.size() << endl;

			if(unitigSequenceExpanded.size() < 100000){
				
            	_graph->getConnectedComponent(u._startNode, component);

				unordered_set<u_int32_t> validNodes;
				for (u_int32_t unitigIndex : component){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						//if(_debug_groundTruthNodeNames.find(nodeName) == _debug_groundTruthNodeNames.end()) continue;
						validNodes.insert(nodeName);
					}
				}
				string outputFilename = _inputDir + "/minimizer_graph_sub_3.gfa";
				GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);

				getchar();
			}

			//exit(1);

			//for(u_int32_t nodeIndex : nodePath){
			//	file_asmResult << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << contigIndex << endl;
			//}
			//file_asmResult.flush();
			//getchar();

			//contigIndex += 1;
			//continue;



			//bool isValid = true;
			/*
			bool isValid = dereplicateContig(contigs, nodePath);

			nodePath = contig._nodePath;
			for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
				nodePath.erase(nodePath.begin());
			}
			for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
				nodePath.pop_back();
			}

			if(isValid && nodePath.size() > _kminmerSize*2){

				//contigs.push_back(contig);
			*/

				/*
				string unitigSequenceLala;
				_toBasespace.createSequence(contig._nodePath, unitigSequenceLala);

				vector<string> bin = {unitigSequenceLala};
				int ret = computeBinStats(clusterDir, bin, filename_binStats);


				if(ret == 0){

					string lastLine = "";
					string line = "";
					ifstream file_binStats(filename_binStats);
					while (std::getline(file_binStats, line)){
						if(line.empty()) continue;
						lastLine = line;
					}
					file_binStats.close();
					
					vector<string>* fields = new vector<string>();
					GfaParser::tokenize(lastLine, fields, ',');
					float completeness = stof((*fields)[0]);
					float contamination = stof((*fields)[1]);
					delete fields;
					cout << "Result: " << completeness << " " << contamination << endl;

					//if(completeness > 0.4){
					//	for(const string& unitigName : hifiasmUnitigNames){
					//		fileHifiasmAll << unitigName << "," << clusterIndex << endl;
					//	}
					//}

					//if(contamination > 0.2){
					//	getchar();
					//}

				}
				*/
				

				/*
				//for(u_int32_t nodeIndex : nodePath){
				//	if(_contigNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _contigNodeNames.end()){
				//		_contigNodeNames.erase(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				//	}
				//}

				u_int64_t size = nodePath.size();
				gzwrite(outputContigFile, (const char*)&size, sizeof(size));
				//gzwrite(_outputContigFile, (const char*)&isCircular, sizeof(isCircular));
				gzwrite(outputContigFile, (const char*)&nodePath[0], size * sizeof(u_int32_t));
				gzwrite(outputContigFile, (const char*)&nodePath_supportingReads[0], size * sizeof(u_int64_t));

				
				string unitigSequence;
				_toBasespace.createSequence(nodePath, unitigSequence);

				cout << "Write contigs: " << nodePath.size() << " " << unitigSequence.size() << endl;

				string header = ">ctg" + to_string(contigIndex) + '\n';
				gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
				unitigSequence +=  '\n';
				gzwrite(outputContigFile_fasta, (const char*)&unitigSequence[0], unitigSequence.size());
				*/
				/*
				const string& filename_binStats = _inputDir + "/binStats.txt";
				vector<string> bin = {unitigSequence};			
				int ret = computeBinStats(_inputDir, bin, filename_binStats);

				if(ret == 0){

					string lastLine = "";
					string line = "";
					ifstream file_binStats(filename_binStats);
					while (std::getline(file_binStats, line)){
						if(line.empty()) continue;
						lastLine = line;
					}
					file_binStats.close();
					
					vector<string>* fields = new vector<string>();
					GfaParser::tokenize(lastLine, fields, ',');
					float completeness = stof((*fields)[0]);
					float contamination = stof((*fields)[1]);
					delete fields;
					cout << "Result: " << completeness << " " << contamination << endl;

					if(_truthInputFilename != ""){

						if(completeness > 0.4){

							unordered_set<string> hifiasmWrittenUnitigs;

							for(u_int32_t nodeIndex : nodePath){

								u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

								if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
									for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
										if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
										fileHifiasmAll << unitigName << "," << contigIndex << endl;
										hifiasmWrittenUnitigs.insert(unitigName);

									}
								}
							}

						}
					}

				}
				*/

				//contigIndex += 1;

			//}

			int validContigs = 0;
			for(Contig& contig : contigs){
				if(contig._nodePath.size() > 0){
					validContigs += 1;
				}
			}
			cout << "\tNb valid contigs: " << validContigs << endl;
			/*
			if(validContigs % 100 == 0){

				for(size_t i=0; i<contigs.size(); i++){
					if(contigs[i]._nodePath.size() == 0) continue;
					string unitigSequenceLala;
					_toBasespace.createSequence(contigs[i]._nodePath, unitigSequenceLala);
					cout << "Contig length: " << unitigSequenceLala.size() << " " << contigs[i]._nodePath.size() << endl;
				}

				for(size_t i=0; i<contigs.size(); i++){
					for(size_t j=i+1; j<contigs.size(); j++){
						if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
						//cout << endl;
						//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
						//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
						u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
						if(nbShared == 0) continue;
						cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
					
						if(nbShared /((float)contigs[i]._nodePath.size()) > 0.2){

							cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
							unordered_set<u_int32_t> sharedElements;
							Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

							cout << endl;
							for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
								if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
									cout << "0";
								}
								else{
									cout << "1";
								}
							}
							cout << endl;

							getchar();
						}
						if(nbShared /((float)contigs[j]._nodePath.size()) > 0.2){
							
							cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
							
							unordered_set<u_int32_t> sharedElements;
							Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

							cout << endl;
							for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
								if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
									cout << "0";
								}
								else{
									cout << "1";
								}
							}
							cout << endl;

							getchar();
						}

					}
				}

			}

			cout << "----" << endl;
			*/
		}

		/*
		cout << contigs.size() << endl;
		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);
		//vector<Contig> contigs_reverse = contigs;
		//std::reverse(contigs_reverse.begin(), contigs_reverse.end());

		for(long i=0; i<contigs.size(); i++){
			for(long j=contigs.size()-1; j>=0; j--){
				if(i == j) continue;

				Contig& contig1 = contigs[i];
				Contig& contig2 = contigs[j];

				if(contig1._nodePath.size() == 0) continue;
				if(contig2._nodePath.size() == 0) continue;

				//cout << i << " " << j << " " << contigs.size()-j-1 << endl;

				u_int64_t nbShared = Utils::computeSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted);
				if(nbShared == 0) continue;

				if(nbShared / ((float)contig2._nodePath.size()) > 0.98){
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
					continue;
				}


				unordered_set<u_int32_t> sharedElements;
				Utils::collectSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted, sharedElements);


				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, false);
				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, true);
				
				if(contig2._nodePath.size() > _kminmerSize*2){

					vector<u_int32_t> nodePath_sorted;
					for(u_int32_t nodeIndex : contig2._nodePath){
						nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}
					std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
					contig2._nodePath_sorted = nodePath_sorted;
				}
				else{
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
				}

			}
		}
		*/


		
		/*
		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.2 || nbShared /((float)contigs[j]._nodePath.size()) > 0.2){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					getchar();
				//}

			}
		}
		*/


		ofstream file_contigToNode(outputFilename_fasta + ".nodes.csv");

		contigIndex = 0;
		for(Contig& contig : contigs){
			if(contig._nodePath.size() > 0){

				string unitigSequenceLala;
				_toBasespace.createSequence(contig._nodePath, unitigSequenceLala);

				/*
				string header = ">ctg" + to_string(contigIndex) + '\n';
				gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
				unitigSequenceLala +=  '\n';
				gzwrite(outputContigFile_fasta, (const char*)&unitigSequenceLala[0], unitigSequenceLala.size());

				contigIndex += 1;

				continue;
				*/

				
				size_t splitSize = 100000;

				if(unitigSequenceLala.size() < splitSize) continue;


				for (size_t i = 0; i < unitigSequenceLala.length(); i += splitSize) {

					string unitigSequence_split = unitigSequenceLala.substr(i, splitSize);
					//cout << str.substr(i, 4) << endl;

					if(unitigSequence_split.size() < splitSize/2) break;

					string header = ">ctg" + to_string(contigIndex) + '\n';
					gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
					unitigSequence_split +=  '\n';
					gzwrite(outputContigFile_fasta, (const char*)&unitigSequence_split[0], unitigSequence_split.size());


					//for(u_int32_t nodeIndex : contig._nodePath){
					//	file_contigToNode << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ";" << contigIndex << endl;
					//}

					contigIndex += 1;

				}
				

			}
		}

		gzclose(outputContigFile);
		gzclose(outputContigFile_fasta);
		fileHifiasmAll.close();
		file_contigToNode.close();

		extractContigKminmers(outputFilename_fasta);
		file_asmResult.close();
		
	}

	bool isContigBinned(const vector<u_int32_t>& nodePath, unordered_set<u_int32_t>& processedNodeNames){

		u_int64_t nbBinnedNodes = 0;

		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
				nbBinnedNodes += 1;
			}
		}

		float sharedRate = nbBinnedNodes / ((float) nodePath.size());

		return sharedRate > 0.3;
	}

	bool isContigAssembled(const vector<u_int32_t>& nodePath, unordered_map<u_int32_t, u_int32_t>& processedNodeNames){

		unordered_set<u_int32_t> distinctContigIndex;
		u_int64_t nbBinnedNodes = 0;

		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
				distinctContigIndex.insert(processedNodeNames[nodeName]);
			}
			else{
				return false;
			}
		}

		//float sharedRate = nbBinnedNodes / ((float) nodePath.size());

		//return sharedRate > 0.98;
		return distinctContigIndex.size() == 1;
	}

	ofstream _fileTestLala;

	void extractContigKminmers (const string& outputFilename_fasta){
		
		_fileTestLala = ofstream(outputFilename_fasta + ".nodes.csv");
		_fileTestLala << "Name,Color" << endl;

		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_inputDir + "/mdbg_nodes.gz");

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly2::extractContigKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		for(auto& it: _contigFeature._nodeNameDuplicate){
			if(it.second == 1){
				_fileTestLala << it.first << "," << "blue" << endl;
			}
			else{
				_fileTestLala << it.first << "," << "red" << endl;
			}
		}

		delete _mdbg;
		_fileTestLala.close();
	}


	void extractContigKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
		}

	}

	ofstream _file_contigToNode;

	void extractContigKminmers2 (const string& outputFilename_fasta){

		_file_contigToNode = ofstream(_inputDir + "/nodeToContig.csv");
		_file_contigToNode << "Name,Color" << endl;
		_contigFeature.loadAbundanceFile_metabat(_filename_abundance);

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly2::extractContigKminmers_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		_file_contigToNode.close();
		//vector<u_int32_t> nodePath = {933376, 1651014, 1762772, 732878};
		//vector<float> lala1;
		//vector<float> lala2;
		//_contigFeature.sequenceToAbundance(nodePath, lala1, lala2);

		//exit(1);
	}


	void extractContigKminmers_read2(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		if(seq.size() < 15000) return;
		//cout << seq.size() << endl;
		//if(readIndex == 20010) getchar(); 

		vector<float> composition;
		_contigFeature.sequenceToComposition(seq, composition);
		_contigFeature._contigCompositions[readIndex] = composition;

		_contigFeature._contigSequences[readIndex] = seq;
		unordered_set<u_int32_t> nodeNames;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//_file_contigToNode << nodeName << "," << readIndex << endl;
			//if(nodeName == 933376) cout << "HAHAHAHA" << endl;

			//if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;

			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			_contigFeature._nodeName_to_contigIndex[nodeName] = readIndex;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
			nodeNames.insert(nodeName);
		}
		
		//std::sort(nodeNames.begin(), nodeNames.end());
		_contigFeature._contigNodes[readIndex] = nodeNames;

	}

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}

	static bool ContigComparator_ByLength_Reverse(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}

	void dereplicateContig2(vector<Contig>& contigs, vector<u_int32_t> nodePath, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){


		if(nodePath.size() < _kminmerSize*2) return;

		vector<u_int32_t> component;
		for(u_int32_t nodeIndex : nodePath){
			component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}
		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	component.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}
		std::sort(component.begin(), component.end());

		bool isValid = true;

		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(contig._nodePath.size() > component.size()){
				return;
			}
			else{
				contig._nodePath.clear();
				contig._nodePath_sorted.clear();
			}
		}
		
		addContig(nodePath, contigs, processedNodeNames, nodeName_to_contigIndex, contigIndex);
	}

	void dereplicateContig(vector<Contig>& contigs, vector<u_int32_t> nodePath, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){

		vector<vector<u_int32_t>> contigsToAdds;
		/*
		vector<u_int32_t> nodePath_sorted;
		for(u_int32_t nodeIndex : nodePath){
			nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}
		std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
		contigResult = {nodePath, nodePath_sorted};

		return true;
		*/

		vector<u_int32_t> component;
		for(u_int32_t nodeIndex : nodePath){
			component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}
		std::sort(component.begin(), component.end());

		bool isValid = true;

		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(nbShared / ((float)component.size()) > 0.98){
			//if(nbShared / ((float)component.size()) > 0.5){
				return;
			}
		}

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(contig._nodePath.size() < nodePath.size()){
				if(nbShared / ((float)contig._nodePath.size()) > 0.98){
				//if(nbShared / ((float)contig._nodePath.size()) > 0.){
					contig._nodePath.clear();
					contig._nodePath_sorted.clear();
					continue;
				}
			}
		}

		//std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength_Reverse);

		bool hasInternalRepeat = false;

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			unordered_set<u_int32_t> sharedElements;
			Utils::collectSharedElements(contig._nodePath_sorted, component, sharedElements);


			//removeOverlap(contig._nodePath, nodePath, sharedElements, false);
			//removeOverlap(contig._nodePath, nodePath, sharedElements, true);
			
			if(contig._nodePath.size() > component.size()){
				removeOverlap(contig._nodePath, nodePath, sharedElements, false);
				
				component.clear();
				for(u_int32_t nodeIndex : nodePath){
					component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
				std::sort(component.begin(), component.end());
				Utils::collectSharedElements(contig._nodePath_sorted, component, sharedElements);

				removeOverlap(contig._nodePath, nodePath, sharedElements, true);

				component.clear();
				for(u_int32_t nodeIndex : nodePath){
					component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
				std::sort(component.begin(), component.end());
				Utils::collectSharedElements(contig._nodePath_sorted, component, sharedElements);

				hasInternalRepeat = removeInternalDuplicate(contig._nodePath, nodePath, sharedElements, contigsToAdds);
			}
			else{
				removeOverlap(nodePath, contig._nodePath, sharedElements, false);

				vector<u_int32_t> nodePath_sorted;
				for(u_int32_t nodeIndex : contig._nodePath){
					nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
				std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
				contig._nodePath_sorted = nodePath_sorted;
				Utils::collectSharedElements(contig._nodePath_sorted, component, sharedElements);

				removeOverlap(nodePath, contig._nodePath, sharedElements, true);

				nodePath_sorted.clear();
				for(u_int32_t nodeIndex : contig._nodePath){
					nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
				std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
				contig._nodePath_sorted = nodePath_sorted;
				Utils::collectSharedElements(contig._nodePath_sorted, component, sharedElements);

				bool hasInternalRepeat2 = removeInternalDuplicate(nodePath, contig._nodePath, sharedElements, contigsToAdds);
				//cout << hasInternalRepeat2 << endl;
				//getchar();
				if(hasInternalRepeat2){
					contig._nodePath.clear();
					contig._nodePath_sorted.clear();
					//cout << "crush" << endl;
					//getchar();
				}
				else{
					if(contig._nodePath.size() > _minimizerSize*2){
						nodePath_sorted.clear();
						for(u_int32_t nodeIndex : contig._nodePath){
							nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
						}
						std::sort(nodePath_sorted.begin(), nodePath_sorted.end());

						contig._nodePath_sorted = nodePath_sorted;
					}
					else{
						contig._nodePath.clear();
						contig._nodePath_sorted.clear();
					}
				}
				
			}

			/*
			if(contig._nodePath.size() > component.size()){
				cout << endl;
				cout << endl;
				cout << endl;
				for(size_t p=0; p<nodePath.size(); p++){
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath[p])) == sharedElements.end()){
						cout << "0";
					}
					else{
						cout << "1";
					}
				}
				cout << endl;

				cout << endl;
				u_int32_t overlapSize_left = getOverlapSize(contig._nodePath, nodePath, sharedElements);
				std::reverse(nodePath.begin(), nodePath.end());
				u_int32_t overlapSize_right = getOverlapSize(contig._nodePath, nodePath, sharedElements);
				std::reverse(nodePath.begin(), nodePath.end());

				cout << nodePath.size() << " " << overlapSize_left << " " << overlapSize_right << endl;
				overlapSize_right = nodePath.size() - overlapSize_right;

				for(size_t p=0; p<nodePath.size(); p++){

					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath[p])) == sharedElements.end()){
						cout << "0";
					}
					else{
						cout << "1";
					}

					if(overlapSize_left != -1 && overlapSize_left == p){
						cout << "LLLL";
					}
					if(overlapSize_right != -1 && overlapSize_right == p){
						cout << "RRRR";
					}
				}
				cout << endl;

			}
			else{
				cout << endl;
				cout << endl;
				cout << endl;
				for(size_t p=0; p<contig._nodePath.size(); p++){
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contig._nodePath[p])) == sharedElements.end()){
						cout << "0";
					}
					else{
						cout << "1";
					}
				}
				cout << endl;

				cout << endl;
				u_int32_t overlapSize_left = getOverlapSize(nodePath, contig._nodePath, sharedElements);
				std::reverse(contig._nodePath.begin(), contig._nodePath.end());
				u_int32_t overlapSize_right = getOverlapSize(nodePath, contig._nodePath, sharedElements);
				std::reverse(contig._nodePath.begin(), contig._nodePath.end());

				cout << contig._nodePath.size() << " " << overlapSize_left << " " << overlapSize_right << endl;

				overlapSize_right = contig._nodePath.size() - overlapSize_right;

				for(size_t p=0; p<contig._nodePath.size(); p++){

					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contig._nodePath[p])) == sharedElements.end()){
						cout << "0";
					}
					else{
						cout << "1";
					}

					if(overlapSize_left != -1 && overlapSize_left == p){
						cout << "LLLL";
					}
					if(overlapSize_right != -1 && overlapSize_right == p){
						cout << "RRRR";
					}
				}
				cout << endl;

			}

			getchar();
			*/

			/*
			if(contig._nodePath.size() < nodePath.size()){
				if(nbShared / ((float)contig._nodePath.size()) > 0.8){
					contig._nodePath.clear();
					contig._nodePath_sorted.clear();
					continue;
				}
			}
			else{
				if(nbShared / ((float)component.size()) > 0.8){
					isValid = false;
					break;
				}
			}
			
			
			unordered_set<u_int32_t> sharedElements;
			Utils::collectSharedElements(component, contig._nodePath_sorted, sharedElements);

			u_int32_t overlapSize_left = getOverlapSize(nodePath, contig._nodePath, sharedElements);
			std::reverse(nodePath.begin(), nodePath.end());
			u_int32_t overlapSize_right = getOverlapSize(nodePath, contig._nodePath, sharedElements);
			std::reverse(nodePath.begin(), nodePath.end());

			//cout << endl;
			//cout << "Existing: " << contig._nodePath.size() << " " << (nbShared/((float)contig._nodePath.size())) << endl;
			//cout << "New: " << nodePath.size() << " " << component.size() << " " << (nbShared/((float)component.size())) << endl;
			//cout << "Nb shared nodes: " << sharedElements.size() << endl;
			//cout << "Overlaps: " << overlapSize_left << " " << overlapSize_right << endl;

			if(overlapSize_left != 0 || overlapSize_right != 0){
				
				if(contig._nodePath.size() < nodePath.size()){
					if(overlapSize_left > overlapSize_right){
						contig._nodePath.erase(contig._nodePath.begin(), contig._nodePath.begin()+overlapSize_left);
					}
					else{
						std::reverse(contig._nodePath.begin(), contig._nodePath.end());
						contig._nodePath.erase(contig._nodePath.begin(), contig._nodePath.begin()+overlapSize_right);
						std::reverse(contig._nodePath.begin(), contig._nodePath.end());
					}
					vector<u_int32_t> nodePath_sorted;
					for(u_int32_t nodeIndex : contig._nodePath){
						nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}
					std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
					contig._nodePath_sorted = nodePath_sorted;
				}
				else{
					if(overlapSize_left > overlapSize_right){
						nodePath.erase(nodePath.begin(), nodePath.begin()+overlapSize_left);
					}
					else{
						std::reverse(nodePath.begin(), nodePath.end());
						nodePath.erase(nodePath.begin(), nodePath.begin()+overlapSize_right);
						std::reverse(nodePath.begin(), nodePath.end());
					}
					component.clear();
					for(u_int32_t nodeIndex : nodePath){
						component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}
					std::sort(component.begin(), component.end());
				}
			}

			
			nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(contig._nodePath.size() < nodePath.size()){
				if(nbShared / ((float)contig._nodePath.size()) > 0.8){
					contig._nodePath.clear();
					contig._nodePath_sorted.clear();
					continue;
				}
			}
			else{
				if(nbShared / ((float)component.size()) > 0.8){
					isValid = false;
					break;
				}
			}
			

			//if((nbShared/((float)contig._nodePath.size())) > 0.95) getchar();
			//if(overlapSize_left > 0 || overlapSize_right > 0) getchar();
			*/
		}

		if(!hasInternalRepeat){
			addContig(nodePath, contigs, processedNodeNames, nodeName_to_contigIndex, contigIndex);
		}

		for(vector<u_int32_t>& nodePath: contigsToAdds){
			dereplicateContig(contigs, nodePath, processedNodeNames, nodeName_to_contigIndex, contigIndex);
			//addContig(nodePath, contigs);
		}

	}

	void removeOverlap(const vector<u_int32_t>& nodePath, vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements, bool right){

		/*
		cout << endl;
		cout << endl;
		cout << endl;
		for(size_t p=0; p<nodePath.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;

		cout << endl;
		for(size_t p=0; p<nodePath_shorter.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;
		*/

		if(right){
			std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
		}
		
		u_int32_t overlapSize = getOverlapSize(nodePath, nodePath_shorter, sharedElements);

		if(right && overlapSize > 0) overlapSize += 1;

		if(overlapSize > 0 && overlapSize < nodePath_shorter.size()){
			nodePath_shorter.erase(nodePath_shorter.begin(), nodePath_shorter.begin() + overlapSize);
		}

		if(right){
			std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
		}

	}

	u_int32_t getOverlapSize(const vector<u_int32_t>& nodePath, const vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements){

		if(sharedElements.size() == 0) return 0;

		size_t uniqueRunLength = 0;
		bool isUnique = true;

		u_int32_t nonUniquePos = -1;

		for(size_t i=0; i<nodePath_shorter.size(); i++){

			//cout << i << ": " << " " << (sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()) << " " << uniqueRunLength << "    " << nonUniquePos << endl;

			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()){
				if(isUnique) uniqueRunLength += 1;
				isUnique = true;
			}
			else{
				nonUniquePos = i;
				isUnique = false;
				uniqueRunLength = 0;
			}

			if(uniqueRunLength > 10) break;

		}

		return nonUniquePos;
		/*
		if(sharedElements.size() == 0) return 0;

		//u_int32_t maxOverlapSize = ; //nodePath1.size() * (sharedElement.size() / ((float) nodePath1.size()));

		u_int32_t maxMiss = sharedElements.size() * 0.03;
		if(maxMiss == 0) maxMiss = 1;
		u_int32_t nbElements = 0;
		u_int32_t nbSharedElements = 0;
		u_int32_t nbMiss = 0;
		//u_int32_t overlapSize = 0;
		//u_int32_t nbUnsharedElements = 0;

		for(size_t i=0; i<sharedElements.size(); i++){
			//cout << i << " " << nbSharedElements << " " << sharedElements.size() << endl;
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath1[i])) != sharedElements.end()){
				nbSharedElements += 1;
				//nbSharedElements += 1;
			}
			else{
				nbMiss += 1;
			}

			nbElements += 1;

			if(nbMiss >= maxMiss){
				float similarity = nbSharedElements / ((float)nbElements);
				if(similarity < 0.97) break;
			}
			//else{
			//	overlapSize += 1;
			//}

			//if(nbUnsharedElements > )

			//float similarity = nbSharedElements / ((float)nbElements);
			//if(similarity < 0.97) break;
		}

		//cout << "overlap sim: " << nbSharedElements / ((float)sharedElements.size()) << endl;
		//if(nbSharedElements / ((float)sharedElements.size()) > 0.9){
		//	return sharedElements.size();
		//}
		if(nbSharedElements >= 1) nbSharedElements -= 1;

		return nbSharedElements;
		//if(overlapSize >= 1) overlapSize -=1;

		//return overlapSize;
		*/
	}
	

	bool removeInternalDuplicate(const vector<u_int32_t>& nodePath, vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements, vector<vector<u_int32_t>>& contigsToAdd){

		/*
		cout << endl;
		cout << endl;
		cout << endl;
		for(size_t p=0; p<nodePath.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;

		cout << endl;
		for(size_t p=0; p<nodePath_shorter.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;
		*/
		pair<u_int32_t, u_int32_t> repeatPos;
		bool isInternalDuplicate = getInternalDuplicateSize(nodePath, nodePath_shorter, sharedElements, repeatPos);


		cout << endl;
		for(size_t p=0; p<nodePath_shorter.size(); p++){
			if(isInternalDuplicate){
				if(p == repeatPos.first) cout << "||";
				if(p == repeatPos.second) cout << "||||";
			}
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;

		if(!isInternalDuplicate) return false;

		vector<u_int32_t> nodePath1;
		vector<u_int32_t> nodePath2;

		for(size_t i=0; i<repeatPos.first && i<nodePath_shorter.size(); i++){
			nodePath1.push_back(nodePath_shorter[i]);
		}

		for(size_t i=repeatPos.second; i<nodePath_shorter.size(); i++){
			nodePath2.push_back(nodePath_shorter[i]);
		}

		/*
		cout << endl;
		cout << endl;
		for(size_t p=0; p<nodePath1.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath1[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;
		for(size_t p=0; p<nodePath2.size(); p++){
			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath2[p])) == sharedElements.end()){
				cout << "0";
			}
			else{
				cout << "1";
			}
		}
		cout << endl;
		*/

		//cout << nodePath1.size() << endl;
		//addContig(nodePath1, contigs);
		//cout << nodePath2.size() << endl;
		//addContig(nodePath2, contigs);
		contigsToAdd.push_back(nodePath1);
		contigsToAdd.push_back(nodePath2);

		cout << "internal repeat" << endl;
		getchar();

		/*
		if(right){
			std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
		}
		
		u_int32_t overlapSize = getOverlapSize(nodePath, nodePath_shorter, sharedElements);

		if(right && overlapSize > 0) overlapSize += 1;

		if(overlapSize > 0 && overlapSize < nodePath_shorter.size()){
			nodePath_shorter.erase(nodePath_shorter.begin(), nodePath_shorter.begin() + overlapSize);
		}

		if(right){
			std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
		}
		*/

		return true;

	}

	void addContig(vector<u_int32_t>& nodePath, vector<Contig>& contigs, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){

		if(nodePath.size() <= _kminmerSize*2) return;

		vector<u_int32_t> nodePath_sorted;
		for(u_int32_t nodeIndex : nodePath){
			nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}

		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}

		std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
		
		contigs.push_back({nodePath, nodePath_sorted});



		//bool isLala = false;
		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			nodeName_to_contigIndex[nodeName] = contigIndex;
			//if(nodeName == 289994){
			//	isLala = true;
			//}
			processedNodeNames.insert(nodeName);
		}

		contigIndex += 1;

		/*
		if(isLala){
			cout << "Size: " << nodePath.size() << endl;
			for(u_int32_t nodeIndex : nodePath){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				cout << nodeName << endl;
			}


			unordered_set<u_int32_t> component;
			_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(289994, false)), component);

			unordered_set<u_int32_t> validNodes;
			for (u_int32_t unitigIndex : component){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					validNodes.insert(nodeName);
				}
			}
			string outputFilename = _inputDir + "/minimizer_graph_sub_2.gfa";
			GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			

			getchar();
		}
		*/


	}

	bool getInternalDuplicateSize(const vector<u_int32_t>& nodePath, const vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements, pair<u_int32_t, u_int32_t>& repeatPos){

		if(sharedElements.size() == 0) return false;

		size_t repeatedRunLength = 0;
		u_int32_t repeatStartPos = -1;

		bool hasInternalRepeat = false;

		for(size_t i=0; i<nodePath_shorter.size(); i++){

			//cout << i << ": " << " " << (sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()) << " " << uniqueRunLength << "    " << nonUniquePos << endl;

			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()){
				repeatStartPos = -1;
				repeatedRunLength = 0;
			}
			else{
				if(repeatStartPos == -1) repeatStartPos = i;
				repeatedRunLength += 1;
			}

			if(repeatedRunLength > 10){
				hasInternalRepeat = true;
				break;
			}
			//if(repeatedRunLength > 10) break;

		}


		if(!hasInternalRepeat) return false;
			

		
		size_t uniqueRunLength = 0;
		bool isUnique = false;
		long i = repeatStartPos;
		u_int32_t nonUniquePos_left = repeatStartPos;

		//Expand left
		while(true){

			if(i < 0) break;

			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()){
				if(isUnique) uniqueRunLength += 1;
				isUnique = true;
			}
			else{
				nonUniquePos_left = i;
				isUnique = false;
				uniqueRunLength = 0;
			}

			if(uniqueRunLength > 10) break;

			i -= 1;
		}

		uniqueRunLength = 0;
		isUnique = false;
		i = repeatStartPos;
		u_int32_t nonUniquePos_right = repeatStartPos;

		//Expand right
		while(true){

			if(i >= nodePath_shorter.size()) break;

			if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(nodePath_shorter[i])) == sharedElements.end()){
				if(isUnique) uniqueRunLength += 1;
				isUnique = true;
			}
			else{
				nonUniquePos_right = i;
				isUnique = false;
				uniqueRunLength = 0;
			}

			if(uniqueRunLength > 10) break;
			
			i += 1;
		}

		repeatPos.first = nonUniquePos_left;
		repeatPos.second = nonUniquePos_right+1;

		//cout << repeatPos.first << " " << repeatPos.second << endl;
		//cout << "internal repeat" << endl;
		//getchar();

		return true;

	}

	/*
	struct ContigFeatures2{
		u_int32_t _unitigIndex;
		vector<float> _composition;
		vector<float> _abundance;
		vector<float> _abundanceVar;
		string _sequenceInit;
		float _length;
		vector<u_int32_t> _nodeNames;
	};

	struct Bin{
		ContigFeatures2 _model;
		vector<ContigFeatures2> _binnedContigs;
	};

	static bool ContigFeatureComparator_ByLength(const ContigFeatures2 &a, const ContigFeatures2 &b){
		return a._length > b._length;
	}

	vector<ContigFeatures2> _contigFeatures2;

	void generateFastaExpanded(const string& outputFilename){

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 


		Assembly assembly(_unitigDatas);

		float prevCutoff = -1;


		unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		for(Unitig& u : _graph->_startingUnitigstest){
			
			cout << processedUnitigs << " " << _graph->_startingUnitigstest.size() << endl;
			processedUnitigs += 1;

			bool isProcessed = false;
			for(u_int32_t nodeIndex : u._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
					isProcessed = true;
					break;
				}
			}
			if(isProcessed) continue;

			
			float cutoff = u._abundance*0.1;

			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff){
				_graph->loadState2(cutoff, u._startNode, _unitigDatas);
				prevCutoff = realCutoff;
			}


			const Unitig& unitigCurrent = _graph->nodeIndex_to_unitig(u._startNode);
			//vector<u_int32_t> nodePathInit = unitigCurrent._nodes;

			string unitigSequenceInit;
			_toBasespace.createSequence(unitigCurrent._nodes, unitigSequenceInit);

			for(u_int32_t nodeIndex : unitigCurrent._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				processedNodeNames.insert(nodeName);
			}


			vector<u_int32_t> nodePath;

			if(unitigSequenceInit.size() < 100000){
				assembly.solveBin2(u._startNode, 0, _graph, 0, 0, false, nodePath);
			}
			else{
				nodePath = unitigCurrent._nodes;
			}


			string unitigSequence;
			_toBasespace.createSequence(nodePath, unitigSequence);

			cout << endl << "\tUnitig: " << u._length << " " <<  unitigSequenceInit.size() << " " << unitigSequence.size() << endl;


			vector<float> composition;
			_contigFeature.sequenceToComposition(unitigSequence, composition);





			vector<u_int32_t> nodeNames;
			for(u_int32_t nodeIndex : nodePath){
				nodeNames.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}

			ContigFeatures2 contigFeature = {u._index, composition, {}, {}, unitigSequenceInit, unitigSequence.size(), nodeNames};
			_contigFeatures2.push_back(contigFeature);

		}

        std::sort(_contigFeatures2.begin(), _contigFeatures2.end(), ContigFeatureComparator_ByLength);

		vector<Bin> bins;

		for(const ContigFeatures2& f : _contigFeatures2){

			cout << bins.size() << endl;

			bool isBinned = false;

			for(Bin& bin : bins){

				if( -log10(_contigFeature.computeCompositionProbability(f._composition, bin._model._composition)) < 5){
					bin._binnedContigs.push_back(f);
					isBinned = true;
					break;
				}
			}

			if(isBinned) continue;

			Bin bin = {f, {f}};
			bins.push_back(bin);
		}
		//for(float t=0.01; t<10; t++){

		//}


		ofstream file_output(clusterDir + "/bin_nodes.csv");
		file_output << "Name,Colour" << endl;

		ofstream fileHifiasmAll(clusterDir + "/bin_hifiasm.csv");
		fileHifiasmAll << "Name,Colour" << endl;

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);
		file_binStats.close();

		int binIndex = 0;
		for(const Bin& bin : bins){
			
			vector<string> binSequences;

			for(const ContigFeatures2& contig : bin._binnedContigs){
				binSequences.push_back(contig._sequenceInit);
			}
			
			int ret = computeBinStats(clusterDir, binSequences, filename_binStats);
			//cout << ret << endl;

			if(ret == 0){

				string lastLine = "";
				string line = "";
				ifstream file_binStats(filename_binStats);
				while (std::getline(file_binStats, line)){
					if(line.empty()) continue;
					lastLine = line;
				}
				file_binStats.close();
				
				vector<string>* fields = new vector<string>();
				GfaParser::tokenize(lastLine, fields, ',');
				float completeness = stof((*fields)[0]);
				float contamination = stof((*fields)[1]);
				delete fields;
				cout << "Result: " << completeness << " " << contamination << endl;

				if(completeness > 0.4){
					for(const ContigFeatures2& contig : bin._binnedContigs){
						for(u_int32_t nodeName : contig._nodeNames){
							file_output << nodeName << "," << binIndex << endl;

							if(!_truthInputFilename.empty()){
							
								unordered_set<string> hifiasmWrittenUnitigs;

								if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
									for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
										if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
										fileHifiasmAll << unitigName << "," << binIndex << endl;
										hifiasmWrittenUnitigs.insert(unitigName);
									}
								}
							}
						}
					}
				}

			}


			binIndex += 1;
		}

		file_output.close();
		fileHifiasmAll.close();
		//gzclose(basespaceUnitigFile);
		exit(1);
	}
	*/

	static bool UnitigComparator_ByLength2(const UnitigLength &a, const UnitigLength &b){
		return a._length > b._length;
	}

	float _minUnitigAbundance;

	void execute_binning2(){


		vector<vector<string>> bins;
		

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;

		u_int64_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;

		unordered_set<u_int32_t> processedNodeNames;
		unordered_set<u_int32_t> processedContigIndex;
		u_int64_t processedUnitigs = 0;


		ofstream file_bin_all(_inputDir + "/bin_all.csv");
		file_bin_all << "Name,Color" << endl;





		vector<float> allCutoffs;

		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			allCutoffs.push_back(saveState._abundanceCutoff_min);
			cout << saveState._abundanceCutoff_min << endl;
		}
		
		std::reverse(allCutoffs.begin(), allCutoffs.end());

		for(float cutoff : allCutoffs){

			_graph->loadState2(cutoff, -1, _unitigDatas);

			vector<UnitigLength> startingUnitigs;

			_minUnitigAbundance = cutoff / 0.2;

			for(const Unitig& unitig : _graph->_unitigs){
				//if(unitig._nbNodes <= _kminmerSize*2) continue;
				//if(unitig._length < 100000) continue;

				if(unitig._abundance < _minUnitigAbundance) continue;

				u_int32_t contigIndex;
				string unitigSequence;
				//_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model);
				_contigFeature.nodepathToContigSequence(unitig._nodes, unitigSequence, contigIndex);
				if(unitigSequence.size() < 50000) continue;

				if(processedContigIndex.find(contigIndex) != processedContigIndex.end()) continue;

				startingUnitigs.push_back({unitig._length, unitig._abundance, unitig._startNode});
			}

			std::sort(startingUnitigs.begin(), startingUnitigs.end(), UnitigComparator_ByLength2);


			for(const UnitigLength& unitigLength : startingUnitigs){

				u_int32_t source_nodeIndex = unitigLength._startNodeIndex;
				const Unitig& unitig = _graph->nodeIndex_to_unitig(source_nodeIndex);

				cout << cutoff << " " << processedUnitigs << " " << startingUnitigs.size() << "     " << unitig._length << endl;
				processedUnitigs += 1;

				//if(isContigBinned(unitig._nodes, processedNodeNames)) continue;
				//if(isContigBinned(unitig._nodes, processedNodeNames)) continue;


				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);


				cout << endl << "Starting unitig: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._length << endl;



				//u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
				//const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];


				for(u_int32_t nodeIndex : unitig._nodes){
					processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}

				binByReadpath(unitig._startNode, processedNodeNames, processedContigIndex, clusterDir, filename_binStats, fileHifiasmAll, fileComponentNodeAll, clusterIndex);


			}
		}




		file_groundTruth.close();
		//file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();
		file_binStats.close();
		file_bin_all.close();

	}


	void generateContigs2(const string& outputFilename, const string& outputFilename_fasta){

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
		if(!fs::exists (path)){
			fs::create_directory(path);
		} 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(_inputDir + "/binning_results_hifiasm.csv");
		fileHifiasmAll << "Name,Colour" << endl;

		gzFile outputContigFile = gzopen(outputFilename.c_str(),"wb");
		gzFile outputContigFile_fasta = gzopen(outputFilename_fasta.c_str(),"wb");

		ofstream file_asmResult = ofstream(_inputDir + "/binning_results.csv");
		file_asmResult << "Name,Colour" << endl;

		//Assembly assembly(_unitigDatas, _contigFeature);

		float prevCutoff = -1;


		unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		//u_int64_t contigIndex = 0;
		u_int64_t __loadState2_index = 0;

		vector<Contig> contigs;
		vector<vector<u_int32_t>> clusters;













		vector<float> allCutoffs;

		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			allCutoffs.push_back(saveState._abundanceCutoff_min);
			cout << saveState._abundanceCutoff_min << endl;
		}
		
		std::reverse(allCutoffs.begin(), allCutoffs.end());



		u_int64_t contigIndexLala = 0;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_contigIndex;

		/*
		vector<u_int32_t> unitigLength_cutoffs = {100000, 50000, 30000, 10000, 5000};
		//vector<vector<UnitigLength>> startingUnitigs;
		//startingUnitigs.resize(unitigLength_cutoffs.size());		

		for(size_t i=0; i<unitigLength_cutoffs.size(); i++){
				
			u_int32_t unitigLength_cutoff_min = unitigLength_cutoffs[i];
			u_int32_t unitigLength_cutoff_max = -1;
			if(i > 0){
				unitigLength_cutoff_max = unitigLength_cutoffs[i-1];
			}

			//vector<UnitigLength>& unitigLengths = _startingUnitigs[i];
			*/

			for(float cutoff : allCutoffs){

				vector<UnitigLength> startingUnitigs;

				_graph->loadState2(cutoff, -1, _unitigDatas);
				_minUnitigAbundance = cutoff / 0.4;


				for(const Unitig& unitig : _graph->_unitigs){
					//cout << unitig._length << " " << unitig._abundance << endl;
					if(unitig._index % 2 == 1) continue;
					/*
					if(unitig._length < unitigLength_cutoff_min) continue;
					if(unitig._length > unitigLength_cutoff_max) continue;
					*/
					if(unitig._abundance < _minUnitigAbundance) continue;
					if(unitig._nbNodes < _kminmerSize*2) continue;
					//if(cutoff == 0 && unitig._length < 10000) continue;
					if(cutoff == 0) continue;

					startingUnitigs.push_back({unitig._length, unitig._abundance, unitig._startNode});
				}

				std::sort(startingUnitigs.begin(), startingUnitigs.end(), UnitigComparator_ByLength2);

				//cout << "Cutoff: " << unitigLength_cutoff_min << "-" << unitigLength_cutoff_max << " " << cutoff << endl;
					
				for(const UnitigLength& unitigLength : startingUnitigs){

					u_int32_t source_nodeIndex = unitigLength._startNodeIndex;
					const Unitig& unitig = _graph->nodeIndex_to_unitig(source_nodeIndex);

					//cout << "Cutoff: " << unitigLength_cutoff_min << "-" << unitigLength_cutoff_max << " " << cutoff << " " << processedUnitigs << " " << startingUnitigs.size() << "     " << unitig._length << endl;
					processedUnitigs += 1;

					/*
					bool isLala = false;
					for(u_int32_t nodeIndex : unitig._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						if(nodeName == 289994){
							isLala = true;
						}
					}

					if(isLala){

						unordered_set<u_int32_t> component;
						_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(289994, false)), component);

						unordered_set<u_int32_t> validNodes;
						for (u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
								validNodes.insert(nodeName);
							}
						}
						string outputFilename = _inputDir + "/minimizer_graph_sub_3.gfa";
						GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			
						cout << component.size() << endl;
						getchar();
					}
					*/

					if(isContigAssembled(unitig._nodes, nodeName_to_contigIndex)) continue;


					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);


					cout << endl << "Starting unitig: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._length << endl;


					for(u_int32_t nodeIndex : unitig._nodes){
						processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}


					vector<u_int32_t> nodePath = unitig._nodes;
					vector<u_int64_t> nodePath_supportingReads;
					//assembly.solveBin2(unitig._startNode, unitig._abundance, _graph, 0, 0, false, nodePath, nodePath_supportingReads, 0);




					string unitigSequenceExpanded;
					_toBasespace.createSequence(nodePath, unitigSequenceExpanded);
					cout << "\tExpanded contigs: " << nodePath.size() << " " << unitigSequenceExpanded.size() << endl;

					dereplicateContig2(contigs, nodePath, processedNodeNames, nodeName_to_contigIndex, contigIndexLala);





					/*
					nodePath = contig._nodePath;
					for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
						nodePath.erase(nodePath.begin());
					}
					for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
						nodePath.pop_back();
					}
					*/

					/*
					//bool isValid = true;
					Contig contig;// = {nodePath, {}};
					bool isValid = dereplicateContig(contigs, nodePath, contig);


					
					if(isValid && nodePath.size() > _kminmerSize*2){

						contigs.push_back(contig);
						contigIndex += 1;

					}
					*/

					int validContigs = 0;
					for(Contig& contig : contigs){
						if(contig._nodePath.size() > 0){
							validContigs += 1;
						}
					}
					cout << "\tNb valid contigs: " << validContigs << endl;
					
					/*
					if(validContigs % 100 == 0){

						for(size_t i=0; i<contigs.size(); i++){
							if(contigs[i]._nodePath.size() == 0) continue;
							string unitigSequenceLala;
							_toBasespace.createSequence(contigs[i]._nodePath, unitigSequenceLala);
							cout << "Contig length: " << unitigSequenceLala.size() << " " << contigs[i]._nodePath.size() << endl;
						}

						for(size_t i=0; i<contigs.size(); i++){
							for(size_t j=i+1; j<contigs.size(); j++){
								if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
								//cout << endl;
								//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
								//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
								u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
								if(nbShared == 0) continue;
								cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
							
								if(nbShared > 0){
								//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05){

									cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
									unordered_set<u_int32_t> sharedElements;
									Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

									cout << endl;
									for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
										if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
											cout << "0";
										}
										else{
											cout << "1";
										}
									}
									cout << endl;

									getchar();
								}
								if(nbShared > 0){
								//if(nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
									
									cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
									
									unordered_set<u_int32_t> sharedElements;
									Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

									cout << endl;
									for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
										if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
											cout << "0";
										}
										else{
											cout << "1";
										}
									}
									cout << endl;

									getchar();
								}

							}
						}

					}
					*/

					cout << "----" << endl;
					
				}







			}


			

		//}








		/*
		for(const Unitig& unitig : _graph->_unitigs){
			//if(unitig._nbNodes <= _kminmerSize*2) continue;
			if(unitig._length < 20000) continue;


			if(unitig._abundance < _minUnitigAbundance) continue;

			startingUnitigs.push_back({unitig._length, unitig._abundance, unitig._startNode});
		}
		*/



		/*
		cout << contigs.size() << endl;
		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);
		//vector<Contig> contigs_reverse = contigs;
		//std::reverse(contigs_reverse.begin(), contigs_reverse.end());

		for(long i=0; i<contigs.size(); i++){
			for(long j=contigs.size()-1; j>=0; j--){
				if(i == j) continue;

				Contig& contig1 = contigs[i];
				Contig& contig2 = contigs[j];

				if(contig1._nodePath.size() == 0) continue;
				if(contig2._nodePath.size() == 0) continue;

				//cout << i << " " << j << " " << contigs.size()-j-1 << endl;

				u_int64_t nbShared = Utils::computeSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted);
				if(nbShared == 0) continue;

				if(nbShared / ((float)contig2._nodePath.size()) > 0.98){
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
					continue;
				}


				unordered_set<u_int32_t> sharedElements;
				Utils::collectSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted, sharedElements);


				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, false);
				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, true);
				
				if(contig2._nodePath.size() > _kminmerSize*2){

					vector<u_int32_t> nodePath_sorted;
					for(u_int32_t nodeIndex : contig2._nodePath){
						nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}
					std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
					contig2._nodePath_sorted = nodePath_sorted;
				}
				else{
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
				}

			}
		}
		*/

		/*
		cout << "check 1" << endl;
		
		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					getchar();
				}

			}
		}
		


		//remove k-1 overlaps

		for(size_t i=0; i<contigs.size(); i++){

			Contig& c = contigs[i];
			if(c._nodePath.size() == 0) continue;

			c._nodePath_sorted.clear();
			for(u_int32_t nodeIndex : c._nodePath){
				c._nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}
			std::sort(c._nodePath_sorted.begin(), c._nodePath_sorted.end());
		}


		cout << "check 2" << endl;

		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

				}

			}
		}



		for(size_t i=0; i<contigs.size(); i++){

			Contig& c1 = contigs[i];
			if(c1._nodePath.size() == 0) continue;

			for(size_t j=i+1; j<contigs.size(); j++){
				
				Contig& c2 = contigs[j];
				if(c2._nodePath.size() == 0) continue;

				unordered_set<u_int32_t> sharedElements;
				Utils::collectSharedElements(c1._nodePath_sorted, c2._nodePath_sorted, sharedElements);
				
				u_int64_t overlapSizeLeft = 0;
				u_int64_t overlapSizeRight = 0;

				//left
				for(size_t p=0; p<_kminmerSize-1; p++){
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(c2._nodePath[p])) != sharedElements.end()){
						overlapSizeLeft += 1;
					}
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(c2._nodePath[c2._nodePath.size() - p - 1])) != sharedElements.end()){
						overlapSizeRight += 1;
					}
				}

				if(overlapSizeLeft == _kminmerSize-1){
					for(size_t p=0; p<_kminmerSize-1; p++){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(c2._nodePath[p]);
						c2._nodePath_sorted.erase(std::remove(c2._nodePath_sorted.begin(), c2._nodePath_sorted.end(), nodeName), c2._nodePath_sorted.end());
					}
					c2._nodePath.erase(c2._nodePath.begin(), c2._nodePath.begin()+_kminmerSize-1);
				}

				if(overlapSizeRight == _kminmerSize-1){
					for(size_t p=0; p<_kminmerSize-1; p++){
						c2._nodePath.pop_back();
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(c2._nodePath[c2._nodePath.size() - p - 1]);
						c2._nodePath_sorted.erase(std::remove(c2._nodePath_sorted.begin(), c2._nodePath_sorted.end(), nodeName), c2._nodePath_sorted.end());
					}
				}

			}
		}




		cout << "check 3" << endl;

		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					getchar();
				}

			}
		}
		*/

		//ofstream file_contigToNode(outputFilename_fasta + ".nodes.csv");
		unordered_map<u_int32_t, u_int32_t> nodeCounts;

		u_int64_t contigIndex = 0;
		for(Contig& contig : contigs){
			if(contig._nodePath.size() > 0){

				for(u_int32_t nodeIndex : contig._nodePath){
					nodeCounts[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
				}


				string unitigSequenceLala;
				_toBasespace.createSequence(contig._nodePath, unitigSequenceLala);

				
				string header = ">ctg" + to_string(contigIndex) + '\n';
				gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
				unitigSequenceLala +=  '\n';
				gzwrite(outputContigFile_fasta, (const char*)&unitigSequenceLala[0], unitigSequenceLala.size());

				contigIndex += 1;

				continue;
				

				
				size_t splitSize = 20000;

				if(unitigSequenceLala.size() < splitSize) continue;


				for (size_t i = 0; i < unitigSequenceLala.length(); i += splitSize) {

					string unitigSequence_split = unitigSequenceLala.substr(i, splitSize);
					//cout << str.substr(i, 4) << endl;

					if(unitigSequence_split.size() < splitSize/2) break;

					string header = ">ctg" + to_string(contigIndex) + '\n';
					gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
					unitigSequence_split +=  '\n';
					gzwrite(outputContigFile_fasta, (const char*)&unitigSequence_split[0], unitigSequence_split.size());


					//for(u_int32_t nodeIndex : contig._nodePath){
					//	file_contigToNode << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ";" << contigIndex << endl;
					//}

					contigIndex += 1;

				}
				

			}
		}


		ofstream fileTestLala(outputFilename_fasta + ".nodes2.csv");
		fileTestLala << "Name,Color" << endl;

		for(auto& it: nodeCounts){
			if(it.second == 1){
				fileTestLala << it.first << "," << "blue" << endl;
			}
			else{
				fileTestLala << it.first << "," << "red" << endl;
			}
		}

		fileTestLala.close();


		gzclose(outputContigFile);
		gzclose(outputContigFile_fasta);
		fileHifiasmAll.close();
		//file_contigToNode.close();

		extractContigKminmers(outputFilename_fasta);
		file_asmResult.close();
		
	}


};





#endif

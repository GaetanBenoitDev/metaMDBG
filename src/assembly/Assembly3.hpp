
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

- singleton erroneous unitigs dans le graphe (si une erreur est isolé, ça créé qq noeuds attaché erroneous, on pourrait ptet supprimé les mini component peut couverte à voir)
*/

//:/mnt/chris-native/chris/Projects/AD_HiFi
//mnt/gpfs/Hackathon/meren_hifi

#ifndef MDBG_METAG_ASSEMBLY3
#define MDBG_METAG_ASSEMBLY3

#define PRINT_DEBUG_PREVRANK

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
//#include "toBasespace/ToBasespaceOnTheFly.hpp"
//#include "contigFeatures/ContigFeature.hpp"










class Assembly3 : public Tool{

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

	string _partitionDir;
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
	MDBG* _mdbgNoFilter;
	GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	u_int64_t _nbReads;
	int _nbCores;

	Assembly3(): Tool (){

	}

	~Assembly3(){

	}

	void execute (){

		loadGraph();
		//generateContigs2(_inputDir + "/contigs.nodepath.gz", _inputDir + "/contigs.fasta.gz");

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
		(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));



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
			_nbCores = result[ARG_NB_CORES].as<int>();
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

		ifstream file_data(_inputDir + "/data.txt");
		file_data.read((char*)&_nbReads, sizeof(_nbReads));
		file_data.close();

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << "Nb reads: " << _nbReads << endl;
		cout << endl;

		//getchar();

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";


		

	}

	unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    vector<bool> _isBubble;
	unordered_map<u_int32_t, vector<u_int64_t>> _nodeName_to_kminmerSequence;
	unordered_set<u_int32_t> _isNodeNameUnsupported;

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		_mdbgNoFilter = new MDBG(_kminmerSize);
		_mdbgNoFilter->load(_inputDir + "/mdbg_nodes_noFilter.gz");

		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		//extractContigKminmers2(_inputDir + "/contigs.fasta.gz");

		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}


		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		//file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		//file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize);
		_graph = graphSimplify;
		
		
		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		_graph->saveUnitigGraph("/home/gats/workspace/run/test_unitig_graph_lala.gfa");
		exit(1);
		//Generate unitigs
		//cout << "Indexing reads" << endl;
		//_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		//_graph->clear(0);
		//_graph->compact(false, _unitigDatas);
		//removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		
		//cout << "done" << endl;
	

		//if(_kminmerSize == 4){
			//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
			//graphSimplify->execute(5, _kminmerSize);
			//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);


			cout << "Cleanning graph 1" << endl;
			_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, true, false, true);
			_isBubble = _graph->_isBubble;
			
			//cout << "Cleanning graph 2" << endl;
			//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, false, false, false, false);

			//!

			cout << "done" << endl;
			//_graph->loadState2(0, -1, _unitigDatas);
		//}

		
		/*
		cout << _graph->_isNodenameRoundabout.size() << endl;
		
		
		ofstream file_correction(_inputDir + "/roundabout.csv");
		file_correction << "Name,Colour" << endl;
		for(u_int32_t nodeName : _graph->_isNodenameRoundabout){
			if(!_isBubble[BiGraph::nodeName_to_nodeIndex(nodeName, true)]){
				file_correction << nodeName << ",red" << endl;
			}
		}
		for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			if(_graph->_isBubble[nodeIndex]){
				//file_correction << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",red" << endl;
			}
		}
		file_correction.close();
		*/
		//cout << _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(12871, true))._nbNodes << " " << _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(12871, true))._length << endl;
		//getchar();
		
		
		//_graph->loadState2(100, -1, _unitigDatas);
		//_graph->saveGraph(_inputDir + "/minimizer_graph_contigs.gfa");
		_graph->loadState2(0, -1, _unitigDatas);
		generateUnitigs();

		_partitionDir = _inputDir + "/" + "partitions";
		fs::path path(_partitionDir);
		if(fs::exists (path)){
			fs::remove_all(path);
		} 
		fs::create_directory(path);

		//_graph->loadState2(0, -1, _unitigDatas);
		_graph->saveGraph(_inputDir + "/minimizer_graph_sub.gfa");
		//_graph->loadState2(7, -1, _unitigDatas);
		//_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		/*
		if(_kminmerSize == 4){
			//cout << "A ENLEVER" << endl;
			_graph->loadState2(10, -1, _unitigDatas);
			_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		}
		else{
			//cout << "A ENLEVER" << endl;
			_graph->loadState2(2.5, -1, _unitigDatas);
			_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		}
		*/

		//cout << Utils::computeSharedReads(_unitigDatas[8228], _unitigDatas[8229]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[6323], _unitigDatas[2778]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[8228], _unitigDatas[2778]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[6323], _unitigDatas[8229]) << endl;
		
		//cout << Utils::computeSharedReads(_unitigDatas[3076], _unitigDatas[4346]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[4506], _unitigDatas[4886]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[3076], _unitigDatas[4886]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[3031], _unitigDatas[8043]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[62], _unitigDatas[1106]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[9895], _unitigDatas[933]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[12321], _unitigDatas[933]) << endl;
		
		//for(auto& it : _evaluation_hifiasmGroundTruth_nodeNamePosition){
		//	if(it.second == 4790){
		//		cout << it.first << " " << _unitigDatas[it.first]._readIndexes.size() << endl;
		//	}
		//}


		//getchar();
		//debug_checkReads();
		//cout << endl << endl;
		//partitionReads();

		extractKminmerSequences();
		correctReads();
		//getchar();

		delete _mdbg;

		//_toBasespace.create(_inputDir);


		//_graph->loadState2(0, -1, _unitigDatas);
		//generateFasta(_inputDir + "/allContigs.fasta.gz");
		//generateFastaExpanded(_inputDir + "/allContigsExpanded.fasta.gz");

	}

	void extractKminmerSequences (){

		cout << "Extracting kminmer sequences" << endl;
		
		
		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");

		while (true) {

			u_int16_t size;
			kminmerFile.read((char*)&size, sizeof(size));

			if(kminmerFile.eof())break;

			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(size);
			kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(u_int64_t));

			u_int32_t nodeName;
			u_int32_t length;
			u_int32_t lengthStart;
			u_int32_t lengthEnd;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&length, sizeof(length));
			kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
			kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

			_nodeName_to_kminmerSequence[nodeName] = minimizerSeq;

		}

		kminmerFile.close();
	}

	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		if(_indexingContigs) readIndex += 2000000000ull;
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

	/*
	void indexReads_contig(const vector<u_int64_t>& minimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		readIndex += 2000000000;
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
	*/

	bool _indexingContigs;

	void removeUnsupportedEdges(const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges, GraphSimplify* graph){

		_indexingContigs = false;
		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly3::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		//_indexingContigs = true;
		//KminmerParser parser2(_inputDir + "/contig_data.txt", _minimizerSize, _kminmerSize, false);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		//auto fp2 = std::bind(&Assembly3::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser2.parse(fp2);
		
		//if(_filename_inputContigs != ""){
		//	KminmerParser parserContig(_filename_inputContigs, _minimizerSize, _kminmerSize, true);
		//	auto fpContig = std::bind(&Assembly3::indexReads_contig, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//	parserContig.parse(fpContig);
		//}

		for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

			vector<u_int32_t> successors;
			graph->getSuccessors(nodeIndex, 0, successors);

			for(u_int32_t successor : successors){
				
				u_int32_t nodeName_succ = BiGraph::nodeIndex_to_nodeName(successor);

				u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[nodeName], _unitigDatas[nodeName_succ]);
				if(nbSharedReads <= 1){
					_isNodeNameUnsupported.insert(nodeName);
					_isNodeNameUnsupported.insert(nodeName_succ);
				}
				//if(Utils::shareAnyRead(_unitigDatas[nodeName], _unitigDatas[nodeName_succ])){
				//	continue;
				//}

				//Unitig& unitig_succ = graph->_unitigs[graph->nodeIndex_to_unitigIndex(successor)];

				//removedEdges.push_back({unitig._endNode, successor});
				
				
			}

		}

		cout << "Nb unsupported nodes: " << _isNodeNameUnsupported.size() << endl;

		return;
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
	


	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers_read(const Read& read){

		u_int64_t readIndex = read._index;
		//ottalSize += strlen(read->seq.s);


		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

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
				//if("utg009434l" == string(read->name.s)){
				//	cout << nodeName << endl;
				//	_hifiasm_startingNodenames.insert(nodeName);
				//}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
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
		
		auto fp = std::bind(&Assembly3::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_truthInputFilename, true, false);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		delete _minimizerParser;
	}


	unordered_map<u_int32_t, ofstream> _readPartitions;

	
	void partitionReads(){

		//cout << "Partitionning reads" << endl;

		u_int32_t cutoffLevel = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			const string filename = _partitionDir + "/part_" + to_string(cutoffLevel) + ".gz";
			_readPartitions[cutoffLevel] = ofstream(filename);
			cutoffLevel += 1;
		}

		_graph->loadState2(0, -1, _unitigDatas);

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&Assembly3::partitionReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);

		for(auto& it : _readPartitions){
			it.second.close();
		}
	}

	void partitionReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){
	
		vector<u_int32_t> nodePath;

		for(const KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				//cout << "XXXXX ";
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodePath.push_back(nodeName);
			//cout << nodeName << " ";
		}
		//cout << endl;

		double n = 0;
		double sum = 0;
		unordered_set<u_int32_t> writtenUnitigs;

		vector<float> readpathAbudance_values;
		for(u_int32_t nodeName : nodePath){

			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, false);
			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
			const Unitig& u = _graph->_unitigs[unitigIndex];

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			for(u_int32_t nodeIndex : u._nodes){
				if(_graph->_isNodenameRoundabout.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _graph->_isNodenameRoundabout.end()) continue;
				if(_graph->_isBubble[nodeIndex]) continue;

				readpathAbudance_values.push_back(u._abundance);
				//cout << u._abundance << " ";
				n += 1;
				sum += u._abundance;
			}
		}
		//cout << endl;
		
		float cutoff = 0;

		if(n > 0){
			float readPathAbundance = sum / n; ////Utils::compute_median_float(readpathAbudance_values);
			//cout << "Read path abundance: " << readPathAbundance << " " << (sum / n) << endl;
			//if(readPathAbundance < 10) getchar();

			//if(readPathAbundance > 50) getchar();
			cutoff = readPathAbundance * 0.2;
			//cout << cutoff << endl;
		}

		//cutoff = 4;
		//if(_kminmerSize != 4) 
		//cutoff = 0;

		u_int32_t cutoffLevel = 0;
		float realCutoff = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			if(saveState._abundanceCutoff_min >= cutoff) break;
			realCutoff = saveState._abundanceCutoff_min;
			cutoffLevel += 1;
		}

		ofstream& file = _readPartitions[cutoffLevel];
		
		u_int32_t size = minimizers.size();
		file.write((const char*)&size, sizeof(size));
		file.write((const char*)&minimizers[0], size * sizeof(u_int64_t));
	}


	//ofstream file_correction;
	//gzFile _outputFile_correctedReads;
	ofstream _file_uncorrectedReads;

	void correctReads(){

		_nbreadsLala = 0;
		_nbUncorrectedReads = 0;

		_graph->loadState2(0, -1, _unitigDatas);
		_graph->_isBubble = _isBubble;

		//const string& outputFilename_correctedReads = _inputDir + "/correctedReads_" + to_string(_kminmerSize) + ".min.gz";
		//_outputFile_correctedReads = gzopen(outputFilename_correctedReads.c_str(),"wb");

		_file_uncorrectedReads = ofstream(_inputDir + "/read_uncorrected.txt");
		
		/*
		u_int32_t cutoffLevel = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){

			cout << cutoffLevel << " " << saveState._abundanceCutoff_min << endl;
			//getchar();

			_graph->loadState2(saveState._abundanceCutoff_min, -1, _unitigDatas);
			//_graph->saveGraph(_inputDir + "minimizer_graph_sub.gfa");
			const string readPartitionFilename = _partitionDir + "/part_" + to_string(cutoffLevel) + ".gz";
			
			KminmerParser parser(readPartitionFilename, _minimizerSize, _kminmerSize, false);
			auto fp = std::bind(&Assembly3::correctReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
			parser.parse(fp);





			//_readPartitions[cutoffLevel] = gzopen(filename.c_str(),"wb");
			cutoffLevel += 1;
		}
		*/



		
		//KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		//auto fp = std::bind(&Assembly3::correctReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parse(fp);
	

		KminmerParserParallel readParser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, _nbCores);
		readParser.parse(ReadCorrectionFunctor(*this));


		//file_correction.close();
		//gzclose(_outputFile_correctedReads);
		_file_uncorrectedReads.close();

		cout << "Nb reads: " << _nbreadsLala << endl;
		cout << "Nb uncorrected reads: " << _nbUncorrectedReads << endl;
		//getchar();
		
	}

	u_int64_t _nbreadsLala;
	u_int64_t _nbUncorrectedReads = 0;
	unordered_map<u_int32_t, vector<u_int64_t>> _nodeIndex_to_kminmerSequence;

	void writeCorrectedRead(const vector<u_int64_t>& readMinimizers, bool success, bool print_read){
		
		#pragma omp critical
		{

			//checkCorrectedSequence(readMinimizers, print_read);
			//if(success){
			u_int32_t size = readMinimizers.size();
			_file_uncorrectedReads.write((const char*)&size, sizeof(size));
			_file_uncorrectedReads.write((const char*)&readMinimizers[0], size*sizeof(u_int64_t));
			//}

			_nbreadsLala += 1;

			if(!success){
				_nbUncorrectedReads += 1;
			}
		
		}
	}

	/*
	void debug_checkReads(){

		cout << "Checking reads" << endl;

		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&Assembly3::debug_checkReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
	}

	void debug_checkReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){
	
		vector<u_int32_t> nodePath;
		bool printRead = false;

		for(const KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				nodePath.push_back(0);
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodePath.push_back(nodeName);

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()) continue;

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName] == 4790){
				printRead = true;
			}
		}

		if(printRead || readIndex == 3427){ //|| readIndex == 2010
			cout << readIndex << endl;
			_debug_readPath[readIndex] = nodePath;
			//cout << readIndex << ": " << endl;
			//for(u_int32_t nodeName : nodePath){
			//	cout << nodeName << " ";
			//}
			//cout << endl;
		}

	}
	*/



	
	void generateUnitigs(){
		//if(_kminmerSize < 10) return;

		const string& outputFilename = _inputDir + "/unitigs.nodepath.gz";

		cout << "Generating contigs" << endl;
		
		//const string& outputFilename = _inputDir + "/contigs.nodepath.gz";
		gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");

		unordered_set<u_int32_t> writtenUnitigs;

		//vector<float> readpathAbudance_values;
		for(const Unitig& u : _graph->_unitigs){
			//if(u._nbNodes < _kminmerSize*2) continue;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			vector<u_int32_t> nodepath = u._nodes;
			if(u._startNode == u._endNode){

				u_int32_t nodeIndex = u._endNode;
				for(size_t i=0; i<_kminmerSize; i++){
					
					vector<u_int32_t> successors;
					_graph->getSuccessors(nodeIndex, 0, successors);

					nodeIndex = successors[0];
					nodepath.push_back(nodeIndex);

				}

			}

			u_int64_t size = nodepath.size();

			//if(size < _kminmerSize*2) continue;

			gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
			gzwrite(outputContigFile_min, (const char*)&nodepath[0], size * sizeof(u_int32_t));
		}
		
		gzclose(outputContigFile_min);
	}
	

	void checkCorrectedSequence(const vector<u_int64_t> minimizers, bool print_read){
		
		//MDBG* mdbg = new MDBG(_kminmerSize);
		//mdbg->load(_inputDir + "/mdbg_nodes_init.gz");

		bool notGood = false;

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		vector<u_int64_t> minimizersPos; 
		vector<u_int64_t> rlePositions;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, -1, false);

		u_int64_t nbFoundMinimizers = 0;
		for(size_t i=0; i<kminmers.size(); i++){
			KmerVec& vec = kminmers[i];
			

			if(_mdbgNoFilter->_dbg_nodes.find(vec) == _mdbgNoFilter->_dbg_nodes.end()){
				//if(i==2){ nbFailed += 1; }
				cout << "Not good: " << i << endl;
				cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << " " << vec._kmers[3] << endl;

				notGood = true;
				continue;	
			}

			//u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;

			nbFoundMinimizers += 1;

		}

		if(print_read){
			cout << "\tNb minimizers: " << kminmers.size() << endl;
			cout << "\tFound minimizers: " << nbFoundMinimizers << endl;
			//if(kminmers.size() != nbFoundMinimizers) getchar();
			if(notGood) getchar();
		}
	}



	class ReadCorrectionFunctor {

		public:

		Assembly3& _assembly3;
		MDBG* _mdbg;
		GraphSimplify* _graph;
		unordered_map<u_int32_t, vector<u_int64_t>>& _nodeName_to_kminmerSequence;
    	size_t _kminmerSize;
		vector<UnitigData>& _unitigDatas;
		unordered_set<u_int32_t>& _isNodeNameUnsupported;

		ReadCorrectionFunctor(Assembly3& assembly3) : _assembly3(assembly3), _nodeName_to_kminmerSequence(assembly3._nodeName_to_kminmerSequence), _unitigDatas(assembly3._unitigDatas), _isNodeNameUnsupported(assembly3._isNodeNameUnsupported){
			_mdbg = _assembly3._mdbg;
			_graph = _assembly3._graph;
			_kminmerSize = _assembly3._kminmerSize;
		}

		ReadCorrectionFunctor(const ReadCorrectionFunctor& copy) : _assembly3(copy._assembly3), _nodeName_to_kminmerSequence(copy._nodeName_to_kminmerSequence), _unitigDatas(copy._unitigDatas), _isNodeNameUnsupported(copy._isNodeNameUnsupported){
			_mdbg = copy._mdbg;
			_graph = copy._graph;
			_kminmerSize = copy._kminmerSize;
		}





		//void correctReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		void operator () (const KminmerList& kminmerList) {
		
			u_int64_t readIndex = kminmerList._readIndex;

			const vector<u_int64_t>& minimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			//u_int32_t readSize = minimizers.size();
			//_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
			//_file_uncorrectedReads.write((const char*)&minimizers[0], readSize*sizeof(u_int64_t));
			//return;

			
			//if(readIndex >= _nbReads){
			//	u_int32_t readSize = minimizers.size();
			//	_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
			//	_file_uncorrectedReads.write((const char*)&minimizers[0], readSize*sizeof(u_int64_t));
			//	return;
			//}
			
			bool print_read = false;
			if(readIndex % 10000 == 0) print_read = true;

			//file_correction = ofstream(_inputDir + "/correction.csv");
			//file_correction << "Name,Color" << endl;

			if(print_read) cout << "-------------" << endl;
			if(print_read) cout << readIndex << endl;
			//if(print_read) cout << minimizers.size() << endl;


			if(print_read){
				cout << "\tRead original: " << endl;
				cout << "\t";
				for(u_int64_t m : minimizers){
					cout << m << " ";
				}
				cout << endl;

				cout << "\tRead nodes: " << endl;
				cout << "\t";
			}

			bool isError = false;
			//_nodeIndex_to_kminmerSequence.clear();
			vector<u_int32_t> nodepath;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& info = kminmersInfos[i];
				const KmerVec& vec = info._vec;

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					if(print_read) cout << "XXXXX ";
					isError = true;
					continue;
				}
				
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

				vector<u_int64_t> kminmerSeq = _nodeName_to_kminmerSequence[nodeName];

				if(_graph->_isBubble[BiGraph::nodeName_to_nodeIndex(nodeName, true)]){
					if(print_read) cout << nodeName << "B ";
					//continue;
				}
				else{
					if(print_read) cout << nodeName << " ";
				}

				if(info._isReversed){
					std::reverse(kminmerSeq.begin(), kminmerSeq.end());
				}

				nodepath.push_back(nodeIndex);
				//_nodeIndex_to_kminmerSequence[nodeIndex] = kminmerSeq;

				//cout << BiGraph::nodeToString(nodeIndex) << ": ";
				//for(u_int64_t m : kminmerSeq){
				//	cout << m << " ";
				//}
				//cout << endl;
			}
			if(print_read) cout << endl;




			vector<u_int32_t> nodePathSolid;
			vector<u_int64_t> readpath;
			applyReadCorrection(nodepath, readpath, print_read, minimizers, kminmersInfos, nodePathSolid);


			if(print_read){
				cout << "\tRead original: " << endl;
				cout << "\t";
				for(u_int64_t m : minimizers){
					cout << m << " ";
				}
				cout << endl;
			}

			if(print_read){
				cout << "\tRead minimizers corrected: " << endl;
				cout << "\t";
				for(u_int64_t minimizer : readpath){
					cout << minimizer << " ";
				}
				cout << endl;
			}

			bool success;

			//cout << nodePathSolid.size() << endl;
			if(nodePathSolid.size() != 2){
				if(print_read) cout << "\tcorrection failed" << endl;
				readpath = minimizers;
				success = false;
				//getchar();
			}
			else{
				success = true;
				extendReadpath(minimizers, kminmersInfos, nodePathSolid, readpath, print_read);
			}
			
			//readpath = minimizers;
			//vector<u_int32_t> contigpath;

			if(print_read){
				cout << "\tRead minimizers corrected extended: " << endl;
				cout << "\t";
				for(u_int64_t minimizer : readpath){
					cout << minimizer << " ";
				}
				cout << endl;
			}



			//if(readIndex == 17) getchar();

			//if(contigpath.size() == 0) return;
			//cout << "decommenter ici dans v final, contig size 0 a discard" << endl;

			//for(u_int32_t nodeIndex : contigpath){
			//	cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			//}
			//cout << endl;
			//getchar();



			//if(readIndex == 6) getchar();
			//getchar();
			
			_assembly3.writeCorrectedRead(readpath, success, print_read);
		}

		void fillCorrectedArea(u_int32_t nodeIndex_source, u_int32_t nodeIndex_dest, const vector<u_int64_t>& readMinimizers, vector<u_int64_t>& readpath, const vector<ReadKminmerComplete>& kminmersInfos, bool isFirstIteration, bool print_read){
			
			if(readpath.size() == 0){

				bool orientation;
				u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndex_source, orientation);
				vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];

				if(orientation){
					if(isFirstIteration){
						for(u_int64_t m : minimizerSeq){
							readpath.push_back(m);
							if(print_read)  cout << "Add: " << m << endl;
						}
					}
					else{
						if(print_read)  cout << "Add: " << minimizerSeq[minimizerSeq.size()-1] << endl;
						readpath.push_back(minimizerSeq[minimizerSeq.size()-1]);
					}
				}
				else{
					if(isFirstIteration){
						std::reverse(minimizerSeq.begin(), minimizerSeq.end());
						for(u_int64_t m : minimizerSeq){
							readpath.push_back(m);
							if(print_read)  cout << "Add: " << m << endl;
						}
					}
					else{
						if(print_read)  cout << "Add: " << minimizerSeq[0] << endl;
						readpath.push_back(minimizerSeq[0]);
					}

				}
			}


			bool orientation;
			u_int32_t nodeName_dest = BiGraph::nodeIndex_to_nodeName(nodeIndex_dest, orientation);
			vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_dest];

			if(orientation){
				if(print_read)  cout << "Add: " << minimizerSeq[minimizerSeq.size()-1] << endl;
				readpath.push_back(minimizerSeq[minimizerSeq.size()-1]);
			}
			else{
				if(print_read)  cout << "Add: " << minimizerSeq[0] << endl;
				readpath.push_back(minimizerSeq[0]);
			}
		}


		void fillUncorrectedArea(u_int32_t position_source, u_int32_t nodeIndex_dest, const vector<u_int64_t>& readMinimizers, vector<u_int64_t>& readpath, const vector<ReadKminmerComplete>& kminmersInfos, bool isFirstIteration, bool print_read){
			
			bool adding = false;


			for(size_t i=position_source; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& info = kminmersInfos[i];
				const KmerVec& vec = info._vec;

				

				if(!adding){
					
					u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
					u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

					if(readpath.size() == 0){
						if(isFirstIteration){
							for(size_t j=i; j<i+_kminmerSize; j++){
								u_int64_t minimizer = readMinimizers[j];
								if(print_read)  cout << "\tFill uncorrected area (start): " << j << " " << minimizer << endl;
								readpath.push_back(minimizer);
							}
						}
					}

					if(print_read)  cout << "Source position: " << nodeName << " " << i << " " << readMinimizers[i] << endl;
					adding = true;
				}
				else if(adding){
					u_int64_t minimizer = readMinimizers[i+_kminmerSize-1];
					readpath.push_back(minimizer);
					if(print_read)  cout << "\tFill uncorrected area: " << i << " " << minimizer << endl;
				}

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

				if(nodeIndex == nodeIndex_dest) return;

			}

		}


		void applyReadCorrection(const vector<u_int32_t>& nodePath, vector<u_int64_t>& readpath, bool print_read, const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, vector<u_int32_t>& nodePathSolid){
			

			bool correctionFailed = false;

			float minSupportingReads = 2;
			u_int32_t maxDepth = nodePath.size()*2;

			readpath.clear();
			


			vector<u_int32_t> nodePath_errorFree;
			for(u_int32_t nodeIndex : nodePath){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

				if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
				
				//if(!_isBubble[nodeIndex]){
					//if(_graph->_isNodenameRoundabout.find(nodeName) != _graph->_isNodenameRoundabout.end()) continue;
				//}

				nodePath_errorFree.push_back(nodeIndex);
			}

			if(print_read){
				cout << "\tRead error free:" << endl;
				cout << "\t";
				for(u_int32_t nodeIndex : nodePath_errorFree){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
				}
				cout << endl;
			}

			if(nodePath_errorFree.size() == 0) return;






			vector<u_int32_t> nodePath_errorFree_fixed;

			for(long i=0; i<((long)nodePath_errorFree.size())-1; i++){
				
				u_int32_t nodeIndex_source = nodePath_errorFree[i];
				u_int32_t nodeIndex_dest = nodePath_errorFree[i+1];

				nodePath_errorFree_fixed.push_back(nodeIndex_source);

				vector<u_int32_t> path;
				//_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, maxDepth, _unitigDatas, unallowedNodeNames, nodeIndex_source, nodePath_anchor);
				_graph->findUniquePath(nodeIndex_source, nodeIndex_dest, path, false, false, maxDepth);

				if(path.size() == 0){
					continue;
				}
				else{
					//isFix = true;
				}

				std::reverse(path.begin(), path.end());

				for(u_int32_t nodeIndex : path){
					nodePath_errorFree_fixed.push_back(nodeIndex);
				}
				
			}
			nodePath_errorFree_fixed.push_back(nodePath_errorFree[nodePath_errorFree.size()-1]);

			if(print_read){
				cout << "\tRead error free (fixed):" << endl;
				cout << "\t";
				for(u_int32_t nodeIndex : nodePath_errorFree_fixed){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
				}
				cout << endl;
				//if(isFix) getchar();
			}



			u_int32_t nodeIndex_source;
			u_int32_t nodeIndexFinal = -1;

			unordered_set<u_int32_t> isNodenameSolid;
			for(u_int32_t nodeIndex : nodePath_errorFree_fixed){
				isNodenameSolid.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}

			vector<u_int32_t> readpath_nodes;
			size_t i=0;
			bool isFirstIteration = true;

			for(size_t i=0; i<kminmersInfos.size()-1; i++){

				const ReadKminmerComplete& info = kminmersInfos[i];
				const KmerVec& vec = info._vec;

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
				
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				if(isNodenameSolid.find(nodeName) == isNodenameSolid.end()) continue;

				u_int32_t nodeIndex_source = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

				//for(long i=0; i<((long)nodePath_connected.size())-1; i++){
				
				//u_int32_t nodeName_source = nodeName; //nodePath_connected[i];
				u_int32_t nodeIndex_dest = -1;
				for(size_t j=i+1; j<kminmersInfos.size(); j++){
					const ReadKminmerComplete& info = kminmersInfos[j];
					const KmerVec& vec = info._vec;
					if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
					u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
					if(isNodenameSolid.find(nodeName) == isNodenameSolid.end()) continue;
					nodeIndex_dest = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);
					break;
				}


				if(print_read)  cout << BiGraph::nodeIndex_to_nodeName(nodeIndex_source) << " -> " << BiGraph::nodeIndex_to_nodeName(nodeIndex_dest) << endl;

				if(nodeIndex_dest == -1) break; //no node index dest found (typically error at the end of the read)

				if(tryFindPathDirect(nodeIndex_source, nodeIndex_dest, readpath_nodes)){
					

					fillCorrectedArea(nodeIndex_source, nodeIndex_dest, readMinimizers, readpath, kminmersInfos, isFirstIteration, print_read);
			
					//fillCorrectedArea();

					//lastSolidNodeIndex = nodeIndex_source;

					nodeIndex_source = nodeIndex_dest; //readpath_nodes[readpath_nodes.size()-1];
					nodeIndexFinal = nodeIndex_dest;


					isFirstIteration = false;
					//continue;
				}
				else{

					vector<u_int32_t> path;
					_graph->findUniquePath(nodeIndex_source, nodeIndex_dest, path, false, true, maxDepth);

					if(path.size() != 0){
						std::reverse(path.begin(), path.end());
						for(u_int32_t nodeIndex : path){
							fillCorrectedArea(nodeIndex_source, nodeIndex, readMinimizers, readpath, kminmersInfos, isFirstIteration, print_read);
							nodeIndex_source = nodeIndex; //readpath_nodes[readpath_nodes.size()-1];
							nodeIndexFinal = nodeIndex;
						}
						isFirstIteration = false;
					}
					else{
						if(print_read) cout << "failed unique path" << endl;
						correctionFailed = true;
						//lastSolidNodeIndex = -1;
						fillUncorrectedArea(i, nodeIndex_dest, readMinimizers, readpath, kminmersInfos, isFirstIteration, print_read);
						nodeIndexFinal = nodeIndex_dest;
						isFirstIteration = false;
						//getchar();
					}


				}


			}
			
			//cout << "node index final: " << nodeIndexFinal << endl;
			nodePathSolid.clear();
			nodePathSolid.push_back(nodePath_errorFree_fixed[0]);
			if(nodeIndexFinal != -1){
				nodePathSolid.push_back(nodeIndexFinal);
			}


			if(print_read){
				cout << "\tRead corrected (minimizers):" << endl;
				cout << "\t";
				for(u_int64_t minimizer : readpath){
					cout << minimizer << " ";
				}
				cout << endl;



			}
	
			if(correctionFailed){
				//cout << "lala" << endl;
				//getchar();
			}
		}



		bool tryFindPathDirect(u_int32_t nodeIndex_source, u_int32_t nodeIndex_dest, vector<u_int32_t>& readpath){

			u_int32_t unitigIndex_source = _graph->nodeIndex_to_unitigIndex(nodeIndex_source);

			vector<u_int32_t> successors;
			_graph->getSuccessors(nodeIndex_source, 0, successors);

			for(u_int32_t nodeIndexSuccessor: successors){

				u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndexSuccessor);
				if(unitigIndex != unitigIndex_source) continue;

				if(nodeIndexSuccessor == nodeIndex_dest) return true;

			}

			return false;
		}

		void extendReadpath(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, const vector<u_int32_t>& nodePathSolid, vector<u_int64_t>& readpath, bool print_read){


			if(nodePathSolid.size() == 0) return;

			u_int32_t nodeIndex_left = nodePathSolid[0];
			u_int32_t nodeIndex_right = nodePathSolid[nodePathSolid.size()-1];
			u_int32_t nodeNameSolidLeft = BiGraph::nodeIndex_to_nodeName(nodeIndex_left);
			u_int32_t nodeNameSolidRight = BiGraph::nodeIndex_to_nodeName(nodeIndex_right);

			bool needExtendLeft = _graph->isEdgeNode(nodeIndex_left) || _graph->isNeighborOfEdgeNode(nodeIndex_left);// || _isNodeNameUnsupported.find(nodeNameSolidLeft) != _isNodeNameUnsupported.end();
			bool needExtendRight = _graph->isEdgeNode(nodeIndex_right) || _graph->isNeighborOfEdgeNode(nodeIndex_right);// || _isNodeNameUnsupported.find(nodeNameSolidRight) != _isNodeNameUnsupported.end();

			if(!needExtendLeft && !needExtendRight) return;

			if(print_read){
				cout << "\tExtending read path" << endl;
				cout << "\t";
				for(u_int32_t nodeIndex : nodePathSolid){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
				}
				cout << endl;
			}


			u_int64_t solidPositionLeft = 0;
			for(long i=0; i < kminmersInfos.size(); i++){
				const KmerVec& vec = kminmersInfos[i]._vec;

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					continue;
				}
				
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

				if(nodeName == nodeNameSolidLeft){
					
					solidPositionLeft = i;
					break;
				}
			}


			u_int64_t solidPositionRight = 0;
			for(long i=kminmersInfos.size()-1; i>=0; i--){
				const KmerVec& vec = kminmersInfos[i]._vec;

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					continue;
				}
				
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

				if(nodeName == nodeNameSolidRight){
					
					solidPositionRight = i;
					break;
				}
			}


			u_int64_t nbExtendLeft = solidPositionLeft;
			u_int64_t nbExtendRight = readMinimizers.size() - (solidPositionRight + _kminmerSize);
			if(print_read){
				cout << "\tLeft solid position: " << nodeNameSolidLeft << " " << solidPositionLeft << endl;
				cout << "\tRight solid position: " << nodeNameSolidRight << " " << solidPositionRight << endl;
			}

			nbExtendLeft += 2;
			nbExtendRight += 2;
			
			if(print_read) cout << "\tExtend: " << nbExtendLeft << " " << nbExtendRight << endl;

			solidPositionLeft -= 1;
			solidPositionRight += _kminmerSize;

			unordered_set<u_int32_t> unallowedNodeNames;
			for(u_int32_t nodeIndex : nodePathSolid){
				unallowedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}

			if(needExtendRight) extendReadpath2(readpath, nodePathSolid, true, nbExtendRight, unallowedNodeNames, solidPositionRight, readMinimizers, print_read);
			if(needExtendLeft) extendReadpath2(readpath, nodePathSolid, false, nbExtendLeft, unallowedNodeNames, solidPositionLeft, readMinimizers, print_read);
			

		}

		struct SuccessorData3{
			u_int32_t _nodeIndex;
			bool _isSupported;
		};

		void extendReadpath2(vector<u_int64_t>& readpath, const vector<u_int32_t>& nodePath, bool forward, size_t nbSuccessors, unordered_set<u_int32_t>& unallowedNodeNames, u_int64_t extendPosition, const vector<u_int64_t>& readMinimizers, bool print_read){

			if(print_read){
				if(forward){
					cout << "\tExtending right" << endl;
				}
				else{
					cout << "\tExtending left" << endl;
				}
			}

			if(nodePath.size() == 0) return;

			u_int32_t nodeIndex = -1;
			if(forward){
				nodeIndex = nodePath[nodePath.size()-1];
			}
			else{
				nodeIndex = nodePath[0];
			}

			//bool orientation;
			//u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
			//cout << nodeName_start << endl;
			//vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];
			//cout << orientation << endl;
			//for(u_int64_t minimizer : minimizerSeq){
			//	cout << minimizer << " ";
			//}
			//cout << endl;

			bool failed = false;

			for(size_t i=0; i<nbSuccessors; i++){

				vector<u_int32_t> successors;
				if(forward){
					_graph->getSuccessors(nodeIndex, 0, successors);
				}
				else{
					_graph->getPredecessors(nodeIndex, 0, successors);
				}


				if(successors.size() == 0){
					failed = true;
					break;
				}
				else if(successors.size() == 1){
					u_int32_t nodeIndexSuccessor = successors[0];
					if(unallowedNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndexSuccessor)) != unallowedNodeNames.end()){
						failed = true;
						break;
					}
					if(_graph->nodeIndex_to_unitigIndex(nodeIndex) != _graph->nodeIndex_to_unitigIndex(nodeIndexSuccessor)){
						failed = true;
						break;
					}
					
					
					bool orientation;
					u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndexSuccessor, orientation);
					vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];
					if(!forward){
						std::reverse(minimizerSeq.begin(), minimizerSeq.end());
					}

					u_int64_t minimizer = -1;
					if(orientation){
						minimizer = minimizerSeq[minimizerSeq.size()-1];
					}
					else{
						minimizer = minimizerSeq[0];
					}


					//cout << "Add: " << nodeName_start << " " <<  minimizer << endl;
					//for(u_int64_t m : minimizerSeq){
					//	cout << m << " ";
					//}
					//cout << endl;

					if(forward){
						readpath.push_back(minimizer);
						extendPosition += 1;
						//nodePath.push_back(nodeIndexSuccessor);
					}
					else{
						//for(u_int64_t m : minimizerSeq){
						//	cout << m <<  " ";
						//}
						//cout << endl;
						readpath.insert(readpath.begin(), minimizer);
						extendPosition -= 1;
						//nodePath.insert(nodePath.begin(), nodeIndexSuccessor);
					}

					nodeIndex = nodeIndexSuccessor;

				}
				else{
					failed = true;
					break;
					
					//nodeIndex = determineBestSupportedSuccessors(nodePath, forward, unallowedNodeNames, successors);
					//if(nodeIndex == -1) return;

					//if(unallowedNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != unallowedNodeNames.end()) return;

					//if(forward){
					//	nodePath.push_back(nodeIndex);
					//}
					//else{
					//	nodePath.insert(nodePath.begin(), nodeIndex);
					//}
					
				}
				
			}
			//vector<u_int32_t> prevNodes;
			//if()

		}

		u_int32_t determineBestSupportedSuccessors(vector<u_int32_t>& nodePath, bool forward, unordered_set<u_int32_t>& unallowedNodeNames, const vector<u_int32_t>& successors){

			#ifdef PRINT_DEBUG_PREVRANK
				cout << "\t-------------------------" << endl;
			#endif
			
			vector<SuccessorData3> possibleSuccessors;
			for(u_int32_t nodeIndex_successor : successors){
				possibleSuccessors.push_back({nodeIndex_successor, true});
			}

			size_t prevRank = 0;

			while(true){


				long prevRankIndex = -1; 

				if(forward){
					prevRankIndex = nodePath.size() - prevRank - 1;
					if(prevRankIndex < 0) return -1;
				}
				else{
					prevRankIndex = prevRank;
					if(prevRankIndex >= nodePath.size()) return -1;
				}

				//cout << "\t" << prevRankIndex << endl;
				u_int32_t nodeIndex_prev = nodePath[prevRankIndex];
				u_int32_t nodeName_prev = BiGraph::nodeIndex_to_nodeName(nodeIndex_prev);

				#ifdef PRINT_DEBUG_PREVRANK
					cout << "\t" << prevRank << " " << nodeName_prev << "    ";
				#endif

				bool isFinished = false;
				for(SuccessorData3& successorData : possibleSuccessors){
					if(!successorData._isSupported) continue;

					u_int32_t nodeName_successor = BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex);

					u_int64_t nbShared = Utils::computeSharedReads(_unitigDatas[nodeName_successor], _unitigDatas[nodeName_prev]);
					
					#ifdef PRINT_DEBUG_PREVRANK
						cout << nodeName_successor << " " << nbShared << "    ";
					#endif

					if(nbShared <= 2){
						successorData._isSupported = false;
					}

				}

				#ifdef PRINT_DEBUG_PREVRANK
					cout << endl;
				#endif

				u_int32_t supporteNodeIndex = -1;
				size_t nbSupportedSuccessors = 0;
				for(SuccessorData3& successorData : possibleSuccessors){
					if(successorData._isSupported){
						nbSupportedSuccessors += 1;
						supporteNodeIndex = successorData._nodeIndex;
					}
				}

				if(nbSupportedSuccessors == 1){
					return supporteNodeIndex;
				}

				prevRank += 1;

				if(prevRank > 10) return -1;
			}

			return -1;
		}
		


	};

};





#endif

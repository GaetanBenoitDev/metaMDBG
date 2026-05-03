


/*

- toujours des probleme de resolution de la dernier bubble dans un contig circulaire
- condition de repetition du repeat solver pas forcement bonne: if(nbUnitigs == _unitigGraph->_unitigs.size()) break;
	- faudrait juste litteralement compter le nb de changement pour etre sur qu'on a fait une modif

- reprise: dans la zymo on a un cas bizarre, une repeaty qui devrait etre compacté (surement des problemes de edges)

- overlap_test_201: on arrive pas a resoudre une perfect repeat, revoir le code pour ça: utiliser des doublet (solve edges)

"reprise: enlever les contig deja considerer comme erreur (derep)"

cout << "Repeat solving: attention aux repeats qui passe plusieurs fois dans le meme cycle (il faut bien utiliser la sequence du contig directement dans ce cas)" << endl;
*/

#ifndef MDBG_METAG_GENERATECONTIGGRAPH
#define MDBG_METAG_GENERATECONTIGGRAPH

#include "../Commons.hpp"
#include "graph/Graph.hpp"
//#include "graph/ProgressiveAbundanceFilter.hpp"
//#include "graph/RepeatSolver.hpp"


class GenerateContigGraph{
    
public:

	Parameters _params;
	string _inputDir;
	int _nbCores;

	struct ContigPosition{
		UnitigType _contigIndex;
		u_int32_t _position;
	};

	struct ReferencePosition{
		string _referenceName;
		u_int32_t _position;
	};

    UnitigGraph2* _unitigGraph;
	vector<vector<MinimizerType>> _unitigName_to_minimizers;
	ankerl::unordered_dense::map<u_int128_t, UnitigType> _kminmer_to_unitigIndex;
	ankerl::unordered_dense::map<UnitigType, vector<ContigPosition>> _unitigName_to_contigIndexes;
	//ankerl::unordered_dense::set<UnitigType> _validContigIndexes;
	ankerl::unordered_dense::set<KmerVec> _solidDoublets;
	ankerl::unordered_dense::set<KmerVec> _solidTriplets;
	ankerl::unordered_dense::set<u_int128_t> _repeatedKminmers;
	ankerl::unordered_dense::set<u_int128_t> _repeatedUnitigNames;

	unordered_set<MinimizerType> _isRepetitiveMinimizers;
	ankerl::unordered_dense::map<UnitigType, vector<ReferencePosition>> _unitigName_to_referenceName;
	ankerl::unordered_dense::map<UnitigType, vector<ReferencePosition>> _unitigName_to_referenceName_init;


	//ofstream _contigGraphFile;

	GenerateContigGraph(const string& inputDir, const int nbCores, const Parameters& params, UnitigGraph2* unitigGraph){
		cout << inputDir << endl;
		_inputDir = inputDir;
		_nbCores = nbCores;
		_params = params;
		_unitigGraph = unitigGraph;
	}
	/*
	GenerateContigGraph(): Tool (){

	}
	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("toBasespace", ""); 
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::ValueFlag<int> arg_maxBubbleLength(parser, "", "Max bubble popping length", {ARG_MAX_BUBBLE_LENGTH}, 50000);
		args::ValueFlag<int> arg_maxTipLength(parser, "", "Max tip clipping length", {ARG_MAX_TIP_LENGTH}, 50000);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);

		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (const std::exception& e)
		{
			cerr << parser;
			cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			cerr << parser;
			exit(0);
		}

		_inputDir = args::get(arg_outputDir);
		_nbCores = args::get(arg_nbCores);

		_maxBubbleLength = args::get(arg_maxBubbleLength);
		_maxTipLength = args::get(arg_maxTipLength);

		_contigGraphDir = _inputDir + "/contigGraph/";
		_params.load(_contigGraphDir + "/parameters.gz");

		openLogFile(_inputDir);
	
	}
	*/

	void setup(){


		//cout << "Indexing invalid contig pas actif et _validContigIndexes desactivé dans mapContigFunctor" << endl;
		//cout << "Indexing valid contig (non erroneous ones), to remove in final version" << endl;
		//indexValidContigs(_inputDir+ "/contigs.fasta.gz");
		
		//_kminmerSize = 21;
		//_params._kminmerSize = _kminmerSize;

		cout << "Loading unitig minimizers" << endl;
		loadUnitigSequences();


		_isRepetitiveMinimizers = Commons::loadRepetitiveMinimizers(_inputDir + "/repetitiveMinimizers.bin", 0);
		saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph_init.gfa", false);

		clearPass();
		indexUnitigs();
		mapReferences();
		_unitigName_to_referenceName_init = _unitigName_to_referenceName;
		
		//getchar();
	}

	/*
    bool execute (){
		
		//_contigGraphFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/contigGraph.tsv");



		//cout << "Loading unitig graph" << endl;
		//UnitigGraph2* unitigGraph2 = new UnitigGraph2(_params._kminmerSize, _params._kminmerLengthMean, _params._kminmerOverlapMean, _params._kminmerLengthMean-_params._kminmerOverlapMean, _params._minimizerSpacingMean, _nbCores);
		//unitigGraph2->load(_contigGraphDir + "/unitigGraph.stats.bin", _contigGraphDir + "/unitigGraph.nodes.bin", _contigGraphDir + "/unitigGraph.nodes.abundances.bin", _contigGraphDir + "/unitigGraph.edges.successors.bin");



		//clearPass();
		//indexUnitigs();
		//mapContigs();


		bool isModification = false;
		int nbPasses = 0;

		while(true){
			
			cout << "\tPass: " << nbPasses << endl;
			u_int64_t nbUnitigs = _unitigGraph->_unitigs.size();



			
			cout << "Cleaning chimeric unitigs" << endl;
			clearPass();
			indexUnitigs();
			mapContigs();
			cleanChimericUnitigs();
			recompactGraph();

			//saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", true);
			//getchar();
			
			cout << "Solving perfect edge repeats" << endl;
			clearPass();
			indexUnitigs();
			mapContigs();
			solvePerfectEdgeRepeats();
			recompactGraph();
			
			cout << "Solving repeats" << endl;
			clearPass();
			indexUnitigs();
			mapContigs();
			solveRepeats();
			recompactGraph();
			
			//saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", true);
			//getchar();




			if(nbUnitigs == _unitigGraph->_unitigs.size()) break;

			nbPasses += 1;
			isModification = true;
		}
		

		//saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", true);
		//getchar();

		return isModification;
	}
	*/

	bool executeSolveRepeats(){

		bool isModification = false;
		int pass = 0;

		while(true){
			cout << "Solving perfect edge repeats " << pass  << endl;
			
			//saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", true);
			//getchar();

			clearPass();
			indexUnitigs();
			mapContigs();
			bool isMod1 = solvePerfectEdgeRepeats();
			recompactGraph();
			
			cout << "Solving repeats" << endl;
			clearPass();
			indexUnitigs();
			mapContigs();
			bool isMod2 = solveRepeats();
			recompactGraph();

			if(isMod1 || isMod2){
				isModification = true;
			}
			else{
				break;
			}

			pass += 1;
		}

		saveGraph("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", true);
		getchar();

		return isModification;
	}

	void saveGraph(const string& filename, const bool doesMapReference){
		cout << "Saving graph" << endl;
		recompactGraph();
		clearPass();
		indexUnitigs();
		mapContigs();
		if(doesMapReference) mapReferences();
		_unitigGraph->save(filename, _inputDir);
		dumpContigNames();


		if(doesMapReference) dumpReferenceNames();
		cout << "done" << endl;
	}

	void recompactGraph(){

		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			_unitigGraph->recompact(node);
		}
	}

	/*
	void indexValidContigs(const string& contigFilename){

		auto fp = std::bind(&GenerateContigGraph::indexValidContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void indexValidContigs_read(const Read& read){

		string contigName = Utils::shortenHeader(read._header);
		contigName.erase(0, 3); //remove "ctg"
		u_int32_t contigIndex = stoull(contigName);

		_validContigIndexes.insert(contigIndex);
	}
	*/

	void loadUnitigSequences(){

		ifstream nodeFile(_inputDir + "/unitigGraph.nodes.bin");

		while(true){

			u_int32_t size;
			nodeFile.read((char*)&size, sizeof(size));
			
			if(nodeFile.eof()) break;

			vector<MinimizerType> minimizers;
			minimizers.resize(size);
			nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));

			UnitigType unitigIndex;
			nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			const UnitigType unitigName = unitigIndex/2;

			while(_unitigName_to_minimizers.size() <= unitigName){
				_unitigName_to_minimizers.push_back({});
			}

			_unitigName_to_minimizers[unitigName] = minimizers;
		}

		nodeFile.close();
	}

	void clearPass(){

		ankerl::unordered_dense::map<UnitigType, vector<ContigPosition>>().swap(_unitigName_to_contigIndexes);
		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmer_to_unitigIndex);
		ankerl::unordered_dense::set<KmerVec>().swap(_solidTriplets);
		ankerl::unordered_dense::set<KmerVec>().swap(_solidDoublets);
	}

	void indexUnitigs(){



        for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

            UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
            if(node == nullptr) continue;

            vector<UnitigType> unitigs;
            if(node->_unitigMerge.size() == 0){
                unitigs.push_back(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));    
            }
            else{
                unitigs = node->_unitigMerge;
            }

			if(node->_isReversed){
				vector<UnitigType> unitigsRC;
				UnitigGraph2::reverseComplementUnitigs(unitigs, unitigsRC);

				unitigs = unitigsRC;
			}

            //u_int32_t size = unitigs.size();
            //u_int8_t isCircular = _unitigGraph2->isCircular(node);
            //float abundance = node->_abundance;
            //u_int32_t nbMinimizers = node->_nbMinimizers;

			vector<MinimizerType> minimizers;
			unitigsToMinimizers(unitigs, minimizers);

	
			//indexUnitig(minimizers, node->_unitigName);
			
			indexUnitig(minimizers, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));

			std::reverse(minimizers.begin(), minimizers.end());
			indexUnitig(minimizers, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true));
			

        }


	}

	void indexUnitig(const vector<MinimizerType>& minimizers, const UnitigType unitigIndex){

		const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _params._kminmerSize);

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t vecHashNorm = kminmer.normalize().hash128();
			if(_repeatedKminmers.find(vecHashNorm) != _repeatedKminmers.end()) continue; //Do not index kminmers from solved repeats

			const u_int128_t vecHash = kminmer.hash128();

			_kminmer_to_unitigIndex[vecHash] = unitigIndex;
			//_kminmer_to_unitigIndex[vecHashNorm] = unitigIndex;
		}

	}

	void unitigsToMinimizers(const vector<UnitigType>& unitigs, vector<MinimizerType>& outputMinimizers){

		outputMinimizers.clear();
		//vector<MinimizerType> prevMinimizers;

		for(size_t i=0; i<unitigs.size(); i++){

			UnitigType unitigIndex = unitigs[i];


			bool isReversed;
			const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

			//cout << "\t" << unitigName << " " << isReversed << endl;
			//cout << "\t" << unitigName << " " << _parent._unitigName_to_minimizers.size() << " " << endl;
			//cout << "\t" << _parent._unitigName_to_minimizers[unitigName].size() << " " << endl;

			vector<MinimizerType> minimizers = _unitigName_to_minimizers[unitigName];

			if(isReversed){
				std::reverse(minimizers.begin(), minimizers.end());
			}

			if(i==0){

				outputMinimizers = minimizers;
			}
			else{

				int overlapSize = _params._kminmerSize - 1;
				
				outputMinimizers.insert(outputMinimizers.end(), minimizers.begin()+overlapSize, minimizers.end());

			}

			//prevMinimizers = minimizers;
		}

	}

	/*
	class FilterFunctor {

		public:

		GenerateContigGraph& _parent;
		bool _hasFilteredDummy;

		FilterFunctor(GenerateContigGraph& parent) : _parent(parent){
			_hasFilteredDummy = false;
		}

		//FilterFunctor(const FilterFunctor& copy) : _parent(copy._parent){
		//}

		~FilterFunctor(){
		}

		void operator () (const float cutoff) {

			cout << "Functor: " << cutoff << endl;

			if(cutoff < 6){
				cout << "Skip: " << cutoff << endl;
				return;
			}

			_parent._minUnitigAbundance = cutoff / 0.5;
			_parent.indexUnitigs();



			KminmerParserParallel parser1(_parent._inputDir + "/contig_data_final.bin", 0, _parent._params._kminmerSize, false, false, _parent._nbCores);
			parser1.parse(MapContigFunctor(_parent));

			


			//_parent._progressiveAbundanceFilter->_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", _parent._inputDir);

			
			if(!_hasFilteredDummy){
				
				cout << "Faudra faire ça une seule fois juste apres avoir load le graphe" << endl;
				
				_hasFilteredDummy = true;

				for(size_t i=0; i<_parent._progressiveAbundanceFilter->_unitigGraph2->_unitigs.size(); i++){

					UnitigGraph2::UnitigNode* node = _parent._progressiveAbundanceFilter->_unitigGraph2->_unitigs[i];
					if(node == nullptr) continue;

					if(_parent._unitigName_to_contigIndexes.find(node->_unitigName) != _parent._unitigName_to_contigIndexes.end()) continue;

					_parent._progressiveAbundanceFilter->_unitigGraph2->removeNode(node);

				}

				for(UnitigGraph2::UnitigNode* node : _parent._progressiveAbundanceFilter->_unitigGraph2->_unitigs){

					if(node == nullptr) continue;

					_parent._progressiveAbundanceFilter->_unitigGraph2->recompact(node);
				}
				
				cout << "check" << endl;



				return;
			}

			

			//KminmerParserParallel parser2(_parent._inputDir + "/contig_data_final.bin", 0, _parent._params._kminmerSize, false, false, _parent._nbCores);
			//parser2.parse(ContigFunctor(_parent));



			ankerl::unordered_dense::map<UnitigType, vector<UnitigType>>().swap(_parent._unitigIndex_to_contigIndexStarts);
			ankerl::unordered_dense::map<UnitigType, vector<UnitigType>>().swap(_parent._unitigIndex_to_contigIndexEnds);
			ankerl::unordered_dense::map<UnitigType, pair<UnitigType,UnitigType>>().swap(_parent._contigIndex_to_endUnitigIndex);



		}


	};
	
	*/


	void mapReferences(){
		
		ReadParserParallel readParser("/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG_1.4/metaMDBG/build/myloasm_ctg1.fasta", true, false, _nbCores);
		//ReadParserParallel readParser("/pasteur/appa/scratch/gbenoit/data/genomes/genomes/201/GCF_000816365.1_ASM81636v1_genomic.fna", true, false, _nbCores);
		//ReadParserParallel readParser("/pasteur/appa/homes/gbenoit/appa/data/nanopore/mock/Zymo_D6331/references.fasta", true, false, _nbCores);
		readParser.parse(MapReferenceFunctor(*this));
	}


	class MapReferenceFunctor {

		public:

		GenerateContigGraph& _parent;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		MapReferenceFunctor(GenerateContigGraph& parent) : _parent(parent){
			_minimizerParser = new MinimizerParser(_parent._params._minimizerSize, _parent._params._minimizerDensity_assembly, _parent._isRepetitiveMinimizers);
		}

		MapReferenceFunctor(const MapReferenceFunctor& copy) : _parent(copy._parent){
			_minimizerParser = new MinimizerParser(_parent._params._minimizerSize, _parent._params._minimizerDensity_assembly, _parent._isRepetitiveMinimizers);
		}

		~MapReferenceFunctor(){
			delete _minimizerParser;
		}
		


		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			if(readIndex % 100000 == 0) Logger::get().debug() << readIndex;


			string seq = read._seq; //.substr(0, 100);

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(seq.c_str(), seq.size(), rleSequence, rlePositions, _parent._params._useHomopolymerCompression);

			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos;
			vector<u_int8_t> minimizerDirections;
			_minimizerParser->parse(rleSequence, minimizers, minimizerPos, minimizerDirections);

			cout << read._header << " " << read._seq.size() << " " << minimizers.size() << " " << _parent._params._kminmerSizeFirst << " " << _parent._params._kminmerSize << endl;

			//minimizers = Commons::purgePalindrome(minimizers, _parent._params._kminmerSizeFirst, _parent._params._kminmerSize);


			
			vector<UnitigType> unitigIndexes;
			
			UnitigType prevUnitigIndex = -1;
			//UnitigType firstUnitigName = -1;
			//UnitigType lastUnitigName = -1;
			
			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _parent._params._kminmerSize);

			for(size_t i=0; i<kminmers.size(); i++){

				const KmerVec& kminmer = kminmers[i];
				const u_int128_t vecHash = kminmer.hash128();

				if(_parent._kminmer_to_unitigIndex.find(vecHash) == _parent._kminmer_to_unitigIndex.end()) continue;

				const UnitigType unitigIndex = _parent._kminmer_to_unitigIndex[vecHash];
				
				//if(contigIndex == 4798) cout << i << "\t" << unitigIndex << endl;
				//const UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(unitigName, isReversed);

				if(unitigIndex != prevUnitigIndex){
					//if(contigIndex == 4798) cout << "utg" << unitigName << "\t" << isReversed << endl;
					prevUnitigIndex = unitigIndex;
					unitigIndexes.push_back(unitigIndex);
					
					//if(readIndex == 18365) cout << "utg" << unitigName << endl;
				}
				//if(firstUnitigName == -1) firstUnitigName = unitigName;
				//lastUnitigName = unitigName;
			}

			
			//if(readIndex == 18365) cout << unitigNames.size() << endl;
			if(unitigIndexes.size() == 0) return;

			const string referenceName = Utils::shortenHeader(read._header);

			indexUnitigs(referenceName, unitigIndexes);

			vector<UnitigType> unitigIndexesRC;
			UnitigGraph2::reverseComplementUnitigs(unitigIndexes, unitigIndexesRC);

			indexUnitigs(referenceName, unitigIndexesRC);

		}

		void indexUnitigs(const string& referenceName, const vector<UnitigType>& unitigIndexes){


			#pragma omp critical(indexUnitigs)
			{
				
				for(size_t i=0; i<unitigIndexes.size(); i++){

					const UnitigType unitigIndex = unitigIndexes[i];
					const u_int32_t contigPosition = i;

					_parent._unitigName_to_referenceName[unitigIndex].push_back({referenceName, contigPosition});
				}

			}

		}

	};

	void mapContigs(){
		
		KminmerParserParallel parser(_inputDir + "/contig_data_final.bin", 0, 0, false, false, _nbCores);
		parser.parseSequences(MapContigFunctor(*this));
	}

	class MapContigFunctor {

		public:

		GenerateContigGraph& _parent;

		MapContigFunctor(GenerateContigGraph& parent) : _parent(parent){
		}

		MapContigFunctor(const MapContigFunctor& copy) : _parent(copy._parent){
		}

		~MapContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			//"reprise: inverse indexUnitig et mapContig, les contigs ne changent jamais donc indexer les contigs, puis traverser le graphe"
			const UnitigType contigIndex = kminmerList._readIndex;

			//if(_parent._validContigIndexes.find(contigIndex) == _parent._validContigIndexes.end()) return;
			//if(_parent._isContigProcessed.find(readIndex) != _parent._isContigProcessed.end()) return;

			const vector<MinimizerType>& minimizers = kminmerList._readMinimizers;

			if(minimizers.size() == 0) return; //Invalid contig index (erroneous etc)

			//const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			//u_int8_t isCircular = kminmerList._isCircular;

			vector<UnitigType> unitigIndexes;
			
			UnitigType prevUnitigIndex = -1;
			//UnitigType firstUnitigName = -1;
			//UnitigType lastUnitigName = -1;
			
			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _parent._params._kminmerSize);

			for(size_t i=0; i<kminmers.size(); i++){

				const KmerVec& kminmer = kminmers[i];
				const u_int128_t vecHash = kminmer.hash128();

				if(_parent._kminmer_to_unitigIndex.find(vecHash) == _parent._kminmer_to_unitigIndex.end()) continue;

				const UnitigType unitigIndex = _parent._kminmer_to_unitigIndex[vecHash];
				
				//if(contigIndex == 4798) cout << i << "\t" << unitigIndex << endl;
				//const UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(unitigName, isReversed);

				if(unitigIndex != prevUnitigIndex){
					//if(contigIndex == 4798) cout << "utg" << unitigName << "\t" << isReversed << endl;
					prevUnitigIndex = unitigIndex;
					unitigIndexes.push_back(unitigIndex);
					
					//if(readIndex == 18365) cout << "utg" << unitigName << endl;
				}
				//if(firstUnitigName == -1) firstUnitigName = unitigName;
				//lastUnitigName = unitigName;


			}

			
			//if(readIndex == 18365) cout << unitigNames.size() << endl;
			if(unitigIndexes.size() == 0) return;


			indexUnitigs(contigIndex, unitigIndexes);

			vector<UnitigType> unitigIndexesRC;
			UnitigGraph2::reverseComplementUnitigs(unitigIndexes, unitigIndexesRC);

			indexUnitigs(contigIndex, unitigIndexesRC);
		}
		
		void indexUnitigs(const UnitigType& contigIndex, const vector<UnitigType>& unitigIndexes){


			vector<KmerVec> doublets = unitigNamesToTriplets(unitigIndexes, 2);
			vector<KmerVec> triplets = unitigNamesToTriplets(unitigIndexes, 3);

			#pragma omp critical(indexUnitigs)
			{
				
				for(size_t i=0; i<unitigIndexes.size(); i++){

					const UnitigType unitigIndex = unitigIndexes[i];
					const u_int32_t contigPosition = i;

					//"reprise: je pense qu'il y a un pb de contigPosition en fonction de l'orientation du contig"
					//if(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) == 37566 || UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) == 38983){
						
					//	cout <<  "utg" << unitigIndex << "\tctg" << contigIndex << "\t" << contigPosition << endl;
					//}

					//if(contigIndex == 4798){
					//	cout << i << "\tutg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << "\t" << unitigIndex << endl;
					//}

					if(_parent._unitigName_to_contigIndexes.find(unitigIndex) == _parent._unitigName_to_contigIndexes.end()){
						_parent._unitigName_to_contigIndexes[unitigIndex].push_back({contigIndex, contigPosition});
					}
					else{
							
						_parent._unitigName_to_contigIndexes[unitigIndex].push_back({contigIndex, contigPosition});
						//auto& it = _parent._unitigName_to_contigIndexes[unitigName];
						//if(std::find(it.begin(), it.end(), readIndex) == it.end()){
						//	it.push_back(readIndex);
						//}
					}
				}

				for(const KmerVec& vec : triplets){
					_parent._solidTriplets.insert(vec);

					//if(contigIndex == 4798) cout << vec.toString() << endl;
				}

				for(const KmerVec& vec : doublets){
					_parent._solidDoublets.insert(vec);

					//if(contigIndex == 4798) cout << vec.toString() << endl;
				}

			}

		}


		vector<KmerVec> unitigNamesToTriplets(const vector<UnitigType>& unitigNames, const size_t& kminmerSize){
			
			vector<KmerVec> kminmers;

			if(unitigNames.size() < kminmerSize) return kminmers;

			for(size_t i=0; i<unitigNames.size()-kminmerSize+1; i++){
				KmerVec vec;
				for(size_t j=0; j<kminmerSize; j++){
					vec._kmers.push_back(unitigNames[i+j]);
				}
				
				kminmers.push_back(vec);
			}

			return kminmers;
			
		}

		/*
		UnitigType determineUnitigIndexStart(const vector<UnitigType>& unitigNames){

			UnitigType unitigIndexStart = UnitigGraph2::unitigName_to_unitigIndex(unitigNames[0], false);
			UnitigType unitigIndexStartRev = UnitigGraph2::unitigName_to_unitigIndex(unitigNames[0], true);

			vector<UnitigType> successors;
			_parent._progressiveAbundanceFilter->_unitigGraph2->getSuccessors(unitigIndexStart, successors);

			for(const auto& unitigIndexSucc : successors){

				if(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) == unitigNames[1]){
					return unitigIndexStart;
				}
			}

			return unitigIndexStartRev;

		}

		UnitigType determineUnitigIndexEnd(const vector<UnitigType>& unitigNames){

			UnitigType unitigIndexEnd = UnitigGraph2::unitigName_to_unitigIndex(unitigNames[unitigNames.size()-1], false);
			UnitigType unitigIndexEndRev = UnitigGraph2::unitigName_to_unitigIndex(unitigNames[unitigNames.size()-1], true);

			vector<UnitigType> predecessors;
			_parent._progressiveAbundanceFilter->_unitigGraph2->getPredecessors(unitigIndexEnd, predecessors);

			for(const auto& unitigIndexPred : predecessors){

				if(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred) == unitigNames[unitigNames.size()-2]){
					return unitigIndexEnd;
				}
			}

			return unitigIndexEndRev;

		}
		*/
	};

	/*
	class ContigFunctor {

		public:

		GenerateContigGraph& _parent;

		ContigFunctor(GenerateContigGraph& parent) : _parent(parent){
		}

		ContigFunctor(const ContigFunctor& copy) : _parent(copy._parent){
		}

		~ContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;

			bool isProcessed = false;
			#pragma omp critical(ContigFunctor)
			{
				if(_parent._isContigProcessed.find(readIndex) != _parent._isContigProcessed.end()) isProcessed = true;
			}
			if(isProcessed) return;

			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			u_int8_t isCircular = kminmerList._isCircular;

			UnitigType lastUnitigIndex = -1;
			vector<u_int32_t> abundances;
			
			u_int32_t minAbundance = -1;
			float sum = 0;
			float n = 0;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const KmerVec& vec = kminmersInfos[i]._vec;
				const u_int128_t& vecHash = vec.hash128();

				if(_parent._kminmer_to_unitigName.find(vecHash) == _parent._kminmer_to_unitigName.end()) continue;

				const UnitigType unitigName = _parent._kminmer_to_unitigName[vecHash];

				u_int32_t abundance = _parent._progressiveAbundanceFilter->_unitigGraph2->_unitigs[unitigName]->_abundance;
				//cout << i << "\t" << abundance << endl;

				abundances.push_back(abundance);

				sum += abundance;
				n += 1;

				if(abundance < minAbundance){
					minAbundance = abundance;
				}
				
			
				//if(readIndex == 8234){
				//	cout << abundance << endl;
				//}
			}



			//float median = Utils::compute_median(abundances);
			//cout << "\t" << median << endl;

			//if(readIndex == 8234){
			//	cout << median << endl;
			//	cout << sum / n << endl;
			//	getchar();
			//}

			if(minAbundance < _parent._minUnitigAbundance){
				
				#pragma omp critical(ContigFunctor)
				{
					cout << "\tProcess: " << readIndex << endl;
					_parent._isContigProcessed.insert(readIndex);
					//cout << "attention _isContigProcessed pas protégé" << endl;
					//_parent.createColorFile();
					//getchar();
					processContig(readIndex);
				}
			
			}

		}

		struct SuccessorQueue{
			UnitigType _unitigIndex;
			u_int32_t _nbMinimizers;
		};

		void processContig(const UnitigType contigIndex){
			
			if(_parent._contigIndex_to_endUnitigIndex.find(contigIndex) == _parent._contigIndex_to_endUnitigIndex.end()) return;

			processContigSub(contigIndex, _parent._contigIndex_to_endUnitigIndex[contigIndex].first);
			processContigSub(contigIndex, _parent._contigIndex_to_endUnitigIndex[contigIndex].second);
		}

		void processContigSub(const UnitigType contigIndex, const UnitigType unitigIndexSource){
			if(contigIndex == 3601){
				_parent._progressiveAbundanceFilter->_unitigGraph2->save("/pasteur/appa/homes/gbenoit/zeus/tmp/assembly_graph.gfa", _parent._inputDir);
				cout << "done" << endl;
				getchar();
			}
			
			//const UnitigType unitigIndexSource = _parent._contigIndex_to_endUnitigIndex[contigIndex];
			
			_parent._contigGraphFile << "ctg" << contigIndex << "\tutg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSource) << "\n";

			std::queue<SuccessorQueue> queue;
			queue.push({unitigIndexSource, 0});
			unordered_set<UnitigType> isVisited;

			//if(contigIndex == 13237) cout << "go" << endl;

			//std::sort(currentPartition._readsToLoad.begin(), currentPartition._readsToLoad.end());
			while (!queue.empty()) {


				const UnitigType unitigIndex = queue.front()._unitigIndex;
				const u_int64_t nbMinimizers = queue.front()._nbMinimizers;
				
				queue.pop();

				u_int64_t distance = nbMinimizers * _parent._params._minimizerSpacingMean;

				if(distance > 100000) continue;
				if(isVisited.find(unitigIndex) != isVisited.end()) continue;

				isVisited.insert(unitigIndex);

				vector<UnitigType> successors;
				_parent._progressiveAbundanceFilter->_unitigGraph2->getSuccessors(unitigIndex, successors);

				//if(contigIndex == 13237){
				//	cout << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << " " << successors.size() << endl;
				//}

				for(const UnitigType unitigIndexSucc : successors){

					const UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc);
					UnitigGraph2::UnitigNode* node = _parent._progressiveAbundanceFilter->_unitigGraph2->_unitigs[unitigNameSucc];
					u_int32_t nbMinimizersSucc = node->_nbMinimizers - (_parent._params._kminmerSize-1);
						
					//if(contigIndex == 13237){
					//	cout << "utg" << unitigNameSucc << endl;
					//}

					queue.push({unitigIndexSucc, nbMinimizers+nbMinimizersSucc});

					//if(node->_unitigMerge.size() == 0){
					if(_parent._unitigIndex_to_contigIndexStarts.find(unitigIndexSucc) != _parent._unitigIndex_to_contigIndexStarts.end()){

						for(const UnitigType contigIndexStart : _parent._unitigIndex_to_contigIndexStarts[unitigIndexSucc]){
							_parent._contigGraphFile << "\tctg" << contigIndexStart << "\tutg" << unitigNameSucc << "\t" << distance << "\n";
						}
					}
					//}
					//else{
					//	for(const UnitigType unitigIndexSucc2 : node->_unitigMerge){
					//		if(_parent._unitigIndex_to_contigIndexStarts.find(unitigIndexSucc2) != _parent._unitigIndex_to_contigIndexStarts.end()){
					//			
					//			for(const UnitigType contigIndexStart : _parent._unitigIndex_to_contigIndexStarts[unitigIndexSucc2]){
					//				_parent._contigGraphFile << "\tctg" << contigIndexStart << "\tutg" << unitigNameSucc << "\t" << distance << "\n";
					//			}
					//		}
					//
					//	}
					//}
				}
				
			}

			//if(contigIndex == 13237){
			//	cout << "gogo" << endl;
			//	getchar();
			//}
		}

	};

	void createColorFile(){


		unordered_map<UnitigType, string> unitigName_to_text;

		for(const auto& it : _unitigIndex_to_contigIndexStarts){

			const UnitigType& unitigIndex = it.first;
			const UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
			const vector<UnitigType>& contigIndexes = it.second;

			string s = "";

			for(const UnitigType contigIndex : contigIndexes){
				s += "S-ctg" + to_string(contigIndex) + " ";
			}

			unitigName_to_text[unitigName] += s;
			//colorFile << "utg" << unitigName << "," << s << endl;
		}

		for(const auto& it : _unitigIndex_to_contigIndexEnds){

			const UnitigType& unitigIndex = it.first;
			const UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
			const vector<UnitigType>& contigIndexes = it.second;

			string s = "";

			for(const UnitigType contigIndex : contigIndexes){
				s += "E-ctg" + to_string(contigIndex) + " ";
			}
			
			unitigName_to_text[unitigName] += s;
			//colorFile << "utg" << unitigName << "," << s << endl;
		}



		ofstream colorFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/contigMap.csv");
		colorFile << "Name,ContigIndex" << endl;

		for(const auto& it : unitigName_to_text){
			colorFile << "utg" << it.first << "," << it.second << endl;
		}

		colorFile.close();

	}
	*/



	bool cleanChimericUnitigs(){

		unordered_set<UnitigType> isUnitigIndexChimeric;

		bool isModification = false;
		clearPass();
		indexUnitigs();
		mapContigs();


		ofstream colorFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/chimericUnitigs.csv");
		colorFile << "Name,Color" << endl;

		vector<UnitigGraph2::UnitigNode*> nodeToRemove;


		#pragma omp parallel for num_threads(_nbCores)
		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){
			
			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			const UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
			const UnitigType unitigIndexRev = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

			if(_unitigName_to_contigIndexes.find(unitigIndex) != _unitigName_to_contigIndexes.end()) continue;

			if(isChimeric(unitigIndex)){
				#pragma omp critical(cleanChimericUnitigs)
				{
					isUnitigIndexChimeric.insert(unitigIndex);
				}
			}

			if(isChimeric(unitigIndexRev)){
				#pragma omp critical(cleanChimericUnitigs)
				{
					isUnitigIndexChimeric.insert(unitigIndexRev);
				}
			}
		}

		//cout << "lul: " << isUnitigIndexChimeric.size() << endl;

		#pragma omp parallel for num_threads(_nbCores)
		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			const UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
			const UnitigType unitigIndexRev = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

			if(_unitigName_to_contigIndexes.find(unitigIndex) != _unitigName_to_contigIndexes.end()) continue;

			if(isChimeric(unitigIndex) && isChimeric(unitigIndexRev)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
			else if(isChimeric(unitigIndex) && !hasContigSuccessorBFS(unitigIndex, isUnitigIndexChimeric)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
			else if(isChimeric(unitigIndexRev) && !hasContigSuccessorBFS(unitigIndexRev, isUnitigIndexChimeric)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
		}

		/*
		#pragma omp parallel for num_threads(_nbCores)
		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			const UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
			const UnitigType unitigIndexRev = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

			if(_unitigName_to_contigIndexes.find(unitigIndex) != _unitigName_to_contigIndexes.end()) continue;

			if(isChimeric(unitigIndex) && isChimeric(unitigIndexRev)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
			else if(isChimeric(unitigIndex) && !hasContigSuccessor(unitigIndex)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
			else if(isChimeric(unitigIndexRev) && !hasContigSuccessor(unitigIndexRev)){
				//colorFile << "utg" << node->_unitigName << "," << "red" << endl;

				#pragma omp critical(cleanChimericUnitigs)
				{
					nodeToRemove.push_back(node);
				}
			}
			//else if(){
			//	colorFile << "utg" << node->_unitigName << "," << "red" << endl;
			//}
			//_parent._progressiveAbundanceFilter->_unitigGraph2->removeNode(node);

		}
		*/
		//cout << nodeToRemove.size() << endl;

		for(UnitigGraph2::UnitigNode* node : nodeToRemove){

			if(node == nullptr) continue;

			colorFile << "utg" << node->_unitigName << "," << "red" << "\n";
			
			_unitigGraph->removeNode(node);
			isModification = true;
		}
		



		colorFile.close();

		
		recompactGraph();

		//cout << "check" << endl;
		//getchar();
		return isModification;
	}

	bool isChimeric(const UnitigType unitigIndex){


		vector<UnitigType> predecessors;
		_unitigGraph->getPredecessors(unitigIndex, predecessors);

		for(const UnitigType unitigIndexPred : predecessors){

			//const UnitigType unitigNamePred = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred);
			//UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[unitigNamePred];

			if(_unitigName_to_contigIndexes.find(unitigIndexPred) == _unitigName_to_contigIndexes.end()) continue;

			unordered_set<UnitigType> contigIndexPreds;
			for(const ContigPosition& contigPosition : _unitigName_to_contigIndexes[unitigIndexPred]){
				contigIndexPreds.insert(contigPosition._contigIndex);
			}

			vector<UnitigType> successors;
			_unitigGraph->getSuccessors(unitigIndexPred, successors);

			for(const UnitigType unitigIndexSucc : successors){
				
				//const UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc);

				if(_unitigName_to_contigIndexes.find(unitigIndexSucc) == _unitigName_to_contigIndexes.end()) continue;

				for(const ContigPosition& contigPosition : _unitigName_to_contigIndexes[unitigIndexSucc]){
					if(contigIndexPreds.find(contigPosition._contigIndex) != contigIndexPreds.end()) return true;
				}


			}
		}

		return false;

	}

	bool hasContigSuccessorBFS(const UnitigType unitigIndexSource, const unordered_set<UnitigType>& isUnitigIndexChimeric){


		unordered_set<UnitigType> isVisited;
		std::queue<UnitigType> queue;

		queue.push(unitigIndexSource);

		while(queue.size() > 0){

			const UnitigType unitigIndex = queue.front();
			queue.pop();

			if(isVisited.find(unitigIndex) != isVisited.end()) continue;
			isVisited.insert(unitigIndex);

			vector<UnitigType> successors;
			_unitigGraph->getSuccessors(unitigIndex, successors);

			for(const UnitigType unitigIndexSucc : successors){

				if(isUnitigIndexChimeric.find(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc)) != isUnitigIndexChimeric.end()) continue;
				if(_unitigName_to_contigIndexes.find(unitigIndexSucc) != _unitigName_to_contigIndexes.end()) return true;

				queue.push(unitigIndexSucc);
			}

		}
		/*
		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		for(const UnitigType unitigIndexSucc : successors){
			
			//const UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc);

			if(_unitigName_to_contigIndexes.find(unitigIndexSucc) != _unitigName_to_contigIndexes.end()) return true;

		}
		*/
		return false;
	}

	ofstream _colorFile;

	bool solvePerfectEdgeRepeats(){

		bool isModification = false;

		_colorFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/chimericUnitigs.csv");
		_colorFile << "Name,Color" << endl;

		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
			UnitigType unitigIndexRev = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

			bool isMod1 = solvePerfectEdgeRepeat(unitigIndex);
			bool isMod2 = solvePerfectEdgeRepeat(unitigIndexRev);

			if(isMod1 || isMod2) isModification = true;
		}

		_colorFile.close();

		return isModification;
	}
	
	
	bool solvePerfectEdgeRepeat(const UnitigType unitigIndex){
		
		//cout << "Solve perfect repeat: " << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;

		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		if(successors.size() < 2) return false;

		//cout << "\ta" << endl;
		vector<UnitigType> bridgedSuccessors;

		for(const UnitigType unitigIndexSucc : successors){
			
			
			KmerVec doublet;
			doublet._kmers = {unitigIndex, unitigIndexSucc};

			if(_solidDoublets.find(doublet) == _solidDoublets.end()) continue;

			bridgedSuccessors.push_back(unitigIndexSucc);
		}

		//cout << "\tb: " << bridgedSuccessors.size() << endl;
		if(bridgedSuccessors.size() != 1) return false;


		for(const UnitigType unitigIndexSucc : successors){
			//cout << "\tc: " << bridgedSuccessors.size() << "\t" << UnitigGraph2::unitigIndexToString(bridgedSuccessors[0]) << endl;
			
			vector<UnitigType> predecessors;
			_unitigGraph->getPredecessors(unitigIndexSucc, predecessors);

			bool isBridged = false;

			for(const UnitigType unitigIndexPred : predecessors){

				//cout << "\t\tPred: " << UnitigGraph2::unitigIndexToString(unitigIndexPred) << endl;

				//vector<UnitigType> successors;
				//_unitigGraph->getSuccessors(unitigIndexPred, successors);

				
				//for(const UnitigType unitigIndexSucc : successors){
				KmerVec doublet;
				doublet._kmers = {unitigIndexPred, unitigIndexSucc};

				//cout << "\t\t\tSucc: " << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << "\t" << (_solidDoublets.find(doublet) != _solidDoublets.end()) << endl;

				if(_solidDoublets.find(doublet) != _solidDoublets.end()){
					isBridged = true;
					break;
				}
				//}

			}
			
			if(!isBridged) return false;
		}


		//cout << "\tc" << endl;


		_unitigGraph->getSuccessors(unitigIndex, successors);

		for(const UnitigType unitigIndexSucc : successors){
			
			if(unitigIndexSucc == bridgedSuccessors[0]) continue;


			//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << ",red" << endl;
			//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) << ",red" << endl;
			_unitigGraph->removeSuccessor(unitigIndex, unitigIndexSucc);
			_unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));

		}

		return true;

		/*
		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		if(successors.size() < 2) return;

		vector<UnitigType> predecessors;

		for(const UnitigType unitigIndexSucc : successors){
			
			_unitigGraph->getPredecessors(unitigIndexSucc, predecessors);

			for(const UnitigType unitigIndexPred : predecessors){
				if(!hasSameSuccessors(unitigIndexPred, successors)) return;
			}

		}

		if(predecessors.size() != successors.size()) return;

		vector<pair<UnitigType, UnitigType>> bridgedEdges;

		for(const UnitigType unitigIndexPred : predecessors){
			vector<UnitigType> bridgedSuccessors = getBridgedSuccessors(unitigIndexPred);

			if(bridgedSuccessors.size() != 1) return;

			bridgedEdges.push_back({unitigIndexPred, bridgedSuccessors[0]});
		}

		
		for(const UnitigType unitigIndexPred : predecessors){
			
			vector<UnitigType> successors;
			_unitigGraph->getSuccessors(unitigIndexPred, successors);

			for(const UnitigType unitigIndexSucc : successors){
				pair<UnitigType, UnitigType> edge = {unitigIndexPred, unitigIndexSucc};
				if(std::find(bridgedEdges.begin(), bridgedEdges.end(), edge) == bridgedEdges.end()){
					
					//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred) << ",red" << endl;
					//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) << ",red" << endl;
					_unitigGraph->removeSuccessor(unitigIndexPred, unitigIndexSucc);
					_unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPred));

				}
			}
		}
		*/
	}

	/*
	void solvePerfectEdgeRepeat(const UnitigType unitigIndex){

		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		if(successors.size() < 2) return;

		vector<UnitigType> predecessors;

		for(const UnitigType unitigIndexSucc : successors){
			
			_unitigGraph->getPredecessors(unitigIndexSucc, predecessors);

			for(const UnitigType unitigIndexPred : predecessors){
				if(!hasSameSuccessors(unitigIndexPred, successors)) return;
			}

		}

		if(predecessors.size() != successors.size()) return;

		vector<pair<UnitigType, UnitigType>> bridgedEdges;

		for(const UnitigType unitigIndexPred : predecessors){
			vector<UnitigType> bridgedSuccessors = getBridgedSuccessors(unitigIndexPred);

			if(bridgedSuccessors.size() != 1) return;

			bridgedEdges.push_back({unitigIndexPred, bridgedSuccessors[0]});
		}

		
		for(const UnitigType unitigIndexPred : predecessors){
			
			vector<UnitigType> successors;
			_unitigGraph->getSuccessors(unitigIndexPred, successors);

			for(const UnitigType unitigIndexSucc : successors){
				pair<UnitigType, UnitigType> edge = {unitigIndexPred, unitigIndexSucc};
				if(std::find(bridgedEdges.begin(), bridgedEdges.end(), edge) == bridgedEdges.end()){
					
					//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred) << ",red" << endl;
					//_colorFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) << ",red" << endl;
					_unitigGraph->removeSuccessor(unitigIndexPred, unitigIndexSucc);
					_unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPred));

				}
			}
		}

	}

	
	bool hasSameSuccessors(const UnitigType unitigIndex, const vector<UnitigType>& checkSuccessors){

		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		if(successors.size() != checkSuccessors.size()) return false;

		for(const UnitigType unitigIndexSucc : successors){
			
			if(std::find(checkSuccessors.begin(), checkSuccessors.end(), unitigIndexSucc) == checkSuccessors.end()) return false;

		}

		return true;
	}

	vector<UnitigType> getBridgedSuccessors(const UnitigType unitigIndex){

		vector<UnitigType> bridgedSuccessors;

		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		for(const UnitigType unitigIndexSucc : successors){
			
			if(haveSameContigIndexes(unitigIndex, unitigIndexSucc)){
				bridgedSuccessors.push_back(unitigIndexSucc);
			}

		}

		return bridgedSuccessors;
	}

	bool haveSameContigIndexes(const UnitigType unitigIndex1, const UnitigType unitigIndex2){
		
		//const UnitigType unitigName1 = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex1);
		//const UnitigType unitigName2 = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex2);

		if(_unitigName_to_contigIndexes.find(unitigIndex1) == _unitigName_to_contigIndexes.end()) return false;
		if(_unitigName_to_contigIndexes.find(unitigIndex2) == _unitigName_to_contigIndexes.end()) return false;

		const vector<ContigPosition>& contigIndexes1 = _unitigName_to_contigIndexes[unitigIndex1];
		const vector<ContigPosition>& contigIndexes2 = _unitigName_to_contigIndexes[unitigIndex2];

		for(const ContigPosition& contigPosition1 : contigIndexes1){
			for(const ContigPosition& contigPosition2 : contigIndexes2){
				if(contigPosition1._contigIndex == contigPosition2._contigIndex) return true;
			}
			//if(std::find(contigIndexes2.begin(), contigIndexes2.end(), contigPosition._contigIndex) != contigIndexes2.end()) return true;
		}

		return false;
	}
	*/



	bool solveRepeats(){

		bool isModification = false;

		_colorFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/chimericUnitigs.csv");
		_colorFile << "Name,Color" << endl;

		for(size_t i=0; i<_unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
			//UnitigType unitigIndexRev = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);

			bool isMod = solveRepeat(unitigIndex);

			if(isMod) isModification = true;
			//solvePerfectEdgeRepeat(unitigIndexRev);
		}

		_colorFile.close();

		return isModification;

	}

	bool solveRepeat(const UnitigType unitigIndex){
		
		unordered_set<UnitigType> bridgedPreds;
		unordered_set<UnitigType> bridgedSuccs;

		const UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
		//cout << "utg" << unitigName << endl;
		vector<UnitigType> predecessors;
		_unitigGraph->getPredecessors(unitigIndex, predecessors);

		vector<UnitigType> successors;
		_unitigGraph->getSuccessors(unitigIndex, successors);

		if(successors.size() <= 1 && predecessors.size() <= 1) return false;

		for(const UnitigType unitigIndexPred : predecessors){
			
			//if(unitigName == 34977) cout << "\tPred: " << UnitigGraph2::unitigIndexToString(unitigIndexPred) << endl;

			for(const UnitigType unitigIndexSucc : successors){
				
				KmerVec triplet;
				triplet._kmers = {unitigIndexPred, unitigIndex, unitigIndexSucc};
				//triplet = triplet.normalize();

				//if(unitigName == 34977){
				//	cout << "\t\t" << UnitigGraph2::unitigIndexToString(triplet._kmers[0]) << "\t" << UnitigGraph2::unitigIndexToString(triplet._kmers[1]) << "\t" << UnitigGraph2::unitigIndexToString(triplet._kmers[2]) << "\t" << (_solidTriplets.find(triplet) != _solidTriplets.end()) << endl;
				//}

				if(_solidTriplets.find(triplet) != _solidTriplets.end()){
					bridgedPreds.insert(unitigIndexPred);
					bridgedSuccs.insert(unitigIndexSucc);
				}

			}
		}

		if(bridgedPreds.size() == predecessors.size() && bridgedSuccs.size() == successors.size()){
			_colorFile << "utg" << unitigName << "," << "red" << "\n";
			
			vector<UnitigGraph2::UnitigNode*> edgeNodes;
			UnitigGraph2::UnitigNode* node = _unitigGraph->_unitigs[unitigName];
			int nbEdgeNodes = max(predecessors.size(), successors.size());

			for(int i=0; i<nbEdgeNodes; i++){
				UnitigGraph2::UnitigNode* edgeNode = createEdgeNode(node);
				edgeNodes.push_back(edgeNode);
			}
			


            vector<UnitigType> unitigs;
            if(node->_unitigMerge.size() == 0){
                unitigs.push_back(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));    
            }
            else{
                unitigs = node->_unitigMerge;
            }

			vector<MinimizerType> minimizers;
			unitigsToMinimizers(unitigs, minimizers);

			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _params._kminmerSize);

			for(const KmerVec& kminmer : kminmers){
				const u_int128_t vecHash = kminmer.normalize().hash128();
				_repeatedKminmers.insert(vecHash);
			}
			
			for(const UnitigType unitigIndex : unitigs){
				_repeatedUnitigNames.insert(UnitigGraph2::unitigIndex_to_unitigName(unitigIndex));
			}

			_unitigGraph->removeNode(node);

			connectEdgesNodes(predecessors, successors, edgeNodes, unitigIndex);
			


			return true;

		}

		return false;
	}
	
	UnitigGraph2::UnitigNode* createEdgeNode(const UnitigGraph2::UnitigNode* parentNode){

		UnitigGraph2::UnitigNode* edgeNode = _unitigGraph->addNode(_unitigGraph->_unitigs.size(), parentNode->_nbMinimizers, true);
		edgeNode->_abundances = parentNode->_abundances;
		edgeNode->_abundance = parentNode->_abundance;
		edgeNode->_isReversed = parentNode->_isReversed;
		edgeNode->_unitigMerge = parentNode->_unitigMerge;

		_unitigName_to_minimizers.push_back(_unitigName_to_minimizers[parentNode->_unitigName]);
		
		return edgeNode;
	}

	void connectEdgesNodes(const vector<UnitigType>& predecessors, const vector<UnitigType>& successors, const vector<UnitigGraph2::UnitigNode*>& edgeNodes, const UnitigType repeatUnitigIndex){

		ankerl::unordered_dense::set<pair<UnitigType, UnitigType>> existingEdges;

		UnitigGraph2::UnitigNode* edgeNode = nullptr;
		UnitigType edgeNodeUnitigIndex = -1;

		for(size_t i=0; i<predecessors.size(); i++){

			const UnitigType unitigIndexSource = predecessors[i];
			//if(UnitigGraph2::unitigIndex_to_unitigName(repeatUnitigIndex) == 34977) cout << "\tPred: " << i << "\t" << UnitigGraph2::unitigIndexToString(unitigIndexSource) << endl;

			if(predecessors.size() >= successors.size()){
				edgeNode = edgeNodes[i];
				edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
			}

			for(size_t j=0; j<successors.size(); j++){
				

				if(predecessors.size() < successors.size()){
					edgeNode = edgeNodes[j];
					edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
				}

				const UnitigType unitigIndexDest = successors[j];

				//if(UnitigGraph2::unitigIndex_to_unitigName(repeatUnitigIndex) == 34977) cout << "\t\tSucc: " << j << "\t" << UnitigGraph2::unitigIndexToString(unitigIndexDest) << endl;

				KmerVec triplet;
				triplet._kmers = {unitigIndexSource, repeatUnitigIndex, unitigIndexDest};
				//triplet = triplet.normalize();

				if(_solidTriplets.find(triplet) != _solidTriplets.end()){


					std::pair<UnitigType, UnitigType> edge = {unitigIndexSource, edgeNodeUnitigIndex};

					if(existingEdges.find(edge) == existingEdges.end()){
						
						//if(UnitigGraph2::unitigIndex_to_unitigName(repeatUnitigIndex) == 34977){
						//	cout << "\t\tConnect: " << UnitigGraph2::unitigIndexToString(unitigIndexSource) << "\t" << UnitigGraph2::unitigIndexToString(edgeNodeUnitigIndex) << endl;
						//}

						_unitigGraph->addSuccessor(unitigIndexSource, edgeNodeUnitigIndex);
						_unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSource));
						existingEdges.insert(edge);
					}

					edge = {edgeNodeUnitigIndex, unitigIndexDest};

					if(existingEdges.find(edge) == existingEdges.end()){

						
						//if(UnitigGraph2::unitigIndex_to_unitigName(repeatUnitigIndex) == 34977){
						//	cout << "\t\tConnect: " << UnitigGraph2::unitigIndexToString(edgeNodeUnitigIndex) << "\t" << UnitigGraph2::unitigIndexToString(unitigIndexDest) << endl;
						//}

						_unitigGraph->addSuccessor(edgeNodeUnitigIndex, unitigIndexDest);
						_unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexDest), UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex));
						existingEdges.insert(edge);
					}
				}
			}	
		}
	}


	void writeUnitigSequences(){


		ofstream outputFile(_inputDir + "/unitigGraph.nodes.bin");

		for(size_t i=0; i<_unitigName_to_minimizers.size(); i++){
			writeUnitigSequence(outputFile, i, _unitigName_to_minimizers[i]);
		}

		outputFile.close();
	}

	void writeUnitigSequence(ofstream& outputFile, const UnitigType unitigName, const vector<MinimizerType>& minimizers){
	
		#pragma omp critical(writeUnitigSequence)
		{

			u_int32_t size = minimizers.size();
			UnitigType unitigIndex = unitigName*2;

			outputFile.write((const char*)&size, sizeof(size));
			outputFile.write((const char*)&minimizers[0], size * sizeof(MinimizerType));
			outputFile.write((const char*)&unitigIndex, sizeof(unitigIndex));

		}

	}

	void dumpContigNames(){

		ofstream contigNameFile("/pasteur/appa/homes/gbenoit/zeus/tmp/contigNames.csv");
		contigNameFile << "Name,ContigIndex" << endl;

		for(const auto& it : _unitigName_to_contigIndexes){

			int total = it.second.size();

			string s = "";
			for(size_t i=0; i<it.second.size() && i<5; i++){
				s += "ctg" + to_string(it.second[i]._contigIndex) + "-" + to_string(it.second[i]._position) + " ";
			}
			if(total > 5){
				s += " (" + to_string(total) + ")"; 
			}

			contigNameFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(it.first) << "," << s << "\n";
		}
		contigNameFile.close();
	}

	void dumpReferenceNames(){

		ofstream contigNameFile("/pasteur/appa/homes/gbenoit/zeus/tmp/referenceNames.csv");
		contigNameFile << "Name,ContigIndex" << endl;

		for(const auto& it : _unitigName_to_referenceName){

			int total = it.second.size();

			string s = "";
			for(size_t i=0; i<it.second.size() && i<5; i++){
				s += it.second[i]._referenceName + "-" + to_string(it.second[i]._position) + " ";
			}
			if(total > 5){
				s += " (" + to_string(total) + ")"; 
			}

			contigNameFile << "utg" << UnitigGraph2::unitigIndex_to_unitigName(it.first) << "," << s << "\n";
		}
		contigNameFile.close();
	}

};	


#endif 


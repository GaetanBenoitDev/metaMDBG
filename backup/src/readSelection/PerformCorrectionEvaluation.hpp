

#ifndef MDBG_METAG_PERFORMCORRECTIONEVALUATION
#define MDBG_METAG_PERFORMCORRECTIONEVALUATION

#include "../Commons.hpp"
#include "../utils/spoa64/include/spoa64/spoa.hpp"


class PerformCorrectionEvaluation : public Tool{
    
public:

	string _filename_exe;
	string _inputDir;
	int _nbCores;
	size_t _minimizerSize;
	size_t _kminmerSize;
	size_t _minimizerDensity;
	//string _contigFilename;
	string _inputFilename;

	vector<string>* _fields;
	string _alignFilename;

	u_int64_t _nbMatchesOriginal;
	u_int64_t _nbMatchesCorrected;
	long double _totalMinimizers;

	PerformCorrectionEvaluation(): Tool (){
	}


	void parseArgs(int argc, char* argv[]){

		_nbMatchesOriginal = 0;
		_nbMatchesCorrected = 0;
		_filename_exe = argv[0];

		args::ArgumentParser parser("lala", "");
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_contigFilename(parser, "contigFilename", "", args::Options::Required);
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
		//_contigFilename = args::get(arg_contigFilename);
		_nbCores = args::get(arg_nbCores);
		_inputFilename = _inputDir + "/input.txt";
		_alignFilename = _inputDir + "/align.paf";
		
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		_logFile << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		_kminmerSize = 2;

		_fields = new vector<string>();
	}

    void execute (){

		cout << "Loading read to contig mapping" << endl;
		loadReadToContigMapping();

		cout << "Loading contigs" << endl;
		loadContigs();

		cout << "Evaluating correction accuracy" << endl;
		evaluateCorrection();

		//closeLogFile();
	}

	struct ContigMatch{
		u_int32_t _contigIndex;
		bool _isReversed;
		u_int32_t _coverage;
	};

	unordered_map<u_int32_t, ContigMatch> _readIndex_to_bestContigMatch;

	void loadReadToContigMapping(){

		if(!fs::exists(_inputDir + "/readIndex_to_contigIndex.bin")){
			cout << "File does not exists, run setup first: " << _inputDir + "/readIndex_to_contigIndex.bin" << endl;
			exit(1);
		}

		ifstream file(_inputDir + "/readIndex_to_contigIndex.bin");

		while (true) {

			u_int32_t readIndex;
			u_int32_t contigIndex;
			bool isReversed;
			u_int32_t contigCoverage;

			file.read((char*)&readIndex, sizeof(readIndex));

			if(file.eof()) break;

			file.read((char*)&contigIndex, sizeof(contigIndex));
			file.read((char*)&isReversed, sizeof(isReversed));
			file.read((char*)&contigCoverage, sizeof(contigCoverage));

			_readIndex_to_bestContigMatch[readIndex] = {contigIndex, isReversed, contigCoverage};
			//cout << readIndex << " " << contigIndex << " " << isReversed << " " << contigCoverage << endl;
		}

		file.close();

	}


	vector<vector<MinimizerType>> _contigs;

	void loadContigs(){
		KminmerParserParallel parser(_inputDir + "/contig_data_init.txt", _minimizerSize, _kminmerSize, false, true, 1);
		parser.parseSequences(LoadContigFunctor(*this));
	}

	class LoadContigFunctor {

		public:

		PerformCorrectionEvaluation& _parent;

		LoadContigFunctor(PerformCorrectionEvaluation& parent) : _parent(parent){
		}

		LoadContigFunctor(const LoadContigFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			_parent._contigs.push_back({kminmerList._readMinimizers});


		}
	};

	bool _print_debug;

	void evaluateCorrection(){

		cout << "Did you copy? : cp " << _inputDir + "/read_data_corrected.txt" << " " << _inputDir + "/read_data_corrected_backup.txt" << endl;
		getchar();

		_print_debug = false;

		//fs::copy(_inputDir + "/read_data_corrected.txt", _inputDir + "/read_data_corrected_backup.txt", fs::copy_options::overwrite_existing);

		KminmerParserParallelPaired parser(_inputDir + "/read_data_init.txt", _inputDir + "/read_data_corrected_backup.txt", true, true, 1);
		//KminmerParserParallel parser2(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, 1);
		parser.parseSequences(EvaluationFunctor(*this));

	}

	class EvaluationFunctor {

		public:

		PerformCorrectionEvaluation& _parent;
		std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngineTrimming;

		EvaluationFunctor(PerformCorrectionEvaluation& parent) : _parent(parent){
		}

		EvaluationFunctor(const EvaluationFunctor& copy) : _parent(copy._parent){
			_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 5, -1, -1, -1, -1, -1);
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 5, -1, -1, -1, -1, -1);
			//_alignmentEngineTrimming = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -1, -1, 0, 0, 0);
		}

		~EvaluationFunctor(){
		}


		void operator () (const KminmerListPaired& kminmerList) {
			//_parent.correctRead(kminmerList, _alignmentEngine, _alignmentEngineTrimming);

			u_int32_t readIndex = kminmerList._readIndex;

			//if(readIndex != 155) return;

			vector<MinimizerType> readMinimizersOriginal = kminmerList._readMinimizers1;
			vector<MinimizerType> readMinimizersCorrected = kminmerList._readMinimizers2;

			int readLength = kminmerList._readMinimizers1pos[kminmerList._readMinimizers1pos.size()-1];
			//if(readLength < MIN_NB_MINIMIZER_REFERENCE || readMinimizersOriginal.size() < 4) return;
			if(readMinimizersOriginal.size() < 20) return;
			//if(readIndex != 3339) return;

			if(_parent._print_debug){
				cout << endl << endl;
				cout << "-------------------" << endl;
				cout << "Read index: " << readIndex << endl;
				cout << "Read size original:  " << readMinimizersOriginal.size() << endl;
				cout << "Read size corrected: " << readMinimizersCorrected.size() << endl;
			}

			//return;

			if(_parent._readIndex_to_bestContigMatch.find(readIndex) == _parent._readIndex_to_bestContigMatch.end()) return;

			const ContigMatch& contig = _parent._readIndex_to_bestContigMatch[readIndex];
			u_int32_t contigIndex = contig._contigIndex;
			bool isContigReversed = contig._isReversed;
			u_int32_t contigCoverage = contig._coverage;
			vector<MinimizerType> contigMinimizers = _parent._contigs[contigIndex];
			if(contigCoverage < 20){
				cout << "Skip low coverage contig" << endl;
				return;
			}

			if(_parent._print_debug){
				cout << "Contig index: " << contigIndex << " " << contigMinimizers.size() << " " << isContigReversed << endl;
				cout << "Contig coverage: " << contigCoverage << endl;
			}

			if(isContigReversed){
				//std::reverse(readMinimizersOriginal.begin(), readMinimizersOriginal.end());
				//std::reverse(readMinimizersCorrected.begin(), readMinimizersCorrected.end());
				std::reverse(contigMinimizers.begin(), contigMinimizers.end());
			}

			int contigStart = 0;
			int contigEnd = 0;
			int contigLength = 0;
			bool isHighNbErrors = false;
			if(_parent._print_debug) cout << endl << "Original alignment:" << endl;
			computeAccuracy(readMinimizersOriginal, contigMinimizers, false, contigLength, contigStart, contigEnd);

			if(readMinimizersCorrected.size() > 0){
				int contigLengthDummy = 0;
				int contigStartDummy = 0;
				int contigEndDummy = 0;
				if(_parent._print_debug) cout << endl << "Corrected alignment:" << endl;
				isHighNbErrors = computeAccuracy(readMinimizersCorrected, contigMinimizers, true, contigLengthDummy, contigStartDummy, contigEndDummy);
			}


			if(_parent._print_debug || isHighNbErrors){


				unordered_set<MinimizerType> subContigMinimizers;
				for(size_t i=contigStart; i<=contigEnd; i++){
					u_int64_t minimizer = contigMinimizers[i];
					subContigMinimizers.insert(minimizer);
				}
				
				cout << "Contig index: " << contigIndex << " " << contigMinimizers.size() << " " << isContigReversed << endl;
				cout << "Contig coverage: " << contigCoverage << endl;

				cout << endl << "\tRead original: " << endl;
				for(size_t i=0; i<readMinimizersOriginal.size(); i++){

					if(subContigMinimizers.find(readMinimizersOriginal[i]) == subContigMinimizers.end()){
						cout << "\t" << readMinimizersOriginal[i] << endl;
					}
					else{
						cout << "\t" << readMinimizersOriginal[i] << " *" << endl;
					}
				}

				cout << endl << "\tRead corrected: " << endl;
				for(size_t i=0; i<readMinimizersCorrected.size(); i++){

					if(subContigMinimizers.find(readMinimizersCorrected[i]) == subContigMinimizers.end()){
						cout << "\t" << readMinimizersCorrected[i] << endl;
					}
					else{
						cout << "\t" << readMinimizersCorrected[i] << " *" << endl;
					}
				}

				cout << endl << "\tContig sequences: " << endl;
				for(size_t i=contigStart; i<=contigEnd; i++){

					u_int64_t minimizer = contigMinimizers[i];
					
					if(std::find(readMinimizersCorrected.begin(), readMinimizersCorrected.end(), minimizer) == readMinimizersCorrected.end()){

						cout << "\t\t" << minimizer << endl;
					}
					else{
						cout << "\t\t" << minimizer << " *"<< endl;
					}
				}
			}

			#pragma omp atomic
			_parent._totalMinimizers += contigLength;

			cout << readIndex << endl;
			cout << "Contig length: " << contigLength << endl;
			cout << "Original nb matches : " << FormatWithCommas(_parent._nbMatchesOriginal) << " " << _parent._nbMatchesOriginal / _parent._totalMinimizers << endl;
			cout << "Corrected nb matches: " << FormatWithCommas(_parent._nbMatchesCorrected) << " " << _parent._nbMatchesCorrected / _parent._totalMinimizers<< endl;
			//cout << _parent._totalMinimizers << endl;
			
			if(_parent._print_debug || isHighNbErrors){
				getchar();
			}
			//if(_parent._print_debug) getchar();
		}


		bool computeAccuracy(const vector<MinimizerType>& readMinimizers, const vector<MinimizerType>& contigMinimizers, bool isCorrectedRead, int& contigLength, int& firstMatchPos, int& lastMatchPos){

			int nbErrors = 0;

			firstMatchPos = -1;
			lastMatchPos = -1;

			contigLength = 0;
			int startPos = -1;
			int endPos = -1;

			vector<MinimizerType> weights(readMinimizers.size(), 1);
			
			spoa64::Graph graph{};

			graph.AddAlignment(spoa64::Alignment(), readMinimizers, readMinimizers.size(), weights);
			spoa64::Alignment alignment = _alignmentEngine->Align(contigMinimizers, contigMinimizers.size(), graph);



			for (const auto& it : alignment) {

				int64_t v1 = it.first;
				int64_t v2 = it.second;
				
				
				if(startPos == -1) startPos = v1;
				endPos = v1;
				//cout << "\t" << v1 << " " << v2 << endl;

				if(v1 == -1){ //insert in

					if(isCorrectedRead){

						#pragma omp atomic
						_parent._nbMatchesCorrected -= 1;
					}
					else{

						#pragma omp atomic
						_parent._nbMatchesOriginal -= 1;
					}

					nbErrors += 1;
					contigLength += 1;
					if(_parent._print_debug) cout << "\tInsertion " << endl;

				}
				else if(v2 == -1){ //insert in

					nbErrors += 1;
					if(isCorrectedRead){

						#pragma omp atomic
						_parent._nbMatchesCorrected -= 1;
					}
					else{

						#pragma omp atomic
						_parent._nbMatchesOriginal -= 1;
					}

					if(_parent._print_debug) cout << "\tDeletion  "  << endl;
				}
				else if(readMinimizers[v1] == contigMinimizers[v2]){

					//cout << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;
					if(firstMatchPos == -1) firstMatchPos = v2;
					lastMatchPos = v2;

					contigLength += 1;

					if(isCorrectedRead){

						#pragma omp atomic
						_parent._nbMatchesCorrected += 1;
					}
					else{

						#pragma omp atomic
						_parent._nbMatchesOriginal += 1;
					}
					if(_parent._print_debug) cout << "\tMatch     " << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;

				}
				else{

					
					if(isCorrectedRead){

						#pragma omp atomic
						_parent._nbMatchesCorrected -= 1;
					}
					else{

						#pragma omp atomic
						_parent._nbMatchesOriginal -= 1;
					}
					

					nbErrors += 1;
					contigLength += 1;
					if(_parent._print_debug) cout << "\tMissmatch " << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;

				}
			}

			int nbStartingMissmatches = startPos;
			int nbEndingMissmatches = readMinimizers.size() - 1 - endPos;

			//cout << nbStartingMissmatches << " " << nbEndingMissmatches << endl;

			//contigLength += nbStartingMissmatches;
			//contigLength += nbEndingMissmatches;

			//cout << contigLength << endl;
			//cout << "lala" << firstMatchPos << " " << lastMatchPos << endl;
			
			if(isCorrectedRead && nbErrors > 5){
				cout << "Nb errors: " << nbErrors << endl;
				//getchar();
				return true;
			}

			return false;
		}

		template<class T>
		std::string FormatWithCommas(T value)
		{
			std::stringstream ss;
			ss.imbue(std::locale(""));
			ss << std::fixed << value;
			return ss.str();
		}

	};
};	


#endif 



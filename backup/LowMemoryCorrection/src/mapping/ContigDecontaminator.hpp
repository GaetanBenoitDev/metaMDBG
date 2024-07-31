

#ifndef MDBG_METAG_ContigDecontaminator
#define MDBG_METAG_ContigDecontaminator

#include "../Commons.hpp"
#include "../contigFeatures/ContigFeature.hpp"

class ContigDecontaminator : public Tool{
    
public:

	string _contigFilename;
	string _outputDir;
	string _kminmerCountFilename;
	ContigFeature _contigFeature;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	size_t _nbDatasets;
	KminmerCountMap _kminmerCounts;
	float _trimming;

	ContigDecontaminator(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("decontaminate", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		args::Positional<std::string> arg_kminmerCountFilename(parser, "kminmerCountFilename", "", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
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

		_contigFilename = args::get(arg_contigs);
		_kminmerCountFilename = args::get(arg_kminmerCountFilename);
		_outputDir = args::get(arg_outputDir);
		_trimming = 0.05;

		if(!fs::exists(_outputDir)){
			fs::create_directories(_outputDir);
		}

	}


    void execute (){

		//loadKminmerCounts();

		cout << "Processing contigs" << endl;
		ReadParserParallel readParser(_contigFilename, true, false, 1, _logFile);
		readParser.parse(DecontaminationFunctor(*this));
	}

	void loadKminmerCounts(){

		cout << "Loading kminmer counts" << endl;

		ifstream inputFile(_kminmerCountFilename);

		inputFile.read((char*)&_minimizerSize, sizeof(_minimizerSize));
		inputFile.read((char*)&_minimizerDensity, sizeof(_minimizerDensity));
		inputFile.read((char*)&_kminmerSize, sizeof(_kminmerSize));
		inputFile.read((char*)&_nbDatasets, sizeof(_nbDatasets));

		cout << _minimizerSize << " " << _minimizerDensity << " " << _kminmerSize << " "<< _nbDatasets << endl;
		//getchar();

		while(true){
			vector<u_int64_t> minimizerSeq;
			vector<u_int32_t> counts;

			minimizerSeq.resize(_kminmerSize);
			inputFile.read((char*)&minimizerSeq[0], _kminmerSize*sizeof(u_int64_t));

			if(inputFile.eof()) break;

			counts.resize(_nbDatasets);
			inputFile.read((char*)&counts[0], _nbDatasets*sizeof(u_int32_t));

			KmerVec vec;
			vec._kmers = minimizerSeq;

			_kminmerCounts[vec] = counts;
			//for(u_int32_t count : counts){
			//	cout << count << " ";
			//}
			//cout << endl;

			//getchar();
		}

		/*
		for(auto& it : _kminmerCounts){

			vector<u_int64_t> minimizerSeq = it.first._kmers;
			vector<u_int32_t> counts = it.second;

			outputFile.write((const char*)&minimizerSeq[0], _minimizerSize*sizeof(uint64_t));
			outputFile.write((const char*)&counts[0], _nbDatasets*sizeof(uint32_t));

			//for(u_int32_t count : counts){
				//cout << count << " ";
			//}
			//cout << endl;

		}
		*/

		inputFile.close();
	}

	class DecontaminationFunctor {

		public:

		ContigDecontaminator& _contigDecontaminator;		
		size_t _nbDatasets;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		DecontaminationFunctor(ContigDecontaminator& contigDecontaminator) : _contigDecontaminator(contigDecontaminator){
			_nbDatasets = _contigDecontaminator._nbDatasets;
			_minimizerParser = new MinimizerParser(_contigDecontaminator._minimizerSize, _contigDecontaminator._minimizerDensity);
		}

		DecontaminationFunctor(const DecontaminationFunctor& copy) : _contigDecontaminator(copy._contigDecontaminator){
			_nbDatasets = _contigDecontaminator._nbDatasets;
			_minimizerParser = new MinimizerParser(_contigDecontaminator._minimizerSize, _contigDecontaminator._minimizerDensity);
		}

		~DecontaminationFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {

			const string& sequence = read._seq;
			cout << sequence.size() << endl;

			/*
			string referenceSequence = sequence.substr(3000000, 200000);
			vector<float> referenceComposition;
			_contigDecontaminator._contigFeature.sequenceToComposition(referenceSequence, referenceComposition);

			for(size_t i=0; i<sequence.size()-200000; i+=200000){

				string substring = sequence.substr(i, 200000);


				vector<float> composition;
				_contigDecontaminator._contigFeature.sequenceToComposition(substring, composition);

				cout << i << "    " << (_contigDecontaminator._contigFeature.cal_tnf_dist(referenceComposition, composition, referenceSequence.size(), substring.size())) << endl;
			}
			*/


			size_t fragLengthInc = 50000;
			size_t fragLength = 200000;

			if(sequence.size() < fragLength *4) return;


			ofstream outputMatrix(_contigDecontaminator._outputDir + "/" + read._header + ".csv");

			string header = ";";
			for(size_t i=0; i<sequence.size()-fragLength; i+=fragLengthInc){
				header += to_string(i) + ";";
			}
			header.pop_back();

			outputMatrix << header << endl;



			for(size_t i=0; i<sequence.size()-fragLength; i+=fragLengthInc){

				string line = to_string(i) + ";";
				

				string referenceSequence = sequence.substr(i, fragLength);
				vector<float> referenceComposition;
				_contigDecontaminator._contigFeature.sequenceToComposition(referenceSequence, referenceComposition);

				vector<float> abundancesMean_reference;
				vector<float> abundancesVar_reference;
				//computeFragmentAbundance(referenceSequence, abundancesMean_reference, abundancesVar_reference);




				for(size_t j=0; j<sequence.size()-fragLength; j+=fragLengthInc){

					string substring = sequence.substr(j, fragLength);


					vector<float> composition;
					_contigDecontaminator._contigFeature.sequenceToComposition(substring, composition);

					vector<float> abundancesMean;
					vector<float> abundancesVar;
					//computeFragmentAbundance(substring, abundancesMean, abundancesVar);

					/*
					cout << endl << endl;
					for(float count : abundancesMean_reference){
						cout << count << " ";
					}
					cout << endl;
					for(float count : abundancesMean){
						cout << count << " ";
					}
					cout << endl;
					*/

					int nnz = 0;
					//float distance = _contigDecontaminator._contigFeature.cal_abd_dist_new(abundancesMean_reference, abundancesVar_reference, abundancesMean, abundancesVar, nnz);
					
					//cout << i << " " << j << "    " << (_contigDecontaminator._contigFeature.cal_tnf_dist(referenceComposition, composition, referenceSequence.size(), substring.size())) << endl;

					float distance = _contigDecontaminator._contigFeature.computeEuclideanDistance(referenceComposition, composition);
					//float distance = _contigDecontaminator._contigFeature.cal_tnf_dist(referenceComposition, composition, referenceSequence.size(), substring.size());
					//cout << i << " " << j << " " << distance << endl;

					//getchar();
					line += to_string(distance) + ";";
				}
				
				line.pop_back();
				outputMatrix << line << endl;
			}

			outputMatrix.close();
		}

		void computeFragmentAbundance(const string& sequence, vector<float>& abundancesMean, vector<float>& abundancesVar){

			abundancesMean.resize(_nbDatasets, 0);
			abundancesVar.resize(_nbDatasets, 0);

			vector<long double> sums(_nbDatasets, 0);
			vector<long double> nbValues(_nbDatasets, 0);
			vector<long double> means(_nbDatasets, 0);
			vector<long double> sumsVar(_nbDatasets, 0);

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(sequence.c_str(), sequence.size(), rleSequence, rlePositions);

			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfos;
			MDBG::getKminmers(_contigDecontaminator._minimizerSize, _contigDecontaminator._kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfos, rlePositions, 0, false);

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				KmerVec vec = kminmers[i];
				const vector<u_int32_t>& counts = _contigDecontaminator._kminmerCounts[vec];

				for(size_t j=0; j<_nbDatasets; j++){

					sums[j] += counts[j];
					nbValues[j] += 1;
				}

			}

			for(size_t i=0; i<_nbDatasets; i++){
				if(nbValues[i] == 0) continue;
				abundancesMean[i] = sums[i] / nbValues[i];
			}

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				KmerVec vec = kminmers[i];
				const vector<u_int32_t>& counts = _contigDecontaminator._kminmerCounts[vec];

				for(size_t j=0; j<_nbDatasets; j++){

					double count = counts[j];
					sumsVar[j] += ((count - abundancesMean[j]) * (count - abundancesMean[j]));

				}

			}

			for(size_t i=0; i<_nbDatasets; i++){
				if(nbValues[i] <= 1) continue;
				abundancesVar[i] = sumsVar[i] / (nbValues[i]-1);
			}


		}

	};



};	


#endif 

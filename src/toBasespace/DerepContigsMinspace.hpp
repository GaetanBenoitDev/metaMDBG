

#ifndef MDBG_METAG_TOBASESPACENOCORRECTION
#define MDBG_METAG_TOBASESPACENOCORRECTION

#include "../Commons.hpp"
#include "OverlapRemover.hpp"
#include "OverlapRemover2.hpp"
#include "RepeatRemover.hpp"

class ToBasespaceNoCorrection : public Tool{
    
public:

	//string _inputFilename;
	string _inputFilenameContig;
	//string _inputFilenameContig_fasta;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	size_t _kminmerSizeFirst;
	//bool _isFirstPass;
	//bool _isOutputFasta;
	int _nbCores;
	bool _hasQuality;

	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;
    size_t _kminmerSizeLast;
    size_t _meanReadLength;

	//string _filename_outputContigs;
	//string _filename_kminmerSequences;
	//MDBG* _mdbg;
	//EncoderRLE _encoderRLE;


	ToBasespaceNoCorrection(): Tool (){


	}

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("toBasespaceFast", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		args::Flag arg_hasQuality(parser, "", "Is quality in read data", {"has-quality"});
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

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
		_inputFilenameContig = args::get(arg_inputContigFilename);
		//_filename_outputContigs = args::get(arg_outputContigFilename);
		//_inputFilename = args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);

		_hasQuality = false;
		if(arg_hasQuality){
			_hasQuality = true;
		}

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzclose(file_parameters);

		_kminmerSize = _kminmerSizeFirst;

		openLogFile(_inputDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//_logFile << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity;
		Logger::get().debug() << "";

	}


    void execute (){
		


		const string& contigFilenameNoOverlaps = _inputFilenameContig + ".nooverlaps";

		{

			//cout << "A remettre overlap remover" << endl;
			OverlapRemover overlapRemover(_inputDir, _inputFilenameContig, contigFilenameNoOverlaps, _kminmerSize);
			//OverlapRemover2 overlapRemover(_inputDir, _inputFilenameContig, contigFilenameNoOverlaps, _kminmerSizeFirst, _nbCores);
			overlapRemover.execute();

		}
		
		const string& contigFilenameNoRepeats = _inputFilenameContig + ".norepeats";

		{
			RepeatRemover repeatRemover(_inputDir, contigFilenameNoOverlaps, contigFilenameNoRepeats, _kminmerSizeFirst+1, _minimizerDensity, _hasQuality, _nbCores);
			repeatRemover.execute();
		}
		//closeLogFile();

	}

};	


#endif 



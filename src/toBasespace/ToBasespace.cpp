

#include "ToBasespace.hpp"


ToBasespace::ToBasespace () : Tool("toBasespace"){

	getParser()->push_back (new OptionOneParam (STR_INPUT, "input contig filename in minimizer space", true));
	getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename in basespace", true));

}


void ToBasespace::execute ()
{

	_inputContigFilename = getInput()->getStr(STR_INPUT);
	_outputContigFilename = getInput()->getStr(STR_OUTPUT);

	cout << endl;
	cout << "Input filename: " << _inputContigFilename << endl;
	cout << "Output filename: " << _outputContigFilename << endl;
	//cout << "Kminmer length: " << _kminmerSize << endl;
	//cout << "Minimizer length: " << _minimizerSize << endl;
	//cout << "Density: " << _minimizerDensity << endl;
	cout << endl;
}

/*
int main (int argc, char* argv[])
{
    try
    {
    	ToBasespace().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}*/
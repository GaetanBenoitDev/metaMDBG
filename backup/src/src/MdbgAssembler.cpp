#include "Bloocoo.hpp"
#include "assembly/Assembly2.hpp"
#include "assembly/Assembly3.hpp"
#include "assembly/GenerateContigs.hpp"
#include "toBasespace/ToBasespace.hpp"
#include "toBasespace/ToBasespaceNoCorrection.hpp"
#include "toBasespace/ToMinspace.hpp"
#include "pipeline/AssemblyPipeline.hpp"
#include "pipeline/BinningPipeline.hpp"
#include "readSelection/ReadSelection.hpp"
#include "contigFeatures/KminmerCounter.hpp"
#include "contigFeatures/KmerCounter.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/Mapping_BinReads.hpp"

void displayHelp(){
	cout << "Usage: ./simkaMin [option]" << endl;

	cout << endl << "[Distance computation options]" << endl;
	cout << "\tsketch      : transform datasets in small sketches of k-mers and their abundance" << endl;
	cout << "\tdistance    : compute Jaccard and Bray-Curtis distances between sketches" << endl;

	cout << endl << "[Distance matrix manipulation options]" << endl;
	cout << "\texport      : export distance matrices stored in binary format" << endl;
	//cout << "\tmatrix-update       : update existing distance matrices" << endl;

	cout << endl << "[Sketch options]" << endl;
	cout << "\tappend      : merge multiple sketch files into a single one" << endl;
	cout << "\tinfo        : list datasets contained in a sketch file" << endl;

	cout << endl;
}

int main (int argc, char* argv[])
{
    try
    {
    	if(argc < 2){
    		displayHelp();
    	}
    	else{

    		//std::vector<char*>  args;
    		vector<char*> argsTemp( argv, argv + argc );
    		argsTemp.erase(argsTemp.begin()+1);
    		//std::transform(argsTemp.begin(), argsTemp.end(), std::back_inserter(vc), convert);

    		char** args = &argsTemp[0];
    		//char* args[];

    		//for(string& arg: argsTemp){

    		//}
    		//rArray = new char*[argc+1];
    		//for(int i=0; i <= argc; i++) {
    		//    rArray[i] = argv[i];
    		//}
    		// use rArray
    		//delete [] rArray;


    		//char* args = new char*[argc-1];
    		//vector<string> test;

    		//for(size_t i=0; i<argc; i++){
    		//	if (i==1) continue;
    		//	args[i] = argv[i];
    		//}

    		argc -= 1;
    		//vector<char*> args(argv);

    		string programName = string(argv[1]);

			if(programName == "asm"){
                AssemblyPipeline().run (argc, args);
    		}
			else if(programName == "bin"){
                BinningPipeline().run (argc, args);
    		}
			else if(programName == "readSelection"){
                ReadSelection().run (argc, args);
    		}
    		else if(programName == "graph"){
                Bloocoo().run (argc, args);
    		}
    		else if(programName == "binPass"){
                Assembly2().run (argc, args);
    		}
    		else if(programName == "contig"){
                GenerateContigs().run (argc, args);
    		}
    		else if(programName == "multik"){
                Assembly3().run (argc, args);
    		}
    		else if(programName == "toBasespace"){
                ToBasespace().run (argc, args);
    		}
    		else if(programName == "toBasespaceFast"){
                ToBasespaceNoCorrection().run (argc, args);
    		}
    		else if(programName == "toMinspace"){
                ToMinspace().run (argc, args);
    		}
    		else if(programName == "count"){
                KminmerCounter().run (argc, args);
    		}
    		else if(programName == "countKmer"){
                KmerCounter().run (argc, args);
    		}
    		else if(programName == "map"){
                Mapping().run (argc, args);
    		}
    		else if(programName == "map-binreads"){
                Mapping_BinReads().run (argc, args);
    		}
    		else{
    			displayHelp();
    		}
    	}
    	//cout << argc << endl;
    	//cout << argv[0] << endl;
    	//cout << argv[1] << endl;
    	//cout << argv[2] << endl;
    	//
    }
    catch (int e)
    {
        cout << "EXCEPTION: " << e << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

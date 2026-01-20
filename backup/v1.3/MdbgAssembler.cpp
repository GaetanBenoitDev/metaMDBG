

#ifndef MDBG_METAG_MDBGASSEMBLER
#define MDBG_METAG_MDBGASSEMBLER


#include "toBasespace/ToBasespace2.hpp"
#include "pipeline/AssemblyPipeline.hpp"
/*
//#include "polish/ContigPolisher.hpp" //!

#include "readSelection/ReadCorrection.hpp"

//#include "readSelection/SetupCorrectionEvaluationSimulation.hpp"
#include "readSelection/SetupCorrectionEvaluation.hpp"

//#include "readSelection/PerformCorrectionEvaluation.hpp"

#include "readSelection/ReadSelection.hpp"

#include "toBasespace/ToMinspace.hpp"
//#include "mapping/MappingContigToGraph.hpp"



//#include "assembly/Circulizer.hpp"
//#include "readSelection/ReadSelectionPaired.hpp"
#include "graph/CreateMdbg.hpp"
//#include "assembly/Assembly2.hpp"
//#include "assembly/Assembly3.hpp"
#include "assembly/GenerateContigs.hpp"
#include "toBasespace/ToBasespace.hpp"
//#include "toBasespace/ToBasespace3.hpp"
#include "toBasespace/ToBasespaceNoCorrection.hpp"
//#include "pipeline/ReassemblePipeline.hpp"
//#include "pipeline/BinningPipeline.hpp"
//#include "contigFeatures/KminmerCounter.hpp"
//#include "contigFeatures/KmerCounter.hpp"
//#include "mapping/Mapping.hpp"
//#include "mapping/ReferencePath.hpp"
//#include "mapping/Mapping_BinReads.hpp"
//#include "mapping/Mapping_ContigToMDBG.hpp"
//#include "mapping/Mapping_BinToMDBG.hpp"
//#include "polish/PurgeDups.hpp" //!
//#include "contigFeatures/KminmerCounter.hpp"
//#include "assembly/Circulizer2.hpp"
//#include "contigFeatures/KmerHicLinker.hpp"
#include "graph/GenerateGfa.hpp"

//#include "toBasespace/Dereplicater.hpp"
//#include "polish/PurgeGraph.hpp"
//#include "mapping/ContigDecontaminator.hpp"
*/


//"reprire: tester remove small contigs humanO1 avec minimap2 au lieu de wfmash"
void displayHelp(string programName){
	cout << " Program: metaMDBG (assembly of long metagenomics reads)" << endl;
	cout << " Version: " << METAMDBG_VERSION << endl;
	cout << " Contact: Gaëtan Benoit (gaetanbenoitdev@gmail.com)" << endl;
	cout << endl;
	cout << " Usage: " + programName + " [command]" << endl;
	cout << endl;
	cout << " command:" << endl;
	cout << " \tasm:    perform read assembly" << endl;
	cout << " \tpolish: polish contigs" << endl;
	//cout << " \tderep    : purge strain duplication" << endl;
	cout << " \tgfa:    generate an assembly graph (.gfa). Require a finished metaMDBG run" << endl;
	//cout << "\treadSelection      : transform readset into its minimizer reprentation" << endl;
	//cout << "\tdgraph   : create minimizer de-bruijn graph" << endl;

	cout << endl;
}


int main (int argc, char* argv[])
{
    try
    {
    	if(argc < 2){
    		displayHelp(argv[0]);
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

if(programName == "toBasespace"){
                ToBasespace2().run (argc, args);
    		}
			else  if(programName == "asm"){
                AssemblyPipeline().run (argc, args);
    		}
			/*
			//if(programName == "polish"){
            //    ContigPolisher().run (argc, args); //!
    		//}
			
			
			//else if(programName == "readSelectionPaired"){
            //    ReadSelectionPaired().run (argc, args);
    		//}
    		if(programName == "readCorrection"){
                ReadCorrection().run (argc, args);
    		}
			
    		else if(programName == "setupCorr"){
                SetupCorrectionEvaluation().run (argc, args);
    		}
			
    		//else if(programName == "evalCorr"){
            //    PerformCorrectionEvaluation().run (argc, args);
    		//}
			else if(programName == "readSelection"){
                ReadSelection().run (argc, args);
    		}

			
			//else if(programName == "extend"){
                //Circulizer().run (argc, args);
    		//}
    		else if(programName == "toMinspace"){
                ToMinspace().run (argc, args);
    		}
			//else if(programName == "map"){
              //  MappingContigToGraph().run (argc, args);
    		//}
	
			
			else if(programName == "graph"){
                CreateMdbg().run (argc, args);
    		}
			else if(programName == "asm"){
                AssemblyPipeline().run (argc, args);
    		}
			//else if(programName == "reasm"){
            //    ReassemblePipeline().run (argc, args);
    		//}
			//else if(programName == "bin"){
                //BinningPipeline().run (argc, args);
    		//}
    		//else if(programName == "binPass"){
            //    Assembly2().run (argc, args);
    		//}
    		else if(programName == "contig"){
                GenerateContigs().run (argc, args);
    		}
    		//else if(programName == "multik"){
            //    Assembly3().run (argc, args);
    		//}
    		//else if(programName == "toBasespaceFast"){
            //    ToBasespaceNoCorrection().run (argc, args);
    		//}
    		//else if(programName == "circ"){
            //    Circulizer2().run (argc, args);
    		//}
    		//else if(programName == "toBasespace_hifi"){
            //    ToBasespace().run (argc, args);
    		//}
    		else if(programName == "toBasespace"){
                ToBasespace2().run (argc, args);
    		}
    		//else if(programName == "derepgraph"){
                //PurgeGraph().run (argc, args);
    		//}
    		//else if(programName == "derep"){
            //    PurgeDups().run (argc, args); //!
    		//}
    		//else if(programName == "derepOld"){
            //    Dereplicater().run (argc, args);
    		//}
    		else if(programName == "toBasespaceFast"){
                ToBasespaceNoCorrection().run (argc, args);
    		}
    		//else if(programName == "count"){
            //    KminmerCounter().run (argc, args);
    		//}
    		//else if(programName == "countKmer"){
            //    KmerCounter().run (argc, args);
    		//}
    		//else if(programName == "mapBin"){
            //    Mapping().run (argc, args);
    		//}
    		//else if(programName == "mapRef"){
            //    ReferencePath().run (argc, args);
    		//}
    		//else if(programName == "map-binreads"){
            //    Mapping_BinReads().run (argc, args);
    		//}
    		//else if(programName == "map-binToMDBG"){
            //    Mapping_BinToMDBG().run (argc, args);
    		//}
    		else if(programName == "gfa"){
                GenerateGfa().run (argc, args);
    		}
    		//else if(programName == "decontaminate"){
                //ContigDecontaminator().run (argc, args);
    		//}
    		//else if(programName == "map-contigToMDBG"){
            //    Mapping_ContigToMDBG().run (argc, args);
    		//}
    		else{
    			displayHelp(argv[0]);
    		}
			*/
			
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


#endif

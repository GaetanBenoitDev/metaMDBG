

#include <Bloocoo.hpp>

using namespace std;

void displayVersion(std::ostream& os){
	
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	os << "* Bloocoo version "<< BLOOCOO_VERSION_MAJOR << "."
	<< BLOOCOO_VERSION_MINOR << "."
	<< BLOOCOO_VERSION_PATCH
	<< "                   *" << endl; //<< " AGPL licence" <<endl;
	os << "* Using gatb-core version "<< System::info().getVersion() <<  "           *" << endl;
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
}


/********************************************************************************/

int main (int argc, char* argv[])
{
	
	
	if(argc > 1 && (   strcmp(argv[1],STR_VERSION)==0 || strcmp(argv[1],"-v")==0    )     ){
		displayVersion(cout);
		return EXIT_FAILURE;
	}
	
	
    // We define a try/catch block in case some method fails
    try
    {
        Bloocoo().run (argc, argv);
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

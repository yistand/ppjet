#include "DetaDphi.h"

#include <string>
#include <cstring>
#include <vector>

int main( int argc, const char** argv ) {

	//const char *defaults[] = {"macroDetaDphi","~/Scratch/pp200Y12_jetunderlying/","FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160314","14","200","0","0","5","200","1"};
	const char *defaults[] = {"macroDetaDphi","~/Scratch/pp200Y12_jetunderlying/","R0.2FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_HadrCorr_160314","10","200","0","0","5","200","1"};

        if ( argc==1 ) {
                argv=defaults;
                argc=sizeof (defaults ) / sizeof (defaults[0] );
        }
        // Throw arguments in a vector
        // ---------------------------
        std::vector<std::string> arguments(argv + 1, argv + argc);

	DetaDphi *etaphi = new DetaDphi();
	etaphi->SetSaveRoot(1);

	etaphi->DoDetaDphi(arguments.at(0),arguments.at(1),std::atof(arguments.at(2).c_str()),std::atof(arguments.at(3).c_str()),std::atoi(arguments.at(4).c_str()),std::atoi(arguments.at(5).c_str()),std::atof(arguments.at(6).c_str()),std::atof(arguments.at(7).c_str()),std::atoi(arguments.at(8).c_str()));

}




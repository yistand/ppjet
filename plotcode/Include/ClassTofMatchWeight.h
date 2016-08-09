// class ClassTofMatchWeight
// 
// 2016.03.29	Li YI
// read weight by tof matching efficiency vs eta from histogram in root file and return weight by input eta.
// read eff using function from fitting result
//
//
#ifndef ROOT_TofMatchWeight
#define ROOT_TofMatchWeight


#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"

#define VZARRAYLENGTH 12

class ClassTofMatchWeight {

	double vzarray4etaeff[VZARRAYLENGTH+1];

	TH1D *hrVseta;
	TFile *f4eta;

	TF1 *feffpt;

	double effcutoff;                       // if efficiency too low, weight will be set as zero to reduce too large correction


	public:
		ClassTofMatchWeight();
		~ClassTofMatchWeight();

		double GetTofMatchWeightVsEta(double eta, double vz=0);
		double GetTofMatchWeightVsPt(double pt);

		double TofMatchPPY12(double eta, double pt, double vz=0);

                // Helper
                double GetEffCutOff();
                void SetEffCutOff(double cutoff);

};
		



#endif

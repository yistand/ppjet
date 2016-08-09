// class ClassTPCWeight
// 
// 2016.03.29	Li YI
// read weight by TPC efficiency vs eta, pt from fitting functions
// pt one including the whole eff value, therefore eta one only has shape dependence
//
//
#ifndef ROOT_TPCEFF
#define ROOT_TPCEFF


#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"

class ClassTPCWeight {

	TF1 *feffpt;
	TF1 *feffeta;

	double effcutoff;			// if efficiency too low, weight will be set as zero to reduce too large correction

	public:
		ClassTPCWeight();
		~ClassTPCWeight();

		double GetTpcWeightVsEta(double eta);
		double GetTpcWeightVsPt(double pt);

		double TpcWeightPPY12(double eta, double pt);

		// Helper
		double GetEffCutOff();
		void SetEffCutOff(double cutoff);

};
		

#endif

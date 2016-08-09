#include "ClassTofMatchWeight.h"

// Define the constructor
ClassTofMatchWeight::ClassTofMatchWeight() {

	effcutoff = 0.2;
	
	// vs eta
	f4eta = new TFile("/home/hep/caines/ly247/ppjet/plotcode/Include/ppTOFeff.root");
	//hrVseta = (TH1D*)f4eta->Get("hrVseta2");
	// not sure how to initialize array simply
	double a[VZARRAYLENGTH+1] = {-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30};
	for(int i = 0; i<VZARRAYLENGTH; i++) {
		vzarray4etaeff[i] = a[i];
	}

	// vs pt
	feffpt=new TF1("feffpt","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);
	feffpt->SetParameters(7.28477e-01, 1.52059e-01, 5.28031); 
}


// Define the destructor.
ClassTofMatchWeight::~ClassTofMatchWeight() {
	// Deallocate the memory that was previously reserved
	if(f4eta) f4eta->Close();
	if(feffpt) feffpt->Delete();
}

double ClassTofMatchWeight::GetTofMatchWeightVsEta(double eta, double vz) {

	int ivz = 0;
	while(vz>vzarray4etaeff[ivz]) ivz++;
	if(ivz>=VZARRAYLENGTH) ivz=VZARRAYLENGTH-1;

	hrVseta = (TH1D*)f4eta->Get(Form("hrVseta_vzdep%d",ivz));

	double eff = hrVseta->GetBinContent(hrVseta->FindBin(eta));

        double weight = 0;
        if(eff>0) weight = 1./eff;
        
        return weight;

}


double ClassTofMatchWeight::GetTofMatchWeightVsPt(double pt) {

	double eff = 0;

	if(pt>4.5) eff = feffpt->Eval(4.5);
	else if(pt<0.15) eff = 0;
	else eff = feffpt->Eval(pt);

        double weight = 0;
        if(eff>0) weight = 1./eff;
        
        return weight;

}

double ClassTofMatchWeight::TofMatchPPY12(double eta, double pt, double vz) {
	
	double weight = GetTofMatchWeightVsPt(pt)*GetTofMatchWeightVsEta(eta, vz);

        if(effcutoff>0 && weight>(1./effcutoff)) {
                weight = 0;
        }

        return weight;

}





// Helper
double ClassTofMatchWeight::GetEffCutOff() { 
        return effcutoff;
}

void ClassTofMatchWeight::SetEffCutOff(double cutoff) {
        effcutoff = cutoff;     
}       
        





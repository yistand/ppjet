#include "ClassTPCWeight.h"

double etafun_tpc(double *x, double *par) ;

// Define the constructor
ClassTPCWeight::ClassTPCWeight() {

	effcutoff = 0.2;

	// vs eta
	feffeta=new TF1("feffeta",etafun_tpc,-1,1,6);
	feffeta->SetParameters(1.00723e+00,4.18242e+00,8.10867e+01,1.00394e+00,4.20742e+00,9.46997e+01);

	// vs pt
	feffpt=new TF1("feffpt","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);
	feffpt->SetParameters(0.874739, 0.156624, 5.67316); 
}


// Define the destructor.
ClassTPCWeight::~ClassTPCWeight() {
	// Deallocate the memory that was previously reserved
	if (feffeta)  feffeta->Delete() ;
	if (feffpt)  feffpt->Delete() ;
}


double ClassTPCWeight::GetTpcWeightVsEta(double eta) {

	double eff = 0;

	if(fabs(eta)<=1) eff = feffeta->Eval(eta);		// do not use particle outside |eta|<1

	double weight = 0; 
	if(eff>0) weight = 1./eff;

	return weight;

}


double ClassTPCWeight::GetTpcWeightVsPt(double pt) {

	double eff = 0;

	if(pt>4.5) eff = feffpt->Eval(4.5);
	else if(pt<0.15) eff = 0;
	else eff = feffpt->Eval(pt);

	double weight = 0; 
	if(eff>0) weight = 1./eff;

	return weight;

}

double ClassTPCWeight::TpcWeightPPY12(double eta, double pt) {
	
	double weight = GetTpcWeightVsPt(pt)*GetTpcWeightVsEta(eta);
	if(effcutoff>0 && weight>(1./effcutoff)) {
		weight = 0;
	}

	return weight;

}


// Helper
double ClassTPCWeight::GetEffCutOff() {
	return effcutoff;
}

void ClassTPCWeight::SetEffCutOff(double cutoff) {
	effcutoff = cutoff;	
}



double etafun_tpc(double *x, double *par) {

	double eta = x[0];
	double y = 0; 
		
	double turnpoint = 0.5;
	double maxeta = 5.3;
	double p0 = par[0];
	double p1 = par[1];
	double p2 = par[2];
	double p3 = par[3];
	double p4 = par[4];
	double p5 = par[5];
	if(eta<-turnpoint&&eta>=-1) {
		y = p0*exp(-pow(p1/(eta+maxeta),p2));
	}
	if(fabs(eta)<=turnpoint) {
		double k = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))-(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/(fabs(turnpoint)*2.);
		double b = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))+(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/2.;
		y = k*eta+b;
	}
	if(eta>turnpoint&&eta<=1) {
		y = p3*exp(-pow(p4/(maxeta-eta),p5));
	}

	return y;

}




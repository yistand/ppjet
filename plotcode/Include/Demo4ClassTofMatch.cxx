#include "ClassTofMatchWeight.h"
#include "ClassTPCWeight.h"
#include <iostream>

using namespace std;


void Demo4ClassTpc (double eta = 0.9, double pt=0.5) {

	ClassTPCWeight *tpcweight = new ClassTPCWeight();	
	
	cout<<"eta = "<<eta<<" pt = "<<pt<<"     tpcw = "<<tpcweight->TpcWeightPPY12(eta,pt)<<endl;

}

void Demo4ClassTofMatch (double eta = 0.9, double pt=0.5) {

	ClassTofMatchWeight *tofweight = new ClassTofMatchWeight();	
	
	cout<<"eta = "<<eta<<" pt = "<<pt<<"     tofw = "<<tofweight->TofMatchPPY12(eta,pt)<<endl;

}


int main() {
	Demo4ClassTofMatch();
	Demo4ClassTpc();
}


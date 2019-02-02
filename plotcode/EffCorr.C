//====================================================================================================
//
//	2018.01.21	Li Yi
//	for PID yield, apply a rough eff. corr. 
//	using TPC eff. from PRC 79 034909 (2009)
//
//====================================================================================================

double tpceff(double pt, string species) {

	//cout<<"tpceff "<<pt<<" "<<species<<endl;

	if(pt<0.2) return 0;
	double tmppt;
	if(pt>1.4) tmppt=1.4;
	else tmppt = pt;

	if(!species.compare("pion")) {
		double p0 = 0.856;
		double p1 = 0.075;
		double p2 = 1.668;
		//cout<<tmppt<<" "<<p0*exp(-pow(p1/tmppt,p2))<<endl;
		return p0*exp(-pow(p1/tmppt,p2));
	}
	else if(!species.compare("kaon")) {
		double p0 = 0.527;
		double p1 = 0.241;
		double p2 = 3.496;
		double p3 = 0.160;
		return p0*exp(-pow(p1/tmppt,p2))+p3*tmppt;
	}
	else if(!species.compare("proton")) {
		double p0 = 0.830;
		double p1 = 0.295;
		double p2 = 7.005;
		double p3 = -0.029;
		//return (p0*exp(-pow(p1/tmppt,p2))+p3*tmppt)*exp(pow(p4/tmppt,p5));
		return p0*exp(-pow(p1/tmppt,p2))+p3*tmppt;
	}
	else {
		cout<<species<<" cannot find tpceff"<<endl;
		return 0;
	}
}

double tofeff(double pt, string species) {

	if(pt<0.2) return 0;
	if(pt<0.25) return 0.57;
	if(pt<0.35) return 0.7;
	if(pt<0.45) return 0.72;
	else return 0.73;

}

TGraphErrors* Dat2TGraph(string datfilename) {
	ifstream datin;
	datin.open(datfilename);
	vector<double> pt, dpt, y, ey, sy;
	double ipt, idpt, iy, iey, isy;
	char c;
	string s;

	TGraphErrors *gr = new TGraphErrors();
	string name =  datfilename;
	cout<<name<<endl;
	name.insert(0,"gr_");
	cout<<name<<endl;
	name.erase(name.end()-4,name.end());
	cout<<name<<endl;
	gr->SetName(TString(name));
	gr->SetTitle(TString(name));

	int ip = 0;

	while(datin >> ipt >> idpt >> iy >> c >> c >> c >> iey >> s) {
		cout<<ipt<<" "<<iy<<" "<<iey<<endl;
		pt.push_back(ipt);
		dpt.push_back(idpt/2.);
		y.push_back(iy);
		ey.push_back(iey);
		sy.push_back(isy);
		gr->SetPoint(ip,ipt,iy);
		gr->SetPointError(ip,idpt/2.,iey);
		ip++;
	}

	return gr;
}

TGraphErrors *Add2TGraph(TGraphErrors *gr1, TGraphErrors *gr2, TString name) {

	double *x1 = gr1->GetX();
	double *y1 = gr1->GetY();
	double *ex1 = gr1->GetEX();
	double *ey1 = gr1->GetEY();

	double *x2 = gr2->GetX();
	double *y2 = gr2->GetY();
	double *ex2 = gr2->GetEX();
	double *ey2 = gr2->GetEY();

	double y, ey;
	TGraphErrors *gr = new TGraphErrors();
	gr->SetName(name);
	gr->SetTitle(name);
	for(int i=0; i<gr1->GetN(); i++) {
		if(fabs(x1[i]-x2[i])<1e-6) {
			y = (y1[i]+y2[i]);
			ey = sqrt(ey1[i]*ey1[i]+ey2[i]*ey2[i]);
			gr->SetPoint(i,x1[i],y);
			gr->SetPointError(i,ex1[i],ey);
		}
	}

	return gr;

}

void EffCorr() {
	TFile *fin = new TFile("yieldpid.root");
	const int Nspecies = 3;
	string Species[Nspecies] = {"pion","kaon","proton"};
	TH1D *hyraw[Nspecies];
	for(int i = 0; i<Nspecies; i++) {
		hyraw[i] = (TH1D*)fin->Get(Form("hYieldspid%d",i));
	}
	TH1D *hpt = (TH1D*)fin->Get("hptave");

	TH1D *hy[Nspecies];
	for(int i = 0; i<Nspecies; i++) {
		hy[i] = (TH1D*)hyraw[i]->Clone(Form("hCorrYieldspid%d",i));
	}
	for(int i = 0; i<hy[0]->GetNbinsX(); i++) {
		double pt = hpt->GetBinContent(i+1);
		double norm = 1./((2*TMath::Pi())*hy[0]->GetBinWidth(i+1)*pt*2.);
		// test if(pt>0.2) {
		if(pt>0.2 && pt<2) {	//test
			//cout<<pt<<" "<<hy[0]->GetBinContent(i+1)<<" "<<tpceff(pt,"pion")<<" "<<tofeff(pt,"pion")<<endl;
			hy[0]->SetBinContent(i+1,norm*hy[0]->GetBinContent(i+1)/(tpceff(pt,Species[0])*tofeff(pt,Species[0])));
			hy[1]->SetBinContent(i+1,norm*hy[1]->GetBinContent(i+1)/(tpceff(pt,Species[1])*tofeff(pt,Species[1])));
			hy[2]->SetBinContent(i+1,norm*hy[2]->GetBinContent(i+1)/(tpceff(pt,Species[2])*tofeff(pt,Species[2])));
			hy[0]->SetBinError(i+1,norm*hy[0]->GetBinError(i+1)/(tpceff(pt,Species[0])*tofeff(pt,Species[0])));
			hy[1]->SetBinError(i+1,norm*hy[1]->GetBinError(i+1)/(tpceff(pt,Species[1])*tofeff(pt,Species[1])));
			hy[2]->SetBinError(i+1,norm*hy[2]->GetBinError(i+1)/(tpceff(pt,Species[2])*tofeff(pt,Species[2])));
		}
		else {
			hy[0]->SetBinContent(i+1,0);
			hy[1]->SetBinContent(i+1,0);
			hy[2]->SetBinContent(i+1,0);
			hy[0]->SetBinError(i+1,0);
			hy[1]->SetBinError(i+1,0);
			hy[2]->SetBinError(i+1,0);
		}
		cout<<hy[0]->GetBinContent(i+1)<<", "<<hy[1]->GetBinContent(i+1)<<", "<<hy[2]->GetBinContent(i+1)<<endl;
		cout<<hy[0]->GetBinError(i+1)<<", "<<hy[1]->GetBinError(i+1)<<", "<<hy[2]->GetBinError(i+1)<<endl;
	}

	string pubdatname[6] = {"pionplus.dat","pionminus.dat","kaonplus.dat","kaonminus.dat","proton.dat","pbar.dat"};
	TGraphErrors* gr[6];
	for(int i = 0; i<6; i++) {
		gr[i] = Dat2TGraph(pubdatname[i]);
	}
	TGraphErrors *gsum[3];
	gsum[0] = Add2TGraph(gr[0],gr[1],"gsum_pion");
	gsum[1] = Add2TGraph(gr[2],gr[3],"gsum_kaon");
	gsum[2] = Add2TGraph(gr[4],gr[5],"gsum_proton");


	TFile *fout = new TFile("corryieldpid.root","RECREATE");
	for(int i = 0; i<Nspecies; i++) {
		hy[i]->Write();
	}
	for(int i = 0; i<6; i++) {
		gr[i]->Write();
	}
	for(int i = 0; i<3; i++) {
		gsum[i]->Write();
	}
	fout->Close();
}

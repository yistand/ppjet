//====================================================================================================
//
//	2017.08.30	Li YI
//	plot unfolding closure test Half MC for train, Half for test 
//
//====================================================================================================

#include "plotcompare30Vs60.C"
#include "plothrecoWunfolderr.C"

void plotcompareTestHalfHalf() {
	const int Nregion = 3;
	const int Nvar = 2;
	const char *nametag[Nregion*Nvar] = {"LeadAreaNtrk","SubAreaNtrk","TranTotNtrk","LeadPtAve","SubPtAve","TranPtAve"};
	const char *nametitle[Nregion*Nvar] = {"Toward #LTdN/d#etad#phi#GT","Away #LTdN/d#etad#phi#GT","Transverse #LTdN/d#etad#phi#GT","Toward #LTp_{T}#GT","Away #LTp_{T}#GT","Transverse #LTp_{T}#GT"};

	TH2D *htrain[Nregion][Nvar];
	TH2D *htraintrue[Nregion][Nvar];
	TH2D *hmeas[Nregion][Nvar];
	TH2D *hreco[Nregion][Nvar];
	TH2D *htrue[Nregion][Nvar];
	
	TProfile *ptrain[Nregion][Nvar];
	TProfile *ptraintrue[Nregion][Nvar];
	TProfile *pmeas[Nregion][Nvar];
	TProfile *preco[Nregion][Nvar];
	TProfile *ptrue[Nregion][Nvar];
	
	TFile *f[Nregion][Nvar];

	for(int i = 0; i<Nregion; i++) {
		for(int j = 0; j<Nvar; j++) {

			f[i][j] = new TFile(Form("TrainTestHalfHalf_%sMBCharged_McPt02.root",nametag[i*Nvar+j]));

			htrain[i][j]=(TH2D*)f[i][j]->Get("htrain");
			htraintrue[i][j]=(TH2D*)f[i][j]->Get("htraintrue");
			hmeas[i][j]=(TH2D*)f[i][j]->Get("hmeas");
			hreco[i][j]=(TH2D*)f[i][j]->Get("hreco");
			htrue[i][j]=(TH2D*)f[i][j]->Get("htrue");

			ptrain[i][j]=(TProfile*)htrain[i][j]->ProfileX("ptrain");
			ptraintrue[i][j]=(TProfile*)htraintrue[i][j]->ProfileX("ptraintrue");
			pmeas[i][j]=(TProfile*)hmeas[i][j]->ProfileX("pmeas");
			preco[i][j]=(TProfile*)hreco[i][j]->ProfileX("preco");
			ptrue[i][j]=(TProfile*)htrue[i][j]->ProfileX("ptrue");
		}
	}

	// Normalized ProfileX to density
	for(int i = 0; i<Nregion; i++) {
		for(int j = 0; j<Nvar; j++) {
			if(strstr(nametag[i*Nvar+j],"Ntrk")||strstr(nametag[i*Nvar+j],"PtSum")) {
				cout<<nametag[i*Nvar+j]<<" scale"<<endl;
				float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
				ptrain[i][j]->Scale(DeDpNorma);
				ptraintrue[i][j]->Scale(DeDpNorma);
				pmeas[i][j]->Scale(DeDpNorma);
				preco[i][j]->Scale(DeDpNorma);
				ptrue[i][j]->Scale(DeDpNorma);
			}
		}
	}

	TCanvas *c[Nregion][Nvar];
	TCanvas *cratio[Nregion][Nvar];
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TProfile *rp_true[Nregion][Nvar];

	TLegend *leg[Nregion][Nvar];
	TLine *line = new TLine(0,1,55,1);

	for(int i = 0; i<Nregion; i++) {
		for(int j = 0; j<Nvar; j++) {
			c[i][j] = new TCanvas();
			ptraintrue[i][j]->GetXaxis()->SetTitle("Leading jet p_{T}");
			ptraintrue[i][j]->GetYaxis()->SetTitle(nametitle[i*Nvar+j]);
			ptraintrue[i][j]->SetLineStyle(2);
			ptraintrue[i][j]->Draw("HIST");
			ptrue[i][j]->Draw("HISTsame");
			preco[i][j]->Draw("HISTsame");
			leg[i][j] = new TLegend(0.54,0.2,0.84,0.5);
			leg[i][j]->AddEntry(ptraintrue[i][j],"Train Truth","l");
			leg[i][j]->AddEntry(ptrue[i][j],"Test Truth","l");
			leg[i][j]->AddEntry(preco[i][j],"Test Unfolded","l");
			leg[i][j]->Draw();


			if(0) {
				cratio[i][j] = new TCanvas();
				rp_true[i][j] = (TProfile*)Ratio2Profile(ptrue[i][j],preco[i][j]);
				rp_true[i][j]->SetName(Form("ratio_ptrue%s",nametag[i*Nvar+j]));
				rp_true[i][j]->GetXaxis()->SetTitle("Leading jet p_{T}");
				rp_true[i][j]->GetYaxis()->SetTitle(Form("Ratio Truth/Unfolded %s",nametitle[i*Nvar+j]));
				rp_true[i][j]->SetMaximum(1.2);
				rp_true[i][j]->SetMinimum(0.8);
				rp_true[i][j]->Draw("HIST");
				line->Draw();
			}
		}
	}




}

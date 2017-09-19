//=====================================================================================================
//
//	2017.05.10	Li YI
//	Plot Transverse density vs collision energy 
//	CDF PhysRevD.92.092009, ALICE JHEP07(2012)116
//	|eta|<0.8, pT>0.5, pTMax (CDF) 5-6 GeV/c, (STAR) 5-7 GeV/c, (ALICE) 4-10 GeV/c
//
//=====================================================================================================



void EnergyDep() {
	// STAR from //maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_WideBin_CDFcutPT05eta08.root 
	// 		maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin_CDFcutPT05eta08
	const int N_star = 1;	
	double ene_star[N_star] = {200};
	//double NtrkDen_star[N_star] = {1.1385269/(2*0.8*2*TMath::Pi()/3.)};		// this (maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_WideBin_CDFcutPT05eta08.root) assume TPC effiency eta dependence is flat and uses the same eff as |eta|<1 for |eta|<0.8. Actually there is a small dependence in eta, which is eff_0.8eta = 1.006 x eff_1eta. While, we don't want to rerun the code, and simply assume there is no TPC eff pT dependence at pT>0.5, use 0.875 from JetMinWorkshop_160526YL.pdf for |eta|<1, and then multiply by 1.006, get 0.8802. The root files used in below calculations is maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin_CDFcutPT05eta08.root. Note: even the efficiency is corrected, the pile up is not removed
	double NoEffCorr_NtrkDen_star[N_star] = {0.99562487/(2*0.8*2*TMath::Pi()/3.)};	
	double NoEffCorr_err_star[N_star] = {0.00082724/(2*0.8*2*TMath::Pi()/3.)};
	double NtrkDen_star[N_star] = {0};	
	double err_star[N_star] = {};
	for(int i = 0; i<N_star; i++) {NtrkDen_star[i]=NoEffCorr_NtrkDen_star[i]/0.8802; err_star[i]=NoEffCorr_err_star[i]/0.8802; }
	// if not TPC eff correction,  Add absolute +- 5% error
	double maxNtrkDen_star[N_star] = {0};
	double minNtrkDen_star[N_star] = {0};
	double syserrmax_star[N_star] = {0};
	double syserrmin_star[N_star] = {0};
	for(int i = 0; i<N_star; i++) { 
		maxNtrkDen_star[i] = NoEffCorr_NtrkDen_star[i]/(0.8802-0.05);
		minNtrkDen_star[i] = NoEffCorr_NtrkDen_star[i]/(0.8802+0.05);
		syserrmax_star[i] = maxNtrkDen_star[i]-NtrkDen_star[i]; 
		syserrmin_star[i] = NtrkDen_star[i]-minNtrkDen_star[i]; 
	}

	// vs leading full jet pT 11-15 GeV/c, |eta|<1, pT>0.5 GeV/c
	double NtrkDen_star_jet12[N_star] = {0.295548};		// from SysErr4Unfolding_TranTotNtrkJPCharged_NFWeight_12JetBinv2_McPtRC02MC05_embedMB_Baye5.root
	double err_star_jet12[N_star] = {0.000143950};
	double syserrmax_star_jet12[N_star] = {0.0451854};
	double syserrmin_star_jet12[N_star] = {0.0299265};
	// vs leading full jet pT 20-25 GeV/c, |eta|<1, pT>0.5 GeV/c
	double NtrkDen_star_jet20[N_star] = {0.22077};		// from SysErr4Unfolding_TranTotNtrkJPCharged_NFWeight_12JetBinv2_McPtRC02MC05_embedMB_Baye5.root
	double err_star_jet20[N_star] = {0.0007};
	double syserrmax_star_jet20[N_star] = {0.0508};
	double syserrmin_star_jet20[N_star] = {0.0204};


	// CDF PhysRevD.92.092009, 300GeV, 900GeV, 1.96TeV, x: leading charged track. |eta|<0.8, pT>0.5 GeV/c
	// 	data from https://hepdata.net/record/ins1388868 Table 3, Table 11, Table 19. pTMax 5-6 GeV
	const int N_cdf = 3;	
	double ene_cdf[N_cdf] = {300, 900, 1960};
	double NtrkDen_cdf[N_cdf] = {0.3031,0.4855, 0.6013};	
	double err_cdf[N_cdf] = {0.0092,0.0028,0.0025};
	double syserr_cdf[N_cdf] = {0.0103,0.0141,0.0211};
	// CDF PhysRevD.65.092002, 1.8 TeV, x: leading charged jet. |eta|<1, pT>0.5 GeV/c 
	// 	data from https://hepdata.net/record/ins564673 Table 4. pT^lead 19-20 GeV, need x1.09 for tracking efficiency
	const int N_cdf_chargejet = 1;	
	double ene_cdf_chargejet[N_cdf_chargejet] = {1800};
	double NtrkDen_cdf_chargejet[N_cdf_chargejet] = {2.34/(2*2*TMath::Pi()/3.)*1.09};	
	double err_cdf_chargejet[N_cdf_chargejet] = {0.166/(2*2*TMath::Pi()/3.)*1.09};
	double syserr_cdf_chargejet[N_cdf_chargejet] = {0.166/(2*2*TMath::Pi()/3.)*1.09};
	
	//// ALICE JHEP07(2012)116, Table 7, Leading charged track, fit to pT 4-10 GeV/c for 900 GeV, 4-25 for 7 TeV.
	//const int N_alice = 2;	
	//double ene_alice[N_alice] = {900, 7000};
	//double NtrkDen_alice[N_alice] = {0.45,0.95};
	//double err_alice[N_alice] = {0.02,0.03};

	// ALICE JHEP07(2012)116, Table 7, Leading charged track
	// 	data from https://hepdata.net/record/ins1080735 Table 3. leading charged pT: 5-6 GeV. |eta|<0.8, pT>0.5
	const int N_alice = 2;	
	double ene_alice[N_alice] = {900, 7000};
	double NtrkDen_alice[N_alice] = {0.4603,0.886};
	double err_alice[N_alice] = {0.023,0.0042};
	double syserr_alice[N_alice] = {0.01,0.04};

	// CMS JHEP 1509(2015) 137,2015, x: leading charged jet. Fig 1. Data from https://hepdata.net/record/ins1385107 Table 1. 17GeV-22GeV
	// |eta|<2, pT>0.5
	const int N_cms_chargejet = 1;
	double ene_cms_chargejet[N_cms_chargejet] = {2760};
	double NtrkDen_cms_chargejet[N_cms_chargejet] = {0.721766};
	double err_cms_chargejet[N_cms_chargejet] = {0.009122};
	double syserr_cms_chargejet[N_cms_chargejet] = {0.023943};
	
	// ATLAS Phys.Rev. D83 (2011) 112001, 2011, x: leading track. pp 900 GeV and 7 TeV
	//	data from https://hepdata.net/record/ins879407. Table 1 and Table 2.  Leading Track pT 5-5.5 GeV/c, |eta|<2.5, pT>0.5
	// ATLAS arxiv1701.05390 or JHEP03(2017)157, x: leading charged particle. pp 13 TeV
	// 	data from https://hepdata.net/record/ins1509919. Table 21. Leading charged particle pT 5-5.5 GeV/c, |eta|<2.5, pT>0.5
	const int N_atlas = 3;
	double ene_atlas[N_atlas] = {900,7000,13000};
	double NtrkDen_atlas[N_atlas] = {0.4412,0.87588,1.06};
	double err_atlas[N_atlas] = {0.01895,0.00264,0.00959};
	double syserr_atlas[N_atlas] = {0.01572,0.03197,0.005922};
	
	// ATLAS Eur.Phys.J. C74 (2014) 2965, 2014, x: leading full jet. pp 7 TeV
	//	data from https://hepdata.net/record/ins1298811 Table 4. Inclusive jet 20-30 GeV/c, |eta|<2.5, pT>0.5
	const int N_atlas_jet = 1;
	double ene_atlas_jet[N_atlas_jet] = {7000};
	double NtrkDen_atlas_jet[N_atlas_jet] = {1.044};
	double err_atlas_jet[N_atlas_jet] = {0.001651};
	double syserr_atlas_jet[N_atlas_jet] = {0.04};
	

	// STAR
	// vs Track pTMax
	TGraphErrors *gr_star = new TGraphErrors();
	gr_star->SetName("star");
	//gr_star->SetTitle("STAR Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_star->SetTitle("STAR Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_star; i++) {
		gr_star->SetPoint(i,ene_star[i],NtrkDen_star[i]);
		gr_star->SetPointError(i,0,err_star[i]);
	}
	gr_star->SetMarkerColor(2);
	gr_star->SetLineColor(2);
	gr_star->SetMarkerStyle(29);
	gr_star->SetMarkerSize(3);
	
	TGraphAsymmErrors *sysgr_star = new TGraphAsymmErrors();
	sysgr_star->SetName("star_sys");
	//sysgr_star->SetTitle("STAR Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_star->SetTitle("STAR Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_star; i++) {
		sysgr_star->SetPoint(i,ene_star[i],NtrkDen_star[i]);
		sysgr_star->SetPointEYlow(i,syserrmin_star[i]);
		sysgr_star->SetPointEYhigh(i,syserrmax_star[i]);
	}
	sysgr_star->SetMarkerColor(gr_star->GetMarkerColor());
	sysgr_star->SetLineColor(gr_star->GetLineColor());
	sysgr_star->SetMarkerStyle(gr_star->GetMarkerStyle());
	sysgr_star->SetMarkerSize(gr_star->GetMarkerSize());

	// vs leading jet pT, |eta|<1
	TGraphErrors *gr_star_jet = new TGraphErrors();
	gr_star_jet->SetName("star_jet");
	//gr_star_jet->SetTitle("STAR Transverse #LTN_{ch}/#delta#eta#delta#phi#GT vs p_{T}^{leading jet}");
	gr_star_jet->SetTitle("STAR Transverse #LTdN_{ch}/d#etad#phi#GT vs p_{T}^{leading jet}");
	for(int i = 0; i<N_star; i++) {
		gr_star_jet->SetPoint(i,ene_star[i],NtrkDen_star_jet20[i]);
		gr_star_jet->SetPointError(i,0,err_star_jet20[i]);
	}
	gr_star_jet->SetMarkerColor(2);
	gr_star_jet->SetLineColor(2);
	gr_star_jet->SetMarkerStyle(29);
	gr_star_jet->SetMarkerSize(3);
	
	TGraphAsymmErrors *sysgr_star_jet = new TGraphAsymmErrors();
	sysgr_star_jet->SetName("star_jet_sys");
	//sysgr_star_jet->SetTitle("STAR Transverse #LTN_{ch}/#delta#eta#delta#phi#GT vs p_{T}^{leading jet}");
	sysgr_star_jet->SetTitle("STAR Transverse #LTdN_{ch}/d#etad#phi#GT vs p_{T}^{leading jet}");
	for(int i = 0; i<N_star; i++) {
		sysgr_star_jet->SetPoint(i,ene_star[i],NtrkDen_star_jet20[i]);
		sysgr_star_jet->SetPointEYlow(i,syserrmin_star_jet20[i]);
		sysgr_star_jet->SetPointEYhigh(i,syserrmax_star_jet20[i]);
	}
	sysgr_star_jet->SetMarkerColor(gr_star_jet->GetMarkerColor());
	sysgr_star_jet->SetLineColor(gr_star_jet->GetLineColor());
	sysgr_star_jet->SetMarkerStyle(gr_star_jet->GetMarkerStyle());
	sysgr_star_jet->SetMarkerSize(gr_star_jet->GetMarkerSize());


	// CDF
	TGraphErrors *gr_cdf = new TGraphErrors();
	gr_cdf->SetName("cdf");
	//gr_cdf->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_cdf->SetTitle("Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cdf; i++) {
		gr_cdf->SetPoint(i,ene_cdf[i],NtrkDen_cdf[i]);
		gr_cdf->SetPointError(i,0,err_cdf[i]);
	}
	gr_cdf->SetMarkerColor(1);
	gr_cdf->SetLineColor(1);
	gr_cdf->SetMarkerStyle(8);
	gr_cdf->SetMarkerSize(2);

	TGraphAsymmErrors *sysgr_cdf = new TGraphAsymmErrors();
	sysgr_cdf->SetName("cdf_sys");
	//sysgr_cdf->SetTitle("CDF Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_cdf->SetTitle("CDF Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cdf; i++) {
		sysgr_cdf->SetPoint(i,ene_cdf[i],NtrkDen_cdf[i]);
		sysgr_cdf->SetPointEYlow(i,syserr_cdf[i]);
		sysgr_cdf->SetPointEYhigh(i,syserr_cdf[i]);
	}
	sysgr_cdf->SetMarkerColor(gr_cdf->GetMarkerColor());
	sysgr_cdf->SetLineColor(gr_cdf->GetLineColor());
	sysgr_cdf->SetMarkerStyle(gr_cdf->GetMarkerStyle());
	sysgr_cdf->SetMarkerSize(gr_cdf->GetMarkerSize());


	TGraphErrors *gr_cdf_chargejet = new TGraphErrors();
	gr_cdf_chargejet->SetName("cdf_chargejet2");
	//gr_cdf_chargejet->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_cdf_chargejet->SetTitle("CDF Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cdf_chargejet; i++) {
		gr_cdf_chargejet->SetPoint(i,ene_cdf_chargejet[i],NtrkDen_cdf_chargejet[i]);
		gr_cdf_chargejet->SetPointError(i,0,err_cdf_chargejet[i]);
	}
	gr_cdf_chargejet->SetMarkerColor(1);
	gr_cdf_chargejet->SetLineColor(1);
	gr_cdf_chargejet->SetMarkerStyle(25);
	gr_cdf_chargejet->SetMarkerSize(1.5);

	TGraphAsymmErrors *sysgr_cdf_chargejet = new TGraphAsymmErrors();
	sysgr_cdf_chargejet->SetName("cdf_chargejet_sys");
	//sysgr_cdf_chargejet->SetTitle("CDF Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_cdf_chargejet->SetTitle("CDF Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cdf_chargejet; i++) {
		sysgr_cdf_chargejet->SetPoint(i,ene_cdf_chargejet[i],NtrkDen_cdf_chargejet[i]);
		sysgr_cdf_chargejet->SetPointEYlow(i,syserr_cdf_chargejet[i]);
		sysgr_cdf_chargejet->SetPointEYhigh(i,syserr_cdf_chargejet[i]);
	}
	sysgr_cdf_chargejet->SetMarkerColor(gr_cdf_chargejet->GetMarkerColor());
	sysgr_cdf_chargejet->SetLineColor(gr_cdf_chargejet->GetLineColor());
	sysgr_cdf_chargejet->SetMarkerStyle(gr_cdf_chargejet->GetMarkerStyle());
	sysgr_cdf_chargejet->SetMarkerSize(gr_cdf_chargejet->GetMarkerSize());



	// ALICE
	TGraphErrors *gr_alice = new TGraphErrors();
	gr_alice->SetName("alice");
	//gr_alice->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_alice->SetTitle("ALICE Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_alice; i++) {
		gr_alice->SetPoint(i,ene_alice[i],NtrkDen_alice[i]);
		gr_alice->SetPointError(i,0,err_alice[i]);
	}
	int color_alice = kAzure+1;
	gr_alice->SetMarkerColor(color_alice);
	gr_alice->SetLineColor(color_alice);
	gr_alice->SetMarkerStyle(8);
	gr_alice->SetMarkerSize(2);


	// CMS
	TGraphErrors *gr_cms_chargejet = new TGraphErrors();
	gr_cms_chargejet->SetName("cms_chargejet");
	//gr_cms_chargejet->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_cms_chargejet->SetTitle("CMS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cms_chargejet; i++) {
		gr_cms_chargejet->SetPoint(i,ene_cms_chargejet[i],NtrkDen_cms_chargejet[i]);
		gr_cms_chargejet->SetPointError(i,0,err_cms_chargejet[i]);
	}
	int color_cms = kOrange+2;
	gr_cms_chargejet->SetMarkerColor(color_cms);
	gr_cms_chargejet->SetLineColor(color_cms);
	gr_cms_chargejet->SetMarkerStyle(25);
	gr_cms_chargejet->SetMarkerSize(2);

	TGraphAsymmErrors *sysgr_cms_chargejet = new TGraphAsymmErrors();
	sysgr_cms_chargejet->SetName("cms_sys");
	//sysgr_cms_chargejet->SetTitle("CMS Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_cms_chargejet->SetTitle("CMS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_cms_chargejet; i++) {
		sysgr_cms_chargejet->SetPoint(i,ene_cms_chargejet[i],NtrkDen_cms_chargejet[i]);
		sysgr_cms_chargejet->SetPointEYlow(i,syserr_cms_chargejet[i]);
		sysgr_cms_chargejet->SetPointEYhigh(i,syserr_cms_chargejet[i]);
	}
	sysgr_cms_chargejet->SetMarkerColor(gr_cms_chargejet->GetMarkerColor());
	sysgr_cms_chargejet->SetLineColor(gr_cms_chargejet->GetLineColor());
	sysgr_cms_chargejet->SetMarkerStyle(gr_cms_chargejet->GetMarkerStyle());
	sysgr_cms_chargejet->SetMarkerSize(gr_cms_chargejet->GetMarkerSize());

	// ATLAS
	TGraphErrors *gr_atlas = new TGraphErrors();
	gr_atlas->SetName("atlas");
	//gr_atlas->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_atlas->SetTitle("ATLAS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_atlas; i++) {
		gr_atlas->SetPoint(i,ene_atlas[i],NtrkDen_atlas[i]);
		gr_atlas->SetPointError(i,0,err_atlas[i]);
	}
	int color_atlas = kGreen+2;
	gr_atlas->SetMarkerColor(color_atlas);
	gr_atlas->SetLineColor(color_atlas);
	gr_atlas->SetMarkerStyle(8);
	gr_atlas->SetMarkerSize(2);

	TGraphAsymmErrors *sysgr_atlas = new TGraphAsymmErrors();
	sysgr_atlas->SetName("atlas_sys");
	//sysgr_atlas->SetTitle("ATLAS Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_atlas->SetTitle("ATLAS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_atlas; i++) {
		sysgr_atlas->SetPoint(i,ene_atlas[i],NtrkDen_atlas[i]);
		sysgr_atlas->SetPointEYlow(i,syserr_atlas[i]);
		sysgr_atlas->SetPointEYhigh(i,syserr_atlas[i]);
	}
	sysgr_atlas->SetMarkerColor(gr_atlas->GetMarkerColor());
	sysgr_atlas->SetLineColor(gr_atlas->GetLineColor());
	sysgr_atlas->SetMarkerStyle(gr_atlas->GetMarkerStyle());
	sysgr_atlas->SetMarkerSize(gr_atlas->GetMarkerSize());

	TGraphErrors *gr_atlas_jet = new TGraphErrors();
	gr_atlas_jet->SetName("atlas_jet");
	//gr_atlas_jet->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	gr_atlas_jet->SetTitle("ATLAS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_atlas_jet; i++) {
		gr_atlas_jet->SetPoint(i,ene_atlas_jet[i],NtrkDen_atlas_jet[i]);
		gr_atlas_jet->SetPointError(i,0,err_atlas_jet[i]);
	}
	gr_atlas_jet->SetMarkerColor(color_atlas);
	gr_atlas_jet->SetLineColor(color_atlas);
	gr_atlas_jet->SetMarkerStyle(25);
	gr_atlas_jet->SetMarkerSize(2);

	TGraphAsymmErrors *sysgr_atlas_jet = new TGraphAsymmErrors();
	sysgr_atlas_jet->SetName("atlas_jet_sys");
	//sysgr_atlas_jet->SetTitle("ATLAS Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	sysgr_atlas_jet->SetTitle("ATLAS Transverse #LTdN_{ch}/d#etad#phi#GT");
	for(int i = 0; i<N_atlas_jet; i++) {
		sysgr_atlas_jet->SetPoint(i,ene_atlas_jet[i],NtrkDen_atlas_jet[i]);
		sysgr_atlas_jet->SetPointEYlow(i,syserr_atlas_jet[i]);
		sysgr_atlas_jet->SetPointEYhigh(i,syserr_atlas_jet[i]);
	}
	sysgr_atlas_jet->SetMarkerColor(gr_atlas_jet->GetMarkerColor());
	sysgr_atlas_jet->SetLineColor(gr_atlas_jet->GetLineColor());
	sysgr_atlas_jet->SetMarkerStyle(gr_atlas_jet->GetMarkerStyle());
	sysgr_atlas_jet->SetMarkerSize(gr_atlas_jet->GetMarkerSize());




	//TH1D *h = new TH1D("h","",100000,90,12000);
	TH1D *h = new TH1D("h","",100000,8e1,3e4);
	h->SetMinimum(0);
	//h->SetMaximum(1.1);
	h->SetMaximum(1.3);
	h->GetYaxis()->SetNdivisions(504);
	h->GetXaxis()->SetTitle("Center-of-Mass Energy (GeV)");
	//h->GetYaxis()->SetTitle("Transverse #LTN_{ch}/#delta#eta#delta#phi#GT");
	h->GetYaxis()->SetTitle("Transverse #LTdN_{ch}/d#etad#phi#GT");
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.2);

	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.14);
	c->GetFrame()->SetLineWidth(2);

	c->SetLogx();

	h->Draw();
	gr_cdf->Draw("psame");
	gr_cdf_chargejet->Draw("psame");
	gr_cms_chargejet->Draw("psame");
	gr_atlas->Draw("psame");
	gr_atlas_jet->Draw("psame");
	gr_alice->Draw("psame");
	//gr_star->Draw("psame");
	gr_star_jet->Draw("psame");
	//sysgr_star->Draw("[]same");
	sysgr_star_jet->Draw("[]same");
	sysgr_cdf->Draw("[]same");
	sysgr_cdf_chargejet->Draw("[]same");
	sysgr_cms_chargejet->Draw("[]same");
	sysgr_atlas->Draw("[]same");
	sysgr_atlas_jet->Draw("[]same");


	//TLegend *leg = new TLegend(0.15,0.55,0.5,0.85);
	TLegend *leg = new TLegend(0.13,0.52,0.52,0.87);
	leg->AddEntry(gr_star,"STAR Full Jet 20-25 GeV/#it{c} |#eta|<1","p");
	leg->AddEntry(gr_cdf,"CDF Track 5-6 GeV/#it{c} |#eta|<0.8","p");
	leg->AddEntry(gr_cdf_chargejet,"CDF Charged Jet 19-20 GeV/#it{c} |#eta|<1","p");
	leg->AddEntry(gr_alice,"ALICE Track 5-6 GeV/#it{c} |#eta|<0.8","p");
	leg->AddEntry(gr_cms_chargejet,"CMS Charged Jet 17-22 GeV/#it{c} |#eta|<2","p");
	leg->AddEntry(gr_atlas,"ATLAS Track 5-5.5 GeV/#it{c} |#eta|<2.5","p");
	leg->AddEntry(gr_atlas_jet,"ATLAS Full Jet 20-30 GeV/#it{c} |#eta|<2.5","p");
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->Draw();


	cout<<"STAR data 200GeV:"<<endl;
	cout<<"MaxTrack pT 5-7 GeV/c bin "<<endl;
	cout<<"  Tran Density = "<<NtrkDen_star[0]<<" +- "<<err_star[0]<<"(stat) + "<<syserrmax_star[0]<<" - "<<syserrmin_star[0]<<"(sys)"<<endl;
	cout<<"LeadJet pT 11-15 GeV/c bin "<<endl;
	cout<<"  Tran Density = "<<NtrkDen_star_jet12[0]<<" +- "<<err_star_jet12[0]<<"(stat) + "<<syserrmax_star_jet12[0]<<" - "<<syserrmin_star_jet12[0]<<"(sys)"<<endl;
	cout<<"LeadJet pT 20-25 GeV/c bin "<<endl;
	cout<<"  Tran Density = "<<NtrkDen_star_jet20[0]<<" +- "<<err_star_jet20[0]<<"(stat) + "<<syserrmax_star_jet20[0]<<" - "<<syserrmin_star_jet20[0]<<"(sys)"<<endl;

}

void runall();	// plot all matrix figures


void PlotUnfoldMatrix(TString filetag="TranTotNtrk") {

//	TString filetag = "TranTotNtrk";//"TranPtAve"; //"LeadPtAve"; //"SubPtAve"; //"SubAreaNtrk"; //"LeadAreaNtrk"; //
	TFile *f = new TFile("Unfolding_"+filetag+"JPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root");

	TH2D *hx = (TH2D*)f->Get("hXRcVsMc");
	TH2D *hy = (TH2D*)f->Get("hYRcVsMc");

	TString xtitle = "leading jet p_{T} (GeV/#it{c})";
	TString ytitle = "Transverse N_{ch}"; //"Transverse #LTp_{T}#GT"; //"Toward #LTp_{T}#GT"; //"Away #LTp_{T}#GT"; //"Away N_{ch}"; // "Toward N_{ch}";// 
	if(filetag.EqualTo("TranPtAve")) ytitle = "Transverse #LTp_{T}#GT";
	if(filetag.EqualTo("LeadPtAve")) ytitle = "Toward #LTp_{T}#GT";
	if(filetag.EqualTo("SubPtAve")) ytitle = "Away #LTp_{T}#GT";
	if(filetag.EqualTo("SubAreaNtrk")) ytitle = "Away N_{ch}";
	if(filetag.EqualTo("LeadAreaNtrk")) ytitle = "Toward N_{ch}";

	hx->GetXaxis()->SetTitle("Generator-level "+xtitle);
	//hx->GetXaxis()->SetTitle("Particle-level "+xtitle);
	hx->GetYaxis()->SetTitle("Detector-level "+xtitle);

	hy->GetXaxis()->SetTitle("Generator-level "+ytitle);
	//hy->GetXaxis()->SetTitle("Particle-level "+ytitle);
	hy->GetYaxis()->SetTitle("Detector-level "+ytitle);


	hx->GetXaxis()->SetTitleSize(0.05);
	hx->GetYaxis()->SetTitleSize(0.05);
	hx->GetXaxis()->SetLabelSize(0.05);
	hx->GetYaxis()->SetLabelSize(0.05);
	hx->GetZaxis()->SetLabelSize(0.05);
	hx->GetZaxis()->SetNdivisions(5);

	hy->GetXaxis()->SetTitleSize(0.05);
	hy->GetYaxis()->SetTitleSize(0.05);
	hy->GetXaxis()->SetLabelSize(0.05);
	hy->GetYaxis()->SetLabelSize(0.05);
	hy->GetZaxis()->SetLabelSize(0.05);
	hy->GetZaxis()->SetNdivisions(5);


	hx->SetMaximum(1);
	hx->SetMinimum(1e-11);

	hy->SetMaximum(1);
	hy->SetMinimum(1e-11); //1e-7);


	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(8);

	TCanvas *cx = new TCanvas();
	cx->SetLogz();
	cx->SetLeftMargin(0.12);
	cx->SetRightMargin(0.13);
	cx->SetBottomMargin(0.12);

	hx->Draw("colz");

	double lx = 0.01;
	double ly = 0.95;
	TLatex *la = new TLatex(lx,ly,"(a)");
	la->SetNDC();
	la->SetTextFont(42);
	la->Draw();


	TCanvas *cy = new TCanvas();
	cy->SetLogz();
	cy->SetLeftMargin(0.12);
	cy->SetRightMargin(0.13);
	cy->SetBottomMargin(0.12);

	double ymax = 20;//35; // 20;
	if(filetag.EqualTo("SubAreaNtrk")) ymax=35;
	if(filetag.EqualTo("LeadAreaNtrk")) ymax=35;
	hy->GetXaxis()->SetRangeUser(-1,ymax);
	hy->GetYaxis()->SetRangeUser(-1,ymax);

	if(filetag.Contains("PtAve")) {
		cy->SetLogx();
		cy->SetLogy();
		hy->SetMinimum(1e-7);
	}

	hy->Draw("colz");

	TLatex *lb = new TLatex(lx,ly,"(b)");
	lb->SetNDC();
	lb->SetTextFont(42);
	lb->Draw();



	if( 1 && filetag.EqualTo("TranTotNtrk")) {

		//cx->SaveAs("/Users/li/Research/Underlying/PaperDraft170405/ResponseMatrixX.pdf");
		//cy->SaveAs("/Users/li/Research/Underlying/PaperDraft170405/ResponseMatrixY.pdf");

		//cx->SaveAs("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/ResponseMatrixX.pdf");
		//cy->SaveAs("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/ResponseMatrixY.pdf");
		cx->SaveAs("fig/ResponseMatrixX.pdf");
		cy->SaveAs("fig/ResponseMatrixY.pdf");

	}


	if( 1 && !filetag.EqualTo("TranTotNtrk")) {
		TString outtag = "TransPtAve"; //"LeadPtAve";//"AwayPtAve"; // "AwayNtrk"; //"LeadNtrk";
		if(filetag.EqualTo("LeadPtAve")) outtag = "LeadPtAve";
		if(filetag.EqualTo("SubPtAve")) outtag = "AwayPtAve";
		if(filetag.EqualTo("SubAreaNtrk")) outtag = "AwayNtrk";
		if(filetag.EqualTo("LeadAreaNtrk")) outtag = "LeadNtrk";
		//cy->SaveAs("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/ResponseM_"+outtag+".pdf");
		cy->SaveAs("fig/ResponseM_"+outtag+".pdf");
	}

}



void runall() {
	PlotUnfoldMatrix("TranTotNtrk");
	PlotUnfoldMatrix("TranPtAve");
	PlotUnfoldMatrix("LeadAreaNtrk");
	PlotUnfoldMatrix("LeadPtAve");
	PlotUnfoldMatrix("SubAreaNtrk");
	PlotUnfoldMatrix("SubPtAve");
}


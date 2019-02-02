//=======================================================================
//
//      2013.06.03      Li Yi
//      Plot 1D graph for paper 2- and 4-particle
//      flow and nonflow
//      plot 1D figures in the same format for \Delta\V{2}, \Delta\V{4} projections, and \delta for running integral
//
//=======================================================================
//
//      2013.06.11  note about \sigma'
//      \sigma'(\Delta\eta) = k (2 - \Delta\eta)
//                          = k (2 - |\eta_\alpha - \eta_\beta|)
//      \sigma'(\eta_\alpha) = Integral( \sigma'(\Delta\eta), (\eta_\beta, -1, 1) ) /2
//                           = ( 3 k - k (\eta_\alpha)^2 ) /2
//      <\sigma'(\eta_\alpha)> = Integral( \sigma'(\eta_\alpha) , (\eta_\alpha, -1, 1) )/2
//                             = 4/3 k
//
//=============================================================================
//
//      2013.11.18  Li YI
//      treat fitting error on D as syst. err as the suggestion from GPC
//
//=============================================================================


#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision


#include "FigFormat1D.h"
#include "FigColor.h"          // 2014.05.02 ly change light color to dark one
#include "Utility.C"          // 2014.05.02 ly change light color to dark one
#include "/home/li/commonmacro/ClassSysErr.C"

const int savefig = 1;

void plotV2vsDeta(int vn, char symbol[]="",char centtag[]="5");             // V{2} vs \Delta\eta
void plotV4(int vn, char symbol[]="", char centtag[]="5") ;               // average \Delta V{4} vs \Delta(\Delta\eta)
// centtag could be 012, 3, 4, 5, 678
// 2014.05.02 ly void plotV2(int vn, double ymax=8., char symbol[]="", char centtag[]="5");                // \Delta \deta{2} vs \Delta(\Delta\eta)
void plotV2(int vn, double ymin = -2, double ymax=8., char symbol[]="", char centtag[]="5");                // \Delta \deta{2} vs \Delta(\Delta\eta)     // 2014.05.02 ly
void plotetaintegral(char centtag[]="5");       // V{2}, <v^2>, V{4} vs eta
void plotdelta(char centtag[]="5");                   // <delta> vs eta-gap
void plotdeltasigma(char centtag[]="5");              // <delta>, <sigma^2> vs eta-gap
void plotcent(double etagap=0.7);                    // centrality dependence
void plotDcent(double etagap=0.7,char symbol[]="", int vn=2);                    // Sqrt(D-sigma') centrality dependence

const char pttext[100] = "0.15<p_{T}<2 GeV/#it{c}";

const int Nocent = 5;
const char centtext[Nocent][20] = {"50-80%","40-50%","30-40%","20-30%","0-20%"};
char *centtag2text(char centtag[]) {
    switch (centtag)
    {
        case "012":
            return centtext[0];
        case "3":
            return centtext[1];
        case "4":
            return centtext[2];
        case "5":
            return centtext[3];
        case "678":
            return centtext[4];
        default:
            cout<<"ERR!! Wrong centtag[] for centtag2text() func."<<endl;
            return "";
    }
};


void plot1D() {

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);


//    plotV4(2,"(c)");
//    plotV2(2,-2,8,"(a)");
//    plotV2(3,-1.6 ,6,"(b)");
    plotetaintegral();
//    plotdelta();
//    plotdeltasigma();
//    plotDcent(0.7,"(a)",2);
//    plotDcent(0.7,"(b)",3);

}




//=======================================================================================
void plotV2vsDeta(int vn, char symbol[],char centtag[]) {
    
    // read the data
    char filename[500] = "fit1cent%seta20dca3_3dataset.130530";
    char histname[100] = "hV%d2";
    char outname[500] = "";
    sprintf(outname,"V%d2vsDeta.%s",vn,Form(filename,centtag));
    TFile *f;
    f = new TFile(Form("%s.root",Form(filename,centtag)));
    TH2D *hV = (TH2D*)f->Get(Form(histname,vn));

    // Project hV onto \Delta\eta = \eta_1 - \eta_2
    const int Neta = 20;
    double deta = 0.1;                  // histogram bin size
    TGraphErrors *gr[2*Neta];
    int color=2;
    for(int i = 0; i<2*Neta; i++) { 
        gr[i] = new TGraphErrors();
        gr[i]->SetTitle(Form("grV%d2vsDeta",vn)); 
        gr[i]->SetMarkerStyle(8);
        gr[i]->SetMarkerSize(2);
        gr[i]->SetMarkerColor(2);
        gr[i]->SetLineColor(2);
    }
    int Npoint[2*Neta] = {0};           // TGraph gr number of points
    for(int i = 0; i<Neta; i++) { 
        for(int j = 0; j<Neta; j++) {
            gr[i+j]->SetPoint(Npoint[i+j],fabs((i-j)*deta), hV->GetBinContent(i+1,j+1));
            gr[i+j]->SetPointError(Npoint[i+j],0, hV->GetBinError(i+1,j+1));
            Npoint[i+j]++;
        }
    }


    // setup enviroment for drawing
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

    TH2D *htmp = new TH2D("htmp","",100,-0.1,2.1,100,38, 52);    // x-range, y-range go here
    htmp->GetXaxis()->SetTitle("| #eta_{#alpha}-#eta_{#beta}|");
    htmp->GetYaxis()->SetTitle(Form("V_{%d}{2}",vn));

    FigFormat1D(c,pad,htmp);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    for(int j = 0; j<Neta*2 ; j++) {
        if(gr[j]->GetN()>0) {
            gr[j]->Draw("psame");
        }
    }

    wLatex(0.05,0.91,symbol,1,12,ZTITLESIZE);


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
//        c->SaveAs(Form("figures/%s.gif",outname));
    }

        
    
    
    
}
//=======================================================================================
void plotV4(int vn, char symbol[],char centtag[]) {

    // read the data
    const int Nofile = 4;
    char filename[Nofile][500] = {"fit1cent%seta10dca3_3dataset.130603","fit1cent%seta10dca2_3dataset.130601","fit1cent%seta10dca3hit15_3dataset.130605","fit1cent%seta10dca3Vz25_3dataset.130603"};
    char histname[100] = "";
    sprintf(histname,"gavedd%d4",vn);
    char outname[500] = "";
    sprintf(outname,"%s.%s",histname,Form(filename[0],centtag));
    TFile *f[Nofile];
    TGraphErrors *gr[Nofile];
    for(int i = 0; i<Nofile; i++) {
        f[i] = new TFile(Form("%s.root",Form(filename[i],centtag)));
        gr[i] = (TGraphErrors*)f[i]->Get(histname);
    }
    ClassSysErr *syserr = new ClassSysErr(Nofile-1, gr[0],gr+1, "m");          // sys. err. from multi source


    // setup enviroment for drawing
    
    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-3,7);    // x-range, y-range go here
//    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-2,8);    // x-range, y-range go here
//    gStyle->SetOptFit(1);                            // This is to plot every centrality, need fit parameter
//    gStyle->SetOptStat(0);                           // This is to plot every centrality                 
    htmp->GetXaxis()->SetTitle("#Delta#eta_{2}-#Delta#eta_{1}");
// show y axis title as \Delta V{m} as GPC 2nd comments asked 2013.11.24 ly
//    htmp->GetYaxis()->SetTitle(Form("#sigma_{%d}'(#Delta#eta_{1})-#sigma_{%d}'(#Delta#eta_{2})",vn,vn,vn));
//    it turns out to be a misunderstanding. It has been asked to add those vn back in title. 2013.12.16 ly
    htmp->GetYaxis()->SetTitle(Form("#DeltaV_{%d}^{#scale[0.8]{#frac{1}{2}}}{4}",vn));
    
    FigFormat1D(c,pad,htmp);
    htmp->GetYaxis()->SetNdivisions(205);           // 2014.05.02. ly upon comments of NPI group

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    wLatex(0.05,0.91,symbol,1,12,ZTITLESIZE);

    wLatex(0.3,0.29,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality

    TLine *line = new TLine(0,0,2,0);
//    graphSystBox(syserr->GetN(),syserr->GetX(),syserr->GetY(),0,gr[0]->GetEY(),syserr->GetEYhigh(),syserr->GetEYlow(),5,1);
// 2014.05.02 ly    graphSystBand(syserr->GetN(),syserr->GetX(),syserr->GetY(),0,gr[0]->GetEY(),syserr->GetEYhigh(),syserr->GetEYlow(),kGray);
    graphSystBand(syserr->GetN(),syserr->GetX(),syserr->GetY(),0,gr[0]->GetEY(),syserr->GetEYhigh(),syserr->GetEYlow(),kGray,20,1,1,1,1001);           // 2014.05.02 ly NPI comments
    line->Draw("same");
    gr[0]->SetMarkerSize(1.5);//2);           // 2014.05.02 ly NPI comments
    gr[0]->SetMarkerColor(1);                   // 2014.05.06 ly fqwang comments
    gr[0]->SetLineColor(1);                   // 2014.05.06 ly fqwang comments
    gr[0]->GetFunction("f0D")->SetLineStyle(7);
    gr[0]->GetFunction("f0D")->SetLineColor(1);     // 2014.05.06 ly fqwang comments
    gr[0]->Draw("psame");


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
//        c->SaveAs(Form("figs/%s.gif",outname));
    }


    // print fit parameters to screen 2013.11.24  ly
    cout<<"\DeltaV{4} = k*x"<<endl;
    cout<<"k = "<<gr[0]->GetFunction("f0D")->GetParameter(0)<<" +- "<<gr[0]->GetFunction("f0D")->GetParError(0)<<"("<<gr[0]->GetFunction("f0D")->GetParError(0)/gr[0]->GetFunction("f0D")->GetParameter(0)<<")"<<endl;
    cout<<"Chi^2/ndf = "<<gr[0]->GetFunction("f0D")->GetChisquare()<<"/"<<gr[0]->GetFunction("f0D")->GetNDF()<<endl;
    cout<<endl;

}



//=======================================================================================


ClassSysErr* SysErrMultSingleS(TGraphErrors *gr0, int msource, TGraphErrors **grm, int ssource, TGraphErrors *grs) {          // for <v> and \delta, there are both the fits and cuts sys. err.
// different fit function only contribute one source sys. err. sum it together with all other sources.
    ClassSysErr *syserrfit, *syserrcut;
    syserrfit = new ClassSysErr(ssource, gr0, grs,"s");
    syserrcut = new ClassSysErr(msource, gr0, grm,"m");
    double *x, *y, *efitl, *efith, *ecut;
    x = gr0->GetX();
    y = gr0->GetY();
    efitl = syserrfit->GetEYlow();
    efith = syserrfit->GetEYhigh();
    ecut = syserrcut->GetEYhigh();          // for 'm' the low and high err. are the same abs
    const int MAXPOINT = 1024;
    if(MAXPOINT<gr0->GetN()) { cout<<"increase MAXPOINT > "<<gr0->GetN()<<" in function SysErrMultSingleS()"<<endl;}
    double xx[MAXPOINT] = {0}, y0[MAXPOINT]  = {0}, eyfitl[MAXPOINT] = {0}, eyfith[MAXPOINT] = {0}, eycut[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
    for(int i = 0; i<gr0->GetN(); i++) {
        xx[i] = x[i];
        y0[i] = y[i];
        eyfitl[i] = efitl[i];
        eyfith[i] = efith[i];
        eycut[i] = ecut[i];
        eyl[i] = -sqrt(pow(efitl[i],2)+pow(ecut[i],2));
        eyh[i] = sqrt(pow(efith[i],2)+pow(ecut[i],2));
    }

    ClassSysErr *syserr = new ClassSysErr(gr0->GetN(), xx, y0, eyl, eyh);
    
    return syserr;

}



ClassSysErr* SysErrMultSingleS(TGraphAsymmErrors *gr0, int msource, TGraphAsymmErrors *grm, int ssource, TGraphAsymmErrors *grs) {          // for <v> and \delta, there are both the fits and cuts sys. err.
// different fit function only contribute one source sys. err. sum it together with all other sources.
    ClassSysErr *syserrfit, *syserrcut;
    syserrfit = new ClassSysErr(ssource, gr0, grs,"s");
    syserrcut = new ClassSysErr(msource, gr0, grm,"m");
    double *x, *y, *efitl, *efith, *ecut;
    x = gr0->GetX();
    y = gr0->GetY();
    efitl = syserrfit->GetEYlow();
    efith = syserrfit->GetEYhigh();
    ecut = syserrcut->GetEYhigh();          // for 'm' the low and high err. are the same abs
    const int MAXPOINT = 1024;
    if(MAXPOINT<gr0->GetN()) { cout<<"increase MAXPOINT > "<<gr0->GetN()<<" in function SysErrMultSingleS()"<<endl;}
    double xx[MAXPOINT] = {0}, y0[MAXPOINT]  = {0}, eyfitl[MAXPOINT] = {0}, eyfith[MAXPOINT] = {0}, eycut[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
    for(int i = 0; i<gr0->GetN(); i++) {
        xx[i] = x[i];
        y0[i] = y[i];
        eyfitl[i] = efitl[i];
        eyfith[i] = efith[i];
        eycut[i] = ecut[i];
        eyl[i] = -sqrt(pow(efitl[i],2)+pow(ecut[i],2));
        eyh[i] = sqrt(pow(efith[i],2)+pow(ecut[i],2));
    }

    ClassSysErr *syserr = new ClassSysErr(gr0->GetN(), xx, y0, eyl, eyh);
    
    return syserr;

}

ClassSysErr *SumTwoSSysErr(TGraphErrors *gr, double *e1l, double *e1h, double *e2l, double *e2h) {

    const int MAXPOINT = 1024;
    if(MAXPOINT<gr->GetN()) { cout<<"increase MAXPOINT > "<<gr->GetN()<<" in function SysErrMultSingleS()"<<endl;}
    double x[MAXPOINT] = {0}, y[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
    double *xx, *yy; 
    xx = gr->GetX();
    yy = gr->GetY();
    for(int i = 0; i<gr->GetN(); i++) {
	    x[i] = xx[i];
	    y[i] = yy[i];
	    eyl[i] = -sqrt( pow(e1l[i],2) + pow(e2l[i],2) );
	    eyh[i] = sqrt( pow(e1h[i],2) + pow(e2h[i],2) );
    }

    ClassSysErr *syserr = new ClassSysErr(gr->GetN(), x, y, eyl, eyh);

	return syserr;

}


ClassSysErr *SumTwoSSysErr(TGraphAsymmErrors *gr, double *e1l, double *e1h, double *e2l, double *e2h) {

    const int MAXPOINT = 1024;
    if(MAXPOINT<gr->GetN()) { cout<<"increase MAXPOINT > "<<gr->GetN()<<" in function SysErrMultSingleS()"<<endl;}
    double x[MAXPOINT] = {0}, y[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
	double *xx, *yy; 
	xx = gr->GetX();
	yy = gr->GetY();
    for(int i = 0; i<gr->GetN(); i++) {
		x[i] = xx[i];
		y[i] = yy[i];
		eyl[i] = -sqrt( pow(e1l[i],2) + pow(e2l[i],2) );
		eyh[i] = sqrt( pow(e1h[i],2) + pow(e2h[i],2) );
	}

    ClassSysErr *syserr = new ClassSysErr(gr->GetN(), x, y, eyl, eyh);

	return syserr;

}


void plotV2(int vn, double ymin, double ymax, char symbol[],char centtag[]) {                // \Delta \deta{2} vs \Delta(\Delta\eta)

    // read the data
    const int Nofile = 4; 
//    char filename[Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530","fit1cent5eta20dca2_3dataset.130506","fit1cent5eta20dca2hit15_3dataset.130531","fit1cent5eta20dca3Vz25_3dataset.130530"};
    char filename[Nofile][1024] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3hit15_3dataset.130328","fit1cent%seta20dca3Vz25_3dataset.130530"};
    const int Nograph = 4;
    char graphtag[100] = "g3_dd%d2_%d";                 // data points  + linear fit
    char functag[100] = "f2D1Dcp%d2%d";                 // 2D fit 
    char histtag[100] = "fitresult_dd%d2";              // fit parameters results
    char outname[500] = "";
    sprintf(outname,"Ddprojfitv%d.%s",vn,Form(filename[0],centtag));
    TFile *f[Nofile];
    TGraphErrors *gr[Nograph][Nofile];
    for(int i = 0; i<Nofile ; i++) {
        f[i] = new TFile(Form("%s.root",Form(filename[i],centtag)));
        for(int j = 0; j<Nograph ; j++) {
            gr[j][i] = (TGraphErrors*)f[i]->Get(Form(graphtag,vn,j+1));
        }
    }
    
    ClassSysErr *syserr[Nograph];
    for(int j = 0; j<Nograph ; j++) {
        syserr[j] = new ClassSysErr(Nofile-1, gr[j][0], gr[j]+1, "m");        // sys. err. from mult-source
    }


    TF1 *hf[Nograph];
    f[0]->cd();
    for(int i = 0; i<Nograph ; i++) {
        hf[i] = (TF1*)f[0]->Get(Form(functag,vn,i+1));
    }

    TH1D *hfitpar = (TH1D*)f[0]->Get(Form(histtag,vn));
    double Chi = hfitpar->GetBinContent(5);
    int ndf = (int*)hfitpar->GetBinContent(6);

    

    // setup enviroment for drawing
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

// 2014.05.02 ly    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-2,ymax);    // x-range, y-range go here
    TH2D *htmp = new TH2D("htmp","",100,0,2,100,ymin,ymax);    // x-range, y-range go here            // 2014.05.02 ly
    htmp->GetXaxis()->SetTitle("#Delta#eta_{2}-#Delta#eta_{1}");
// show y axis title as \Delta V{m} as suggested by GPC 2nd comments 2013.11.24 ly
//    htmp->GetYaxis()->SetTitle(Form("D_{%d}(#Delta#eta_{1})-D_{%d}(#Delta#eta_{2})",vn,vn,vn));
//    it turns out to be a misunderstanding. It has been asked to add those vn back in title. 2013.12.16 ly
    htmp->GetYaxis()->SetTitle(Form("#DeltaV_{%d}{2}",vn));
    
    FigFormat1D(c,pad,htmp);
    htmp->GetYaxis()->SetNdivisions(205);           // 2014.05.02. ly upon comments of NPI group

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    TLine *line = new TLine(0,0,2,0);
// 2014.05.02 ly    line->Draw("same");
    for(int j = 0; j<Nograph ; j++) {
//        graphSystBox(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,gr[j][0]->GetEY(),syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),kGray,1);
// 2014.05.02 ly        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,gr[j][0]->GetEY(),syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),kGray);
        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,gr[j][0]->GetEY(),syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),kGray+1,0,1,1,1,1001);      // 2014.05.02 ly NPI comments
        line->Draw("same");                                 // 2014.05.02 ly NPI comments
        gr[j][0]->SetMarkerSize(1.5); // 1);               // 2014.05.02 ly NPI comments
        gr[j][0]->SetMarkerColor(FigColor(j+1));        // 2014.05.02 ly NPI comments
        gr[j][0]->SetLineColor(FigColor(j+1));          // 2014.05.02 ly NPI comments
        gr[j][0]->GetFunction("f1D")->SetLineStyle(7);          // 2014.05.02 ly NPI comments
//        gr[j][0]->GetFunction("f1D")->SetLineWidth(2);          // 2014.05.02 ly NPI comments
        gr[j][0]->GetFunction("f1D")->SetLineColor(FigColor(j+1));          // 2014.05.02 ly NPI comments
        hf[j]->SetMarkerColor(FigColor(j+1));           // 2014.05.02 ly NPI comments
        hf[j]->SetLineColor(FigColor(j+1));             // 2014.05.02 ly NPI comments
        hf[j]->SetLineWidth(3);             // 2014.05.02 ly NPI comments
        gr[j][0]->Draw("psame");
        hf[j]->Draw("same");
    }

    wLatex(0.3,0.8,Form("#chi^{2}/ndf = %.1f/%d",Chi,ndf));

    wLatex(0.05,0.91,symbol,1,12,ZTITLESIZE);

    wLatex(0.3,0.29,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality


    double ytxt = 0.25; // 0.3;     // 2014.05.02 ly NPI comments
    double markersize = 1.5; //1;           // 2014.05.02 ly NPI comments
    double tsize = 0.05;//0.04;     //2014.05.02 ly upon NPI group comments
    int Neta = 20;
    double etabin = 2./Neta;
    for(int j = 0; j<Nograph ; j++) {
// 2014.05.02 ly        keySymbol(0.75,ytxt,Form("#Delta#eta_{1}=%3.1f",(j+1)*etabin),j+1,20+j,tsize,markersize);
        keySymbol(0.7,ytxt,Form("#Delta#eta_{1}=%3.1f",(j+1)*etabin),FigColor(j+1),20+j,tsize,markersize);    //2014.05.02 ly upon NPI group comments
        ytxt+=0.07;
    }

    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
//        c->SaveAs(Form("figures/%s.gif",outname));
    }


        
    
// print fit result to screen 2013.11.24 ly	
// assume fit1
    const int Nofitpar = 4;
    char fitparname[Nofitpar][10] = {"a","b","A","#sigma"};
    cout<<"\Delta\V{2} = a*exp(-x/b) + A*exp(-x^2/(2*\sigma^2))"<<endl;
    for(int i = 0; i<4; i++) {
	    cout<<fitparname[i]<<" = "<<hfitpar->GetBinContent(i+1)<<" +- "<<hfitpar->GetBinError(i+1)<<endl;
    }
    cout<<"#Chi^2/ndf = "<<hfitpar->GetBinContent(5)<<"/"<<hfitpar->GetBinContent(6)<<endl;
    cout<<endl;



}



//=======================================================================================
double IntegralMean(int N, double *y) {
    if(N<=0) return 0;
    double sum = 0; 
    for(int i = 0; i<N ; i++) {
        sum+=y[i];
    }
    return sum/N;
}

double IntegralError(int N, double *ey) {
    if(N<=0) return 0;
    double sum = 0; 
    for(int i = 0; i<N ; i++) {
        sum+=pow(ey[i],2);
    }
    return sqrt(sum)/N;
}


double Average(int N, double *y) {
    if(N<=0) return 0;
    double sum = 0; 
    for(int i = 0; i<N ; i++) {
        sum+=y[i];
    }
    return sum/N;
}

double ratio4sigma2v2(double vv2, double V4) {              // return sqrt( (<v^2>-V{4})/(<v^2>+V{4} )
    double ratio = 0;
    if(vv2+V4 > 0 && vv2 - V4 > 0) {
        ratio = sqrt( (vv2 - V4) /(vv2+V4) );
    }
    return ratio;
}

double StatErrRatio4sigma2v2(double vv2, double V4, double evv2, double eV4) {      // return the error propagation for  sqrt( (<v^2>-V{4})/(<v^2>+V{4}) )
    double error = 0;
    double value  = ratio4sigma2v2(vv2, V4);
    double dfdx = 2*V4/pow(vv2+V4,2);
    double dfdy = -2*vv2/pow(vv2+V4,2);
    error = pow(dfdx*evv2,2)+pow(dfdy*eV4,2);
    error = sqrt(error);
    error = error/(2.*value);
    return error;
}

void plotetaintegral(char centtag[]) {       // V{2}, <v^2>, V{4} vs eta

    // read the data for V{2}, <v^2>

    const int V2Nofile = 7;
//    char V_22filename[V2Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530","fit1cent5eta20dca2_3dataset.130506","fit1cent5eta20dca3Vz25_3dataset.130530","fit1cent5eta20dca2hit15_3dataset.130531","fit0cent5eta20dca3_3dataset.130530","fit8cent5eta20dca3_3dataset.130530"};              // fit function for v2
    char V_22filename[V2Nofile][1024] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};              // fit function for v2
//    char V_32filename[V2Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530","fit1cent5eta20dca2_3dataset.130506","fit1cent5eta20dca3Vz25_3dataset.130530","fit0cent5eta20dca3_3dataset.130530","fit3cent5eta20dca3_3dataset.130530","fit1cent5eta20dca2hit15_3dataset.130531"};              // fit function for v3, note the order of the files are different that v2, due to the last root file doesn't fit v3 well.
    char V_32filename[V2Nofile][1024] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit2cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};              // fit function for v3
    // read the data for V{4}
    const int V4Nofile = 4;
    char V4filename[V4Nofile][500] = {"fit1cent%seta10dca3_3dataset.130603","fit1cent%seta10dca2_3dataset.130601","fit1cent%seta10dca3Vz25_3dataset.130603","fit1cent%seta10dca3hit15_3dataset.130605"};             // files are in the same order as V2{2} for later fluct/flow ratio calculation

    char outname2[500] = "";
    char outname3[500] = "";
    sprintf(outname2,"V2etadep.%s.%s",Form(V_22filename[0],centtag),Form(V4filename[0],centtag));
    sprintf(outname3,"V3etadep.%s.%s",Form(V_32filename[0],centtag),Form(V4filename[0],centtag));

    const int NoGraph = 6+2;                                 // V2{2}, <v2^2>, V2{4}, V3{2}, <v3^2>, V3{4}, <d2>, <d3>
    int Nofile[NoGraph];                                  
    Nofile[0] = V2Nofile;
    Nofile[1] = V2Nofile;
    Nofile[3] = V2Nofile;
    Nofile[4] = V2Nofile;
    Nofile[2] = V4Nofile;
    Nofile[5] = V4Nofile;
    Nofile[6] = V2Nofile;
    Nofile[7] = V2Nofile;
//    for(int i = 0; i<NoGraph ;i++) { Nofile[i] = 4; }       // 4 for different events, tracks cuts
//    Nofile[1] = 3;                                  // <v2^2> only  consider different fit function here they are large contributor                                      //V2Nofile;                                          // <v2^2> has fit diff. func. err.
//    Nofile[2] = V4Nofile;                                          //  V2{4} need err. from \sigma' with 1sigma variance, it is calculated from the default file \Delta V{4} fit parameters. This effect will be plot as band later
//    Nofile[4] = 3;                                          //  <v3^2> fit : exp, exp+linear, exp+gaus. 
    TH1D *hV[NoGraph][V2Nofile+V4Nofile];                  // read in as hist  , array size is about 2x larger than real needed
    TGraphErrors *gV[NoGraph][V2Nofile+V4Nofile];          // process as graph  , array size is  about 2x larger than real needed
    int graphcolor[NoGraph] = {634,419,603,634,419,603};                   // should keep the same with plotlego.C            (634 for kRed+2, 603 for kBlue+3, 419 for kGreen+3)
//    int graphbandcolor[NoGraph] = {2,920,4,2,400,4};                   // should keep the same with plotlego.C             (920 for kGray. 400 for kYellow)
// 2014.05.02 ly    int graphbandcolor[NoGraph] = {625,920,593,625,393,593};                   // should keep the same with plotlego.C             (625 for kRed-7, 920 for kGray. 400 for kYellow, 393 for kYellow-7, 593 for kBlue-7)
    int graphbandcolor[NoGraph] = {920,920,920,920,920,920};                   // should keep the same with plotlego.C             (625 for kRed-7, 921 for kGray+1. 400 for kYellow, 798 for kOrange-3, 593 for kBlue-7)                 // 2014.05.02 ly
// 2014.05.06 ly    int graphbandcolor[NoGraph] = {625,921,593,625,797,593};                   // should keep the same with plotlego.C             (625 for kRed-7, 921 for kGray+1. 400 for kYellow, 798 for kOrange-3, 593 for kBlue-7)                 // 2014.05.02 ly
// 2014.05.02 ly    Int_t sigmacolor = (Int_t) (kCyan+2);
    Int_t sigmacolor = (Int_t) (kCyan+1);               // 2014.05.02 ly
//    int graphstyle[NoGraph] = {4,8,23,4,8,23};
    int graphstyle[NoGraph] = {21,8,23,21,8,23};				// open circle (4) is hard to see, use solid square (21) 
//-------------------------------------------------2013.11.25 ly-------------------------------------------------------------------------
// GPC 2nd comments asked not to write V_n{m}, instead V{m}. ly
//    char graphtitle[NoGraph][200] = {"V_{2}{2}","#LT v_{2}^{2} #GT","","V_{3}{2}","#LT v_{3}^{2} #GT","V_{3}{4}","#delta_{2}","#delta_{3}"};
//    sprintf(graphtitle[2],Form("V_{2}{4}#color[%d]{+#sigma'}",sigmacolor));
//    it turns out to be a misunderstanding. It has been asked to add those vn back in title. 2013.12.16 ly
//    char graphtitle[NoGraph][200] = {"V{2}","#LT v^{2} #GT","","V{2}","#LT v^{2} #GT","V{4}","#delta","#delta"};     // 2013.12.16 ly
//    sprintf(graphtitle[2],Form("V{4}#color[%d]{+#sigma'}",sigmacolor));       // 2013.12.16 ly
//    char graphtitle[NoGraph][200] = {"V_{2}{2}","#LT #font[12]{v}_{2}^{2} #GT","","V_{3}{2}","#LT #font[12]{v}_{3}^{2} #GT","V_{3}^{#scale[0.8]{#frac{1}{2}}}{4}","#delta_{2}","#delta_{3}"};
    char graphtitle[NoGraph][200] = {"V_{2}{2}","#LT v_{2}^{2} #GT","","V_{3}{2}","#LT v_{3}^{2} #GT","V_{3}^{#scale[0.8]{#frac{1}{2}}}{4}","#delta_{2}","#delta_{3}"};
    sprintf(graphtitle[2],Form("V_{2}^{#scale[0.8]{#frac{1}{2}}}{4}#color[%d]{+#sigma'}",sigmacolor));
//-------------------------------------------------2013.11.25 ly-------------------------------------------------------------------------


    // read in hist.
    TFile *fV2[V2Nofile][2], *fV4[V4Nofile];            // v2 and v3 have different fV2 due to the fit functions
    double scalev3 = 1;
    //TH1D *hV[NoGraph][V2Nofile+V4Nofile];                  // read in as hist  , array size is about 2x larger than real needed
    TH1D *hVtmp;                // Since we are gonna change histograms, it is better to make copy before hand.
    for(int i = 0; i<V2Nofile ; i++) {
        fV2[i][0] = new TFile(Form("%s.root",Form(V_22filename[i],centtag)));
        fV2[i][0]->cd();
//        hV[0][i] = (TH1D*)fV2[i][0]->Get("hproj_V22");
//        hV[1][i] = (TH1D*)fV2[i][0]->Get("hproj_v22sq");
//        hV[6][i] = (TH1D*)fV2[i][0]->Get("hproj_delta22");            // delta histogram
        hVtmp = (TH1D*)fV2[i][0]->Get("hproj_V22");
        hV[0][i] = (TH1D*)hVtmp->Clone("hproj_V22cp");
        hVtmp = (TH1D*)fV2[i][0]->Get("hproj_v22sq");
        hV[1][i] = (TH1D*)hVtmp->Clone("hproj_v22sqcp");
        hVtmp = (TH1D*)fV2[i][0]->Get("hproj_delta22");            // delta histogram
        hV[6][i] = (TH1D*)hVtmp->Clone("hproj_delta22cp");
        fV2[i][1] = new TFile(Form("%s.root",Form(V_32filename[i],centtag)));
        fV2[i][1]->cd();
//        hV[3][i] = (TH1D*)fV2[i][1]->Get("hproj_V32");
//        hV[4][i] = (TH1D*)fV2[i][1]->Get("hproj_v32sq");
//        hV[7][i] = (TH1D*)fV2[i][0]->Get("hproj_delta32");            // delta histogram
        hVtmp = (TH1D*)fV2[i][1]->Get("hproj_V32");
        hV[3][i] = (TH1D*)hVtmp->Clone("hproj_V32cp");
        hVtmp = (TH1D*)fV2[i][1]->Get("hproj_v32sq");
        hV[4][i] = (TH1D*)hVtmp->Clone("hproj_v32sqcp");
        hVtmp = (TH1D*)fV2[i][1]->Get("hproj_delta32");            // delta histogram
        hV[7][i] = (TH1D*)hVtmp->Clone("hproj_delta32cp");
//        for(int ib = 0; ib<hV[3][i]->GetNbinsX(); ib++) {             // no need to scale v3 on v2, they are ploted separately on two graphs
//            hV[3][i]->SetBinContent(ib+1,hV[3][i]->GetBinContent(ib+1)*scalev3);
//        }

    }
    for(int i = 0; i<V4Nofile ; i++) {
        fV4[i] = new TFile(Form("%s.root",Form(V4filename[i],centtag)));
        fV4[i]->cd();
//        hV[2][i] = (TH1D*)fV4[i]->Get("hproj_V24");
//        hV[5][i] = (TH1D*)fV4[i]->Get("hproj_V34");
        hVtmp = (TH1D*)fV4[i]->Get("hproj_V24");
        hV[2][i] = (TH1D*)hVtmp->Clone("hproj_V24cp");
        hVtmp = (TH1D*)fV4[i]->Get("hproj_V34");
        hV[5][i] = (TH1D*)hVtmp->Clone("hproj_V34cp");
    }

//----------------------------------------2013. 11. 25------------------------------------------------------------------------------------------------    
// delta fitting error is treated as sys. err. as asked for by GPC 2nd comments 2013.11.25 ly
//    // <v^2> stat. err. need to include the one from delta fitting error.
//    for(int i = 0; i<V2Nofile ; i++) {
//        double tmp1 = 0, tmp2 = 0;
//        for(int j = 0; j<hV[1][i]->GetNbinsX(); j++) {
//            tmp1 = hV[1][i]->GetBinError(j+1);                              // <v2^2>
//            tmp2 = hV[6][i]->GetBinError(j+1);                              // delta_2
//            hV[1][i]->SetBinError(j+1,sqrt(pow(tmp1,2)+pow(tmp2,2)));       
//            tmp1 = hV[4][i]->GetBinError(j+1);                              // <v3^2>
//            tmp2 = hV[7][i]->GetBinError(j+1);                              // delta_3
//            hV[4][i]->SetBinError(j+1,sqrt(pow(tmp1,2)+pow(tmp2,2)));
//        }
//    }
//----------------------------------------end------------------------------------------------------------------------------------------------    


    // convert to graph
    for(int j = 0; j<NoGraph ; j++) {
        for(int i = 0 ; i<Nofile[j] ; i++) {
            if(hV[j][i]) {
                hV[j][i]->SetMarkerColor(graphcolor[j]);
                hV[j][i]->SetMarkerStyle(graphstyle[j]);
                hV[j][i]->SetMarkerSize(1.5);
                gV[j][i] = new TGraphErrors(hV[j][i]);
                gV[j][i]->SetMarkerColor(graphcolor[j]);
                gV[j][i]->SetMarkerStyle(graphstyle[j]);
                gV[j][i]->SetMarkerSize(1.5);
            }
        }
    }

    //============================== for V2{4} err. from \sigma'. need a graph gV for the value V2{4} + \sigma'{4}
    TGraphErrors *gavedd24 = (TGraphErrors*)fV4[0]->Get("gavedd24");            // used for \sigma' in V{4}
    double k = gavedd24->GetFunction("f0D")->GetParameter(0);
    double ek = gavedd24->GetFunction("f0D")->GetParError(0);
    TF1 *sigmaprime = new TF1("sigmaprime","[0]*(3-x^2)/2.",-1,1);                 // \sigma'
    double *px = gV[2][0]->GetX(), *py = gV[2][0]->GetY();
    gV[2][V4Nofile] = new TGraphErrors();                       // if memeory error message here, please check for V2{4}
    gV[2][V4Nofile+1] = new TGraphErrors();
    for(int i = 0; i<gV[2][0]->GetN(); i++) {            
        sigmaprime->SetParameter(0,k+ek);               // 1-sigma variance
        gV[2][V4Nofile]->SetPoint(i,px[i],py[i] + sigmaprime->Eval(px[i]));     // V2{4} + \sigma'
        sigmaprime->SetParameter(0,k-ek);               // 1-sigma variance
        gV[2][V4Nofile+1]->SetPoint(i,px[i],py[i] + sigmaprime->Eval(px[i]));     // V2{4} + \sigma'
    }

    //==================================   Sys. Err. for All Graphs     =====================
    ClassSysErr *syserr[NoGraph+1];                 // the last one is for V2{4} + \sigma'2 above
///    const int MAXPOINT = 50;
///    if(MAXPOINT<gV[0][0]->GetN()) {cout<<"increase MAXPOINT > "<<gV[0][0]->GetN()<<" in function plotetaintegral() "<<endl;}
///    double av22[4][MAXPOINT] = {{0}}, av32[4][MAXPOINT] = {{0}};        // for v2^2 and v3^2 sys. err. they have two contributor: cuts and fits.          for tmp use
///    double ad22[4][MAXPOINT] = {{0}}, ad32[4][MAXPOINT] = {{0}};        // for d2 and d3 sys. err. they have two contributor: cuts and fits.          for tmp use
///    double *ax, *ay, *aeyl, *aeyh;               // for tmp use
//    for(int j = 0; j<NoGraph ; j++) {
//        syserr[j] = new ClassSysErr(Nofile[j]-1, gV[j][0],gV[j]+1);
//    }
    syserr[0] = new ClassSysErr(3,gV[0][0],gV[0]+1,"m");          // V2{2} 3 other different cuts
    syserr[2] = new ClassSysErr(3,gV[2][0],gV[2]+1,"m");          // V2{4} 3 other different cuts
    syserr[3] = new ClassSysErr(3,gV[3][0],gV[3]+1,"m");          // V3{2} 3 other different cuts
    syserr[5] = new ClassSysErr(3,gV[5][0],gV[5]+1,"m");          // V3{4} 3 other different cuts

//------------------------------------------2013.11.25----------------------------------------------------------------------------------------------    
// <d> sys. err. now include two sources: fitting error from parameters ans different fit function as asked by GPC 2nd comments 2013.11.25 ly
//    syserr[6] = new ClassSysErr(3,gV[6][0],gV[6]+4,"s");				// <d2> only different fit func. err.
//    syserr[7] = new ClassSysErr(3,gV[7][0],gV[7]+4,"s");				// <d3> only different fit func. err.
    ClassSysErr *syserrtmp1[2], *syserrtmp2[2];
    const int MAXPOINT = 1024;
    double tmpgey[MAXPOINT], *tmppy;
    for(int ig = 6; ig<=7; ig++) {							// <d2> and <d3>
	    int igtmp = ig-6;
	    tmppy = gV[ig][0]->GetEY();
	    for(int ip = 0; ip<gV[ig][0]->GetN(); ip++) {
		    tmpgey[ip] = -tmppy[ip];						// ClassSysErr: low err is defined as negative
//		    cout<<igtmp<<"\t"<<ip<<"\t"<<tmpgey[ip]<<"\t"<<tmppy[ip]<<endl;
	    }
	    syserrtmp1[igtmp] = new ClassSysErr(gV[ig][0]->GetN(),gV[ig][0]->GetX(),gV[ig][0]->GetY(),tmpgey,gV[ig][0]->GetEY());																// <d> fitting parameters errors also treated as sys. err. as asked by GPC 2nd comments
	    syserrtmp2[igtmp] = new ClassSysErr(3,gV[ig][0],gV[ig]+4,"s");				// <d> sys. err from  different fit func. err.
	    syserr[ig] = (ClassSysErr*) SumTwoSSysErr(gV[ig][0],syserrtmp1[igtmp]->GetEYlow(), syserrtmp1[igtmp]->GetEYhigh(), syserrtmp2[igtmp]->GetEYlow(),syserrtmp2[igtmp]->GetEYhigh());		// <d> sys. err. = sqrt( <d> fit parameters err. ^2 + <d> different fit function err. ^2 )

//	    // print to screen for debug
//	    cout<<"d"<<ig-4<<"{2} histogram err. bar"<<endl;
//	    for(int i = 0; i<hV[ig][0]->GetNbinsX(); i++) {
//	    	cout<<hV[ig][0]->GetBinError(i+1)<<endl;
//	    }
//	    cout<<endl;
//    	    cout<<"<d"<<ig-4<<"> sys. err. from fit pars. "<<endl;
//    	    syserrtmp1[igtmp]->Print();
//    	    cout<<endl;
//    	    cout<<"<d"<<ig-4<<"> sys. err. from diff. fit func. "<<endl;
//    	    syserrtmp2[igtmp]->Print();
//    	    cout<<endl;
//	    cout<<"<d"<<ig-4<<"> sys. err. sum"<< endl;
//	    syserr[ig]->Print();

}
//--------------------------------------------end--------------------------------------------------------------------------------------------    



    syserr[1] = (ClassSysErr*) SumTwoSSysErr(gV[1][0],syserr[0]->GetEYlow(),syserr[0]->GetEYhigh(),syserr[6]->GetEYlow(),syserr[6]->GetEYhigh());							// <v2^2> sys. err. = sqrt( V2{2} cuts err ^2 + <d2> fit func. err ^2  )
    syserr[4] = (ClassSysErr*) SumTwoSSysErr(gV[4][0],syserr[3]->GetEYlow(),syserr[3]->GetEYhigh(),syserr[7]->GetEYlow(),syserr[7]->GetEYhigh());							// <v2^2> sys. err. = sqrt( V2{2} cuts err ^2 + <d2> fit func. err ^2  )

//-------------------------------------------------------------------------------------------------------------------------------------------
//    // sys. err for <v2^2> and <v3^2> have two sources, consider separately and sum together
//    // sys. err for <d2> and <d3> have two sources, consider separately and sum together
//    syserr[1] = (ClassSysErr*) SysErrMultSingleS(gV[1][0],3,gV[1]+1,2,gV[1]+4);
//    syserr[4] = (ClassSysErr*) SysErrMultSingleS(gV[4][0],3,gV[4]+1,2,gV[4]+4);
//    syserr[6] = (ClassSysErr*) SysErrMultSingleS(gV[6][0],3,gV[6]+1,2,gV[6]+4);
//    syserr[7] = (ClassSysErr*) SysErrMultSingleS(gV[7][0],3,gV[7]+1,2,gV[7]+4);
///    ClassSysErr *syserrfit[4];                      // v2^2 and v3^2 need imediate step for different fit 's'
///    // so do d2 and d3
///    ClassSysErr *syserrcut[4];                      // v2^2 and v3^2 need imediate step for different cut 'm'
///    // so do d2 and d3
///    // <v2^2>
///    syserrfit[0] = new ClassSysErr(2,gV[1][0],gV[1]+4,"s");          // v2^2 2 other different fits
///    ax = syserrfit[0]->GetX();
///    ay = syserrfit[0]->GetY(); 
///    aeyl = syserrfit[0]->GetEYlow(); 
///    aeyh = syserrfit[0]->GetEYhigh(); 
///    for(int i = 0; i<syserrfit[0]->GetN(); i++) {
///        av22[0][i] = ay[i];
///        av22[1][i] = aeyl[i];     // itself is negative
///        av22[2][i] = aeyh[i];
///    }
///    syserrcut[0] = new ClassSysErr(3,gV[1][0],gV[1]+1,"m");          // v2^2 3 other different cuts
///    aeyh = syserrcut[0]->GetEYhigh();        // low and high err are same magnititude for 'm' option
///    for(int i = 0; i<syserrcut[0]->GetN(); i++) {
///        av22[3][i] = aeyh[i];
///    }
///    for(int i = 0; i<syserrcut[0]->GetN(); i++) {               // diff. fit and cut are two sources
///        av22[1][i] = -sqrt(pow(av22[1][i],2) + pow(av22[3][i],2));      // low err.
///        av22[2][i] = sqrt(pow(av22[2][i],2) + pow(av22[3][i],2));       // high err.
///    }
///    syserr[1] = new ClassSysErr(gV[1][0]->GetN(), gV[1][0]->GetX(), av22[0], av22[1], av22[2] );
///
///
///    // <v3^2>
///    syserrfit[1] = new ClassSysErr(2,gV[4][0],gV[4]+4,"s");          // v3^2 2 other different fits
///    ax = syserrfit[1]->GetX();
///    ay = syserrfit[1]->GetY(); 
///    aeyl = syserrfit[1]->GetEYlow(); 
///    aeyh = syserrfit[1]->GetEYhigh(); 
///    for(int i = 0; i<syserrfit[1]->GetN(); i++) {
///        av32[0][i] = ay[i];
///        av32[1][i] = aeyl[i];     // itself is negative
///        av32[2][i] = aeyh[i];
///    }
///    syserrcut[1] = new ClassSysErr(3,gV[4][0],gV[4]+1,"m");          // v3^2 3 other different cuts
///
///    aeyh = syserrcut[1]->GetEYhigh();        // low and high err are same magnititude for 'm' option
///    for(int i = 0; i<syserrcut[1]->GetN(); i++) {
///        av32[3][i] = aeyh[i];
///    }
///    for(int i = 0; i<syserrcut[1]->GetN(); i++) {               // diff. fit and cut are two sources
///        av32[1][i] = -sqrt(pow(av22[1][i],2) + pow(av22[3][i],2));           // low err.
///        av32[2][i] = sqrt(pow(av22[2][i],2) + pow(av22[3][i],2));           // high err.
///    }
///    syserr[5] = new ClassSysErr(gV[4][0]->GetN(), gV[4][0]->GetX(), av32[0], av32[1], av32[2] );
///
///
///    // <d2>
///    syserrfit[2] = new ClassSysErr(2,gV[6][0],gV[6]+4,"s");          // d2 2 other different fits
///    ax = syserrfit[2]->GetX();
///    ay = syserrfit[2]->GetY(); 
///    aeyl = syserrfit[2]->GetEYlow(); 
///    aeyh = syserrfit[2]->GetEYhigh(); 
///    for(int i = 0; i<syserrfit[2]->GetN(); i++) {
///        ad22[0][i] = ay[i];
///        ad22[1][i] = aeyl[i];     // itself is negative
///        ad22[2][i] = aeyh[i];
///    }
///    syserrcut[2] = new ClassSysErr(3,gV[6][0],gV[6]+1,"m");          // d2 3 other different cuts
///    aeyh = syserrcut[2]->GetEYhigh();        // low and high err are same magnititude for 'm' option
///    for(int i = 0; i<syserrcut[2]->GetN(); i++) {
///        ad22[3][i] = aeyh[i];
///    }
///    for(int i = 0; i<syserrcut[2]->GetN(); i++) {               // diff. fit and cut are two sources
///        ad22[1][i] = -sqrt(pow(ad22[1][i],2) + pow(ad22[3][i],2));      // low err.
///        ad22[2][i] = sqrt(pow(ad22[2][i],2) + pow(ad22[3][i],2));       // high err.
///    }
///    syserr[6] = new ClassSysErr(gV[6][0]->GetN(), gV[6][0]->GetX(), ad22[0], ad22[1], ad22[2] );
///
///
///    // <d3>
///    syserrfit[3] = new ClassSysErr(2,gV[7][0],gV[7]+4,"s");          // d3 2 other different fits
///    ax = syserrfit[3]->GetX();
///    ay = syserrfit[3]->GetY(); 
///    aeyl = syserrfit[3]->GetEYlow(); 
///    aeyh = syserrfit[3]->GetEYhigh(); 
///    for(int i = 0; i<syserrfit[3]->GetN(); i++) {
///        ad32[0][i] = ay[i];
///        ad32[1][i] = aeyl[i];     // itself is negative
///        ad32[2][i] = aeyh[i];
///    }
///    syserrcut[3] = new ClassSysErr(3,gV[7][0],gV[7]+1,"m");          // d3 3 other different cuts
///
///    aeyh = syserrcut[3]->GetEYhigh();        // low and high err are same magnititude for 'm' option
///    for(int i = 0; i<syserrcut[3]->GetN(); i++) {
///        ad32[3][i] = aeyh[i];
///    }
///    for(int i = 0; i<syserrcut[3]->GetN(); i++) {               // diff. fit and cut are two sources
///        ad32[1][i] = -sqrt(pow(ad22[1][i],2) + pow(ad22[3][i],2));           // low err.
///        ad32[2][i] = sqrt(pow(ad22[2][i],2) + pow(ad22[3][i],2));           // high err.
///    }
///    syserr[7] = new ClassSysErr(gV[7][0]->GetN(), gV[7][0]->GetX(), ad32[0], ad32[1], ad32[2] );
///
///
//---------------------------------------end----------------------------------------------------------------------------------------------------

    // V2{4}
    syserr[NoGraph] = new ClassSysErr(2,gV[2][0],gV[2]+3,"s");          // V2{4} + \sigma'_{2}   . \sigma' is larger than the difference between other sys. err. for V2{4}

    
    // setup enviroment for drawing
    
    //============  draw v2
    gStyle->SetOptFit(0);                            // This is to plot every centrality, need fit parameter
    gStyle->SetOptStat(0);                           // This is to plot every centrality                 

    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

    double ymin = 25, ymax = 49;
//    double ymin = 10, ymax = 20;
    TH2D *htmp = new TH2D("htmp","",100,-1,1,100,ymin,ymax);    // x-range, y-range go here
    htmp->GetXaxis()->SetTitle("#eta");
    
    FigFormat1D(c,pad,htmp);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

//    graphSystBand(syserr[NoGraph]->GetN(),syserr[NoGraph]->GetX(),syserr[NoGraph]->GetY(),0,0,syserr[NoGraph]->GetEYhigh(),syserr[NoGraph]->GetEYlow(),sigmacolor,0,0,1,1,3545);          // \sigma'2 + d{4}
// 2014.05.02 ly    graphSystBand(syserr[NoGraph]->GetN(),syserr[NoGraph]->GetX(),syserr[NoGraph]->GetY(),0,0,syserr[NoGraph]->GetEYhigh(),syserr[NoGraph]->GetEYlow(),sigmacolor,0,0,1,1,3005);          // \sigma'2 + d{4}
    graphSystBand(syserr[NoGraph]->GetN(),syserr[NoGraph]->GetX(),syserr[NoGraph]->GetY(),0,0,syserr[NoGraph]->GetEYhigh(),syserr[NoGraph]->GetEYlow(),sigmacolor,0,0,1,1,1001);          // \sigma'2 + d{4}                 // 2014.05.02 ly
    for(int j = 0; j<3; j++) {          // v_{2} only
//        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),graphbandcolor[j],0);
        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),graphbandcolor[j],0,1,1,1,1001);					// fill style is hard to see, 1001 for fill solid 
        hV[j][0]->Draw("epx0same");				// x0 for no err. bar on X-axis

    }

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    wLatex(0.05,0.91,"(c)",1,12,0.08);              // figure's subcaption

    wLatex(0.65,0.83,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality
    wLatex(0.5,0.27,pttext,1,12,0.06);            // This is to plot per centrality

    for(int j = 0; j<3 ; j++) {
// ly 2014.05.06        wLatex(0.03, (hV[j][0]->GetFunction("fpol0")->GetParameter(0)-ymin)/(ymax-ymin)*0.8+0.12,graphtitle[j],graphcolor[j],12,0.06);
        double textx = 0.1;         // 2014.05.06 ly
        if(j==2) textx = 0.06;          // 2014.05.06 ly
        wLatex(textx, (hV[j][0]->GetFunction("fpol0")->GetParameter(0)-ymin)/(ymax-ymin)*0.8+0.12,graphtitle[j],graphcolor[j],12,0.06);            // ly 2014.05.06    fqwang's comments
    }
//    wLatex(0.03, (hV[2][0]->GetFunction("fpol0")->GetParameter(0)*scalev3-ymin)/(ymax-ymin),graphtitle[2],graphcolor[2],12,0.06);               // V3{2}
    

    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname2));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname2));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname2));
//        c->SaveAs(Form("figures/%s.gif",outname2));
    }


    // =========== draw v3


    TCanvas *c1 = new TCanvas("c1","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad1 = new TPad("pad1","pad1",PADXL,PADYL,PADXH,PADYH);

    ymin = 1.9;
    ymax = 5.5;
    TH2D *htmp3 = new TH2D("htmp3","",100,-1,1,100,ymin,ymax);    // x-range, y-range go here
    htmp3->GetXaxis()->SetTitle("#eta");
    
    FigFormat1D(c1,pad1,htmp3);

    c1->cd();
    pad1->Draw();
    pad1->cd();
    htmp3->Draw();

    for(int j = 3; j<3+2; j++) {          // v_{3} only
//        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),graphbandcolor[j],0);
        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),graphbandcolor[j],0,1,1,1,1001);			// fill style is hard to see, change to solid fill 1001
//        hV[j][0]->Draw("psame");
        hV[j][0]->Draw("epx0same");
    }

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    wLatex(0.05,0.91,"(d)",1,12,0.08);              // figure's subcaption

    wLatex(0.65,0.83,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality
    wLatex(0.5,0.27,pttext,1,12,0.06);            // This is to plot per centrality


    for(int j = 3; j<3+2 ; j++) {
// 2014.05.06 ly        wLatex(0.03, (hV[j][0]->GetFunction("fpol0")->GetParameter(0)-ymin)/(ymax-ymin)*0.8+0.12,graphtitle[j],graphcolor[j],12,0.06);
        wLatex(0.1, (hV[j][0]->GetFunction("fpol0")->GetParameter(0)-ymin)/(ymax-ymin)*0.8+0.12,graphtitle[j],graphcolor[j],12,0.06);            // 2014.05.06 ly fqwang's comment
    }
//    wLatex(0.03, (hV[2][0]->GetFunction("fpol0")->GetParameter(0)*scalev3-ymin)/(ymax-ymin),graphtitle[2],graphcolor[2],12,0.06);               // V3{2}
    

    if(savefig) {
        c1->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname3));
        c1->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname3));
        c1->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname3));
//        c1->SaveAs(Form("figures/%s.gif",outname3));
    }

    //=========================== get values for V{2}, V{4}, <v^2>, d sys. err. & stat. err.
    double Mean[NoGraph] = {0}, Err[NoGraph] = {0}, SlErr[NoGraph] = {0}, ShErr[NoGraph] = {0};        // Mean,            stat. err,              sys. err low,          sys. err. high
    char graphtextname[NoGraph][100] = {"V_{2}{2}","<v_{2}^{2}>","V_{2}{4}","V_{3}{2}","<v_{3}^{2}>","V_{3}{4}","D_{2}","D_{3}"};
    for(int i = 0; i<NoGraph ; i++) {
        Mean[i] = Average(gV[i][0]->GetN(), gV[i][0]->GetY());
        Err[i] = Average(gV[i][0]->GetN(), gV[i][0]->GetEY());
        SlErr[i] = Average(syserr[i]->GetN(), syserr[i]->GetEYlow());
        ShErr[i] = Average(syserr[i]->GetN(), syserr[i]->GetEYhigh());

//-----------------2013.11.25 ly \delta fit parameters err. treated as sys. err. as askded by GPC 2nd comment
//        std:cout<<std::setprecision(4)<<graphtextname[i]<<" = "<<Mean[i]<<" +- "<<Err[i]<<"("<<Err[i]/Mean[i]*100<<"% stat.) + "<<ShErr[i]<<"("<<ShErr[i]/Mean[i]*100<<"% sys.) - "<<fabs(SlErr[i])<<"("<<fabs(SlErr[i]/Mean[i]*100)<<"% sys.)"<<endl; 
	if(i<6) {
        	std:cout<<std::setprecision(4)<<graphtextname[i]<<" = "<<Mean[i]<<" +- "<<Err[i]<<"("<<Err[i]/Mean[i]*100<<"% stat.) + "<<ShErr[i]<<"("<<ShErr[i]/Mean[i]*100<<"% sys.) - "<<fabs(SlErr[i])<<"("<<fabs(SlErr[i]/Mean[i]*100)<<"% sys.)"<<endl; 
	}
	else {
        	std:cout<<std::setprecision(4)<<graphtextname[i]<<" = "<<Mean[i]<<" (fit function no stat. err.) + "<<ShErr[i]<<"("<<ShErr[i]/Mean[i]*100<<"% sys.) - "<<fabs(SlErr[i])<<"("<<fabs(SlErr[i]/Mean[i]*100)<<"% sys.)"<<endl; 
	}
//---------------------------------------------------------------------------------------------------
    }
    // for <v2> and <v3>
    std:cout<<std::setprecision(4)<<"sqrt{<v_{2}^{2}>}"<<" = "<<sqrt(Mean[1])<<" +- "<<Err[1]/(2*sqrt(Mean[1]))<<"("<<Err[1]/(2*Mean[1])*100<<"% stat.) + "<<ShErr[1]/(2*sqrt(Mean[1]))<<"("<<ShErr[1]/(2*Mean[1])*100<<"% sys.) - "<<fabs(SlErr[1]/(2*sqrt(Mean[1])))<<"("<<fabs(SlErr[1]/(2*Mean[1])*100)<<"% sys.)"<<endl; 
    std:cout<<std::setprecision(4)<<"sqrt{<v_{3}^{2}>}"<<" = "<<sqrt(Mean[4])<<" +- "<<Err[4]/(2*sqrt(Mean[4]))<<"("<<Err[4]/(2*Mean[4])*100<<"% stat.) + "<<ShErr[4]/(2*sqrt(Mean[4]))<<"("<<ShErr[4]/(2*Mean[4])*100<<"% sys.) - "<<fabs(SlErr[4]/(2*sqrt(Mean[4])))<<"("<<fabs(SlErr[4]/(2*Mean[4])*100)<<"% sys.)"<<endl; 



    // =================== get value and err. for fluct/flow = sigma_{2}/v_{2}
    double vv2[V2Nofile], V4[V4Nofile];
    double evv2[V2Nofile], eV4[V4Nofile];
    for(int i = 0; i<V2Nofile ; i++) { 
        vv2[i] = Average(gV[1][i]->GetN(), gV[1][i]->GetY());
        evv2[i] = Average(gV[1][i]->GetN(), gV[1][i]->GetEY());
    }
    for(int i = 0; i<V4Nofile ; i++) { 
        V4[i] = Average(gV[2][i]->GetN(), gV[2][i]->GetY());
        eV4[i] = Average(gV[2][i]->GetN(), gV[2][i]->GetEY());
    }

    // sys. err. 
    // get mean
    double sigma2v2[V2Nofile+V4Nofile];         // array size is larger than real used.
    for(int i = 0; i<TMath::Min(V2Nofile,V4Nofile); i++) {
        sigma2v2[i] = ratio4sigma2v2(vv2[i],V4[i]+4./3.*k);   // k is fit slope for \sigma'{4}, ek is the error for \sigma'{4} slope. <\sigma'> = 4/3 k. \sigma' is treated as a same number in sys. err. It doesn't contribute to sys. err.
    }
    for(int i = TMath::Min(V2Nofile,V4Nofile); i<TMath::Max(V2Nofile,V4Nofile); i++) {
        sigma2v2[i] = ratio4sigma2v2(vv2[i],V4[0]+4./3.*k);             // rest of V2 root file is different fit function with the same cuts set for the first file as V4
    }
    // get err from diff. fits or cuts.
    double syserr4sigma2v2[V2Nofile+V4Nofile];         // array size is larger than real used.
    for(int i = 1; i<TMath::Min(V2Nofile,V4Nofile); i++) {
        syserr4sigma2v2[i] = sigma2v2[i]-sigma2v2[0];				// diff. cut
    }
//    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)] = TMath::MaxElement(TMath::Max(V2Nofile,V4Nofile)-TMath::Min(V2Nofile,V4Nofile), sigma2v2+TMath::Min(V2Nofile,V4Nofile)) - sigma2v2[0];      // high err from diff. fit functions
//    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)+1] = TMath::MinElement(TMath::Max(V2Nofile,V4Nofile)-TMath::Min(V2Nofile,V4Nofile), sigma2v2+TMath::Min(V2Nofile,V4Nofile)) - sigma2v2[0];      // low err from diff. fit functions
//-------------------------------------------------2013.11.25 ly-----------------------------------------------------------------------------------
// \sigma' fitting error needs to be added in the sys. err. it was treated as stat. err. 2013. 11. 25 ly
//    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)] = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,ShErr[6],0);            // high diff fit func. errs are from \delta and propogate it to the ratio
//    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)+1] = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,fabs(SlErr[6]),0);         // low diff fit func err from \delta prop. to the ratio
    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)] = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,ShErr[6],4./3.*ek);            // high sys. err.: diff fit func. errs and fitting parameters errors (already sqrt sum in syserr[6], which is <d_2>) are from \delta and propogate it to the ratio, \sigma' fitting error propogate to ratio
    syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)+1] = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,fabs(SlErr[6]),4./3.*ek);         // low sys. err.: diff fit func. errs and fitting parameters errors (already sqrt sum in syserr[6], which is <d_2>) are from \delta and propogate it to the ratio, \sigma' fitting error propogate to ratio
    // calculate sys. err for low and high
//-------------------------------------------------end-------------------------------------------------------------------------------------
    double hsyserrRatio = 0, lsyserrRatio = 0;
    double tmpsyserrRatio = 0;
    for(int i = 0; i<TMath::Min(V2Nofile,V4Nofile); i++) {
        tmpsyserrRatio+=pow(syserr4sigma2v2[i],2);				// Diff cuts. are diff. source err.
    }
    hsyserrRatio = sqrt( tmpsyserrRatio + pow(syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)],2) );
    lsyserrRatio = sqrt( tmpsyserrRatio + pow(syserr4sigma2v2[TMath::Min(V2Nofile,V4Nofile)+1],2) );

//    hsyserrRatio = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,ShErr[1],ShErr[2]);     // sys. err. propogate
//    lsyserrRatio = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,SlErr[1],SlErr[2]);
//    cout<<"(sigma/<v>) = "<<sigma2v2[0]*100<<"% + "<<(TMath::MaxElement(TMath::Max(V2Nofile,V4Nofile),sigma2v2)-sigma2v2[0])*100<<"%(sys.) - "<<(sigma2v2[0] - TMath::MinElement(TMath::Max(V2Nofile,V4Nofile),sigma2v2))*100<<"%(sys.)";



//=========== stat. err. =====================
//-------------------------------------------------2013.11.25 ly-----------------------------------------------------------------------------------
// \sigma' fitting error needs to be added in the sys. err. it was treated as stat. err. 2013. 11. 25 ly
//    double hstaterrRatio = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,evv2[0],sqrt(pow(eV4[0],2)+pow(4./3.*(ek),2)));            // k is \sigma'{4}, ek is the error for \sigma'{4}. we are calculating ( <v^2> - (V{4} + \sigma') ) / ( <v^2> + (V{4} + \sigma') )
////    double lstaterrRatio = StatErrRatio4sigma2v2(vv2[0],V4[0],evv2[0],eV4[0]);            // k is \sigma'{4}, ek is the error for \sigma'{4}
//    double lstaterrRatio = hstaterrRatio;
    double hstaterrRatio = StatErrRatio4sigma2v2(vv2[0],V4[0]+4./3.*k,evv2[0],fabs(eV4[0]));            // k is \sigma'{4}, ek is the error for \sigma'{4}, which is treated as sys.err. as asked by GPC 2nd comments 2013.11.25 ly. we are calculating ( <v^2> - (V{4} + \sigma') ) / ( <v^2> + (V{4} + \sigma') )
    double lstaterrRatio = hstaterrRatio;
//-------------------------------------------------end-------------------------------------------------------------------------------------
//    cout<<" + "<<hstaterrRatio*100<<"%(stat.)"<<" - "<<lstaterrRatio*100<<"%(stat.)"<<endl;
    cout<<"sigma' = "<<k*4./3<<" +- "<<ek*4./3<<"(stat.)"<<endl;
    cout<<"(sigma/<v>) = "<<sigma2v2[0]*100<<"% +- "<<hstaterrRatio*100<<"%(stat.)"<<" + "<<hsyserrRatio*100<<"%(sys.) - "<<lsyserrRatio*100<<"%(sys.)"<<endl;
    cout<<"sqrt((V2-V4)/(V2+V4)) = sqrt(("<<Mean[0]<<"-"<<Mean[2]<<")/("<<Mean[0]<<"+"<<Mean[2]<<")) = "<<ratio4sigma2v2(Mean[0],Mean[2])*100<<"% +- "<<StatErrRatio4sigma2v2(Mean[0],Mean[2],Err[0],Err[2])*100<<"%(stat.)"<<endl;
    cout<<"sqrt((vv2-V4)/(vv+V4)) = sqrt(("<<vv2[0]<<"-"<<Mean[2]<<")/("<<vv2[0]<<"+"<<Mean[2]<<")) = "<<ratio4sigma2v2(vv2[0],Mean[2])*100<<"% +- "<<StatErrRatio4sigma2v2(vv2[0],Mean[2],evv2[0],Err[2])*  100<<"%(stat.)"<<endl;


//    for( int i = 0; i<V2Nofile ; i++)   {
//        fV2[i][0]->Close();
//        fV2[i][1]->Close();
//    }
//    for( int i = 0; i<V4Nofile ; i++)   {
//        fV4[i]->Close();
//    }
}



//=======================================================================================
void plotdelta(char centtag[]) {                   // <delta> vs eta-gap

    // read the data
//    const int Nofile = 5;
//    char filename[Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530","fit0cent5eta20dca3_3dataset.130530","fit8cent5eta20dca3_3dataset.130530","fit3cent5eta20dca3_3dataset.130530","fit2cent5eta20dca3_3dataset.130530"};         
//    char filename[Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530","fit1cent5eta20dca2_3dataset.130506","fit1cent5eta20dca2hit15_3dataset.130531","fit1cent5eta20dca3Vz25_3dataset.130530"};
    const int Nofile =6;                                                       // v3 only
    char filename[Nofile][1024] = {"fit1cent%seta20dca3_3dataset.130530","fit0cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530","fit2cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};     // v3 only fit err.
//    const int Nofile = 6;
//    char filename[Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530",                  "fit1cent5eta20dca2_3dataset.130506","fit1cent5eta20dca3Vz25_3dataset.130530",                  "fit1cent5eta20dca3hit15_3dataset.130328","fit0cent5eta20dca3_3dataset.130530",                 "fit3cent5eta20dca3_3dataset.130530"};                       // v3 only  fit and cuts error
//    char filename[Nofile][1024] = {"fit1cent5eta20dca3_3dataset.130530",                  "fit1","fit1cent5eta20dca3Vz25_3dataset.130530",                  "fit1cent5eta20dca3hit15_3dataset.130328","fit0cent5eta20dca3_3dataset.130530",                 "fit3cent5eta20dca3_3dataset.130530"};                       // v3 only  fit and cuts error
    const int Nograph = 1;
    int NofileUsed[Nograph] = {Nofile};
    int vn[Nograph] = {3};                      // v3
//    const int Nograph = 2;
//    int vn[Nograph] = {2, 3};                 // v2
//    int NofileUsed[Nograph] = {3,5};           // the fit2 is gaus only, bad Chi for v2 
                                               // the fit3 is exp only, also not good Chi for v2
                                               // v3 diff. func. sys. err. will not be ploted.
    char graphtag[100] = "gd%d2";             // vs eta-gap
//    char graphtag[100] = "gdVsDeltaeta%d2";     // vs \Delta\eta
    char outname[500] = "";
    sprintf(outname,"deltaetagap.%s",Form(filename[0],centtag));
    TGraphAsymmErrors *gr[Nograph][Nofile];
    TFile *f[Nofile];
    for(int i = 0; i<Nofile; i++) {
        f[i] = new TFile(Form("%s.root",Form(filename[i],centtag)));
        for(int j = 0; j<Nograph ; j++) {
            gr[j][i] = (TGraphAsymmErrors*)f[i]->Get(Form(graphtag,vn[j]));
        }
    }
    char gbandtag[100] = "gerrband_cumu2_d%d2";             // stat. err. from fit par.
//    char gbandtag[100] = "gerrband_cumu2VsDeltaeta_d%d2";             // stat. err. from fit par.
    TGraphAsymmErrors *gband[Nograph];
    for(int j = 0; j<Nograph ; j++) {
        gband[j] = (TGraphAsymmErrors*)f[0]->Get(Form(gbandtag,vn[j]));
    }

//--------------------------------------2013.11.25 ly-------------------------------------------------------------------------
// fitting parameters err. is treated as sys. err. as asked by GPC 2nd comments 2013.11.25 ly    
//    ClassSysErr *syserr[Nograph];
//    for(int j = 0; j<Nograph ; j++) {
//        syserr[j] = new ClassSysErr(Nofile-1, gr[j][0],gr[j]+1,"s");          // sys. err. diff fit func.
////        syserr[j] = (ClassSysErr*) SysErrMultSingleS(gr[j][0],3,gr[j]+1,2,gr[j]+4);       // fits and cuts sys. err. sum
//    }
    ClassSysErr *syserr[Nograph];
    ClassSysErr *syserrtmp1[Nograph], *syserrtmp2[Nograph];
    const int MAXPOINT = 1024;
    double tmpgey[MAXPOINT];
    for(int j = 0; j<Nograph ; j++) {
	for(int ip = 0; ip<gr[j][0]->GetN(); ip++) {
		tmpgey[ip] = - gr[j][0]->GetErrorYlow(ip);					// ClassSysErr: low err is defined as negative
	}
	syserrtmp1[j] = new ClassSysErr(gr[j][0]->GetN(),gr[j][0]->GetX(),gr[j][0]->GetY(),tmpgey,gr[j][0]->GetEYhigh());			// sys. err. fit par err.
    syserrtmp2[j] = new ClassSysErr(Nofile-1, gr[j][0],gr[j]+1,"s");          // sys. err. diff fit func.
	syserr[j] = (ClassSysErr*) SumTwoSSysErr(gr[j][0],syserrtmp1[j]->GetEYlow(),syserrtmp1[j]->GetEYhigh(), syserrtmp2[j]->GetEYlow(), syserrtmp2[j]->GetEYhigh());
    }
//------------------------------------------end-------------------------------------------------------------------------



    // setup enviroment for drawing
    gStyle->SetOptStat(0);
    
    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

//    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-1,3.5);    // x-range, y-range go here       this is range of v2
// 2014.05.02 ly    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-0.3,2.8);    // x-range, y-range go here       this is range of v2
    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-0.3,2.3);    // x-range, y-range go here       this is range of v2                 // 2014.05.02 ly upon NPI comments
//    htmp->GetXaxis()->SetTitle("#Delta#eta-gap");
    htmp->GetXaxis()->SetTitle("|#Delta#eta|>x");

// as asked by 2nd GPC comments, show V{m} instead of Vn{m}		2013.11.25 ly    
//    htmp->GetYaxis()->SetTitle("#bar{D}_{3} = #sigma_{3}' + #delta_{3}");            // only v3
//    it turns out to be a misunderstanding. It has been asked to add those vn back in title. 2013.12.16 ly
//    htmp->GetYaxis()->SetTitle("#bar{D} = #sigma' + #delta");            // only v3   // 2013.12.16 ly
//    htmp->GetYaxis()->SetTitle("#bar{D}_{3} = #sigma_{3}' + #delta_{3}");            // only v3
    htmp->GetYaxis()->SetTitle("#bar{D}_{3}(|#Delta#eta|)");            // only v3
    
    FigFormat1D(c,pad,htmp);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    wLatex(0.05,0.91,"(b)",1,12,ZTITLESIZE);

    wLatex(0.64,0.83,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality
    wLatex(0.47,0.735,pttext,1,12,0.07);            // This is to plot per centrality

    wLatex(0.75,0.67,"|#eta|<1",1,12,0.07);            // This is to plot per centrality

    TLine *line = new TLine(0,0,2,0);
// 2014.05.02 ly    line->Draw("same");
//    for(int j = 0; j<1 ; j++) {             // only for v2
//2014.05.06 ly    int bandcolor[Nograph] = {kOrange+1};           // kYellow};        // 2014.05.02 ly
    int bandcolor[Nograph] = {kGray};                               // 2014.05.06 ly
    for(int j = 0; j<Nograph; j++) {          
        graphSystBand(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),bandcolor[j],20,0,1,1,1001);
//        graphSystLine(syserr[j]->GetN(),syserr[j]->GetX(),syserr[j]->GetY(),0,0,syserr[j]->GetEYhigh(),syserr[j]->GetEYlow(),kBlue,0,1,1,1,2,2);
    }
    line->Draw("same");             // 2014.05.02 ly

//    int bandcolor[Nograph] = {kGray,kYellow};
    for(int j = 0; j<Nograph ; j++) {
        int ic = bandcolor[j];
        gband[j]->SetFillColor(ic);
//        gband[j]->SetFillStyle(3001);
//        gband[j]->Draw("fsame");                  // fitting err. -> stat. err. ploted as error bar on gr[][]. drawn as "pesame" below.
        gr[j][0]->SetMarkerStyle(22);                   // keep consistent with plotcent()
        gr[j][0]->SetMarkerSize(2);
//        gr[j][0]->RemovePoint(gr[j][0]->GetN()-1);			// 2013.11.25 ly no idea why not draw the last point, so remove this sentence

//        gr[j][0]->Draw("pesame");					// 2013.11.25 ly. fit par. err. treated as sys. err as 2nd GPC comments 
        gr[j][0]->SetLineWidth(4);
//        gr[j][0]->Draw("cXsame");					// 2013.11.25 ly no stat. err. for fit function, see above                     // 2013.12.16 ly 3rd GPC comment draw as point not line
        gr[j][0]->SetMarkerColor(1);                // 2014.05.06 ly fqwang's comment
        gr[j][0]->Draw("pXsame");					// 2013.11.25 ly no stat. err. for fit function, see above             // 2013.12.16 ly 3rd GPC comment draw as point not line

//        for(int i = 1;i<Nofile; i++) {
//            gr[j][i]->SetLineColor(i);
//            gr[j][i]->SetMarkerColor(i);
//            gr[j][i]->Draw("pXsame");
//        }
    }


    double ytxt = 0.8;
//    int markerstyle[Nograph] = {8,21};          // solid dots, solid square
    int markerstyle[Nograph] = {21};          // solid dots, solid square
    double markersize = 1.5;
    double tsize = 0.05;
//    for(int j = 0; j<Nograph ; j++) {
//        keySymbol(0.65,ytxt,Form("#LT#delta_{%d}#GT+#LT#sigma_{%d}'#GT",vn[j],vn[j]),vn[j]-1,markerstyle[j],tsize, markersize);
//        ytxt-=0.07;
//    }


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
//        c->SaveAs(Form("figures/%s.gif",outname));
    }
}

//=======================================================================================
TGraphAsymmErrors *TF2TGraphAsymmErrors(TF1 *f1, double xmin, double xmax,int Npar, double *par, double *err) {      // create a TGraphErrors using TF1 *f1 with parameter par, error err. The value is from par. One sigma variance for errors. Graph range from xmin to xmax
    const int gNpoint = 100;
    double x[gNpoint], y[gNpoint], eyl[gNpoint], eyh[gNpoint];

    for(int j = 0; j<Npar; j++) {
        f1->SetParameter(j,par[j]);
    }
    for(int i = 0; i<gNpoint; i++) {
        x[i] = (xmax-xmin)/gNpoint*(i+0.5) + xmin;
        y[i] = f1->Eval(x[i]);
    }
    for(int j = 0; j<Npar; j++) {
        f1->SetParameter(j,par[j]+err[j]);
    }
    for(int i = 0; i<gNpoint; i++) {
        x[i] = (xmax-xmin)/gNpoint*(i+0.5) + xmin;
        eyh[i] = f1->Eval(x[i])-y[i];
    }
    for(int j = 0; j<Npar; j++) {
        f1->SetParameter(j,par[j]-err[j]);
    }
    for(int i = 0; i<gNpoint; i++) {
        x[i] = (xmax-xmin)/gNpoint*(i+0.5) + xmin;
        eyl[i] = y[i]-f1->Eval(x[i]);
    }
    TGraphAsymmErrors *g = new TGraphAsymmErrors(gNpoint,x,y,0,0,eyl,eyh);

    return g;
}

void plotdeltasigma(char centtag[]) {              // <delta>, <sigma^2> vs eta-gap

    // read the data
//    const int dNofile = 3;
//    char dfilename[dNofile][100] = {"fit1cent5eta20dca3_3dataset.130530","fit0cent5eta20dca3_3dataset.130530","fit8cent5eta20dca3_3dataset.130530"};         // read delta{2}
//    TFile *df[dNofile];
    const int dNofile = 7;
    char dfilename[dNofile][100] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530", "fit4cent%seta20dca3_3dataset.130530"};         // read delta{2}
    TFile *df[dNofile];
    TGraphAsymmErrors *gd22[dNofile];
    for(int i = 0 ;i<dNofile; i++) {
        df[i] = new TFile(Form("%s.root",Form(dfilename[i],centtag)));
        gd22[i] = (TGraphAsymmErrors*)df[i]->Get("gd22");                 // vs eta-gap
//        gd22[i] = (TGraphAsymmErrors*)df[i]->Get("gdVsDeltaeta22");         // vs \Delta\eta
    }
    df[0]->cd();    
    TGraph *gd22errband = (TGraph*)df[0]->Get("gerrband_cumu2_d22");      // vs eta-gap
//    TGraph *gd22errband = (TGraph*)df[0]->Get("gerrband_cumu2VsDeltaeta_d22");  // vs \Delta\eta
    gd22[0]->SetMarkerSize(2);
    gd22[0]->SetMarkerColor(1);
    gd22[0]->SetMarkerStyle(8);
    for(int i = 1; i<dNofile; i++) {
        gd22[i]->SetLineColor(kBlue);
        gd22[i]->SetLineStyle(2);
    }
    gd22errband->SetFillColor(kGray);
    gd22errband->SetFillStyle(3001);
    


//--------------------------------------2013.11.25 ly-------------------------------------------------------------------------
// fitting parameters err. is treated as sys. err. as asked by GPC 2nd comments 2013.11.25 ly    
    ClassSysErr *syserrtmp1,*syserrtmp2;
    const int MAXPOINT = 1024;
    double tmpgey[MAXPOINT];
    for(int ip = 0; ip<gd22[0]->GetN(); ip++) {
	    tmpgey[ip] = - gd22[0]->GetErrorYlow(ip);						// ClassSysErr: low err is defined as negative
    }
    syserrtmp1 = new ClassSysErr(gd22[0]->GetN(), gd22[0]->GetX(), gd22[0]->GetY(), tmpgey, gd22[0]->GetEYhigh());			// sys. err. fit par err.
    syserrtmp2 = new ClassSysErr(3,gd22[0],gd22+4,"s");					// sys. err from diff. fit func.
//    syserrtmp2 = (ClassSysErr*) SysErrMultSingleS(gd22[0],3,gd22+1,3,gd22+4);      //  sys. err. from diff. fit func. and diff. cuts		// 2013.11.25 ly D function only have fit func. err. . Diff cuts err. should not be taken into D err.
//    ClassSysErr *syserr = (ClassSysErr*) SysErrMultSingleS(gd22[0],3,gd22+1,3,gd22+4);        // 2013.11.25 ly fitting error treated as sys. err. as advised by GPC
////    ClassSysErr *syserr = (ClassSysErr*) SysErrMultSingleS(gd22[0],2,gd22+1,2,gd22+3);
    ClassSysErr *syserr = (ClassSysErr*) SumTwoSSysErr(gd22[0],syserrtmp1->GetEYlow(),syserrtmp1->GetEYhigh(),syserrtmp2->GetEYlow(),syserrtmp2->GetEYhigh());					// sqrt( diff fit func.^2 + fit par. err.^2 )
//--------------------------------------end-------------------------------------------------------------------------


    char sfilename[100] = "fit1cent%seta10dca3_3dataset.130603";         // read sigma'{4}
    TFile *sf = new TFile(Form("%s.root",Form(sfilename,centtag)));
    TGraphErrors *gave4 = (TGraphErrors*)sf->Get("gavedd24");
    TF1 *fave4 = (TF1*)gave4->GetFunction("f0D");                       // obtain fit func.
    double par0 = fave4->GetParameter(0);
    double err0 = fave4->GetParError(0);
    TF1 *fsigma = new TF1("fsigma","[0]*(2-x)/2",0,2);                // vs eta-gap
//    TF1 *fsigma = new TF1("fsigma","[0]*(2-x)",0,2);                  // vs \Delta\eta
    double xmin, xmax;
    fsigma->GetRange(xmin,xmax);
    TGraphAsymmErrors *gs = (TGraphAsymmErrors*)TF2TGraphAsymmErrors(fsigma,xmin,xmax,1,&par0,&err0);      // sigma' from Func. with one sigma variance
// 2014.05.06 ly    Int_t sigmacolor = (Int_t)(kCyan+2);
    Int_t sigmacolor = (Int_t)(kGray+1);            // 2014.05.02 ly
    gs->SetLineColor(sigmacolor);
    gs->SetFillColor(sigmacolor);
//    gs->SetFillStyle(3545);
// 2014.05.02 ly    gs->SetFillStyle(3005);
    gs->SetFillStyle(3356);         // 2014.05.06 ly
    gStyle->SetHatchesLineWidth(4);         // 2014.05.06 ly
    gStyle->SetHatchesSpacing(3);           // 2014.05.06 ly

    char outname[500] = "";
    sprintf(outname,"delta_sigma.%s.%s",Form(dfilename[0],centtag),Form(sfilename,centtag));

    // setup enviroment for drawing
    
    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

// 2014.05.02 ly    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-0.7,7);    // x-range, y-range go here
    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-0.6,5.5);    // x-range, y-range go here                // 2014.05.02 ly upon NPI comments
//    TH2D *htmp = new TH2D("htmp","",100,0,2,100,-2,8);    // x-range, y-range go here
//    htmp->GetXaxis()->SetTitle("#Delta#eta-gap");
    htmp->GetXaxis()->SetTitle("|#Delta#eta|>x");
// as 2nd GPC comments asked not to draw V_n{m}, instead only draw V{m}  2013.11.25 ly
//    htmp->GetYaxis()->SetTitle("#bar{D}_{2} = #sigma_{2}' + #delta_{2}");
//    it turns out to be a misunderstanding. It has been asked to add those vn back in title. 2013.12.16 ly
//    htmp->GetYaxis()->SetTitle("#bar{D} = #sigma' + #delta");     // 2013.12.16 ly
//    htmp->GetYaxis()->SetTitle("#bar{D}_{2} = #sigma_{2}' + #delta_{2}");
    htmp->GetYaxis()->SetTitle("#bar{D}_{2}(|#Delta#eta|)");
    
    FigFormat1D(c,pad,htmp);

    gStyle->SetOptStat(0);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "#times10^{-4}",1,12,0.06);

    TLine *line = new TLine(0,0,2,0);
    line->Draw("same");
//    gd22errband->Draw("f");               // fitting err. is stat. err. drawn as bar on gd22[0] as below
// 2014.05.02 ly    graphSystBand(syserr->GetN(),syserr->GetX(),syserr->GetY(),0,0,syserr->GetEYhigh(),syserr->GetEYlow(),kGray,0);                   // keep the same with plotcent()
    graphSystBand(syserr->GetN(),syserr->GetX(),syserr->GetY(),0,0,syserr->GetEYhigh(),syserr->GetEYlow(),kGray,0,1,1,1,1001);                   // keep the same with plotcent()                // 2014.05.02 ly
    gs->Draw("e3");
//    gd22[0]->RemovePoint(gd22[0]->GetN()-1);				// 2013.11.25 ly no idea why remove last point, so delete this sentence
//    gd22[0]->Draw("pesame");
    gd22[0]->SetLineWidth(4);
//    gd22[0]->Draw("cXsame");              // 2013.12.16 ly
    gd22[0]->Draw("pXsame");                // 3rd GPC comment, draw as point not line. 2013.12.16 ly
//    gd22[1]->Draw("lX");
//    gd22[2]->Draw("lX");
//    for(int i = 1; i<dNofile; i++) {
//        gd22[i]->SetMarkerColor(i);
//        gd22[i]->Draw("pXsame");
//    }

// 2014.05.02 ly    wLatex(0.3,0.35,"#sigma'",sigmacolor,12,0.08);         // vs eta-gap
    wLatex(0.33,0.37,"#sigma_{2}'",1,12,0.08);         // vs eta-gap            // 2014.05.02 ly
//    wLatex(0.3,0.48,"#sigma'",sigmacolor,12,0.08);           // vs \Delta\eta

    wLatex(0.05,0.91,"(a)",1,12,ZTITLESIZE);

    wLatex(0.64,0.83,centtag2text(centtag),1,12,ZTITLESIZE);            // This is to plot per centrality
    wLatex(0.47,0.735,pttext,1,12,0.07);            // This is to plot per centrality

    wLatex(0.75,0.67,"|#eta|<1",1,12,0.07);            // This is to plot per centrality


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
//        c->SaveAs(Form("figures/%s.gif",outname));
    }
}


//=============================================================================================
void readYeY4Graph(TGraphErrors *gr,double x, double &y, double &ey) {
    y = gr->Eval(x);
    double *gx = gr->GetX();
    double *gy = gr->GetY();
    double *gey = gr->GetEY();
    double tmpx = 0, tmpy = 0, tmpey = 0;
    double tmplx = 0, tmply = 0, tmpley = 0;
    for(int i = 0; i<gr->GetN(); i++) {
        tmplx = tmpx;
        tmply = tmpy;
        tmpley = tmpey;
        tmpx = gx[i];
        tmpy = gy[i];
        tmpey = gey[i];
        if(tmpx>x && i>0) {
            tmpy = (tmpy-tmply)/(tmpx-tmplx)*(x-tmplx)+tmply;
            tmpey = tmpey*(x-tmplx)/(tmpx-tmplx) + tmpley*(tmpx-tmplx)/(tmpx-tmplx);
//            tmpy = 0.5*(tmpy+tmply);
//            tmpey = 0.5*(tmpey+tmpley);
            break;
        }
    }
    if(fabs(y-tmpy)<tmpey) {ey = tmpey;}
    else { cout<<"Please CHECK function readYeY4Graph(): |y-tmpy| = "<<fabs(y-tmpy)<<" > tmey = "<< tmpey <<endl; ey = 0;}

}


void readYeY4Graph(TGraphAsymmErrors *gr,double x, double &y, double &eyl, double &eyh) {
    y = gr->Eval(x);
    double *gx = gr->GetX();
    double *gy = gr->GetY();
    double *geyl = gr->GetEYlow();
    double *geyh = gr->GetEYhigh();
    double tmpx = 0, tmpy = 0, tmpeyl = 0, tmpeyh = 0;
    double tmplx = 0, tmply = 0, tmpleyl = 0, tmpleyh = 0;
    for(int i = 0; i<gr->GetN(); i++) {
        tmplx = tmpx;
        tmply = tmpy;
        tmpleyl = tmpeyl;
        tmpleyh = tmpeyh;
        tmpx = gx[i];
        tmpy = gy[i];
        tmpeyl = geyl[i];
        tmpeyh = geyh[i];
        if(tmpx>x && i>0) {
            tmpy = (tmpy-tmply)/(tmpx-tmplx)*(x-tmplx)+tmply;
            tmpeyl = tmpeyl*(x-tmplx)/(tmpx-tmplx) + tmpleyl*(tmpx-tmplx)/(tmpx-tmplx);
            tmpeyh = tmpeyh*(x-tmplx)/(tmpx-tmplx) + tmpleyh*(tmpx-tmplx)/(tmpx-tmplx);
//            tmpy = 0.5*(tmpy+tmply);
//            tmpeyl = 0.5*(tmpeyl+tmpleyl);
//            tmpeyh = 0.5*(tmpeyh+tmpleyh);
            break;
        }
    }
    if(fabs(y-tmpy)<0.5*(tmpeyl+tmpeyh)) {eyl = tmpeyl;  eyh = tmpeyh;}
    else { cout<<"Please CHECK function readYeY4Graph(): |y-tmpy| = "<<fabs(y-tmpy)<<" > 0.5(tmpeyl+tmpeyh) = "<< 0.5*(tmpeyl+tmpeyh) <<endl; eyl = 0; eyh = 0;}


}

TGraphErrors* EtaGapDepfrom2Dhist(TH2D *hV) {
    const int ncut = 20;                      // 20 points eta-gap
    double etamax = 1;                  
    double deta_cut[ncut] = {0}, v_cut[ncut] = {0}, ev_cut[ncut] = {0};
    for(int icut = 0; icut<ncut; icut++) {
        deta_cut[icut] = 2*etamax/ncut*icut;
        int ncnt = 0;
        double sumv = 0, esumv = 0;
        for(int i = 1; i<=hV->GetNbinsX(); i++) {
            float eta_1 = -etamax+(i-0.5)*2*etamax/hV->GetNbinsX();
            for(int j = 1; j<=hV->GetNbinsY() ; j++) {
                float eta_2 = -etamax+(j-0.5)*2*etamax/hV->GetNbinsY();
                if(fabs(eta_1-eta_2)<deta_cut[icut]) continue;
                sumv+=hV->GetBinContent(i,j);
                esumv+=pow(hV->GetBinError(i,j),2);
                ncnt++;
            }
        }
        if(!ncnt) continue;
        v_cut[icut] = sumv/ncnt;
        ev_cut[icut] = sqrt(esumv)/ncnt;
    }

    TGraphErrors* gvv = new TGraphErrors(ncut,deta_cut,v_cut,0,ev_cut);

    return gvv;
}

void plotcent(double etagap) {                    // centrality dependence

    // read the data
    const int Nocent = 5;
    char centtag[Nocent][100] = {"678","5","4","3","012"};
    int centloc[Nocent] = {10,25,35,45,65};         // %
    double centlocv2[Nocent] = {0};         // %               v2
    double centlocv3[Nocent] = {0};         // %               v3
    double shift = 0.5;
    for(int i = 0 ;i <Nocent  ; i++ ) {
        centlocv2[i] = centloc[i]+shift;
        centlocv3[i] = centloc[i]-shift;
    }
    // read v2
//    const int v2Nofile = 3;
//    char v2filetag[v2Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530",                     "fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530"};
//    const int v2Nofile = 6;
//    char v2filetag[v2Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530"};
    const int v2Nofile = 7;
    char v2filetag[v2Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};
    TFile *v2f[Nocent][v2Nofile];
    TH2D *hV22[Nocent][v2Nofile];			// V_{2}{2} vs eta_a, eta_b
    TGraphErrors *igV_22[Nocent][v2Nofile];			// it is eta-gap dep. from V_{2}{2} 2D hist.
    TGraphAsymmErrors *igd22[Nocent][v2Nofile];         // it is eta-gap dep. read in per cent per file
    TGraphErrors *igv22sq[Nocent][v2Nofile];              // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v22[Nocent][v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v22sq[Nocent][v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphErrors *gV_22[v2Nofile];			// vs cent
    TGraphAsymmErrors *gd22[v2Nofile] ;                  // vs cent
    TGraphAsymmErrors *gd22sq[v2Nofile] ;                  // vs cent
    TGraphErrors *gv22sq[v2Nofile] ;                       // vs cent
    TGraphErrors *gv22[v2Nofile] ;                       // vs cent
    TGraphAsymmErrors *gd2v22[v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *gd2v22sq[v2Nofile];        // it is eta-gap dep. read in per cent per file

    double v2scale = 0.5;                               // draw v2 as v2/2

    for(int i = 0; i<v2Nofile ; i++) {
        gd22[i] = new TGraphAsymmErrors();
        gd22sq[i] = new TGraphAsymmErrors();
        gd2v22[i] = new TGraphAsymmErrors();
        gd2v22sq[i] = new TGraphAsymmErrors();
        gv22[i] = new TGraphErrors();
        gv22sq[i] = new TGraphErrors();
        gV_22[i] = new TGraphErrors();
        for(int j = 0; j<Nocent; j++) {
            v2f[j][i] = new TFile(Form("%s.root",Form(v2filetag[i],centtag[j])));
            double y = 0, ey = 0, eyl = 0, eyh = 0;
            igd22[j][i] = (TGraphAsymmErrors*)v2f[j][i]->Get("gd22");
            readYeY4Graph(igd22[j][i],etagap,y,eyl,eyh);
            gd22sq[i]->SetPoint(j,centlocv2[j],y);
            gd22sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gd22[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
            igd2v22[j][i] = (TGraphAsymmErrors*)v2f[j][i]->Get("gd2v22");
            readYeY4Graph(igd2v22[j][i],etagap,y,eyl,eyh);
            gd2v22sq[i]->SetPoint(j,centlocv2[j],y);
            gd2v22sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd2v22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gd2v22[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
            igv22sq[j][i] = (TGraphErrors*)v2f[j][i]->Get("gv22");
            readYeY4Graph(igv22sq[j][i],etagap,y,ey);
            gv22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gv22[i]->SetPointError(j,0,ey/(2*sqrt(y)));
            gv22sq[i]->SetPoint(j,centlocv2[j],y);
            gv22sq[i]->SetPointError(j,0,ey);
            hV22[j][i] = (TH2D*)v2f[j][i]->Get("hV22");
            igV_22[j][i] = (TGraphErrors*) EtaGapDepfrom2Dhist(hV22[j][i]);
            readYeY4Graph(igV_22[j][i],etagap,y,ey);
            gV_22[i]->SetPoint(j,centlocv2[j],y);
            gV_22[i]->SetPointError(j,0,ey);
        }
    }

//    ClassSysErr *v2syserr = new ClassSysErr(3,gv22[0],gv22+1,"m");
//    ClassSysErr *d2syserr = new ClassSysErr(2,gd22[0],gd22+4,"s");
//    ClassSysErr *d2v2syserr = new ClassSysErr(v2Nofile-1,gd2v22[0],gd2v22+1);          // sigma/v
//    ClassSysErr *d2v2syserrsq = new ClassSysErr(v2Nofile-1,gd2v22sq[0],gd2v22sq+1);          // (sigma/v)^2
///    ClassSysErr *v2syserr = (ClassSysErr*)SysErrMultSingleS(gv22[0],3,gv22+1,2,gv22+4);
///    ClassSysErr *d2syserr = (ClassSysErr*)SysErrMultSingleS(gd22[0],3,gd22+1,2,gd22+4);
    ClassSysErr *d2syserr = new ClassSysErr(v2Nofile-1-3,gd22[0],gd22+4,"s");              // delta only fit function sys. err.
    ClassSysErr *V_22syserr = new ClassSysErr(3,gV_22[0],gV_22+1,"m");		// V_{2}{2} only cut sys. err.
    // <v^2> has two source of sys. err. 
    double V_22syserrl[Nocent] = {0},v2cutsyserrl[Nocent] = {0}, *tmpV22sel;
    double V_22syserrh[Nocent] = {0},v2cutsyserrh[Nocent] = {0}, *tmpV22seh;
    double V22[Nocent] = {0}, *tmpV22;
    double *vv2;
    vv2 = gv22sq[0]->GetY();
    // v2 cuts sys. err
    tmpV22 = gv22sq[0]->GetY();
    tmpV22sel = V_22syserr->GetEYlow();
    tmpV22seh = V_22syserr->GetEYhigh();
    for(int i = 0; i<V_22syserr->GetN();i++) {
	    V22[i] = tmpV22[i];							// v22
	    V_22syserrl[i] = tmpV22sel[i];					// error for V22
	    V_22syserrh[i] = tmpV22seh[i];
	    if(vv2[i]>=0) {
	    	v2cutsyserrl[i] = tmpV22sel[i]/(2*sqrt(vv2[i]));			// error for sqrt(V22)
	    	v2cutsyserrh[i] = tmpV22seh[i]/(2*sqrt(vv2[i]));
	    }
	    else {
		cout<<"vv2["<<i<<"] < 0!!"<<endl;
	    	v2cutsyserrl[i] = 0;
	    	v2cutsyserrh[i] = 0;
	    }
    }
    // v2 fits sys. err.
    double v2fitsyserrl[Nocent] = {0}, v2fitsyserrh[Nocent] = {0};
    double *d2syserrl, *d2syserrh;
    d2syserrl = d2syserr->GetEYlow();               // delta{2} sys. low
    d2syserrh = d2syserr->GetEYhigh();              // delta{2} sys. high
    for(int i = 0; i<Nocent ; i++) {
        v2fitsyserrl[i] = d2syserrl[i]/(2*sqrt(vv2[i]));
        v2fitsyserrh[i] = d2syserrh[i]/(2*sqrt(vv2[i]));
    }
    ClassSysErr *v2syserr = (ClassSysErr*) SumTwoSSysErr(gv22[0],v2cutsyserrl,v2cutsyserrh,v2fitsyserrl,v2fitsyserrh);

    // delta{2}/v^2 and its sqrt also have two kind of sys. err.: fits and cuts
    ClassSysErr *d2v2syserrfit = new ClassSysErr(v2Nofile-1-3,gd2v22[0],gd2v22+4,"s");					// sqrt(delta/v^2) from diff. fit. func.
    ClassSysErr *d2v2syserrsqfit = new ClassSysErr(v2Nofile-1-3,gd2v22sq[0],gd2v22sq+4,"s");					// delta{2}/v^2 from diff. fit. func.
    double delta[Nocent] = {0}, *tmpdelta;
    double d2v2syserrcutl[Nocent] = {0}, d2v2syserrcuth[Nocent] = {0};						// sqrt(delta/<v^2>)
    double d2v2syserrsqcutl[Nocent] = {0}, d2v2syserrsqcuth[Nocent] = {0};						// delta{2}/<v^2>	
    tmpdelta = gd22sq[0]->GetY();						// delta{2}
    for(int i = 0; i<Nocent; i++) {
	    delta[i] = tmpdelta[i];
	    if(delta[i]>=0 && V22[i]>=0) {
	    	d2v2syserrsqcuth[i] = V_22syserrh[i]*delta[i]/pow(vv2[i],2);		// delta{2}/<v^2>
	    	d2v2syserrcuth[i] = V_22syserrh[i]*sqrt(delta[i])*pow(vv2[i],-3./2.)/2.;		// sqrt(delta{2}/v^2)
	    	d2v2syserrsqcutl[i] = V_22syserrl[i]*delta[i]/pow(vv2[i],2);		// delta{2}/<v^2>
	    	d2v2syserrcutl[i] = V_22syserrl[i]*sqrt(delta[i])*pow(vv2[i],-3./2.)/2.;		// sqrt(delta{2}/v^2)
	    }
	    else {
		cout<<"delta or V22["<<i<<"] < 0!!"<<endl;
	    	d2v2syserrsqcuth[i] = 0;
	    	d2v2syserrcuth[i] = 0;
	    	d2v2syserrsqcutl[i] = 0;
	    	d2v2syserrcutl[i] = 0;
	    }
    }
    ClassSysErr *d2v2syserr = (ClassSysErr*) SumTwoSSysErr(gd2v22[0],d2v2syserrcutl,d2v2syserrcuth,d2v2syserrfit->GetEYlow(),d2v2syserrfit->GetEYhigh());
    ClassSysErr *d2v2syserrsq = (ClassSysErr*) SumTwoSSysErr(gd2v22sq[0],d2v2syserrsqcutl,d2v2syserrsqcuth,d2v2syserrsqfit->GetEYlow(),d2v2syserrsqfit->GetEYhigh());
///    ClassSysErr *d2v2syserr = (ClassSysErr*)SysErrMultSingleS(gd2v22[0],3,gd2v22+1,2,gd2v22+4);
///    ClassSysErr *d2v2syserrsq = (ClassSysErr*)SysErrMultSingleS(gd2v22sq[0],3,gd2v22sq+1,2,gd2v22sq+4);

     // read v3
//    const int v3Nofile = 3;
//    char v3filetag[v3Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit0cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530"};           // 0: exp + linear, 1: exp + gaus, 2: gaus, 3: exp, 8: exp + 4th gaus
//    const int v3Nofile = 6;
//    char v3filetag[v3Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530"};           // 0: exp + linear, 1: exp + gaus, 2: gaus, 3: exp, 8: exp + 4th gaus
    const int v3Nofile = 7;
    char v3filetag[v3Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};           // 0: exp + linear, 1: exp + gaus, 2: gaus, 3: exp, 8: exp + 4th gaus, 4: gaus + linear
    TFile *v3f[Nocent][v3Nofile];
    TH2D *hV32[Nocent][v3Nofile];			// V_{3}{2} vs eta_a, eta_b
    TGraphErrors *igV_32[Nocent][v3Nofile];			// it is eta-gap dep. from V_{3}{2} 2D hist.
    TGraphAsymmErrors *igd32[Nocent][v3Nofile];         // it is eta-gap dep. read in per cent per file
    TGraphErrors *igv32sq[Nocent][v3Nofile];              // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v32[Nocent][v3Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v32sq[Nocent][v3Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphErrors *gV_32[v3Nofile];			// vs cent
    TGraphAsymmErrors *gd32[v3Nofile] ;                  // vs cent
    TGraphAsymmErrors *gd32sq[v3Nofile] ;                  // vs cent
    TGraphErrors *gv32sq[v3Nofile] ;                       // vs cent
    TGraphErrors *gv32[v3Nofile] ;                       // vs cent
    TGraphAsymmErrors *gd2v32[v3Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *gd2v32sq[v3Nofile];        // it is eta-gap dep. read in per cent per file

    double v3scale = 0.5;                               // draw v3 as v3/2

    for(int i = 0; i<v3Nofile ; i++) {
        gd32[i] = new TGraphAsymmErrors();
        gd32sq[i] = new TGraphAsymmErrors();
        gd2v32[i] = new TGraphAsymmErrors();
        gd2v32sq[i] = new TGraphAsymmErrors();
        gv32[i] = new TGraphErrors();
        gv32sq[i] = new TGraphErrors();
        gV_32[i] = new TGraphErrors();
        for(int j = 0; j<Nocent; j++) {
            v3f[j][i] = new TFile(Form("%s.root",Form(v3filetag[i],centtag[j])));
            double y = 0, ey = 0, eyl = 0, eyh = 0;
            igd32[j][i] = (TGraphAsymmErrors*)v3f[j][i]->Get("gd32");
            readYeY4Graph(igd32[j][i],etagap,y,eyl,eyh);
            gd32sq[i]->SetPoint(j,centlocv3[j],y);
            gd32sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd32[i]->SetPoint(j,centlocv3[j],sqrt(y));
            gd32[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
            igd2v32[j][i] = (TGraphAsymmErrors*)v3f[j][i]->Get("gd2v32");
            readYeY4Graph(igd2v32[j][i],etagap,y,eyl,eyh);
            gd2v32sq[i]->SetPoint(j,centlocv3[j],y);
            gd2v32sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd2v32[i]->SetPoint(j,centlocv3[j],sqrt(y));
            gd2v32[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
            igv32sq[j][i] = (TGraphErrors*)v3f[j][i]->Get("gv32");
            readYeY4Graph(igv32sq[j][i],etagap,y,ey);
            gv32[i]->SetPoint(j,centlocv3[j],sqrt(y));
            gv32[i]->SetPointError(j,0,ey/(2*sqrt(y)));
            gv32sq[i]->SetPoint(j,centlocv3[j],y);
            gv32sq[i]->SetPointError(j,0,ey);
            hV32[j][i] = (TH2D*)v3f[j][i]->Get("hV32");
            igV_32[j][i] = (TGraphErrors*) EtaGapDepfrom2Dhist(hV32[j][i]);
            readYeY4Graph(igV_32[j][i],etagap,y,ey);
            gV_32[i]->SetPoint(j,centlocv3[j],y);
            gV_32[i]->SetPointError(j,0,ey);
        }
    }

//    ClassSysErr *v3syserr = new ClassSysErr(3,gv32[0],gv32+1,"m");
//    ClassSysErr *d3syserr = new ClassSysErr(2,gd32[0],gd32+4,"s");
//    ClassSysErr *d2v3syserr = new ClassSysErr(v3Nofile-1,gd2v32[0],gd2v32+1);          // sigma/v
//    ClassSysErr *d2v3syserrsq = new ClassSysErr(v3Nofile-1,gd2v32sq[0],gd2v32sq+1);          // (sigma/v)^2
///    ClassSysErr *v3syserr = (ClassSysErr*)SysErrMultSingleS(gv32[0],3,gv32+1,2,gv32+4);
///    ClassSysErr *d3syserr = (ClassSysErr*)SysErrMultSingleS(gd32[0],3,gd32+1,2,gd32+4);
    ClassSysErr *d3syserr = new ClassSysErr(v3Nofile-1-3,gd32[0],gd32+4,"s");              // delta only fit function sys. err.
    ClassSysErr *V_32syserr = new ClassSysErr(3,gV_32[0],gV_32+1,"m");		// V_{3}{2} only cut sys. err.
    // <v^2> has two source of sys. err. 
    double V_32syserrl[Nocent] = {0},v3cutsyserrl[Nocent] = {0}, *tmpv32sel;
    double V_32syserrh[Nocent] = {0},v3cutsyserrh[Nocent] = {0}, *tmpv32seh;
    double v32[Nocent] = {0}, *tmpv32;
    double *vv3;
    vv3 = gv32sq[0]->GetY();
    // v3 cuts sys. err
    tmpv32 = gv32sq[0]->GetY();
    tmpv32sel = V_32syserr->GetEYlow();
    tmpv32seh = V_32syserr->GetEYhigh();
    for(int i = 0; i<V_32syserr->GetN();i++) {
	    v32[i] = tmpv32[i];							// v32
	    V_32syserrl[i] = tmpv32sel[i];					// error for v32
	    V_32syserrh[i] = tmpv32seh[i];
	    if(vv3[i]>=0) {
	    	v3cutsyserrl[i] = tmpv32sel[i]/(2*sqrt(vv3[i]));			// error for sqrt(v32)
	    	v3cutsyserrh[i] = tmpv32seh[i]/(2*sqrt(vv3[i]));
	    }
	    else {
		cout<<"vv3["<<i<<"] < 0!!"<<endl;
	    	v3cutsyserrl[i] = 0;
	    	v3cutsyserrh[i] = 0;
	    }
    }
    // v3 fits sys. err.
    double v3fitsyserrl[Nocent] = {0}, v3fitsyserrh[Nocent] = {0};
    double *d3syserrl, *d3syserrh;
    d3syserrl = d3syserr->GetEYlow();               // delta{2} sys. low
    d3syserrh = d3syserr->GetEYhigh();              // delta{2} sys. high
    for(int i = 0; i<Nocent ; i++) {
        v3fitsyserrl[i] = d3syserrl[i]/(2*sqrt(vv3[i]));
        v3fitsyserrh[i] = d3syserrh[i]/(2*sqrt(vv3[i]));
    }
    ClassSysErr *v3syserr = (ClassSysErr*) SumTwoSSysErr(gv32[0],v3cutsyserrl,v3cutsyserrh,v3fitsyserrl,v3fitsyserrh);

    // delta{2}/v^2 and its sqrt also have two kind of sys. err.: fits and cuts
    ClassSysErr *d2v3syserrfit = new ClassSysErr(v3Nofile-1-3,gd2v32[0],gd2v32+4,"s");					// sqrt(delta/v^2) from diff. fit. func.
    ClassSysErr *d2v3syserrsqfit = new ClassSysErr(v3Nofile-1-3,gd2v32sq[0],gd2v32sq+4,"s");					// delta{2}/v^2 from diff. fit. func.
    double delta[Nocent] = {0}, *tmpdelta;
    double d2v3syserrcutl[Nocent] = {0}, d2v3syserrcuth[Nocent] = {0};						// sqrt(delta/<v^2>)
    double d2v3syserrsqcutl[Nocent] = {0}, d2v3syserrsqcuth[Nocent] = {0};						// delta{2}/<v^2>	
    tmpdelta = gd32sq[0]->GetY();						// delta{2}
    for(int i = 0; i<Nocent; i++) {
	    delta[i] = tmpdelta[i];
	    if(delta[i]>=0 && v32[i]>=0) {
	    	d2v3syserrsqcuth[i] = V_32syserrh[i]*delta[i]/pow(vv3[i],2);		// delta{2}/<v^2>
	    	d2v3syserrcuth[i] = V_32syserrh[i]*sqrt(delta[i])*pow(vv3[i],-3./2.)/2.;		// sqrt(delta{2}/v^2)
	    	d2v3syserrsqcutl[i] = V_32syserrl[i]*delta[i]/pow(vv3[i],2);		// delta{2}/<v^2>
	    	d2v3syserrcutl[i] = V_32syserrl[i]*sqrt(delta[i])*pow(vv3[i],-3./2.)/2.;		// sqrt(delta{2}/v^2)
	    }
	    else {
		cout<<"delta or v32["<<i<<"] < 0!!"<<endl;
	    	d2v3syserrsqcuth[i] = 0;
	    	d2v3syserrcuth[i] = 0;
	    	d2v3syserrsqcutl[i] = 0;
	    	d2v3syserrcutl[i] = 0;
	    }
    }
    ClassSysErr *d2v3syserr = (ClassSysErr*) SumTwoSSysErr(gd2v32[0],d2v3syserrcutl,d2v3syserrcuth,d2v3syserrfit->GetEYlow(),d2v3syserrfit->GetEYhigh());
    ClassSysErr *d2v3syserrsq = (ClassSysErr*) SumTwoSSysErr(gd2v32sq[0],d2v3syserrsqcutl,d2v3syserrsqcuth,d2v3syserrsqfit->GetEYlow(),d2v3syserrsqfit->GetEYhigh());
///    ClassSysErr *d2v3syserr = (ClassSysErr*)SysErrMultSingleS(gd2v32[0],3,gd2v32+1,2,gd2v32+4);
///    ClassSysErr *d2v3syserrsq = (ClassSysErr*)SysErrMultSingleS(gd2v32sq[0],3,gd2v32sq+1,2,gd2v32sq+4);
//
//
//
//
//    TFile *v3f[Nocent][v3Nofile];
//    TGraphAsymmErrors *igd32[Nocent][v3Nofile];         // it is eta-gap dep. read in per cent per file
//    TGraphErrors *igv32[Nocent][v3Nofile];              // it is eta-gap dep. read in per cent per file
//    TGraphAsymmErrors *igd2v32[Nocent][v3Nofile];         // it is eta-gap dep. read in per cent per file
//    TGraphAsymmErrors *igd2v32sq[Nocent][v3Nofile];         // it is eta-gap dep. read in per cent per file
//    TGraphAsymmErrors *gd32[v3Nofile] = new TGraphAsymmErrors();                  // vs cent
//    TGraphAsymmErrors *gd2v32[v3Nofile] = new TGraphAsymmErrors();                  // vs cent
//    TGraphAsymmErrors *gd2v32sq[v3Nofile] = new TGraphAsymmErrors();                  // vs cent
//    TGraphErrors *gv32[v3Nofile] = new TGraphErrors();                       // vs cent
//
//    for(int i = 0; i<v3Nofile ; i++) {
//        gd32[i] = new TGraphAsymmErrors();
//        gd2v32[i] = new TGraphAsymmErrors();
//        gd2v32sq[i] = new TGraphAsymmErrors();
//        gv32[i] = new TGraphErrors();
//        for(int j = 0; j<Nocent; j++) {
//            v3f[j][i] = new TFile(Form("%s.root",Form(v3filetag[i],centtag[j])));
//            double y = 0, ey = 0, eyl = 0, eyh = 0;
//            igd32[j][i] = (TGraphAsymmErrors*)v3f[j][i]->Get("gd32");
//            readYeY4Graph(igd32[j][i],etagap,y,eyl,eyh);
//            gd32[i]->SetPoint(j,centlocv3[j],sqrt(y));
//            gd32[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
//            igd2v32[j][i] = (TGraphAsymmErrors*)v3f[j][i]->Get("gd2v32");
//            readYeY4Graph(igd2v32[j][i],etagap,y,eyl,eyh);
//            gd2v32sq[i]->SetPoint(j,centlocv3[j],y);
//            gd2v32sq[i]->SetPointError(j,0,0,eyl,eyh);
//            gd2v32[i]->SetPoint(j,centlocv3[j],sqrt(y));
//            gd2v32[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
//            igv32[j][i] = (TGraphErrors*)v3f[j][i]->Get("gv32");
//            readYeY4Graph(igv32[j][i],etagap,y,ey);
//            gv32[i]->SetPoint(j,centlocv3[j],sqrt(y));
//            gv32[i]->SetPointError(j,0,ey/(2*sqrt(y)));
//        }
//    }
//
////    ClassSysErr *v3syserr = new ClassSysErr(v3Nofile-1,gv32[0],gv32+1);
////    ClassSysErr *d3syserr = new ClassSysErr(v3Nofile-1,gd32[0],gd32+1);
////    ClassSysErr *d2v3syserr = new ClassSysErr(v3Nofile-1,gd2v32[0],gd2v32+1);
////    ClassSysErr *d2v3syserrsq = new ClassSysErr(v3Nofile-1,gd2v32sq[0],gd2v32sq+1);
//    ClassSysErr *v3syserr = (ClassSysErr*)SysErrMultSingleS(gv32[0],3,gv32+1,2,gv32+4);
//    ClassSysErr *d3syserr = (ClassSysErr*)SysErrMultSingleS(gd32[0],3,gd32+1,2,gd32+4);
////    ClassSysErr *d3syserr = new ClassSysErr(2,gd32[0],gd32+4,"s");         // delta only fit function sys. err.
//    ClassSysErr *d2v3syserr = (ClassSysErr*)SysErrMultSingleS(gd2v32[0],3,gd2v32+1,2,gd2v32+4);
//    ClassSysErr *d2v3syserrsq = (ClassSysErr*)SysErrMultSingleS(gd2v32sq[0],3,gd2v32sq+1,2,gd2v32sq+4);


    char outname[500] = "";
    sprintf(outname,"vdeltavscent.%s.%s",Form(v2filetag[0],"all"),Form(v3filetag[0],"all"));
    
    // setup enviroment for drawing

    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

    TH2D *htmp = new TH2D("htmp","",8,0,80,100,0,5);    // x-range, y-range go here
    htmp->GetXaxis()->SetTitle("centrality %");
//    htmp->GetXaxis()->CenterTitle(0);
    htmp->GetYaxis()->SetTitle("#sqrt{#LT v^{2} #GT} and #sqrt{#bar{D}} (%)");
    
    FigFormat1D(c,pad,htmp);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();

    wLatex(0.25, 0.94, "|#Delta#eta|>0.7",1,12,0.06);

    wLatex(0.05,0.91,"(c)",1,12,0.08);              // figure's subcaption

    int v2color=1, v3color=2;
    int v2marker=4, d2marker=8;
    int v3marker=26, d3marker=22;

    v2syserr->ScaleMeanYOnly(v2scale);
    graphSystBox(d2syserr->GetN(),d2syserr->GetX(),d2syserr->GetY(),0,0,d2syserr->GetEYhigh(),d2syserr->GetEYlow(),kGray+1,100,0);
    graphSystBox(v2syserr->GetN(),v2syserr->GetX(),v2syserr->GetY(),0,0,v2syserr->GetEYhigh(),v2syserr->GetEYlow(),kGray,100,0);

    graphSystBox(d3syserr->GetN(),d3syserr->GetX(),d3syserr->GetY(),0,0,d3syserr->GetEYhigh(),d3syserr->GetEYlow(),kYellow,100,0);
    graphSystBox(v3syserr->GetN(),v3syserr->GetX(),v3syserr->GetY(),0,0,v3syserr->GetEYhigh(),v3syserr->GetEYlow(),kOrange,100,0);

    // need to scale v2 on the plot
    TGraphErrors* gv22cp = new TGraphErrors();
    double *gv220x = gv22[0]->GetX();
    double *gv220y = gv22[0]->GetY();
    double *gv220ey = gv22[0]->GetEY();
    for(i = 0; i<gv22[0]->GetN() ; i++) {
        gv22cp->SetPoint(i,gv220x[i],gv220y[i]*v2scale);
        gv22cp->SetPointError(i,0,gv220ey[i]);
    }
    gv22cp->SetMarkerStyle(v2marker);
    gv22cp->SetMarkerSize(2);
    gv22cp->SetMarkerColor(v2color);
    gv22cp->SetLineColor(v2color);
    gd22[0]->SetMarkerStyle(d2marker);
    gd22[0]->SetMarkerSize(2);
    gd22[0]->SetMarkerColor(v2color);
    gd22[0]->SetLineColor(v2color);
    gv32[0]->SetMarkerStyle(v3marker);
    gv32[0]->SetMarkerSize(2);
    gv32[0]->SetMarkerColor(v3color);
    gv32[0]->SetLineColor(v3color);
    gd32[0]->SetMarkerStyle(d3marker);
    gd32[0]->SetMarkerSize(2);
    gd32[0]->SetMarkerColor(v3color);
    gd32[0]->SetLineColor(v3color);

    gd22[0]->Draw("epsame");
    gv22cp->Draw("epsame");
    gd32[0]->Draw("epsame");
    gv32[0]->Draw("epsame");

    double xtxt = 0.3;
    double ytxt = 0.85, ystep = 0.08;
    double tsize = 0.05;
    double markersize = 1.5;
    keySymbol(xtxt,ytxt,"#sqrt{#LT v_{2}^{2} #GT}/2",v2color,v2marker,tsize,markersize);
    keySymbol(xtxt,ytxt-ystep,"#sqrt{#bar{D}_{2}}",v2color,d2marker,tsize,markersize);
    keySymbol(xtxt,ytxt-ystep*2,"#sqrt{#LT v_{3}^{2} #GT}",v3color,v3marker,tsize,markersize);
    keySymbol(xtxt,ytxt-ystep*3,"#sqrt{#bar{D}_{3}}",v3color,d3marker,tsize,markersize);


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
    }



    cout<<"|#Delta#eta|>"<<etagap<<" sqrt(d_{2}): "<<endl;
    double *tmppy, *tmpelpy, *tmpehpy, *tmpeslpy, *tmpeshpy;
    tmppy = gd22[0]->GetY();
    tmpelpy = gd22[0]->GetEYlow();
    tmpehpy = gd22[0]->GetEYhigh();
    tmpeslpy = d2syserr->GetEYlow();
    tmpeshpy = d2syserr->GetEYhigh();
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<tmppy[i]<<" + "<<tmpehpy[i]<<" - "<<tmpelpy[i]<<"(stat.) + "<<tmpeshpy[i]<<" "<<tmpeslpy[i]<<"(sys.)"<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<" sqrt(d_{2})/<v_{2}>: "<<endl;
    double *tmppy, *tmpelpy, *tmpehpy, *tmpeslpy, *tmpeshpy;
    tmppy = gd2v22[0]->GetY();
    tmpelpy = gd2v22[0]->GetEYlow();
    tmpehpy = gd2v22[0]->GetEYhigh();
    tmpeslpy = d2v2syserr->GetEYlow();
    tmpeshpy = d2v2syserr->GetEYhigh();
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<tmppy[i]*100<<"% + "<<tmpehpy[i]*100<<"% - "<<tmpelpy[i]*100<<"%(stat.) + "<<tmpeshpy[i]*100<<"% "<<tmpeslpy[i]*100<<"%(sys.)"<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<" d_{2}/<v_{2}>: "<<endl;
    double *tmppy, *tmpelpy, *tmpehpy, *tmpeslpy, *tmpeshpy;
    tmppy = gd2v22sq[0]->GetY();
    tmpelpy = gd2v22sq[0]->GetEYlow();
    tmpehpy = gd2v22sq[0]->GetEYhigh();
    tmpeslpy = d2v2syserrsq->GetEYlow();
    tmpeshpy = d2v2syserrsq->GetEYhigh();
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<tmppy[i]*100<<"% + "<<tmpehpy[i]*100<<"% - "<<tmpelpy[i]*100<<"%(stat.) + "<<tmpeshpy[i]*100<<"% "<<tmpeslpy[i]*100<<"%(sys.)"<<endl;
    }

}





//------------------------------------- 2013.12.02 ly--------------------------------------
// D fit parameters errors are treated as sys. err. instead of sys. err. as asked by 2nd GPC comments.
void plotDcent(double etagap, char symbol[], int vn) {                    // centrality dependence for \delta_{vn} = D - \sigma'

    // read the data
    const int Nocent = 5;
    char centtag[Nocent][100] = {"678","5","4","3","012"};
    int centloc[Nocent] = {10,25,35,45,65};         // %
    double centlocv2[Nocent] = {0};         // %               v2
    double shift = 0.;
    for(int i = 0 ;i <Nocent  ; i++ ) {
        centlocv2[i] = centloc[i]+shift;
    }
    // read v2
    const int v2Nofile = 7;
    char v22filetag[v2Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit8cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};
    const int v3Nofile = 7;					// if not the same with v2Nofile, need to change v2Nofile to be largest
    char v32filetag[v3Nofile][200] = {"fit1cent%seta20dca3_3dataset.130530","fit1cent%seta20dca2_3dataset.130506","fit1cent%seta20dca3Vz25_3dataset.130530","fit1cent%seta20dca3hit15_3dataset.130328","fit0cent%seta20dca3_3dataset.130530","fit3cent%seta20dca3_3dataset.130530","fit4cent%seta20dca3_3dataset.130530"};           // 0: exp + linear, 1: exp + gaus, 2: gaus, 3: exp, 8: exp + 4th gaus, 4: gaus + linear
    char v2filetag[v2Nofile][200];
    if(vn==2) {
        for(int i = 0; i<v2Nofile; i++ ) {
            sprintf(v2filetag[i],"%s",v22filetag[i]);
        }
    }
    else if(vn==3) {
        for(int i = 0; i<v2Nofile; i++ ) {
            sprintf(v2filetag[i],"%s",v32filetag[i]);
        }
    }
    else { cout<<"ERR!! Wrong vn ="<<vn<<" input in function plotDcent()"<<endl; return;}

    TFile *v2f[Nocent][v2Nofile];
    TH2D *hV22[Nocent][v2Nofile];			// V_{vn}{2} vs eta_a, eta_b
    TGraphErrors *igV_22[Nocent][v2Nofile];			// it is eta-gap dep. from V_{vn}{2} 2D hist.
    TGraphAsymmErrors *igd22[Nocent][v2Nofile];         // it is eta-gap dep. read in per cent per file
    TGraphErrors *igv22sq[Nocent][v2Nofile];              // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v22[Nocent][v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *igd2v22sq[Nocent][v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphErrors *gV_22[v2Nofile];			// vs cent
    TGraphAsymmErrors *gd22[v2Nofile] ;                  // vs cent
    TGraphAsymmErrors *gd22sq[v2Nofile] ;                  // vs cent
    TGraphErrors *gv22sq[v2Nofile] ;                       // vs cent
    TGraphErrors *gv22[v2Nofile] ;                       // vs cent
    TGraphAsymmErrors *gd2v22[v2Nofile];        // it is eta-gap dep. read in per cent per file
    TGraphAsymmErrors *gd2v22sq[v2Nofile];        // it is eta-gap dep. read in per cent per file

    TGraphErrors *gDmsigmaprime[v2Nofile];               // \delta = D - sigma' vs cent
    TGraphErrors *gsqrtDmsigmaprime[v2Nofile];               // Sqrt(\delta = D - sigma') vs cent

    TGraphErrors *gD2v2 = new TGraphErrors();			// D/(V{2} - D)	vs centrality
    TGraphErrors *gsqrtD2v2 = new TGraphErrors();			// sqrt( D/(V{2} - D) )	vs centrality
    
    // read the data for sigma'
    const int v4Nofile = 4;
    char v4filetag[v4Nofile][500] = {"fit1cent%seta10dca3_3dataset.130603","fit1cent%seta10dca2_3dataset.130601","fit1cent%seta10dca3Vz25_3dataset.130603","fit1cent%seta10dca3hit15_3dataset.130605"};             // files are in the same order as V2{2} for later sigma calculation
    TFile *v4f[Nocent][v4Nofile];
    TGraphErrors *igavedd24[Nocent][v4Nofile];            // used for \sigma' in V{4}
    TGraphErrors *gsigmaprime[v4Nofile];

    double v2scale = 0.5;                               // draw v2 as v2/2
    if(vn==3) v2scale = 1;				// draw v3 as v3

    // read from v2 file
    for(int i = 0; i<v2Nofile ; i++) {
        gd22[i] = new TGraphAsymmErrors();
        gd22sq[i] = new TGraphAsymmErrors();
        gd2v22[i] = new TGraphAsymmErrors();
        gd2v22sq[i] = new TGraphAsymmErrors();
        gv22[i] = new TGraphErrors();
        gv22sq[i] = new TGraphErrors();
        gV_22[i] = new TGraphErrors();
        gDmsigmaprime[i] = new TGraphErrors();					// D - \sigma'
        gsqrtDmsigmaprime[i] = new TGraphErrors();				// sqrt(D-\sigma')

        for(int j = 0; j<Nocent; j++) {
            v2f[j][i] = new TFile(Form("%s.root",Form(v2filetag[i],centtag[j])));
            double y = 0, ey = 0, eyl = 0, eyh = 0;
            igd22[j][i] = (TGraphAsymmErrors*)v2f[j][i]->Get(Form("gd%d2",vn));
            readYeY4Graph(igd22[j][i],etagap,y,eyl,eyh);
            gd22sq[i]->SetPoint(j,centlocv2[j],y);
            gd22sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gd22[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
//            if(i==0) {cout<<"gd22 = "<<sqrt(y)<<"+-"<<eyh/(2*sqrt(y)) <<"\tgd22sq = "<<y<<"+-"<<eyh<<endl;}
            igd2v22[j][i] = (TGraphAsymmErrors*)v2f[j][i]->Get(Form("gd2v%d2",vn));
            readYeY4Graph(igd2v22[j][i],etagap,y,eyl,eyh);
            gd2v22sq[i]->SetPoint(j,centlocv2[j],y);
            gd2v22sq[i]->SetPointError(j,0,0,eyl,eyh);
            gd2v22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gd2v22[i]->SetPointError(j,0,0,eyl/(2*sqrt(y)),eyh/(2*sqrt(y)));
            igv22sq[j][i] = (TGraphErrors*)v2f[j][i]->Get(Form("gv%d2",vn));
            readYeY4Graph(igv22sq[j][i],etagap,y,ey);
            gv22[i]->SetPoint(j,centlocv2[j],sqrt(y));
            gv22[i]->SetPointError(j,0,ey/(2*sqrt(y)));
            gv22sq[i]->SetPoint(j,centlocv2[j],y);
            gv22sq[i]->SetPointError(j,0,ey);
            hV22[j][i] = (TH2D*)v2f[j][i]->Get(Form("hV%d2",vn));
            igV_22[j][i] = (TGraphErrors*) EtaGapDepfrom2Dhist(hV22[j][i]);
            readYeY4Graph(igV_22[j][i],etagap,y,ey);
            gV_22[i]->SetPoint(j,centlocv2[j],y);
            gV_22[i]->SetPointError(j,0,ey);
        }
    }

    // read from V4 file
    TF1 *fsigma = new TF1("fsigma","x*(2-[0])/2",0,2);                // sigma' vs eta-gap, [0] for eta-gap. x for slope k
    fsigma->SetParameter(0,etagap);        // eta-gap 
    for(int i = 0; i<v4Nofile; i++) {
        gsigmaprime[i] = new TGraphErrors();
        for(int j = 0; j<Nocent; j++) {
            v4f[j][i] = new TFile(Form("%s.root",Form(v4filetag[i],centtag[j])));
            igavedd24[j][i] = (TGraphErrors*)v4f[j][i]->Get(Form("gavedd%d4",vn));
            double k = igavedd24[j][i]->GetFunction("f0D")->GetParameter(0);
            double ek = igavedd24[j][i]->GetFunction("f0D")->GetParError(0);
            double y = fsigma->Eval(k);
            double ey = fabs(fsigma->Eval(k+ek) - fsigma->Eval(k));
            gsigmaprime[i]->SetPoint(j,centlocv2[j],y);
            gsigmaprime[i]->SetPointError(j,0,ey);
        }
    }


    // D sys. err. from par. fit err. + diff. fit func. err. as asked by GPC 2nd round comments			2013.12.01 ly
    // function D doesn't have stat. err. 
    double Delta[Nocent] = {0}, *tmpDelta;
    tmpDelta = gd22sq[0]->GetY();
    for(int i = 0; i<gd22sq[0]->GetN(); i++) {
        Delta[i] = tmpDelta[i];
    }
    ClassSysErr *D2funsyserr = new ClassSysErr(v2Nofile-1-3,gd22sq[0],gd22sq+4,"s");		// D sys. err. from diff. fit func.
    ClassSysErr *D2syserr = (ClassSysErr*)SumTwoSSysErr(gd22sq[0],gd22sq[0]->GetEYlow(),gd22sq[0]->GetEYhigh(),D2funsyserr->GetEYlow(),D2funsyserr->GetEYhigh());						// D has two sources of sys. err. : fit par. + diff. func.
    double D2syserrl[Nocent], D2syserrh[Nocent], *tmpD2syserrl, *tmpD2syserrh;
    tmpD2syserrl = D2syserr->GetEYlow();               // D{2} sys. low
    tmpD2syserrh = D2syserr->GetEYhigh();              // D{2} sys. high
    for(int i = 0; i<D2syserr->GetN() ; i++) {
        D2syserrl[i] = tmpD2syserrl[i];
        D2syserrh[i] = tmpD2syserrh[i];
    }
    // sqrt(D)
    double sqrtDelta[Nocent] = {0};
    double sqrtD2syserrl[Nocent] = {0}, sqrtD2syserrh[Nocent] = {0};
    for(int i = 0; i<Nocent ; i++) {
	    if(Delta[i]>0) {
	    	sqrtDelta[i] = sqrt(Delta[i]);
	    	sqrtD2syserrl[i] = D2syserrl[i]/(2.*sqrt(Delta[i]));
	    	sqrtD2syserrh[i] = D2syserrh[i]/(2.*sqrt(Delta[i]));
	    }    
    }
    ClassSysErr *sqrtD2syserr = new ClassSysErr(Nocent, centlocv2, sqrtDelta, sqrtD2syserrl, sqrtD2syserrh);		// sqrt(D) sys. err. is from D


    // sigma' sys. err.
    ClassSysErr *sigmaprimesyserr;        // sigma' only has fit par. err. as sys. err. 
    double sigmaprime[Nocent] = {0}, *tmpsigmaprime;
    tmpsigmaprime = gsigmaprime[0]->GetY();
    for( int i = 0; i<gsigmaprime[0]->GetN(); i++ ) {
        sigmaprime[i] = tmpsigmaprime[i];
    }
    double sigmaprimeerrl[Nocent], *tmpsigmaprimeerr;
    tmpsigmaprimeerr = gsigmaprime[0]->GetEY();
    for( int ip = 0; ip<gsigmaprime[0]->GetN(); ip++ ) {
	    sigmaprimeerrl[ip] = -fabs(tmpsigmaprimeerr[ip]);					// ClassSysErr: low err is defined as negative
    }
    sigmaprimesyserr = new ClassSysErr(gsigmaprime[0]->GetN(), gsigmaprime[0]->GetX(), gsigmaprime[0]->GetY(), sigmaprimeerrl, gsigmaprime[0]->GetEY());
    double sigmaprimesyserrl[Nocent], sigmaprimesyserrh[Nocent], *tmpsigmaprimesyserrl, *tmpsigmaprimesyserrh; 
    tmpsigmaprimesyserrl = sigmaprimesyserr->GetEYlow();
    tmpsigmaprimesyserrh = sigmaprimesyserr->GetEYhigh();
    for(int i = 0; i<Nocent ; i++) {
        sigmaprimesyserrl[i] = tmpsigmaprimesyserrl[i];
        sigmaprimesyserrh[i] = tmpsigmaprimesyserrh[i];
    }



    // <v^2>
    // stat. err. only from V{2}
    double V_22[Nocent] = {0}, *tmpV_22;					// V_2{2}
    double eV_22[Nocent] = {0}, *tmpeV_22;					// V_2{2} stat. err.
    double vv2[Nocent] = {0}, *tmpvv2;						// v^2
    double evv2[Nocent] = {0};							// v^2 stat. err.
    double sqrtvv2[Nocent] = {0}, esqrtvv2[Nocent] = {0};			// sqrt(v^2)
    tmpvv2 = gv22sq[0]->GetY();							// v^2, <v^2> graph error in root file contains err. from D. which is treated as sys. err. as GPC 2nd comment 2013.12.01 ly
    tmpV_22 = gV_22[0]->GetY();							// V_2{2}
    tmpeV_22 = gV_22[0]->GetEY();						// V_2{2} stat. err.
    for(int i = 0; i<gV_22[0]->GetN();i++) {
	    V_22[i] = tmpV_22[i];							// V_2{2}
	    eV_22[i] = tmpeV_22[i];							// V_2{2} stat. err.
	    vv2[i] = tmpvv2[i];
	    evv2[i] = eV_22[i];							// v^2 stat. err. is from V{2}
	    if(vv2[i]>=0) {
	    	sqrtvv2[i] = sqrt(vv2[i]);						// sqrt(<v^2>) 
	    	esqrtvv2[i] = evv2[i]/(2*sqrt(vv2[i]));				// stat. error for sqrt(V22-D) from V22
	    }
    }


    // v^2 sys. err. : v^2 = V{2} - D. so from V{2} and D
    ClassSysErr *V_22syserr = new ClassSysErr(3,gV_22[0],gV_22+1,"m");		// V_{2}{2} only cut sys. err.
//    for(int i = 0; i<Nocent; i++) {cout<<"d2syserr_Y["<<i<<"] = "<<d2sqsyserr->GetY()[i]<<"+-"<<d2sqsyserr->GetEYhigh()[i]<<"\td2syserr_E="<<d2syserr->GetEYhigh()[i]<<endl;}
    // sqrt( <v^2> ) has two source of sys. err. : from V{2} and from \delta
    // V{2} sys. err. is from cuts
    double V_22syserrl[Nocent] = {0},v2cutsyserrl[Nocent] = {0}, *tmpV22sel;
    double V_22syserrh[Nocent] = {0},v2cutsyserrh[Nocent] = {0}, *tmpV22seh;
    tmpV22sel = V_22syserr->GetEYlow();
    tmpV22seh = V_22syserr->GetEYhigh();
    for(int i = 0; i<V_22syserr->GetN();i++) {
	    V_22syserrl[i] = tmpV22sel[i];					// v22 error from V22
	    V_22syserrh[i] = tmpV22seh[i];
	    if(vv2[i]>=0) {
	    	v2cutsyserrl[i] = tmpV22sel[i]/(2*sqrt(vv2[i]));			// sys. error for sqrt(V22-D) from V22
	    	v2cutsyserrh[i] = tmpV22seh[i]/(2*sqrt(vv2[i]));
	    }
	    else {
		cout<<"vv2["<<i<<"] < 0!!"<<endl;
	    	v2cutsyserrl[i] = 0;
	    	v2cutsyserrh[i] = 0;
	    }
    }
    // sqrt( <v^2> ) sys. err. from D : fit par. err. + diff. fit func.			2013.11.29	ly
    double v2fitsyserrl[Nocent] = {0}, v2fitsyserrh[Nocent] = {0};
    for(int i = 0; i<Nocent ; i++) {
        v2fitsyserrl[i] = D2syserrl[i]/(2*sqrt(vv2[i]));
        v2fitsyserrh[i] = D2syserrh[i]/(2*sqrt(vv2[i]));
    }
    ClassSysErr *v2syserr = (ClassSysErr*) SumTwoSSysErr(gv22[0],v2cutsyserrl,v2cutsyserrh,v2fitsyserrl,v2fitsyserrh);
    double v2syserrl[Nocent] = {0}, v2syserrh[Nocent] = {0}, *tmpv2syserrl, *tmpv2syserrh;
    tmpv2syserrl = v2syserr->GetEYlow();
    tmpv2syserrh = v2syserr->GetEYhigh();
    for(int i = 0; i<v2syserr->GetN(); i++) {
        v2syserrl[i] = tmpv2syserrl[i];
        v2syserrh[i] = tmpv2syserrh[i];
    }



    // \delta = D-\sigma' and sqrt(\delta) for each file
    for(int  i = 0; i<v2Nofile; i++) {
        double *tmpDarray, *tmpsigmaarray;
        tmpDarray = gd22sq[i]->GetY();
        int iv4 = i>(v4Nofile-1)?0:i;                   // v2Nofile has diff fit than v4Nofile
        
        tmpsigmaarray = gsigmaprime[iv4]->GetY();
        for(int j = 0; j<Nocent ; j++ ) {
            double yDmsigmaprime = tmpDarray[j]-tmpsigmaarray[j];
            gDmsigmaprime[i]->SetPoint(j,centlocv2[j],yDmsigmaprime);
            gDmsigmaprime[i]->SetPointError(j,0,0);							// 2013.11.29  no stat. err. for function D-\sigma'
//            if(i==0)             cout<<"Dmsigmaprime["<<j<<"] = "<<yDmsigmaprime<<endl;
            if(yDmsigmaprime>=0) {
                gsqrtDmsigmaprime[i]->SetPoint(j,centlocv2[j],sqrt(yDmsigmaprime));
                gsqrtDmsigmaprime[i]->SetPointError(j,0,0);							//  2013.11.29 no stat. err. for function sqrt(D-\sigma) 
            }
        }
    }


    // \delta = D - sigma' sys. err.  from D and sigma'
    // \delta doesn't have stat. err.
    double Dmsigmaprime[Nocent] = {0}, *tmpDmsigmaprime;
    tmpDmsigmaprime = gDmsigmaprime[0]->GetY();
    for(int i = 0 ; i<gDmsigmaprime[0]->GetN(); i++) {
        Dmsigmaprime[i] = tmpDmsigmaprime[i];
    }
    ClassSysErr *Dmsigmaprimesyserr = (ClassSysErr*) SumTwoSSysErr(gDmsigmaprime[0],D2syserr->GetEYlow(),D2syserr->GetEYhigh(),sigmaprimesyserr->GetEYlow(),sigmaprimesyserr->GetEYhigh());			// D-\sigma' sys. err. from D + \sigma'					// 2013.11.30		ly



    // sqrt(D-\sigma')
    double *Dmsigmaprimesysl=Dmsigmaprimesyserr->GetEYlow(), *Dmsigmaprimesysh=Dmsigmaprimesyserr->GetEYhigh();
    double sqrtDmsigmaprime[Nocent];
    double sqrtDmsigmaprimesysl[Nocent], sqrtDmsigmaprimesysh[Nocent];
//    double *tmpsqrtDmsigmaprimesysl, *tmpsqrtDmsigmaprimesysh;
    for(int i = 0; i<Nocent; i++) {
        if(Dmsigmaprime[i]>=0) {
            sqrtDmsigmaprime[i] = sqrt(Dmsigmaprime[i]);
            sqrtDmsigmaprimesysl[i] = Dmsigmaprimesysl[i]/(2.*sqrt(Dmsigmaprime[i]));					// this is a negative number
            sqrtDmsigmaprimesysh[i] = Dmsigmaprimesysh[i]/(2.*sqrt(Dmsigmaprime[i]));
//            cout<<"Dmsigmaprime["<<i<<"] = "<<Dmsigmaprime[i]<<"\tDmsigmaprimesys_Eh="<<Dmsigmaprimesysh[i]<<"\tsqrtDmsigmaprime_Eh ="<<sqrtDmsigmaprimesysh[i]<<endl;
        }
        else {
            sqrtDmsigmaprimesysl[i] = 0;
            sqrtDmsigmaprimesysh[i] = 0;
        }
    }
//    tmpsqrtDmsigmaprimesysl = sqrtDmsigmaprimesysl;
//    tmpsqrtDmsigmaprimesysh = sqrtDmsigmaprimesysh;
//    ClassSysErr *sqrtDmsigmaprimesyserr = new ClassSysErr(Nocent, centlocv2, gsqrtDmsigmaprime[0]->GetY(), tmpsqrtDmsigmaprimesysl,tmpsqrtDmsigmaprimesysh);
    ClassSysErr *sqrtDmsigmaprimesyserr = new ClassSysErr(Nocent, centlocv2, gsqrtDmsigmaprime[0]->GetY(), sqrtDmsigmaprimesysl,sqrtDmsigmaprimesysh);



    // D{2}/v^2 = D / ( V{2}-D )    and     sqrt( D/(V{2}-D) )

    // D/v^2 stat. err. only from V{2}
    double D2v2[Nocent] = {0}, eD2v2[Nocent] = {0};						// D/(V{2}-D)
    double sqrtD2v2[Nocent] = {0}, esqrtD2v2[Nocent] = {0};					// sqrt( D/(V{2}-D) )
    for(int i = 0; i<V_22syserr->GetN();i++) {
	    if((V_22[i]-Delta[i])!=0) {
	    	D2v2[i] = Delta[i]/(V_22[i]-Delta[i]);
		    eD2v2[i] = eV_22[i]*fabs(-Delta[i]/pow(V_22[i]-Delta[i],2));			// D/(V{2}-D) stat. err.
		    gD2v2->SetPoint(i, centlocv2[i], D2v2[i]);
		    gD2v2->SetPointError(i, 0, eD2v2[i]);
		    if(D2v2[i]>0) {
		    	sqrtD2v2[i] = sqrt(D2v2[i]);						// sqrt( D/(V{2}-D) )
		    	esqrtD2v2[i] = eD2v2[i]/(2.*sqrt(D2v2[i]));				// sqrt( D/(V{2}-D) )  stat. err.
		    	gsqrtD2v2->SetPoint(i, centlocv2[i], sqrtD2v2[i]);
		    	gsqrtD2v2->SetPointError(i, 0, esqrtD2v2[i]);
		    }
	    }
    }

    // D{2}/v^2  has two kind of sys. err.: from D and V{2} 
    // from V{2} sys. err.
    double D2v2syserrV22l[Nocent] = {0}, D2v2syserrV22h[Nocent] = {0};						// sqrt(delta/<v^2>) sys. err. from V2{2}
    for(int i = 0; i<V_22syserr->GetN();i++) {
	    D2v2syserrV22l[i] = V_22syserrl[i]*(-Delta[i]/pow(V_22[i]-Delta[i],2));			// sys. error for D/(V{2}-D) from V{2}
	    D2v2syserrV22h[i] = V_22syserrh[i]*(-Delta[i]/pow(V_22[i]-Delta[i],2));		
    }
    // from D sys. err			2013.11.30	ly
    double D2v2syserrDl[Nocent] = {0}, D2v2syserrDh[Nocent] = {0};						// sqrt(delta/<v^2>) sys. err. from D
    for(int i = 0; i<V_22syserr->GetN();i++) {
	    D2v2syserrDl[i] = D2syserrl[i]*(V_22[i]/pow(V_22[i]-Delta[i],2));			// error for D/(V{2}-D) from D
	    D2v2syserrDh[i] = D2syserrh[i]*(V_22[i]/pow(V_22[i]-Delta[i],2));		
    }
    ClassSysErr *D2v2syserr = (ClassSysErr*) SumTwoSSysErr(gD2v2,D2v2syserrV22l,D2v2syserrV22h,D2v2syserrDl,D2v2syserrDh);

    // sqrt( D/(V{2}-D) )  sys. err. from D/(V{2}-D)
    double D2v2syserrl[Nocent] = {0}, D2v2syserrh[Nocent] = {0}, *tmpD2v2syserrl, *tmpD2v2syserrh;
    double sqrtD2v2syserrl[Nocent] = {0}, sqrtD2v2syserrh[Nocent] = {0};
    tmpD2v2syserrl = D2v2syserr->GetEYlow();
    tmpD2v2syserrh = D2v2syserr->GetEYhigh();
    for(int i = 0; i<Nocent; i++) {
        D2v2syserrl[i] = tmpD2v2syserrl[i];
        D2v2syserrh[i] = tmpD2v2syserrh[i];
	    if(D2v2[i]>0) {
	    	sqrtD2v2syserrl[i] = D2v2syserrl[i]/(2.*sqrt(D2v2[i]));					// ClassSysErr low err is negative number
	    	sqrtD2v2syserrh[i] = D2v2syserrh[i]/(2.*sqrt(D2v2[i]));
	    }
    }
    ClassSysErr *sqrtD2v2syserr = new ClassSysErr(Nocent, centlocv2, sqrtD2v2, sqrtD2v2syserrl, sqrtD2v2syserrh);



    char outname[500] = "";
    sprintf(outname,"v%dDmsigmaprimevscent.%s",vn,Form(v2filetag[0],"all"));
    
    // setup enviroment for drawing

    TCanvas *c = new TCanvas("c","",CANVASWIDTH,CANVASHEIGTH);

    TPad *pad = new TPad("pad","pad",PADXL,PADYL,PADXH,PADYH);

    double yaxismin = -0; //-0.2;
    double yaxismax = 5;
//   if(vn==3) { yaxismin = -0.8; yaxismax = 3.1; }                       // 2014.05.02 ly
    if(vn==3) { yaxismin = 0; yaxismax = 3.7; }                       // 2014.05.02 ly
    TH2D *htmp = new TH2D("htmp","",8,0,80,100,yaxismin,yaxismax);    // x-range, y-range go here
    htmp->GetXaxis()->SetTitle("Centrality %");
//    htmp->GetYaxis()->SetTitle("#sqrt{#LT #font[12]{v}^{2} #GT} and #sqrt{#bar{D}} (%)");
    htmp->GetYaxis()->SetTitle("#sqrt{#LT v^{2} #GT} and #sqrt{#bar{D}} (%)");

    TLine *line = new TLine(-0.1,0,80,0);      // keep same xmin and xmax as htmp
    
    FigFormat1D(c,pad,htmp);

    c->cd();
    pad->Draw();
    pad->cd();
    htmp->Draw();
    line->Draw("same");

    //wLatex(0.7, 0.835, "|#Delta#eta|>0.7",1,12,0.06);
    //wLatex(0.52, 0.355, pttext,1,12,0.06);
    //wLatex(0.456, 0.28, "Au+Au 200GeV",1,12,0.08);
    if(vn==3) {
        wLatex(0.63, 0.85, "Au+Au 200GeV",1,12,0.05);
        wLatex(0.59, 0.78, pttext,1,12,0.05);
        wLatex(0.75, 0.73, "|#Delta#eta|>0.7",1,12,0.05);
    }
    else {
        wLatex(0.75, 0.37, "|#Delta#eta|>0.7",1,12,0.05);
        wLatex(0.59, 0.305, pttext,1,12,0.05);
        wLatex(0.63, 0.26, "Au+Au 200GeV",1,12,0.05);
    }

    wLatex(0.05,0.91,symbol,1,12,0.08);              // figure's subcaption

    int v2color=1, v3color=1;//634;//   2;                  // 634: kRed+2
    int v2band=921, v3band = 921; // 801; // 800;       // 920: kGray, 800: kOrange
    int d2band=1, d3band = 1;// 613; //400;       // 921: kGray+1, 400: kYellow, 613: kMagenta-1
    int v2marker=4, d2marker=8;
    int v3marker=26, d3marker=22;
    int vcolor = v2color, vmarker = v2marker, dmarker = d2marker, vband = v2band, dband = d2band;
    if(vn==3) {
        vcolor=v3color;
	    vmarker=v3marker;
	    dmarker=d3marker;
        vband = v3band;
        dband = d3band;
    }
    int Dmsigmaprimecolor=1;//4;            // 2014.05.06 ly
    int Dmsigmaprimebandcolor=920;//4;            // 2014.05.06 ly
    int Dmsigmaprimemarker=30;             // sqrt(D-sigma)
    float Dmsigmaprimemarkersize=2.5;             // sqrt(D-sigma)

    int dxwidth=150;
    v2syserr->ScaleMeanYOnly(v2scale);
    graphSystBox(v2syserr->GetN(),v2syserr->GetX(),v2syserr->GetY(),0,0,v2syserr->GetEYhigh(),v2syserr->GetEYlow(),vband,dxwidth,0);
// 2014.05.06 ly   graphSystBox(sqrtD2syserr->GetN(),sqrtD2syserr->GetX(),sqrtD2syserr->GetY(),0,0,sqrtD2syserr->GetEYhigh(),sqrtD2syserr->GetEYlow(),dband,dxwidth,0);
    if(vn==2) {
// 2015.05.06 ly        graphSystBox(sqrtDmsigmaprimesyserr->GetN(),sqrtDmsigmaprimesyserr->GetX(),sqrtDmsigmaprimesyserr->GetY(),0,0,sqrtDmsigmaprimesyserr->GetEYhigh(),sqrtDmsigmaprimesyserr->GetEYlow(),-1,dxwidth,0,1,Dmsigmaprimecolor);
        graphSystBox(sqrtDmsigmaprimesyserr->GetN(),sqrtDmsigmaprimesyserr->GetX(),sqrtDmsigmaprimesyserr->GetY(),0,0,sqrtDmsigmaprimesyserr->GetEYhigh(),sqrtDmsigmaprimesyserr->GetEYlow(),Dmsigmaprimebandcolor,dxwidth,0);       // 2014.05.06 ly
    }
    graphSystBox(sqrtD2syserr->GetN(),sqrtD2syserr->GetX(),sqrtD2syserr->GetY(),0,0,sqrtD2syserr->GetEYhigh(),sqrtD2syserr->GetEYlow(),-1,dxwidth,0,1,dband);           // 2014.05.06 ly


    // need to scale v2 on the plot
    TGraphErrors* gv22cp = new TGraphErrors();
    for(i = 0; i<gv22[0]->GetN() ; i++) {
        gv22cp->SetPoint(i,centlocv2[i],sqrtvv2[i]*v2scale);
        gv22cp->SetPointError(i,0,esqrtvv2[i]);
    }
    gv22cp->SetMarkerStyle(vmarker);
    gv22cp->SetMarkerSize(2);
    gv22cp->SetMarkerColor(vcolor);
    gv22cp->SetLineColor(vcolor);
    gsqrtDmsigmaprime[0]->SetMarkerStyle(Dmsigmaprimemarker);
    gsqrtDmsigmaprime[0]->SetMarkerSize(Dmsigmaprimemarkersize);
    gsqrtDmsigmaprime[0]->SetMarkerColor(Dmsigmaprimecolor);
    gsqrtDmsigmaprime[0]->SetLineColor(Dmsigmaprimecolor);
    gd22[0]->SetMarkerStyle(dmarker);
    gd22[0]->SetMarkerSize(2);
    gd22[0]->SetMarkerColor(vcolor);
    gd22[0]->SetLineColor(vcolor);


    if(vn==2) {
        gsqrtDmsigmaprime[0]->Draw("epsame");
    }
    gd22[0]->Draw("pXsame");
    gv22cp->Draw("epsame");

    //double xtxt = 0.31, xstep = 0.22;
    //double ytxt = 0.82, ystep = 0.095;
    double xtxt = 0.29, xstep = 0.22;
    double ytxt = 0.82, ystep = 0.095;
    double tsize = 0.05;
    double markersize = 1.5;
    if(vn==2) {
//        keySymbol(xtxt,ytxt,Form("#sqrt{#LT #font[12]{v}_{%d}^{2} #GT}/2",vn),vcolor,vmarker,tsize,markersize);
        keySymbol(xtxt,ytxt,Form("#sqrt{#LT v_{%d}^{2} #GT}/2",vn),vcolor,vmarker,tsize,markersize);
        keySymbol(xtxt+xstep*1,ytxt-ystep*0,"#sqrt{ #delta_{2}}",Dmsigmaprimecolor,Dmsigmaprimemarker,tsize,Dmsigmaprimemarkersize);
    	keySymbol(xtxt,ytxt-ystep,Form("#sqrt{ #bar{D}_{%d}}",vn),vcolor,dmarker,tsize,markersize);
    }
    else {
//        keySymbol(xtxt,ytxt,Form("#sqrt{#LT #font[12]{v}_{%d}^{2} #GT}",vn),vcolor,vmarker,tsize,markersize);
        keySymbol(xtxt,ytxt,Form("#sqrt{#LT v_{%d}^{2} #GT}",vn),vcolor,vmarker,tsize,markersize);
    	keySymbol(xtxt+xstep*0.8,ytxt,Form("#sqrt{ #bar{D}_{%d}}",vn),vcolor,dmarker,tsize,markersize);
    }


    if(savefig) {
        c->SaveAs(Form("paperfig_collreview_July4/%s.eps",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.gif",outname));
        c->SaveAs(Form("paperfig_collreview_July4/%s.pdf",outname));
    }




    cout<<"|#Delta#eta|>"<<etagap<<Form(" sigma': ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<sigmaprime[i]<<" + "<<sigmaprimesyserrh[i]<<" - "<<sigmaprimesyserrl[i]<<"(sys.) "<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" sqrt(v_{%d}^2): ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<sqrtvv2[i]<<" +- "<<esqrtvv2[i]<<"(stat.) + "<<v2syserrh[i]<<" - "<<v2syserrl[i]<<"(sys.) "<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" D_{%d}: ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<Delta[i]<<" + "<<D2syserrh[i]<<" - "<<D2syserrl[i]<<"(sys.) "<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" sqrt(D_{%d}): ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<sqrtDelta[i]<<" + "<<sqrtD2syserrh[i]<<" + "<<sqrtD2syserrl[i]<<"(sys.) "<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" D_{%d}-sigma': ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<Dmsigmaprime[i]<<" + "<<Dmsigmaprimesysh[i]<<" - "<<Dmsigmaprimesysl[i]<<"(sys.)"<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" sqrt(D_{%d}-sigma'): ",vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<sqrtDmsigmaprime[i]<<" + "<<sqrtDmsigmaprimesysh[i]<<" - "<<sqrtDmsigmaprimesysl[i]<<"(sys.)"<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" D_{%d}/<v_{%d}>: ",vn,vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<D2v2[i]*100<<"% +- "<<eD2v2[i]*100<<"% (stat.) + "<<D2v2syserrh[i]*100<<"% "<<D2v2syserrl[i]*100<<"%(sys.)"<<endl;
    }

    cout<<"|#Delta#eta|>"<<etagap<<Form(" sqrt(D_{%d}/<v_{%d}^2>): ",vn,vn)<<endl;
    for(int i = 0 ;i <Nocent ; i++ ) {
        cout<<"cent "<<centloc[i]<<" "<<sqrtD2v2[i]*100<<"% +- "<<esqrtD2v2[i]*100<<"% (stat.) + "<<sqrtD2v2syserrh[i]*100<<"% "<<sqrtD2v2syserrl[i]*100<<"%(sys.)"<<endl;
    }



}


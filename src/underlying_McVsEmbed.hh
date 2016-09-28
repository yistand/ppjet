#ifndef __JETPT_MCVSEMBED_HH
#define __JETPT_MCVSEMBED_HH


#include "AjParameters.hh"
#include "JetAnalyzer.hh"



#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <string>


#define MAXARRAYLENGTH  5000


class underlying_McVsEmbed
{
	private:
		// Parameters
		int mVerbose; 		// how debug information to print out (0--no, 10 alot);
		// Underlying event particle: charged(1) or neutral(0) or all(2, or whatever val not equal 1 or 0)
		int mUnderlyingParticleCharge;
		// Jet: charged(1) or neutral(0) or all(2, or whatever val not equal 1 or 0)
		int mJetCharge;
		// Underlying event phi size: default 60 degree
		float mTranPhiSize;



		int nEventPassed;
		int nEventOutlierMcpT;		// how many events were cut off due to high Mc leading jet pt
		int nEventOutlierRcpT;		// how many events were cut off due to high Rc leading jet pt

		// Match to trigger or not
		bool mNeedToMatchTrig;

		// Neutral/Total Pt of Jet fraction cut
		bool mNeutralJetFracCut;


		// OutlierMcpTCut flag
		bool DoOutlierMcpTCut;
		double mOutlierMcpTCut;
		bool DoOutlierRcpTCut;
		double mOutlierRcpTCut;

		// result purpose
		bool flagIsTrigger;
		bool flagGoodEtaMcJet;
		bool flagGoodEtaRcJet;

		bool flagMatch2Lead;
		bool flagMatch2Sub;
		bool flagMatch2McJet;		// any Mc jet, including 1st, 2nd ... hardest jet

		bool flagMatch2LeadGood;
		bool flagMatch2SubGood;	

		bool trigmatch;		



		// For output 
		TFile *fout;
		TString OutFileName;
		// Tree and its variables
		TTree* ResultTree;
			
		// Matched:
		TH1D* McMatchedLeadJetPt;
		TH1D* McMatchedLeadOrSubJetPt;
		TH2D* CovMcMatchedLeadJetVsRcJet;
		TH2D* CovMcMatchedLeadOrSubJetVsRcJet;

		// No matching requirement from Mc to Rc
		TH1D* McLeadJetPt;
		TH1D* McSubLeadJetPt;
		TH1D* RcLeadJetPt;
		TH1D* RcSubLeadJetPt;
		

		// event header
		int eventid;
		int runid;

		double Mcrefmult; 
		double Mcvz;

		double Rcrefmult; 
		double Rcvz;

		// event weight for that pt bin
		//double weight;

		// Jet-finding from fastjet
		float Mcj1pt, Mcjaspt, Mcj2pt;
		float Mcj1phi, Mcjasphi, Mcj2phi;
		float Mcj1eta, Mcjaseta, Mcj2eta;
		float Mcj1area,Mcjasarea, Mcj2area;
		float Mcj1area_err, Mcjasarea_err, Mcj2area_err;
		float Mcj1neutralfrac;

		float Mcrho, Mcrhoerr;

		float Mcj3pt, Mcj4pt;
		float Mcj3phi, Mcj4phi;
		float Mcj3eta, Mcj4eta;



		float Rcj1pt, Rcjaspt, Rcj2pt;
		float Rcj1phi, Rcjasphi, Rcj2phi;
		float Rcj1eta, Rcjaseta, Rcj2eta;
		float Rcj1area,Rcjasarea, Rcj2area;
		float Rcj1area_err, Rcjasarea_err, Rcj2area_err;
		float Rcj1neutralfrac;

		float Rcrho, Rcrhoerr;

		float Rcj3pt, Rcj4pt;
		float Rcj3phi, Rcj4phi;
		float Rcj3eta, Rcj4eta;

		// any Mc level jet matched with Rcj1: the n-th hardest matched, pt, phi, eta of the matched MC jet
		int MatchedNthMcj;
		float MatchedMcjpt, MatchedMcjphi, MatchedMcjeta; 	
		// if any Mc level jet matched with Rcj1: also record the leading Mc jet info in the event
		float LeadInMatchedMcjpt, LeadInMatchedMcjphi, LeadInMatchedMcjeta; 	
		
		// any Rc level jet matched with Mcj1 including even outside jet eta region but constituents cut still applied
		int MatchedNthRcj;
		float MatchedRcjpt, MatchedRcjphi, MatchedRcjeta; 	


		// underlying event info
		float McLeadAreaPt;
		float McSubAreaPt;
		float McTranMaxPt;
		float McTranMinPt;
		float McTranPt;
		
		int McLeadAreaNtrk;
		int McSubAreaNtrk;
		int McTranMaxNtrk;
		int McTranMinNtrk;
		
		float McTrkTranMaxPt[MAXARRAYLENGTH];		// for track in tranmax
		float McTrkTranMaxPhi[MAXARRAYLENGTH];		// for track in tranmax
		float McTrkTranMaxEta[MAXARRAYLENGTH];		// for track in tranmax
		float McTrkTranMinPt[MAXARRAYLENGTH];		// for track in tranmin
		float McTrkTranMinPhi[MAXARRAYLENGTH];		// for track in tranmax
		float McTrkTranMinEta[MAXARRAYLENGTH];		// for track in tranmax
		float McTrkLeadAreaPt[MAXARRAYLENGTH];			// for track in lead
		float McTrkLeadAreaPhi[MAXARRAYLENGTH];		// for track in lead
		float McTrkLeadAreaEta[MAXARRAYLENGTH];		// for track in lead
		float McTrkSubAreaPt[MAXARRAYLENGTH];		// for track in sublead
		float McTrkSubAreaPhi[MAXARRAYLENGTH];		// for track in sublead
		float McTrkSubAreaEta[MAXARRAYLENGTH];		// for track in sublead
 



		float RcLeadAreaPt;
		float RcSubAreaPt;
		float RcTranMaxPt;
		float RcTranMinPt;
		float RcTranPt;
		
		int RcLeadAreaNtrk;
		int RcSubAreaNtrk;
		int RcTranMaxNtrk;
		int RcTranMinNtrk;
		
		float RcTrkTranMaxPt[MAXARRAYLENGTH];		// for track in tranmax
		float RcTrkTranMaxPhi[MAXARRAYLENGTH];		// for track in tranmax
		float RcTrkTranMaxEta[MAXARRAYLENGTH];		// for track in tranmax
		float RcTrkTranMinPt[MAXARRAYLENGTH];		// for track in tranmin
		float RcTrkTranMinPhi[MAXARRAYLENGTH];		// for track in tranmax
		float RcTrkTranMinEta[MAXARRAYLENGTH];		// for track in tranmax
		float RcTrkLeadAreaPt[MAXARRAYLENGTH];			// for track in lead
		float RcTrkLeadAreaPhi[MAXARRAYLENGTH];		// for track in lead
		float RcTrkLeadAreaEta[MAXARRAYLENGTH];		// for track in lead
		float RcTrkSubAreaPt[MAXARRAYLENGTH];		// for track in sublead
		float RcTrkSubAreaPhi[MAXARRAYLENGTH];		// for track in sublead
		float RcTrkSubAreaEta[MAXARRAYLENGTH];		// for track in sublead
 








		// for jet finding
		// --------------------
		double R;              ///< Resolution parameter ("jet radius")
		double max_rap;        ///< jet rapidity acceptance
		double ghost_maxrap;   ///< for ghosted area, should be >= max_rap + 2*R

		double max_const_rap;  ///< constituent rapidity cut
		double min_const_pt;   ///< constituent pt cut

		double dPhiCut;        ///< opening angle for dijet requirement. Accept only  |&phi;1 - &phi;2 - &pi;| < &Delta;&phi;.

		fastjet::JetDefinition jet_def;       ///< jet definition
		fastjet::JetDefinition other_jet_def; ///< jet definition with a different radius

		fastjet::Selector select_const_rap;   ///< constituent rapidity selector
		fastjet::Selector select_const_ptmin;   ///< constituent p<SUB>T</SUB>  selector

		fastjet::Selector Mcsconst;                ///< compound selector for constituents for jet finding
		fastjet::Selector Rcsconst;                ///< compound selector for constituents for jet finding

		fastjet::Selector McsUconst;                ///< compound selector for constituents for underlying events 
		fastjet::Selector RcsUconst;                ///< compound selector for constituents for underlying events 


		// Relevant jet candidates
		fastjet::Selector select_jet_rap;        ///< jet rapidity selector
		fastjet::Selector select_jet_pt_min;     ///< jet p<SUB>T</SUB> selector
		fastjet::Selector select_jet_pt_max;     ///< jet p<SUB>T</SUB> selector
		fastjet::Selector sjet;                  ///< compound jet selector
		fastjet::Selector NoGhosts;              ///< compound jet selector

		fastjet::GhostedAreaSpec area_spec;      ///< ghosted area specification
		fastjet::AreaDefinition area_def;        ///< jet area definition

		JetAnalyzer* McnosjetJA;                      ///< JetAnalyzer objectMC particle level jet
		JetAnalyzer* McJA_bkgsub;                      ///< JetAnalyzer object with background subtraction

		JetAnalyzer* RcnosjetJA;                      ///< JetAnalyzer objectRc (reconstructed) detector level jet
		JetAnalyzer* RcJA_bkgsub;                      ///< JetAnalyzer object with background subtraction

		std::vector<fastjet::PseudoJet> McJconstituents;     ///< constituents for jet finding
		std::vector<fastjet::PseudoJet> RcJconstituents;     ///< constituents for jet finding
		std::vector<fastjet::PseudoJet> McUconstituents;     ///< constituents for underlying events
		std::vector<fastjet::PseudoJet> RcUconstituents;     ///< constituents for underlying events 

		std::vector<fastjet::PseudoJet> McnosjetJAResult;  ///< Unaltered clustering result 	no jet eta selection
		std::vector<fastjet::PseudoJet> RcnosjetJAResult;  ///< Unaltered clustering result  	no jet eta selection
		std::vector<fastjet::PseudoJet> McJAResult;  ///< Unaltered clustering result		with jet eta selection 
		std::vector<fastjet::PseudoJet> RcJAResult;  ///< Unaltered clustering result 		with jet eta selection

		std::vector<fastjet::PseudoJet> McJAResult_bkgsub;  ///< Unaltered clustering result with background subtraction
		std::vector<fastjet::PseudoJet> RcJAResult_bkgsub;  ///< Unaltered clustering result with background subtraction



		std::vector<fastjet::PseudoJet> McnosjetDiJets;    ///< Dijet result 
		std::vector<fastjet::PseudoJet> RcDiJets;    ///< Dijet result 


	public:
		underlying_McVsEmbed (	double R = 0.6,
				double max_const_rap = 1.0,
				double min_const_pt = 0.2,  
				std::string jetalgo = "antikt",
				TString name = "JetMcEmbedoutput.root" 		
				);

		~underlying_McVsEmbed();

		int LoopUnderlying (float RefPhi, std::vector<fastjet::PseudoJet> Uconstituents, int &LeadAreaNtrk, int &SubAreaNtrk, int &TranMaxNtrk, int &TranMinNtrk, float &LeadAreaPt, float &SubAreaPt, float &TranMaxPt, float &TranMinPt, float *TrkLeadAreaPt, float *TrkLeadAreaPhi, float *TrkLeadAreaEta, float *TrkSubAreaPt, float *TrkSubAreaPhi, float *TrkSubAreaEta, float *TrkTranMaxPt, float *TrkTranMaxPhi, float *TrkTranMaxEta, float *TrkTranMinPt, float *TrkTranMinPhi, float *TrkTranMinEta);

		int Init();

		int Make(const std::vector<fastjet::PseudoJet>& Mcparticles,
			 const std::vector<fastjet::PseudoJet>& Rcparticles,
			 int ineventid, int inrunid,
			 double inMcrefmult, double inMcvz,
			 double inRcrefmult, double inRcvz,
			 const std::vector<std::pair<float,float> > &ToMatch,
			 Bool_t is_trigger
			 //double weightbyXsec
			);

		int Finish();


		// Getters and Setters
		// -------------------
		// Whether need to match jet found by fastjet with the location fired the trigger
		void SetToMatchJetTrigger(bool val) {mNeedToMatchTrig = val; }
		bool GetToMatchJetTrigger() {return mNeedToMatchTrig; }

		// Wether apply Neutral/Total Pt of Jet fraction cut
		void SetNetraulJetFracCut(bool val) {mNeutralJetFracCut = val; };

		// Underlying event particle: charged(1) or neutral(0) or all(2)
		void SetUnderlyingParticleCharge(int val);
		int GetUnderlyingParticleCharge() { return mUnderlyingParticleCharge; }
		
		// Jet: charged(1) or neutral(0) or all(2)
		void SetJetCharge(int val);
		int GetJetCharge() { return mJetCharge; }

		// Underlying phi size
		void SetTransversePhiSize(float val) { if(val<180&&val>0) {mTranPhiSize = val;} else {mTranPhiSize = 60;} }
		int GetTransversePhiSize() { return mTranPhiSize; }


		int SetVerbose(int val) {mVerbose=val;}

		/// Get jet radius
		inline double GetR ( )                   { return R; }
		/// Set jet radius
		inline void   SetR ( const double newv ) { R=newv;   }


		/// Get jet rapidity acceptance
		inline double GetMax_rap ( )                   { return max_rap; }
		/// Set jet rapidity acceptance
		inline void   SetMax_rap ( const double newv ) { max_rap=newv;   }

		/// Get ghosted area rapidity cut, should be >= max_rap + 2*R
		inline double GetGhost_maxrap ( )                   { return ghost_maxrap; }
		/// Set ghosted area rapidity cut, should be >= max_rap + 2*R
		inline void   SetGhost_maxrap ( const double newv ) { ghost_maxrap=newv;   }

		/// Get dijet opening angle
		inline double GetDPhiCut ( )                   { return dPhiCut; }
		/// Set dijet opening angle
		inline void   SetDPhiCut ( const double newv ) { dPhiCut=newv;   }


		// Objects will be handed by _reference_! Obviates need for setter
		/// Handle to jet definition
		inline fastjet::JetDefinition& GetJet_def () { return jet_def; }
		/// Handle to selector for constituents
		inline fastjet::Selector& GetMcConsSelector () { return Mcsconst; }
		inline fastjet::Selector& GetRcConsSelector () { return Rcsconst; }

		/// Handle to selector for jet candidates
		inline fastjet::Selector& GetJetSelector () { return sjet; }

		/// Handle to ghosted area specification
		inline fastjet::GhostedAreaSpec& GetArea_spec () { return area_spec; }
		/// Handle to jet area definition
		inline fastjet::AreaDefinition& GetArea_def () { return area_def; }


		/// Handle to JetAnalyzer
		inline JetAnalyzer* GetMcJA() {return McnosjetJA; }
		inline JetAnalyzer* GetRcJA() {return RcnosjetJA; }

		/// Handle to unaltered clustering result with high pT constituents
		inline std::vector<fastjet::PseudoJet> GetMcJAResult() {return McnosjetJAResult; }
		inline std::vector<fastjet::PseudoJet> GetRcJAResult() {return RcJAResult; }

		/// Handle to constituents
		inline std::vector<fastjet::PseudoJet> GetMcConstituents() {return McJconstituents; }
		inline std::vector<fastjet::PseudoJet> GetRcConstituents() {return RcJconstituents; }


		// Delete outlier
		//		2016.07.25	Li Yi
		//		pt2_3 has one entry with extremely high Mc pt. This may cause trouble when merging
		//		because low pT bin will have huge weight, therefore we delete this bin
		void SetOutlierMcpTCut(double val) {DoOutlierMcpTCut = true; mOutlierMcpTCut = val;}
		bool IsOutlierMcpTCutApplied() {return DoOutlierMcpTCut;}
		double GetOutlierMcpTCut() {return mOutlierMcpTCut; }
		int GetNEventOutlierMcpTCut() {return nEventOutlierMcpT;}
		void SetOutlierRcpTCut(double val) {DoOutlierRcpTCut = true; mOutlierRcpTCut = val;}
		bool IsOutlierRcpTCutApplied() {return DoOutlierRcpTCut;}
		double GetOutlierRcpTCut() {return mOutlierRcpTCut; }
		int GetNEventOutlierRcpTCut() {return nEventOutlierRcpT;}
};


#endif // __JETPT_MCVSEMBED_HH

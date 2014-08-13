#include<iostream>
#include <TH1.h>
#include "TStyle.h"
#include <TCanvas.h>
#include<Riostream.h>
#include "TROOT.h"
#include<TGraph.h>
#include "TGraphErrors.h"
using namespace std;

void RMSvsnPV(){

gStyle->SetOptStat(0000);

int maxEvents=40000;

double pTmin=300;
double pTmax=4000;

char inputfile[100];

sprintf(inputfile,"QCD300to470csa14Puppiv3.root");

const unsigned int nPVmin=15;
const unsigned int nPVmax=50;


double x_err=2.5;
double binsize=5;
const unsigned double setMinMeanVsnpuTr=-30;
const unsigned double setMaxMeanVsnpuTr=60;

const unsigned double setMinSigmaVsnpuTr=0;
const unsigned double setMaxSigmaVsnpuTr=60;

const unsigned double setMinMeanVsnpuPr=-30;
const unsigned double setMaxMeanVsnpuPr=60;

const unsigned double setMinSigmaVsnpuPr=0;
const unsigned double setMaxSigmaVsnpuPr=60;

const unsigned double setMinMeanVsnpuSd=-30;
const unsigned double setMaxMeanVsnpuSd=60;

const unsigned double setMinSigmaVsnpuSd=0;
const unsigned double setMaxSigmaVsnpuSd=60;


const unsigned double labelsize=0.035;
const unsigned double titlesize=0.045;
const unsigned double linewidth=2;
const unsigned double markerstyle=21;
const unsigned double markersize=1.5;
const unsigned int nPVnumOfbins=7;//should be a divisible number of (nPVmax-nPVmin)

const unsigned int NumtrPar=4;
const unsigned int NumprPar=4;
const unsigned int NumsdPar=4;

  double rangeMin=-100;
  double rangeMax=100;

   int NoOfbins=100;


//////////////////////////////////response
TH1F *JetRsAKpf[nPVnumOfbins];
TH1F *JetRsAKchs[nPVnumOfbins];
TH1F *JetRsAKpuppi[nPVnumOfbins];


TH1F *JetRsAKpfTr_nPV_h[NumtrPar][nPVnumOfbins];
TH1F *JetRsAKchsTr_nPV_h[NumtrPar][nPVnumOfbins];
TH1F *JetRsAKpuppiTr_nPV_h[NumtrPar][nPVnumOfbins];


TH1F *JetRsAKpfPr_nPV_h[NumprPar][nPVnumOfbins];
TH1F *JetRsAKchsPr_nPV_h[NumprPar][nPVnumOfbins];
TH1F *JetRsAKpuppiPr_nPV_h[NumprPar][nPVnumOfbins];


TH1F *JetRsAKpfSd_nPV_h[NumsdPar][nPVnumOfbins];
TH1F *JetRsAKchsSd_nPV_h[NumsdPar][nPVnumOfbins];
TH1F *JetRsAKpuppiSd_nPV_h[NumsdPar][nPVnumOfbins];
/////////////////////////////////////////response

char histnamepf1[100];
char histnamepf2[100];
char histnamepf3[100];

char histnamechs1[100];
char histnamechs2[100];
char histnamechs3[100];

char histnamepuppi1[100];
char histnamepuppi2[100];
char histnamepuppi3[100];

char histnamepf[100];
char histnamechs[100];
char histnamepuppi[100];








for (int ii=0;ii<nPVnumOfbins;ii++) {//loop over nPV bins

sprintf(histnamepf,"ResponsePFnPVbins_%i",ii);
sprintf(histnamechs,"ResponseCHSnPVbins_%i",ii);
sprintf(histnamepuppi,"ResponsePUPPInPVbins_%i",ii);
JetRsAKpf[ii]=new TH1F(histnamepf,histnamepf,NoOfbins,rangeMin,rangeMax);
JetRsAKchs[ii]=new TH1F(histnamechs,histnamechs,NoOfbins,rangeMin,rangeMax);
JetRsAKpuppi[ii]=new TH1F(histnamepuppi,histnamepuppi,NoOfbins,rangeMin,rangeMax);


for(int jj=0;jj<NumtrPar;jj++){//loop over trim pars

sprintf(histnamepf1,"ResponseTrPFnPVbins_%i_%i",ii,jj);
sprintf(histnamechs1,"ResponseTrCHSnPVbins_%i_%i",ii,jj);
sprintf(histnamepuppi1,"ResponseTrPUPPInPVbins_%i_%i",ii,jj);

JetRsAKpfTr_nPV_h[jj][ii]= new TH1F(histnamepf1,histnamepf1,NoOfbins,rangeMin,rangeMax);
JetRsAKchsTr_nPV_h[jj][ii]= new TH1F(histnamechs1,histnamechs1,NoOfbins,rangeMin,rangeMax);
JetRsAKpuppiTr_nPV_h[jj][ii]= new TH1F(histnamepuppi1,histnamepuppi1,NoOfbins,rangeMin,rangeMax);

}//loop over Trim pars


for(int jj=0;jj<NumprPar;jj++){//loop over prun pars

sprintf(histnamepf2,"ResponsePrPFnPVbins_%i_%i",ii,jj);
sprintf(histnamechs2,"ResponsePrCHSnPVbins_%i_%i",ii,jj);
sprintf(histnamepuppi2,"ResponsePrPUPPInPVbins_%i_%i",ii,jj);

JetRsAKpfPr_nPV_h[jj][ii]= new TH1F(histnamepf2,histnamepf2,NoOfbins,rangeMin,rangeMax);
JetRsAKchsPr_nPV_h[jj][ii]= new TH1F(histnamechs2,histnamechs2,NoOfbins,rangeMin,rangeMax);
JetRsAKpuppiPr_nPV_h[jj][ii]= new TH1F(histnamepuppi2,histnamepuppi2,NoOfbins,rangeMin,rangeMax);

}//loop over prun pars

for(int jj=0;jj<NumsdPar;jj++){//loop over sd pars

sprintf(histnamepf3,"ResponseSdPFnPVbins_%i_%i",ii,jj);
sprintf(histnamechs3,"ResponseSdCHSnPVbins_%i_%i",ii,jj);
sprintf(histnamepuppi3,"ResponseSdPUPPInPVbins_%i_%i",ii,jj);

JetRsAKpfSd_nPV_h[jj][ii]= new TH1F(histnamepf3,histnamepf3,NoOfbins,rangeMin,rangeMax);
JetRsAKchsSd_nPV_h[jj][ii]= new TH1F(histnamechs3,histnamechs3,NoOfbins,rangeMin,rangeMax);
JetRsAKpuppiSd_nPV_h[jj][ii]= new TH1F(histnamepuppi3,histnamepuppi3,NoOfbins,rangeMin,rangeMax);


}//loop over sd pars

}//loop over nPV bins


    TFile *f1 = new TFile(inputfile,"READ","My root file1");
    TTree *trpf = (TTree*)f1->Get("pf");
    TTree *trchs = (TTree*)f1->Get("chs");
    TTree *trpuppi = (TTree*)f1->Get("puppi");
    TTree *trsk = (TTree*)f1->Get("softkiller");
    TTree *trgen = (TTree*)f1->Get("gen");

       int   npv;
  //////////////////////////////////////PF Tree
   trpf->SetBranchAddress("npv", &npv);
   vector<float>   *pf_pttrim_Rtrim_020_Ptfrac_005;
   vector<float>   *pf_mtrim_Rtrim_020_Ptfrac_005;
   vector<float>   *pf_pttrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *pf_mtrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *pf_pttrim_Rtrim_010_Ptfrac_003;
   vector<float>   *pf_mtrim_Rtrim_010_Ptfrac_003;
   vector<float>   *pf_pttrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *pf_mtrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *pf_pttrim_Rtrim_020_Ptfrac_003;
   vector<float>   *pf_mtrim_Rtrim_020_Ptfrac_003;
   vector<float>   *pf_pttrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *pf_mtrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *pf_pttrim_Rtrim_030_Ptfrac_003;
   vector<float>   *pf_mtrim_Rtrim_030_Ptfrac_003;
   vector<float>   *pf_pttrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *pf_mtrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *pf_ptconst;
   vector<float>   *pf_mconst;
   vector<float>   *pf_ptpruned_zcut_010_R_cut_050;
   vector<float>   *pf_mpruned_zcut_010_R_cut_050;
   vector<float>   *pf_ptprunedsafe_zcut_010_R_cut_050;
   vector<float>   *pf_mprunedsafe_zcut_010_R_cut_050;
   vector<float>   *pf_ptpruned_zcut_005_R_cut_050;
   vector<float>   *pf_mpruned_zcut_005_R_cut_050;
   vector<float>   *pf_ptprunedsafe_zcut_005_R_cut_050;
   vector<float>   *pf_mprunedsafe_zcut_005_R_cut_050;
   vector<float>   *pf_ptpruned_zcut_005_R_cut_075;
   vector<float>   *pf_mpruned_zcut_005_R_cut_075;
   vector<float>   *pf_ptprunedsafe_zcut_005_R_cut_075;
   vector<float>   *pf_mprunedsafe_zcut_005_R_cut_075;
   vector<float>   *pf_ptpruned_zcut_010_R_cut_075;
   vector<float>   *pf_mpruned_zcut_010_R_cut_075;
   vector<float>   *pf_ptprunedsafe_zcut_010_R_cut_075;
   vector<float>   *pf_mprunedsafe_zcut_010_R_cut_075;
   vector<float>   *pf_ptsoftdrop_beta20;
   vector<float>   *pf_msoftdrop_beta20;
   vector<float>   *pf_ptsoftdropsafe_beta20;
   vector<float>   *pf_msoftdropsafe_beta20;
   vector<float>   *pf_ptsoftdrop_beta00;
   vector<float>   *pf_msoftdrop_beta00;
   vector<float>   *pf_ptsoftdropsafe_beta00;
   vector<float>   *pf_msoftdropsafe_beta00;
   vector<float>   *pf_ptsoftdrop_beta10;
   vector<float>   *pf_msoftdrop_beta10;
   vector<float>   *pf_ptsoftdropsafe_beta10;
   vector<float>   *pf_msoftdropsafe_beta10;
   vector<float>   *pf_ptsoftdrop_betam1;
   vector<float>   *pf_msoftdrop_betam1;
   vector<float>   *pf_ptsoftdropsafe_betam1;
   vector<float>   *pf_msoftdropsafe_betam1;

trpf->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005",&pf_pttrim_Rtrim_020_Ptfrac_005);
trpf->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005",&pf_mtrim_Rtrim_020_Ptfrac_005);
trpf->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005",&pf_pttrimsafe_Rtrim_020_Ptfrac_005);
trpf->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005",&pf_mtrimsafe_Rtrim_020_Ptfrac_005);
trpf->SetBranchAddress("pttrim_Rtrim_010_Ptfrac_003",&pf_pttrim_Rtrim_010_Ptfrac_003);
trpf->SetBranchAddress("mtrim_Rtrim_010_Ptfrac_003",&pf_mtrim_Rtrim_010_Ptfrac_003);
trpf->SetBranchAddress("pttrimsafe_Rtrim_010_Ptfrac_003",&pf_pttrimsafe_Rtrim_010_Ptfrac_003);
trpf->SetBranchAddress("mtrimsafe_Rtrim_010_Ptfrac_003",&pf_mtrimsafe_Rtrim_010_Ptfrac_003);
trpf->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_003",&pf_pttrim_Rtrim_020_Ptfrac_003);
trpf->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_003",&pf_mtrim_Rtrim_020_Ptfrac_003);
trpf->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_003",&pf_pttrimsafe_Rtrim_020_Ptfrac_003);
trpf->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_003",&pf_mtrimsafe_Rtrim_020_Ptfrac_003);
trpf->SetBranchAddress("pttrim_Rtrim_030_Ptfrac_003",&pf_pttrim_Rtrim_030_Ptfrac_003);
trpf->SetBranchAddress("mtrim_Rtrim_030_Ptfrac_003",&pf_mtrim_Rtrim_030_Ptfrac_003);
trpf->SetBranchAddress("pttrimsafe_Rtrim_030_Ptfrac_003",&pf_pttrimsafe_Rtrim_030_Ptfrac_003);
trpf->SetBranchAddress("mtrimsafe_Rtrim_030_Ptfrac_003",&pf_mtrimsafe_Rtrim_030_Ptfrac_003);
trpf->SetBranchAddress("ptconst",&pf_ptconst);
trpf->SetBranchAddress("mconst",&pf_mconst);
trpf->SetBranchAddress("ptpruned_zcut_010_R_cut_050",&pf_ptpruned_zcut_010_R_cut_050);
trpf->SetBranchAddress("mpruned_zcut_010_R_cut_050",&pf_mpruned_zcut_010_R_cut_050);
trpf->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_050",&pf_ptprunedsafe_zcut_010_R_cut_050);
trpf->SetBranchAddress("mprunedsafe_zcut_010_R_cut_050",&pf_mprunedsafe_zcut_010_R_cut_050);
trpf->SetBranchAddress("ptpruned_zcut_005_R_cut_050",&pf_ptpruned_zcut_005_R_cut_050);
trpf->SetBranchAddress("mpruned_zcut_005_R_cut_050",&pf_mpruned_zcut_005_R_cut_050);
trpf->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_050",&pf_ptprunedsafe_zcut_005_R_cut_050);
trpf->SetBranchAddress("mprunedsafe_zcut_005_R_cut_050",&pf_mprunedsafe_zcut_005_R_cut_050);
trpf->SetBranchAddress("ptpruned_zcut_005_R_cut_075",&pf_ptpruned_zcut_005_R_cut_075);
trpf->SetBranchAddress("mpruned_zcut_005_R_cut_075",&pf_mpruned_zcut_005_R_cut_075);
trpf->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_075",&pf_ptprunedsafe_zcut_005_R_cut_075);
trpf->SetBranchAddress("mprunedsafe_zcut_005_R_cut_075",&pf_mprunedsafe_zcut_005_R_cut_075);
trpf->SetBranchAddress("ptpruned_zcut_010_R_cut_075",&pf_ptpruned_zcut_010_R_cut_075);
trpf->SetBranchAddress("mpruned_zcut_010_R_cut_075",&pf_mpruned_zcut_010_R_cut_075);
trpf->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_075",&pf_ptprunedsafe_zcut_010_R_cut_075);
trpf->SetBranchAddress("mprunedsafe_zcut_010_R_cut_075",&pf_mprunedsafe_zcut_010_R_cut_075);
trpf->SetBranchAddress("ptsoftdrop_beta20",&pf_ptsoftdrop_beta20);
trpf->SetBranchAddress("msoftdrop_beta20",&pf_msoftdrop_beta20);
trpf->SetBranchAddress("ptsoftdropsafe_beta20",&pf_ptsoftdropsafe_beta20);
trpf->SetBranchAddress("msoftdropsafe_beta20",&pf_msoftdropsafe_beta20);
trpf->SetBranchAddress("ptsoftdrop_beta00",&pf_ptsoftdrop_beta00);
trpf->SetBranchAddress("msoftdrop_beta00",&pf_msoftdrop_beta00);
trpf->SetBranchAddress("ptsoftdropsafe_beta00",&pf_ptsoftdropsafe_beta00);
trpf->SetBranchAddress("msoftdropsafe_beta00",&pf_msoftdropsafe_beta00);
trpf->SetBranchAddress("ptsoftdrop_beta10",&pf_ptsoftdrop_beta10);
trpf->SetBranchAddress("msoftdrop_beta10",&pf_msoftdrop_beta10);
trpf->SetBranchAddress("ptsoftdropsafe_beta10",&pf_ptsoftdropsafe_beta10);
trpf->SetBranchAddress("msoftdropsafe_beta10",&pf_msoftdropsafe_beta10);
trpf->SetBranchAddress("ptsoftdrop_beta-10",&pf_ptsoftdrop_betam1);
trpf->SetBranchAddress("msoftdrop_beta-10",&pf_msoftdrop_betam1);
trpf->SetBranchAddress("ptsoftdropsafe_beta-10",&pf_ptsoftdropsafe_betam1);
trpf->SetBranchAddress("msoftdropsafe_beta-10",&pf_msoftdropsafe_betam1);

 /////////////////////////////////////////PF Tree



////////////////////////////////////chs

   vector<float>   *chs_pttrim_Rtrim_020_Ptfrac_005;
   vector<float>   *chs_mtrim_Rtrim_020_Ptfrac_005;
   vector<float>   *chs_pttrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *chs_mtrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *chs_pttrim_Rtrim_010_Ptfrac_003;
   vector<float>   *chs_mtrim_Rtrim_010_Ptfrac_003;
   vector<float>   *chs_pttrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *chs_mtrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *chs_pttrim_Rtrim_020_Ptfrac_003;
   vector<float>   *chs_mtrim_Rtrim_020_Ptfrac_003;
   vector<float>   *chs_pttrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *chs_mtrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *chs_pttrim_Rtrim_030_Ptfrac_003;
   vector<float>   *chs_mtrim_Rtrim_030_Ptfrac_003;
   vector<float>   *chs_pttrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *chs_mtrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *chs_ptconst;
   vector<float>   *chs_mconst;
   vector<float>   *chs_ptpruned_zcut_010_R_cut_050;
   vector<float>   *chs_mpruned_zcut_010_R_cut_050;
   vector<float>   *chs_ptprunedsafe_zcut_010_R_cut_050;
   vector<float>   *chs_mprunedsafe_zcut_010_R_cut_050;
   vector<float>   *chs_ptpruned_zcut_005_R_cut_050;
   vector<float>   *chs_mpruned_zcut_005_R_cut_050;
   vector<float>   *chs_ptprunedsafe_zcut_005_R_cut_050;
   vector<float>   *chs_mprunedsafe_zcut_005_R_cut_050;
   vector<float>   *chs_ptpruned_zcut_005_R_cut_075;
   vector<float>   *chs_mpruned_zcut_005_R_cut_075;
   vector<float>   *chs_ptprunedsafe_zcut_005_R_cut_075;
   vector<float>   *chs_mprunedsafe_zcut_005_R_cut_075;
   vector<float>   *chs_ptpruned_zcut_010_R_cut_075;
   vector<float>   *chs_mpruned_zcut_010_R_cut_075;
   vector<float>   *chs_ptprunedsafe_zcut_010_R_cut_075;
   vector<float>   *chs_mprunedsafe_zcut_010_R_cut_075;
   vector<float>   *chs_ptsoftdrop_beta20;
   vector<float>   *chs_msoftdrop_beta20;
   vector<float>   *chs_ptsoftdropsafe_beta20;
   vector<float>   *chs_msoftdropsafe_beta20;
   vector<float>   *chs_ptsoftdrop_beta00;
   vector<float>   *chs_msoftdrop_beta00;
   vector<float>   *chs_ptsoftdropsafe_beta00;
   vector<float>   *chs_msoftdropsafe_beta00;
   vector<float>   *chs_ptsoftdrop_beta10;
   vector<float>   *chs_msoftdrop_beta10;
   vector<float>   *chs_ptsoftdropsafe_beta10;
   vector<float>   *chs_msoftdropsafe_beta10;
   vector<float>   *chs_ptsoftdrop_betam1;
   vector<float>   *chs_msoftdrop_betam1;
   vector<float>   *chs_ptsoftdropsafe_betam1;
   vector<float>   *chs_msoftdropsafe_betam1;
trchs->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005",&chs_pttrim_Rtrim_020_Ptfrac_005);
trchs->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005",&chs_mtrim_Rtrim_020_Ptfrac_005);
trchs->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005",&chs_pttrimsafe_Rtrim_020_Ptfrac_005);
trchs->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005",&chs_mtrimsafe_Rtrim_020_Ptfrac_005);
trchs->SetBranchAddress("pttrim_Rtrim_010_Ptfrac_003",&chs_pttrim_Rtrim_010_Ptfrac_003);
trchs->SetBranchAddress("mtrim_Rtrim_010_Ptfrac_003",&chs_mtrim_Rtrim_010_Ptfrac_003);
trchs->SetBranchAddress("pttrimsafe_Rtrim_010_Ptfrac_003",&chs_pttrimsafe_Rtrim_010_Ptfrac_003);
trchs->SetBranchAddress("mtrimsafe_Rtrim_010_Ptfrac_003",&chs_mtrimsafe_Rtrim_010_Ptfrac_003);
trchs->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_003",&chs_pttrim_Rtrim_020_Ptfrac_003);
trchs->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_003",&chs_mtrim_Rtrim_020_Ptfrac_003);
trchs->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_003",&chs_pttrimsafe_Rtrim_020_Ptfrac_003);
trchs->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_003",&chs_mtrimsafe_Rtrim_020_Ptfrac_003);
trchs->SetBranchAddress("pttrim_Rtrim_030_Ptfrac_003",&chs_pttrim_Rtrim_030_Ptfrac_003);
trchs->SetBranchAddress("mtrim_Rtrim_030_Ptfrac_003",&chs_mtrim_Rtrim_030_Ptfrac_003);
trchs->SetBranchAddress("pttrimsafe_Rtrim_030_Ptfrac_003",&chs_pttrimsafe_Rtrim_030_Ptfrac_003);
trchs->SetBranchAddress("mtrimsafe_Rtrim_030_Ptfrac_003",&chs_mtrimsafe_Rtrim_030_Ptfrac_003);
trchs->SetBranchAddress("ptconst",&chs_ptconst);
trchs->SetBranchAddress("mconst",&chs_mconst);
trchs->SetBranchAddress("ptpruned_zcut_010_R_cut_050",&chs_ptpruned_zcut_010_R_cut_050);
trchs->SetBranchAddress("mpruned_zcut_010_R_cut_050",&chs_mpruned_zcut_010_R_cut_050);
trchs->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_050",&chs_ptprunedsafe_zcut_010_R_cut_050);
trchs->SetBranchAddress("mprunedsafe_zcut_010_R_cut_050",&chs_mprunedsafe_zcut_010_R_cut_050);
trchs->SetBranchAddress("ptpruned_zcut_005_R_cut_050",&chs_ptpruned_zcut_005_R_cut_050);
trchs->SetBranchAddress("mpruned_zcut_005_R_cut_050",&chs_mpruned_zcut_005_R_cut_050);
trchs->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_050",&chs_ptprunedsafe_zcut_005_R_cut_050);
trchs->SetBranchAddress("mprunedsafe_zcut_005_R_cut_050",&chs_mprunedsafe_zcut_005_R_cut_050);
trchs->SetBranchAddress("ptpruned_zcut_005_R_cut_075",&chs_ptpruned_zcut_005_R_cut_075);
trchs->SetBranchAddress("mpruned_zcut_005_R_cut_075",&chs_mpruned_zcut_005_R_cut_075);
trchs->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_075",&chs_ptprunedsafe_zcut_005_R_cut_075);
trchs->SetBranchAddress("mprunedsafe_zcut_005_R_cut_075",&chs_mprunedsafe_zcut_005_R_cut_075);
trchs->SetBranchAddress("ptpruned_zcut_010_R_cut_075",&chs_ptpruned_zcut_010_R_cut_075);
trchs->SetBranchAddress("mpruned_zcut_010_R_cut_075",&chs_mpruned_zcut_010_R_cut_075);
trchs->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_075",&chs_ptprunedsafe_zcut_010_R_cut_075);
trchs->SetBranchAddress("mprunedsafe_zcut_010_R_cut_075",&chs_mprunedsafe_zcut_010_R_cut_075);
trchs->SetBranchAddress("ptsoftdrop_beta20",&chs_ptsoftdrop_beta20);
trchs->SetBranchAddress("msoftdrop_beta20",&chs_msoftdrop_beta20);
trchs->SetBranchAddress("ptsoftdropsafe_beta20",&chs_ptsoftdropsafe_beta20);
trchs->SetBranchAddress("msoftdropsafe_beta20",&chs_msoftdropsafe_beta20);
trchs->SetBranchAddress("ptsoftdrop_beta00",&chs_ptsoftdrop_beta00);
trchs->SetBranchAddress("msoftdrop_beta00",&chs_msoftdrop_beta00);
trchs->SetBranchAddress("ptsoftdropsafe_beta00",&chs_ptsoftdropsafe_beta00);
trchs->SetBranchAddress("msoftdropsafe_beta00",&chs_msoftdropsafe_beta00);
trchs->SetBranchAddress("ptsoftdrop_beta10",&chs_ptsoftdrop_beta10);
trchs->SetBranchAddress("msoftdrop_beta10",&chs_msoftdrop_beta10);
trchs->SetBranchAddress("ptsoftdropsafe_beta10",&chs_ptsoftdropsafe_beta10);
trchs->SetBranchAddress("msoftdropsafe_beta10",&chs_msoftdropsafe_beta10);
trchs->SetBranchAddress("ptsoftdrop_beta-10",&chs_ptsoftdrop_betam1);
trchs->SetBranchAddress("msoftdrop_beta-10",&chs_msoftdrop_betam1);
trchs->SetBranchAddress("ptsoftdropsafe_beta-10",&chs_ptsoftdropsafe_betam1);
trchs->SetBranchAddress("msoftdropsafe_beta-10",&chs_msoftdropsafe_betam1);


////////////////////////////////////chs



//////////////////////////////puppi

   vector<float>   *puppi_pttrim_Rtrim_020_Ptfrac_005;
   vector<float>   *puppi_mtrim_Rtrim_020_Ptfrac_005;
   vector<float>   *puppi_pttrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *puppi_mtrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *puppi_pttrim_Rtrim_010_Ptfrac_003;
   vector<float>   *puppi_mtrim_Rtrim_010_Ptfrac_003;
   vector<float>   *puppi_pttrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *puppi_mtrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *puppi_pttrim_Rtrim_020_Ptfrac_003;
   vector<float>   *puppi_mtrim_Rtrim_020_Ptfrac_003;
   vector<float>   *puppi_pttrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *puppi_mtrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *puppi_pttrim_Rtrim_030_Ptfrac_003;
   vector<float>   *puppi_mtrim_Rtrim_030_Ptfrac_003;
   vector<float>   *puppi_pttrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *puppi_mtrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *puppi_ptconst;
   vector<float>   *puppi_mconst;
   vector<float>   *puppi_ptpruned_zcut_010_R_cut_050;
   vector<float>   *puppi_mpruned_zcut_010_R_cut_050;
   vector<float>   *puppi_ptprunedsafe_zcut_010_R_cut_050;
   vector<float>   *puppi_mprunedsafe_zcut_010_R_cut_050;
   vector<float>   *puppi_ptpruned_zcut_005_R_cut_050;
   vector<float>   *puppi_mpruned_zcut_005_R_cut_050;
   vector<float>   *puppi_ptprunedsafe_zcut_005_R_cut_050;
   vector<float>   *puppi_mprunedsafe_zcut_005_R_cut_050;
   vector<float>   *puppi_ptpruned_zcut_005_R_cut_075;
   vector<float>   *puppi_mpruned_zcut_005_R_cut_075;
   vector<float>   *puppi_ptprunedsafe_zcut_005_R_cut_075;
   vector<float>   *puppi_mprunedsafe_zcut_005_R_cut_075;
   vector<float>   *puppi_ptpruned_zcut_010_R_cut_075;
   vector<float>   *puppi_mpruned_zcut_010_R_cut_075;
   vector<float>   *puppi_ptprunedsafe_zcut_010_R_cut_075;
   vector<float>   *puppi_mprunedsafe_zcut_010_R_cut_075;
   vector<float>   *puppi_ptsoftdrop_beta20;
   vector<float>   *puppi_msoftdrop_beta20;
   vector<float>   *puppi_ptsoftdropsafe_beta20;
   vector<float>   *puppi_msoftdropsafe_beta20;
   vector<float>   *puppi_ptsoftdrop_beta00;
   vector<float>   *puppi_msoftdrop_beta00;
   vector<float>   *puppi_ptsoftdropsafe_beta00;
   vector<float>   *puppi_msoftdropsafe_beta00;
   vector<float>   *puppi_ptsoftdrop_beta10;
   vector<float>   *puppi_msoftdrop_beta10;
   vector<float>   *puppi_ptsoftdropsafe_beta10;
   vector<float>   *puppi_msoftdropsafe_beta10;
   vector<float>   *puppi_ptsoftdrop_betam1;
   vector<float>   *puppi_msoftdrop_betam1;
   vector<float>   *puppi_ptsoftdropsafe_betam1;
   vector<float>   *puppi_msoftdropsafe_betam1;
trpuppi->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005",&puppi_pttrim_Rtrim_020_Ptfrac_005);
trpuppi->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005",&puppi_mtrim_Rtrim_020_Ptfrac_005);
trpuppi->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005",&puppi_pttrimsafe_Rtrim_020_Ptfrac_005);
trpuppi->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005",&puppi_mtrimsafe_Rtrim_020_Ptfrac_005);
trpuppi->SetBranchAddress("pttrim_Rtrim_010_Ptfrac_003",&puppi_pttrim_Rtrim_010_Ptfrac_003);
trpuppi->SetBranchAddress("mtrim_Rtrim_010_Ptfrac_003",&puppi_mtrim_Rtrim_010_Ptfrac_003);
trpuppi->SetBranchAddress("pttrimsafe_Rtrim_010_Ptfrac_003",&puppi_pttrimsafe_Rtrim_010_Ptfrac_003);
trpuppi->SetBranchAddress("mtrimsafe_Rtrim_010_Ptfrac_003",&puppi_mtrimsafe_Rtrim_010_Ptfrac_003);
trpuppi->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_003",&puppi_pttrim_Rtrim_020_Ptfrac_003);
trpuppi->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_003",&puppi_mtrim_Rtrim_020_Ptfrac_003);
trpuppi->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_003",&puppi_pttrimsafe_Rtrim_020_Ptfrac_003);
trpuppi->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_003",&puppi_mtrimsafe_Rtrim_020_Ptfrac_003);
trpuppi->SetBranchAddress("pttrim_Rtrim_030_Ptfrac_003",&puppi_pttrim_Rtrim_030_Ptfrac_003);
trpuppi->SetBranchAddress("mtrim_Rtrim_030_Ptfrac_003",&puppi_mtrim_Rtrim_030_Ptfrac_003);
trpuppi->SetBranchAddress("pttrimsafe_Rtrim_030_Ptfrac_003",&puppi_pttrimsafe_Rtrim_030_Ptfrac_003);
trpuppi->SetBranchAddress("mtrimsafe_Rtrim_030_Ptfrac_003",&puppi_mtrimsafe_Rtrim_030_Ptfrac_003);
trpuppi->SetBranchAddress("ptconst",&puppi_ptconst);
trpuppi->SetBranchAddress("mconst",&puppi_mconst);
trpuppi->SetBranchAddress("ptpruned_zcut_010_R_cut_050",&puppi_ptpruned_zcut_010_R_cut_050);
trpuppi->SetBranchAddress("mpruned_zcut_010_R_cut_050",&puppi_mpruned_zcut_010_R_cut_050);
trpuppi->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_050",&puppi_ptprunedsafe_zcut_010_R_cut_050);
trpuppi->SetBranchAddress("mprunedsafe_zcut_010_R_cut_050",&puppi_mprunedsafe_zcut_010_R_cut_050);
trpuppi->SetBranchAddress("ptpruned_zcut_005_R_cut_050",&puppi_ptpruned_zcut_005_R_cut_050);
trpuppi->SetBranchAddress("mpruned_zcut_005_R_cut_050",&puppi_mpruned_zcut_005_R_cut_050);
trpuppi->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_050",&puppi_ptprunedsafe_zcut_005_R_cut_050);
trpuppi->SetBranchAddress("mprunedsafe_zcut_005_R_cut_050",&puppi_mprunedsafe_zcut_005_R_cut_050);
trpuppi->SetBranchAddress("ptpruned_zcut_005_R_cut_075",&puppi_ptpruned_zcut_005_R_cut_075);
trpuppi->SetBranchAddress("mpruned_zcut_005_R_cut_075",&puppi_mpruned_zcut_005_R_cut_075);
trpuppi->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_075",&puppi_ptprunedsafe_zcut_005_R_cut_075);
trpuppi->SetBranchAddress("mprunedsafe_zcut_005_R_cut_075",&puppi_mprunedsafe_zcut_005_R_cut_075);
trpuppi->SetBranchAddress("ptpruned_zcut_010_R_cut_075",&puppi_ptpruned_zcut_010_R_cut_075);
trpuppi->SetBranchAddress("mpruned_zcut_010_R_cut_075",&puppi_mpruned_zcut_010_R_cut_075);
trpuppi->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_075",&puppi_ptprunedsafe_zcut_010_R_cut_075);
trpuppi->SetBranchAddress("mprunedsafe_zcut_010_R_cut_075",&puppi_mprunedsafe_zcut_010_R_cut_075);
trpuppi->SetBranchAddress("ptsoftdrop_beta20",&puppi_ptsoftdrop_beta20);
trpuppi->SetBranchAddress("msoftdrop_beta20",&puppi_msoftdrop_beta20);
trpuppi->SetBranchAddress("ptsoftdropsafe_beta20",&puppi_ptsoftdropsafe_beta20);
trpuppi->SetBranchAddress("msoftdropsafe_beta20",&puppi_msoftdropsafe_beta20);
trpuppi->SetBranchAddress("ptsoftdrop_beta00",&puppi_ptsoftdrop_beta00);
trpuppi->SetBranchAddress("msoftdrop_beta00",&puppi_msoftdrop_beta00);
trpuppi->SetBranchAddress("ptsoftdropsafe_beta00",&puppi_ptsoftdropsafe_beta00);
trpuppi->SetBranchAddress("msoftdropsafe_beta00",&puppi_msoftdropsafe_beta00);
trpuppi->SetBranchAddress("ptsoftdrop_beta10",&puppi_ptsoftdrop_beta10);
trpuppi->SetBranchAddress("msoftdrop_beta10",&puppi_msoftdrop_beta10);
trpuppi->SetBranchAddress("ptsoftdropsafe_beta10",&puppi_ptsoftdropsafe_beta10);
trpuppi->SetBranchAddress("msoftdropsafe_beta10",&puppi_msoftdropsafe_beta10);
trpuppi->SetBranchAddress("ptsoftdrop_beta-10",&puppi_ptsoftdrop_betam1);
trpuppi->SetBranchAddress("msoftdrop_beta-10",&puppi_msoftdrop_betam1);
trpuppi->SetBranchAddress("ptsoftdropsafe_beta-10",&puppi_ptsoftdropsafe_betam1);
trpuppi->SetBranchAddress("msoftdropsafe_beta-10",&puppi_msoftdropsafe_betam1);
////////////////////////////////puppi

/////////////////////////////////gen
   vector<float>   *gen_pttrim_Rtrim_020_Ptfrac_005;
   vector<float>   *gen_mtrim_Rtrim_020_Ptfrac_005;
   vector<float>   *gen_pttrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *gen_mtrimsafe_Rtrim_020_Ptfrac_005;
   vector<float>   *gen_pttrim_Rtrim_010_Ptfrac_003;
   vector<float>   *gen_mtrim_Rtrim_010_Ptfrac_003;
   vector<float>   *gen_pttrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *gen_mtrimsafe_Rtrim_010_Ptfrac_003;
   vector<float>   *gen_pttrim_Rtrim_020_Ptfrac_003;
   vector<float>   *gen_mtrim_Rtrim_020_Ptfrac_003;
   vector<float>   *gen_pttrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *gen_mtrimsafe_Rtrim_020_Ptfrac_003;
   vector<float>   *gen_pttrim_Rtrim_030_Ptfrac_003;
   vector<float>   *gen_mtrim_Rtrim_030_Ptfrac_003;
   vector<float>   *gen_pttrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *gen_mtrimsafe_Rtrim_030_Ptfrac_003;
   vector<float>   *gen_ptconst;
   vector<float>   *gen_mconst;
   vector<float>   *gen_ptpruned_zcut_010_R_cut_050;
   vector<float>   *gen_mpruned_zcut_010_R_cut_050;
   vector<float>   *gen_ptprunedsafe_zcut_010_R_cut_050;
   vector<float>   *gen_mprunedsafe_zcut_010_R_cut_050;
   vector<float>   *gen_ptpruned_zcut_005_R_cut_050;
   vector<float>   *gen_mpruned_zcut_005_R_cut_050;
   vector<float>   *gen_ptprunedsafe_zcut_005_R_cut_050;
   vector<float>   *gen_mprunedsafe_zcut_005_R_cut_050;
   vector<float>   *gen_ptpruned_zcut_005_R_cut_075;
   vector<float>   *gen_mpruned_zcut_005_R_cut_075;
   vector<float>   *gen_ptprunedsafe_zcut_005_R_cut_075;
   vector<float>   *gen_mprunedsafe_zcut_005_R_cut_075;
   vector<float>   *gen_ptpruned_zcut_010_R_cut_075;
   vector<float>   *gen_mpruned_zcut_010_R_cut_075;
   vector<float>   *gen_ptprunedsafe_zcut_010_R_cut_075;
   vector<float>   *gen_mprunedsafe_zcut_010_R_cut_075;
   vector<float>   *gen_ptsoftdrop_beta20;
   vector<float>   *gen_msoftdrop_beta20;
   vector<float>   *gen_ptsoftdropsafe_beta20;
   vector<float>   *gen_msoftdropsafe_beta20;
   vector<float>   *gen_ptsoftdrop_beta00;
   vector<float>   *gen_msoftdrop_beta00;
   vector<float>   *gen_ptsoftdropsafe_beta00;
   vector<float>   *gen_msoftdropsafe_beta00;
   vector<float>   *gen_ptsoftdrop_beta10;
   vector<float>   *gen_msoftdrop_beta10;
   vector<float>   *gen_ptsoftdropsafe_beta10;
   vector<float>   *gen_msoftdropsafe_beta10;
   vector<float>   *gen_ptsoftdrop_betam1;
   vector<float>   *gen_msoftdrop_betam1;
   vector<float>   *gen_ptsoftdropsafe_betam1;
   vector<float>   *gen_msoftdropsafe_betam1;

trgen->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_005",&gen_pttrim_Rtrim_020_Ptfrac_005);
trgen->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_005",&gen_mtrim_Rtrim_020_Ptfrac_005);
trgen->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_005",&gen_pttrimsafe_Rtrim_020_Ptfrac_005);
trgen->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_005",&gen_mtrimsafe_Rtrim_020_Ptfrac_005);
trgen->SetBranchAddress("pttrim_Rtrim_010_Ptfrac_003",&gen_pttrim_Rtrim_010_Ptfrac_003);
trgen->SetBranchAddress("mtrim_Rtrim_010_Ptfrac_003",&gen_mtrim_Rtrim_010_Ptfrac_003);
trgen->SetBranchAddress("pttrimsafe_Rtrim_010_Ptfrac_003",&gen_pttrimsafe_Rtrim_010_Ptfrac_003);
trgen->SetBranchAddress("mtrimsafe_Rtrim_010_Ptfrac_003",&gen_mtrimsafe_Rtrim_010_Ptfrac_003);
trgen->SetBranchAddress("pttrim_Rtrim_020_Ptfrac_003",&gen_pttrim_Rtrim_020_Ptfrac_003);
trgen->SetBranchAddress("mtrim_Rtrim_020_Ptfrac_003",&gen_mtrim_Rtrim_020_Ptfrac_003);
trgen->SetBranchAddress("pttrimsafe_Rtrim_020_Ptfrac_003",&gen_pttrimsafe_Rtrim_020_Ptfrac_003);
trgen->SetBranchAddress("mtrimsafe_Rtrim_020_Ptfrac_003",&gen_mtrimsafe_Rtrim_020_Ptfrac_003);
trgen->SetBranchAddress("pttrim_Rtrim_030_Ptfrac_003",&gen_pttrim_Rtrim_030_Ptfrac_003);
trgen->SetBranchAddress("mtrim_Rtrim_030_Ptfrac_003",&gen_mtrim_Rtrim_030_Ptfrac_003);
trgen->SetBranchAddress("pttrimsafe_Rtrim_030_Ptfrac_003",&gen_pttrimsafe_Rtrim_030_Ptfrac_003);
trgen->SetBranchAddress("mtrimsafe_Rtrim_030_Ptfrac_003",&gen_mtrimsafe_Rtrim_030_Ptfrac_003);
trgen->SetBranchAddress("ptconst",&gen_ptconst);
trgen->SetBranchAddress("mconst",&gen_mconst);
trgen->SetBranchAddress("ptpruned_zcut_010_R_cut_050",&gen_ptpruned_zcut_010_R_cut_050);
trgen->SetBranchAddress("mpruned_zcut_010_R_cut_050",&gen_mpruned_zcut_010_R_cut_050);
trgen->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_050",&gen_ptprunedsafe_zcut_010_R_cut_050);
trgen->SetBranchAddress("mprunedsafe_zcut_010_R_cut_050",&gen_mprunedsafe_zcut_010_R_cut_050);
trgen->SetBranchAddress("ptpruned_zcut_005_R_cut_050",&gen_ptpruned_zcut_005_R_cut_050);
trgen->SetBranchAddress("mpruned_zcut_005_R_cut_050",&gen_mpruned_zcut_005_R_cut_050);
trgen->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_050",&gen_ptprunedsafe_zcut_005_R_cut_050);
trgen->SetBranchAddress("mprunedsafe_zcut_005_R_cut_050",&gen_mprunedsafe_zcut_005_R_cut_050);
trgen->SetBranchAddress("ptpruned_zcut_005_R_cut_075",&gen_ptpruned_zcut_005_R_cut_075);
trgen->SetBranchAddress("mpruned_zcut_005_R_cut_075",&gen_mpruned_zcut_005_R_cut_075);
trgen->SetBranchAddress("ptprunedsafe_zcut_005_R_cut_075",&gen_ptprunedsafe_zcut_005_R_cut_075);
trgen->SetBranchAddress("mprunedsafe_zcut_005_R_cut_075",&gen_mprunedsafe_zcut_005_R_cut_075);
trgen->SetBranchAddress("ptpruned_zcut_010_R_cut_075",&gen_ptpruned_zcut_010_R_cut_075);
trgen->SetBranchAddress("mpruned_zcut_010_R_cut_075",&gen_mpruned_zcut_010_R_cut_075);
trgen->SetBranchAddress("ptprunedsafe_zcut_010_R_cut_075",&gen_ptprunedsafe_zcut_010_R_cut_075);
trgen->SetBranchAddress("mprunedsafe_zcut_010_R_cut_075",&gen_mprunedsafe_zcut_010_R_cut_075);
trgen->SetBranchAddress("ptsoftdrop_beta20",&gen_ptsoftdrop_beta20);
trgen->SetBranchAddress("msoftdrop_beta20",&gen_msoftdrop_beta20);
trgen->SetBranchAddress("ptsoftdropsafe_beta20",&gen_ptsoftdropsafe_beta20);
trgen->SetBranchAddress("msoftdropsafe_beta20",&gen_msoftdropsafe_beta20);
trgen->SetBranchAddress("ptsoftdrop_beta00",&gen_ptsoftdrop_beta00);
trgen->SetBranchAddress("msoftdrop_beta00",&gen_msoftdrop_beta00);
trgen->SetBranchAddress("ptsoftdropsafe_beta00",&gen_ptsoftdropsafe_beta00);
trgen->SetBranchAddress("msoftdropsafe_beta00",&gen_msoftdropsafe_beta00);
trgen->SetBranchAddress("ptsoftdrop_beta10",&gen_ptsoftdrop_beta10);
trgen->SetBranchAddress("msoftdrop_beta10",&gen_msoftdrop_beta10);
trgen->SetBranchAddress("ptsoftdropsafe_beta10",&gen_ptsoftdropsafe_beta10);
trgen->SetBranchAddress("msoftdropsafe_beta10",&gen_msoftdropsafe_beta10);
trgen->SetBranchAddress("ptsoftdrop_beta-10",&gen_ptsoftdrop_betam1);
trgen->SetBranchAddress("msoftdrop_beta-10",&gen_msoftdrop_betam1);
trgen->SetBranchAddress("ptsoftdropsafe_beta-10",&gen_ptsoftdropsafe_betam1);
trgen->SetBranchAddress("msoftdropsafe_beta-10",&gen_msoftdropsafe_betam1);



int n = (int) trpf->GetEntries();
    

    cout<<"Total Entries = "<<n<<endl;
    cout<<"Code will run on "<<maxEvents<<" Entries"<<endl;
    for (int i = 0; i < n; ++i) {//loop over entries 
if(i != 12256 && i !=17340 && i !=36631 && i != 42571 && i != 42779){//bad events filtered
cout<<"Event = "<<i<<endl;
     if(i==maxEvents){
     break;
                     }

    trpf->GetEntry(i);
    trchs->GetEntry(i);
    trpuppi->GetEntry(i);
    trsk->GetEntry(i);
    trgen->GetEntry(i);

     const unsigned int  nPVbinsize=(nPVmax-nPVmin)/nPVnumOfbins;
            int nPV=npv;


 for(const unsigned int ii2=0;ii2<nPVnumOfbins;ii2++){//loop over the nPV bins

   const unsigned int nPVminLimit=nPVmin+ii2*nPVbinsize;
   const unsigned int nPVmaxLimit=nPVminLimit+nPVbinsize;

  //cout<<"bins = "<<nPVminLimit<<endl;

if(nPV >= nPVminLimit && nPV < nPVmaxLimit){//nPV bins condition
    //JetRsAKpfTr_nPV_h[0][ii2]->Fill(pf_mtrimsafe_Rtrim_020_Ptfrac_005->at(0)-gen_mtrim_Rtrim_020_Ptfrac_005->at(0));
    

////////////////////trimming
if(pf_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) > pTmin && pf_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) < pTmax){
      JetRsAKpfTr_nPV_h[0][ii2]->Fill(pf_mtrimsafe_Rtrim_020_Ptfrac_005->at(0)-gen_mtrim_Rtrim_020_Ptfrac_005->at(0));
          }

if(pf_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) > pTmin&& pf_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) < pTmax){
      JetRsAKpfTr_nPV_h[1][ii2]->Fill(pf_mtrimsafe_Rtrim_010_Ptfrac_003->at(0)-gen_mtrim_Rtrim_010_Ptfrac_003->at(0));
 }

if(pf_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) > pTmin&& pf_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) < pTmax){
      JetRsAKpfTr_nPV_h[2][ii2]->Fill(pf_mtrimsafe_Rtrim_020_Ptfrac_003->at(0)-gen_mtrim_Rtrim_020_Ptfrac_003->at(0));
 }

if(pf_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) > pTmin && pf_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) < pTmax){
      JetRsAKpfTr_nPV_h[3][ii2]->Fill(pf_mtrimsafe_Rtrim_030_Ptfrac_003->at(0)-gen_mtrim_Rtrim_030_Ptfrac_003->at(0));
 }

  if(chs_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) > pTmin && chs_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) < pTmax){
     JetRsAKchsTr_nPV_h [0][ii2]->Fill(chs_mtrimsafe_Rtrim_020_Ptfrac_005->at(0)-gen_mtrim_Rtrim_020_Ptfrac_005->at(0));
 }

if(chs_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) > pTmin&& chs_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) < pTmax){
      JetRsAKchsTr_nPV_h[1][ii2]->Fill(chs_mtrimsafe_Rtrim_010_Ptfrac_003->at(0)-gen_mtrim_Rtrim_010_Ptfrac_003->at(0));
 }

if(chs_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) > pTmin&& chs_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) < pTmax){
      JetRsAKchsTr_nPV_h[2][ii2]->Fill(chs_mtrimsafe_Rtrim_020_Ptfrac_003->at(0)-gen_mtrim_Rtrim_020_Ptfrac_003->at(0));
 }

if(chs_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) > pTmin && chs_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) < pTmax){
      JetRsAKchsTr_nPV_h[3][ii2]->Fill(chs_mtrimsafe_Rtrim_030_Ptfrac_003->at(0)-gen_mtrim_Rtrim_030_Ptfrac_003->at(0));
 }

if(puppi_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) > pTmin && puppi_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_005->at(0) < pTmax){
      JetRsAKpuppiTr_nPV_h[0][ii2]->Fill(puppi_mtrim_Rtrim_020_Ptfrac_005->at(0)-gen_mtrim_Rtrim_020_Ptfrac_005->at(0));
 }
if(puppi_pttrim_Rtrim_010_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) > pTmin && puppi_pttrim_Rtrim_010_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) < pTmax){
     JetRsAKpuppiTr_nPV_h [1][ii2]->Fill(puppi_mtrim_Rtrim_010_Ptfrac_003->at(0)-gen_mtrim_Rtrim_010_Ptfrac_003->at(0));
 }


if(puppi_pttrim_Rtrim_020_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) > pTmin && puppi_pttrim_Rtrim_020_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) < pTmax){
      JetRsAKpuppiTr_nPV_h[2][ii2]->Fill(puppi_mtrim_Rtrim_020_Ptfrac_003->at(0)-gen_mtrim_Rtrim_020_Ptfrac_003->at(0));
 }
if(puppi_pttrim_Rtrim_030_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) > pTmin && puppi_pttrim_Rtrim_030_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) < pTmax){
      JetRsAKpuppiTr_nPV_h[3][ii2]->Fill(puppi_mtrim_Rtrim_030_Ptfrac_003->at(0)-gen_mtrim_Rtrim_030_Ptfrac_003->at(0));
 }




/////////////////////trimming
  
////////////////////////Prunning below
if(pf_ptprunedsafe_zcut_010_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_050->at(0) > pTmin && pf_ptprunedsafe_zcut_010_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_050->at(0) < pTmax){
       JetRsAKpfPr_nPV_h[0][ii2]->Fill(pf_mprunedsafe_zcut_010_R_cut_050->at(0)-gen_mpruned_zcut_010_R_cut_050->at(0));
 }

if(pf_ptprunedsafe_zcut_005_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_050->at(0) > pTmin && pf_ptprunedsafe_zcut_005_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_050->at(0) < pTmax){
       JetRsAKpfPr_nPV_h[1][ii2]->Fill(pf_mprunedsafe_zcut_005_R_cut_050->at(0)-gen_mpruned_zcut_005_R_cut_050->at(0));
 }
if(pf_ptprunedsafe_zcut_005_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_075->at(0) > pTmin && pf_ptprunedsafe_zcut_005_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_075->at(0) < pTmax){
       JetRsAKpfPr_nPV_h[2][ii2]->Fill(pf_mprunedsafe_zcut_005_R_cut_075->at(0)-gen_mpruned_zcut_005_R_cut_075->at(0));
 }
if(pf_ptprunedsafe_zcut_010_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_075->at(0) > pTmin && pf_ptprunedsafe_zcut_010_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_075->at(0) < pTmax){
      JetRsAKpfPr_nPV_h [3][ii2]->Fill(pf_mprunedsafe_zcut_010_R_cut_075->at(0)-gen_mpruned_zcut_010_R_cut_075->at(0));
 }


if(chs_ptprunedsafe_zcut_010_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_050->at(0) > pTmin && chs_ptprunedsafe_zcut_010_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_050->at(0) < pTmax){
       JetRsAKchsPr_nPV_h[0][ii2]->Fill(chs_mprunedsafe_zcut_010_R_cut_050->at(0)-gen_mpruned_zcut_010_R_cut_050->at(0));
 }
if(chs_ptprunedsafe_zcut_005_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_050->at(0) > pTmin && chs_ptprunedsafe_zcut_005_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_050->at(0) < pTmax){
      JetRsAKchsPr_nPV_h[1][ii2]->Fill(chs_mprunedsafe_zcut_005_R_cut_050->at(0)-gen_mpruned_zcut_005_R_cut_050->at(0));
 }
if(chs_ptprunedsafe_zcut_005_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_075->at(0) > pTmin && chs_ptprunedsafe_zcut_005_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_075->at(0) < pTmax){
      JetRsAKchsPr_nPV_h[2][ii2]->Fill(chs_mprunedsafe_zcut_005_R_cut_075->at(0)-gen_mpruned_zcut_005_R_cut_075->at(0));
 }
if(chs_ptprunedsafe_zcut_010_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_075->at(0) > pTmin && chs_ptprunedsafe_zcut_010_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_075->at(0) < pTmax){
       JetRsAKchsPr_nPV_h[3][ii2]->Fill(chs_mprunedsafe_zcut_010_R_cut_075->at(0)-gen_mpruned_zcut_010_R_cut_075->at(0));
 }
if(puppi_ptpruned_zcut_010_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_050->at(0) > pTmin && puppi_ptpruned_zcut_010_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_050->at(0) < pTmax){
       JetRsAKpuppiPr_nPV_h[0][ii2]->Fill(puppi_mpruned_zcut_010_R_cut_050->at(0)-gen_mpruned_zcut_010_R_cut_050->at(0));
 }
if(puppi_ptpruned_zcut_005_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_050->at(0) > pTmin && puppi_ptpruned_zcut_005_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_050->at(0) < pTmax){
      JetRsAKpuppiPr_nPV_h[1][ii2]->Fill(puppi_mpruned_zcut_005_R_cut_050->at(0)-gen_mpruned_zcut_005_R_cut_050->at(0));
 }
if(puppi_ptpruned_zcut_005_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_075->at(0) > pTmin && puppi_ptpruned_zcut_005_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_075->at(0) < pTmax){
       JetRsAKpuppiPr_nPV_h[2][ii2]->Fill(puppi_mpruned_zcut_005_R_cut_075->at(0)-gen_mpruned_zcut_005_R_cut_075->at(0));
 }
if(puppi_ptpruned_zcut_010_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_075->at(0) > pTmin && puppi_ptpruned_zcut_010_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_075->at(0) < pTmax){
       JetRsAKpuppiPr_nPV_h[3][ii2]->Fill(puppi_mpruned_zcut_010_R_cut_075->at(0)-gen_mpruned_zcut_010_R_cut_075->at(0));
 }


//////////////////////////////softdrop below
if(pf_ptsoftdropsafe_beta20->at(0) > pTmin && gen_ptsoftdrop_beta20->at(0) > pTmin && pf_ptsoftdropsafe_beta20->at(0) < pTmax && gen_ptsoftdrop_beta20->at(0) < pTmax){
       JetRsAKpfSd_nPV_h[0][ii2]->Fill(pf_msoftdropsafe_beta20->at(0)-gen_msoftdrop_beta20->at(0));
 }
if(pf_ptsoftdropsafe_beta00->at(0) > pTmin && gen_ptsoftdrop_beta00->at(0) > pTmin && pf_ptsoftdropsafe_beta00->at(0) < pTmax && gen_ptsoftdrop_beta00->at(0) < pTmax){
      JetRsAKpfSd_nPV_h[1][ii2]->Fill(pf_msoftdropsafe_beta00->at(0)-gen_msoftdrop_beta00->at(0));
 }
if(pf_ptsoftdropsafe_beta10->at(0) > pTmin && gen_ptsoftdrop_beta10->at(0) > pTmin && pf_ptsoftdropsafe_beta10->at(0) < pTmax && gen_ptsoftdrop_beta10->at(0) < pTmax){
       JetRsAKpfSd_nPV_h[2][ii2]->Fill(pf_msoftdropsafe_beta10->at(0)-gen_msoftdrop_beta10->at(0));
 }
if(pf_ptsoftdropsafe_betam1->at(0) > pTmin && gen_ptsoftdrop_betam1->at(0) > pTmin && pf_ptsoftdropsafe_betam1->at(0) < pTmax && gen_ptsoftdrop_betam1->at(0) < pTmax){
       JetRsAKpfSd_nPV_h[3][ii2]->Fill(pf_msoftdropsafe_betam1->at(0)-gen_msoftdrop_betam1->at(0));
 }
if(chs_ptsoftdropsafe_beta20->at(0) > pTmin && gen_ptsoftdrop_beta20->at(0) > pTmin && chs_ptsoftdropsafe_beta20->at(0) < pTmax && gen_ptsoftdrop_beta20->at(0) < pTmax){
      JetRsAKchsSd_nPV_h[0][ii2]->Fill(chs_msoftdropsafe_beta20->at(0)-gen_msoftdrop_beta20->at(0));
 }
if(chs_ptsoftdropsafe_beta00->at(0) > pTmin && gen_ptsoftdrop_beta00->at(0) > pTmin && chs_ptsoftdropsafe_beta00->at(0) < pTmax && gen_ptsoftdrop_beta00->at(0) < pTmax){
       JetRsAKchsSd_nPV_h[1][ii2]->Fill(chs_msoftdropsafe_beta00->at(0)-gen_msoftdrop_beta00->at(0));
}
if(chs_ptsoftdropsafe_beta10->at(0) > pTmin && gen_ptsoftdrop_beta10->at(0) > pTmin && chs_ptsoftdropsafe_beta10->at(0) < pTmax && gen_ptsoftdrop_beta10->at(0) < pTmax){
       JetRsAKchsSd_nPV_h[2][ii2]->Fill(chs_msoftdropsafe_beta10->at(0)-gen_msoftdrop_beta10->at(0));
}
if(chs_ptsoftdropsafe_betam1->at(0) > pTmin && gen_ptsoftdrop_betam1->at(0) > pTmin && chs_ptsoftdropsafe_betam1->at(0) < pTmax && gen_ptsoftdrop_betam1->at(0) < pTmax){
       JetRsAKchsSd_nPV_h[3][ii2]->Fill(chs_msoftdropsafe_betam1->at(0)-gen_msoftdrop_betam1->at(0));
}
if(puppi_ptsoftdrop_beta20->at(0) > pTmin && gen_ptsoftdrop_beta20->at(0) > pTmin && puppi_ptsoftdrop_beta20->at(0) < pTmax && gen_ptsoftdrop_beta20->at(0) < pTmax){
       JetRsAKpuppiSd_nPV_h[0][ii2]->Fill(puppi_msoftdrop_beta20->at(0)-gen_msoftdrop_beta20->at(0));
}
if(puppi_ptsoftdrop_beta00->at(0) > pTmin && gen_ptsoftdrop_beta00->at(0) > pTmin && puppi_ptsoftdrop_beta00->at(0) < pTmax && gen_ptsoftdrop_beta00->at(0) < pTmax){
       JetRsAKpuppiSd_nPV_h[1][ii2]->Fill(puppi_msoftdrop_beta00->at(0)-gen_msoftdrop_beta00->at(0));
}
if(puppi_ptsoftdrop_beta10->at(0) > pTmin && gen_ptsoftdrop_beta10->at(0) > pTmin && puppi_ptsoftdrop_beta10->at(0) < pTmax && gen_ptsoftdrop_beta10->at(0) < pTmax){
       JetRsAKpuppiSd_nPV_h[2][ii2]->Fill(puppi_msoftdrop_beta10->at(0)-gen_msoftdrop_beta10->at(0));
}
if(puppi_ptsoftdrop_betam1->at(0) > pTmin && gen_ptsoftdrop_betam1->at(0) > pTmin && puppi_ptsoftdrop_betam1->at(0) < pTmax && gen_ptsoftdrop_betam1->at(0) < pTmax){
       JetRsAKpuppiSd_nPV_h[3][ii2]->Fill(puppi_msoftdrop_betam1->at(0)-gen_msoftdrop_betam1->at(0));
}



}//nPV bins condition

}//loop over the nPV bins


}//bad events filtered
}//loop over entries

////////////////////////trimmed below
double pf_tr_meanX_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_meanY_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_meanYerr_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_meanXerr_nPV[NumtrPar][nPVnumOfbins];


double pf_tr_sigmaX_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_sigmaY_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_sigmaYerr_nPV[NumtrPar][nPVnumOfbins];
double pf_tr_sigmaXerr_nPV[NumtrPar][nPVnumOfbins];


double chs_tr_meanX_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_meanY_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_meanYerr_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_meanXerr_nPV[NumtrPar][nPVnumOfbins];

double chs_tr_sigmaX_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_sigmaY_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_sigmaYerr_nPV[NumtrPar][nPVnumOfbins];
double chs_tr_sigmaXerr_nPV[NumtrPar][nPVnumOfbins];


double puppi_tr_meanX_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_meanY_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_meanYerr_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_meanXerr_nPV[NumtrPar][nPVnumOfbins];


double puppi_tr_sigmaX_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_sigmaY_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_sigmaYerr_nPV[NumtrPar][nPVnumOfbins];
double puppi_tr_sigmaXerr_nPV[NumtrPar][nPVnumOfbins];

////////////////////////////pruned below
double pf_pr_meanX_nPV[NumprPar][nPVnumOfbins];
double pf_pr_meanY_nPV[NumprPar][nPVnumOfbins];
double pf_pr_meanYerr_nPV[NumprPar][nPVnumOfbins];
double pf_pr_meanXerr_nPV[NumprPar][nPVnumOfbins];


double pf_pr_sigmaX_nPV[NumprPar][nPVnumOfbins];
double pf_pr_sigmaY_nPV[NumprPar][nPVnumOfbins];
double pf_pr_sigmaYerr_nPV[NumprPar][nPVnumOfbins];
double pf_pr_sigmaXerr_nPV[NumprPar][nPVnumOfbins];


double chs_pr_meanX_nPV[NumprPar][nPVnumOfbins];
double chs_pr_meanY_nPV[NumprPar][nPVnumOfbins];
double chs_pr_meanYerr_nPV[NumprPar][nPVnumOfbins];
double chs_pr_meanXerr_nPV[NumprPar][nPVnumOfbins];


double chs_pr_sigmaX_nPV[NumprPar][nPVnumOfbins];
double chs_pr_sigmaY_nPV[NumprPar][nPVnumOfbins];
double chs_pr_sigmaYerr_nPV[NumprPar][nPVnumOfbins];
double chs_pr_sigmaXerr_nPV[NumprPar][nPVnumOfbins];


double puppi_pr_meanX_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_meanY_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_meanYerr_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_meanXerr_nPV[NumprPar][nPVnumOfbins];


double puppi_pr_sigmaX_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_sigmaY_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_sigmaYerr_nPV[NumprPar][nPVnumOfbins];
double puppi_pr_sigmaXerr_nPV[NumprPar][nPVnumOfbins];

///////////////////////////softdrop

double pf_sd_meanX_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_meanY_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_meanYerr_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_meanXerr_nPV[NumsdPar][nPVnumOfbins];

double pf_sd_sigmaX_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_sigmaY_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_sigmaYerr_nPV[NumsdPar][nPVnumOfbins];
double pf_sd_sigmaXerr_nPV[NumsdPar][nPVnumOfbins];


double chs_sd_meanX_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_meanY_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_meanYerr_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_meanXerr_nPV[NumsdPar][nPVnumOfbins];


double chs_sd_sigmaX_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_sigmaY_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_sigmaYerr_nPV[NumsdPar][nPVnumOfbins];
double chs_sd_sigmaXerr_nPV[NumsdPar][nPVnumOfbins];


double puppi_sd_meanX_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_meanY_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_meanYerr_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_meanXerr_nPV[NumsdPar][nPVnumOfbins];


double puppi_sd_sigmaX_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_sigmaY_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_sigmaYerr_nPV[NumsdPar][nPVnumOfbins];
double puppi_sd_sigmaXerr_nPV[NumsdPar][nPVnumOfbins];


/////////////////////////////////////////////////////////////

for(const unsigned int ii3=0;ii3<nPVnumOfbins;ii3++){//loop over tr nPV bins
for(const unsigned int jj3=0;jj3<NumtrPar;jj3++){//loop over pars
////////tr
JetRsAKpfTr_nPV_h[jj3][ii3]->Scale(1/JetRsAKpfTr_nPV_h[jj3][ii3]->Integral());
JetRsAKchsTr_nPV_h[jj3][ii3]->Scale(1/JetRsAKchsTr_nPV_h[jj3][ii3]->Integral());
JetRsAKpuppiTr_nPV_h[jj3][ii3]->Scale(1/JetRsAKpuppiTr_nPV_h[jj3][ii3]->Integral());

pf_tr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_tr_meanY_nPV[jj3][ii3]=JetRsAKpfTr_nPV_h[jj3][ii3]->GetMean();
pf_tr_meanYerr_nPV[jj3][ii3]=JetRsAKpfTr_nPV_h[jj3][ii3]->GetMeanError();
pf_tr_meanXerr_nPV[jj3][ii3]=x_err;


pf_tr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_tr_sigmaY_nPV[jj3][ii3]=JetRsAKpfTr_nPV_h[jj3][ii3]->GetRMS();
pf_tr_sigmaYerr_nPV[jj3][ii3]=JetRsAKpfTr_nPV_h[jj3][ii3]->GetRMSError();
pf_tr_sigmaXerr_nPV[jj3][ii3]=x_err;


chs_tr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_tr_meanY_nPV[jj3][ii3]=JetRsAKchsTr_nPV_h[jj3][ii3]->GetMean();
chs_tr_meanYerr_nPV[jj3][ii3]=JetRsAKchsTr_nPV_h[jj3][ii3]->GetMeanError();
chs_tr_meanXerr_nPV[jj3][ii3]=x_err;

chs_tr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_tr_sigmaY_nPV[jj3][ii3]=JetRsAKchsTr_nPV_h[jj3][ii3]->GetRMS();
chs_tr_sigmaYerr_nPV[jj3][ii3]=JetRsAKchsTr_nPV_h[jj3][ii3]->GetRMSError();
chs_tr_sigmaXerr_nPV[jj3][ii3]=x_err;

puppi_tr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_tr_meanY_nPV[jj3][ii3]=JetRsAKpuppiTr_nPV_h[jj3][ii3]->GetMean();
puppi_tr_meanYerr_nPV[jj3][ii3]=JetRsAKpuppiTr_nPV_h[jj3][ii3]->GetMeanError();
puppi_tr_meanXerr_nPV[jj3][ii3]=x_err;

puppi_tr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_tr_sigmaY_nPV[jj3][ii3]=JetRsAKpuppiTr_nPV_h[jj3][ii3]->GetRMS();
puppi_tr_sigmaYerr_nPV[jj3][ii3]=JetRsAKpuppiTr_nPV_h[jj3][ii3]->GetRMSError();
puppi_tr_sigmaXerr_nPV[jj3][ii3]=x_err;

//tr

}//loop over pars
}//loop over tr nPV bins

for(const unsigned int ii3=0;ii3<nPVnumOfbins;ii3++){//loop over pr nPV bins
for(const unsigned int jj3=0;jj3<NumprPar;jj3++){//loop over pars
////////pr
JetRsAKpfPr_nPV_h[jj3][ii3]->Scale(1/JetRsAKpfPr_nPV_h[jj3][ii3]->Integral());
JetRsAKchsPr_nPV_h[jj3][ii3]->Scale(1/JetRsAKchsPr_nPV_h[jj3][ii3]->Integral());
JetRsAKpuppiPr_nPV_h[jj3][ii3]->Scale(1/JetRsAKpuppiPr_nPV_h[jj3][ii3]->Integral());

pf_pr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_pr_meanY_nPV[jj3][ii3]=JetRsAKpfPr_nPV_h[jj3][ii3]->GetMean();
pf_pr_meanYerr_nPV[jj3][ii3]=JetRsAKpfPr_nPV_h[jj3][ii3]->GetMeanError();
pf_pr_meanXerr_nPV[jj3][ii3]=x_err;

pf_pr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_pr_sigmaY_nPV[jj3][ii3]=JetRsAKpfPr_nPV_h[jj3][ii3]->GetRMS();
pf_pr_sigmaYerr_nPV[jj3][ii3]=JetRsAKpfPr_nPV_h[jj3][ii3]->GetRMSError();
pf_pr_sigmaXerr_nPV[jj3][ii3]=x_err;

chs_pr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_pr_meanY_nPV[jj3][ii3]=JetRsAKchsPr_nPV_h[jj3][ii3]->GetMean();
chs_pr_meanYerr_nPV[jj3][ii3]=JetRsAKchsPr_nPV_h[jj3][ii3]->GetMeanError();
chs_pr_meanXerr_nPV[jj3][ii3]=x_err;

chs_pr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_pr_sigmaY_nPV[jj3][ii3]=JetRsAKchsPr_nPV_h[jj3][ii3]->GetRMS();
chs_pr_sigmaYerr_nPV[jj3][ii3]=JetRsAKchsPr_nPV_h[jj3][ii3]->GetRMSError();
chs_pr_sigmaXerr_nPV[jj3][ii3]=x_err;

puppi_pr_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_pr_meanY_nPV[jj3][ii3]=JetRsAKpuppiPr_nPV_h[jj3][ii3]->GetMean();
puppi_pr_meanYerr_nPV[jj3][ii3]=JetRsAKpuppiPr_nPV_h[jj3][ii3]->GetMeanError();
puppi_pr_meanXerr_nPV[jj3][ii3]=x_err;

puppi_pr_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_pr_sigmaY_nPV[jj3][ii3]=JetRsAKpuppiPr_nPV_h[jj3][ii3]->GetRMS();
puppi_pr_sigmaYerr_nPV[jj3][ii3]=JetRsAKpuppiPr_nPV_h[jj3][ii3]->GetRMSError();
puppi_pr_sigmaXerr_nPV[jj3][ii3]=x_err;
//pr

}//loop over pars
}//loop over pr nPV bins


for(const unsigned int ii3=0;ii3<nPVnumOfbins;ii3++){//loop over sd nPV bins
for(const unsigned int jj3=0;jj3<NumsdPar;jj3++){//loop over pars
////////sd
JetRsAKpfSd_nPV_h[jj3][ii3]->Scale(1/JetRsAKpfSd_nPV_h[jj3][ii3]->Integral());
JetRsAKchsSd_nPV_h[jj3][ii3]->Scale(1/JetRsAKchsSd_nPV_h[jj3][ii3]->Integral());
JetRsAKpuppiSd_nPV_h[jj3][ii3]->Scale(1/JetRsAKpuppiSd_nPV_h[jj3][ii3]->Integral());



pf_sd_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_sd_meanY_nPV[jj3][ii3]=JetRsAKpfSd_nPV_h[jj3][ii3]->GetMean();
pf_sd_meanYerr_nPV[jj3][ii3]=JetRsAKpfSd_nPV_h[jj3][ii3]->GetMeanError();
pf_sd_meanXerr_nPV[jj3][ii3]=x_err;

pf_sd_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
pf_sd_sigmaY_nPV[jj3][ii3]=JetRsAKpfSd_nPV_h[jj3][ii3]->GetRMS();
pf_sd_sigmaYerr_nPV[jj3][ii3]=JetRsAKpfSd_nPV_h[jj3][ii3]->GetRMSError();
pf_sd_sigmaXerr_nPV[jj3][ii3]=x_err;

chs_sd_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_sd_meanY_nPV[jj3][ii3]=JetRsAKchsSd_nPV_h[jj3][ii3]->GetMean();
chs_sd_meanYerr_nPV[jj3][ii3]=JetRsAKchsSd_nPV_h[jj3][ii3]->GetMeanError();
chs_sd_meanXerr_nPV[jj3][ii3]=x_err;

chs_sd_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
chs_sd_sigmaY_nPV[jj3][ii3]=JetRsAKchsSd_nPV_h[jj3][ii3]->GetRMS();
chs_sd_sigmaYerr_nPV[jj3][ii3]=JetRsAKchsSd_nPV_h[jj3][ii3]->GetRMSError();
chs_sd_sigmaXerr_nPV[jj3][ii3]=x_err;

puppi_sd_meanX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_sd_meanY_nPV[jj3][ii3]=JetRsAKpuppiSd_nPV_h[jj3][ii3]->GetMean();
puppi_sd_meanYerr_nPV[jj3][ii3]=JetRsAKpuppiSd_nPV_h[jj3][ii3]->GetMeanError();
puppi_sd_meanXerr_nPV[jj3][ii3]=x_err;

puppi_sd_sigmaX_nPV[jj3][ii3]=nPVmin+ii3*nPVbinsize +(binsize/2);
puppi_sd_sigmaY_nPV[jj3][ii3]=JetRsAKpuppiSd_nPV_h[jj3][ii3]->GetRMS();
puppi_sd_sigmaYerr_nPV[jj3][ii3]=JetRsAKpuppiSd_nPV_h[jj3][ii3]->GetRMSError();
puppi_sd_sigmaXerr_nPV[jj3][ii3]=x_err;
//sd

}//loop over pars
}//loop over sd nPV bins





TGraphErrors *pf_graph_trMeanvsNpv[NumtrPar];
TGraphErrors *pf_graph_trSigmavsNpv[NumtrPar];

TGraphErrors *chs_graph_trMeanvsNpv[NumtrPar];
TGraphErrors *chs_graph_trSigmavsNpv[NumtrPar];

TGraphErrors *puppi_graph_trMeanvsNpv[NumtrPar];
TGraphErrors *puppi_graph_trSigmavsNpv[NumtrPar];


TGraphErrors *pf_graph_prMeanvsNpv[NumprPar];
TGraphErrors *pf_graph_prSigmavsNpv[NumprPar];

TGraphErrors *chs_graph_prMeanvsNpv[NumprPar];
TGraphErrors *chs_graph_prSigmavsNpv[NumprPar];

TGraphErrors *puppi_graph_prMeanvsNpv[NumprPar];
TGraphErrors *puppi_graph_prSigmavsNpv[NumprPar];

TGraphErrors *pf_graph_sdMeanvsNpv[NumsdPar];
TGraphErrors *pf_graph_sdSigmavsNpv[NumsdPar];

TGraphErrors *chs_graph_sdMeanvsNpv[NumsdPar];
TGraphErrors *chs_graph_sdSigmavsNpv[NumsdPar];

TGraphErrors *puppi_graph_sdMeanvsNpv[NumsdPar];
TGraphErrors *puppi_graph_sdSigmavsNpv[NumsdPar];


for(int j1=0;j1<NumtrPar;j1++){
pf_graph_trMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_tr_meanX_nPV[j1],pf_tr_meanY_nPV[j1],pf_tr_meanXerr_nPV[j1],pf_tr_meanYerr_nPV[j1]);
pf_graph_trMeanvsNpv[j1]->SetLineColor(j1+1);
pf_graph_trMeanvsNpv[j1]->SetLineWidth(linewidth);
pf_graph_trMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_trMeanvsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_trMeanvsNpv[j1]->SetMarkerSize(markersize);

pf_graph_trSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_tr_sigmaX_nPV[j1],pf_tr_sigmaY_nPV[j1],pf_tr_sigmaXerr_nPV[j1],pf_tr_sigmaYerr_nPV[j1]);
pf_graph_trSigmavsNpv[j1]->SetLineColor(j1+1);
pf_graph_trSigmavsNpv[j1]->SetLineWidth(linewidth);
pf_graph_trSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_trSigmavsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_trSigmavsNpv[j1]->SetMarkerSize(markersize);

chs_graph_trMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_tr_meanX_nPV[j1],chs_tr_meanY_nPV[j1],chs_tr_meanXerr_nPV[j1],chs_tr_meanYerr_nPV[j1]);
chs_graph_trMeanvsNpv[j1]->SetLineColor(j1+1);
chs_graph_trMeanvsNpv[j1]->SetLineWidth(linewidth);
chs_graph_trMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_trMeanvsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_trMeanvsNpv[j1]->SetMarkerSize(markersize);

chs_graph_trSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_tr_sigmaX_nPV[j1],chs_tr_sigmaY_nPV[j1],chs_tr_sigmaXerr_nPV[j1],chs_tr_sigmaYerr_nPV[j1]);
chs_graph_trSigmavsNpv[j1]->SetLineColor(j1+1);
chs_graph_trSigmavsNpv[j1]->SetLineWidth(linewidth);
chs_graph_trSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_trSigmavsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_trSigmavsNpv[j1]->SetMarkerSize(markersize);

puppi_graph_trMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_tr_meanX_nPV[j1],puppi_tr_meanY_nPV[j1],puppi_tr_meanXerr_nPV[j1],puppi_tr_meanYerr_nPV[j1]);
puppi_graph_trMeanvsNpv[j1]->SetLineColor(j1+1);
puppi_graph_trMeanvsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_trMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_trMeanvsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_trMeanvsNpv[j1]->SetMarkerSize(markersize);

puppi_graph_trSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_tr_sigmaX_nPV[j1],puppi_tr_sigmaY_nPV[j1],puppi_tr_sigmaXerr_nPV[j1],puppi_tr_sigmaYerr_nPV[j1]);
puppi_graph_trSigmavsNpv[j1]->SetLineColor(j1+1);
puppi_graph_trSigmavsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_trSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_trSigmavsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_trSigmavsNpv[j1]->SetMarkerSize(markersize);

}

for(int j1=0;j1<NumprPar;j1++){
pf_graph_prMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_pr_meanX_nPV[j1],pf_pr_meanY_nPV[j1],pf_pr_meanXerr_nPV[j1],pf_pr_meanYerr_nPV[j1]);
pf_graph_prMeanvsNpv[j1]->SetLineColor(j1+1);
pf_graph_prMeanvsNpv[j1]->SetLineWidth(linewidth);
pf_graph_prMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_prMeanvsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_prMeanvsNpv[j1]->SetMarkerSize(markersize);

pf_graph_prSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_pr_sigmaX_nPV[j1],pf_pr_sigmaY_nPV[j1],pf_pr_sigmaXerr_nPV[j1],pf_pr_sigmaYerr_nPV[j1]);
pf_graph_prSigmavsNpv[j1]->SetLineColor(j1+1);
pf_graph_prSigmavsNpv[j1]->SetLineWidth(linewidth);
pf_graph_prSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_prSigmavsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_prSigmavsNpv[j1]->SetMarkerSize(markersize);

chs_graph_prMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_pr_meanX_nPV[j1],chs_pr_meanY_nPV[j1],chs_pr_meanXerr_nPV[j1],chs_pr_meanYerr_nPV[j1]);
chs_graph_prMeanvsNpv[j1]->SetLineColor(j1+1);
chs_graph_prMeanvsNpv[j1]->SetLineWidth(linewidth);
chs_graph_prMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_prMeanvsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_prMeanvsNpv[j1]->SetMarkerSize(markersize);

chs_graph_prSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_pr_sigmaX_nPV[j1],chs_pr_sigmaY_nPV[j1],chs_pr_sigmaXerr_nPV[j1],chs_pr_sigmaYerr_nPV[j1]);
chs_graph_prSigmavsNpv[j1]->SetLineColor(j1+1);
chs_graph_prSigmavsNpv[j1]->SetLineWidth(linewidth);
chs_graph_prSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_prSigmavsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_prSigmavsNpv[j1]->SetMarkerSize(markersize);

puppi_graph_prMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_pr_meanX_nPV[j1],puppi_pr_meanY_nPV[j1],puppi_pr_meanXerr_nPV[j1],puppi_pr_meanYerr_nPV[j1]);
puppi_graph_prMeanvsNpv[j1]->SetLineColor(j1+1);
puppi_graph_prMeanvsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_prMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_prMeanvsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_prMeanvsNpv[j1]->SetMarkerSize(markersize);

puppi_graph_prSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_pr_sigmaX_nPV[j1],puppi_pr_sigmaY_nPV[j1],puppi_pr_sigmaXerr_nPV[j1],puppi_pr_sigmaYerr_nPV[j1]);
puppi_graph_prSigmavsNpv[j1]->SetLineColor(j1+1);
puppi_graph_prSigmavsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_prSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_prSigmavsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_prSigmavsNpv[j1]->SetMarkerSize(markersize);
}


for(int j1=0;j1<NumsdPar;j1++){
pf_graph_sdMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_sd_meanX_nPV[j1],pf_sd_meanY_nPV[j1],pf_sd_meanXerr_nPV[j1],pf_sd_meanYerr_nPV[j1]);
pf_graph_sdMeanvsNpv[j1]->SetLineColor(j1+1);
pf_graph_sdMeanvsNpv[j1]->SetLineWidth(linewidth);
pf_graph_sdMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_sdMeanvsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_sdMeanvsNpv[j1]->SetMarkerSize(markersize);
pf_graph_sdSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,pf_sd_sigmaX_nPV[j1],pf_sd_sigmaY_nPV[j1],pf_sd_sigmaXerr_nPV[j1],pf_sd_sigmaYerr_nPV[j1]);
pf_graph_sdSigmavsNpv[j1]->SetLineColor(j1+1);
pf_graph_sdSigmavsNpv[j1]->SetLineWidth(linewidth);
pf_graph_sdSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
pf_graph_sdSigmavsNpv[j1]->SetMarkerColor(j1+1);
pf_graph_sdSigmavsNpv[j1]->SetMarkerSize(markersize);
chs_graph_sdMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_sd_meanX_nPV[j1],chs_sd_meanY_nPV[j1],chs_sd_meanXerr_nPV[j1],chs_sd_meanYerr_nPV[j1]);
chs_graph_sdMeanvsNpv[j1]->SetLineColor(j1+1);
chs_graph_sdMeanvsNpv[j1]->SetLineWidth(linewidth);
chs_graph_sdMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_sdMeanvsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_sdMeanvsNpv[j1]->SetMarkerSize(markersize);
chs_graph_sdSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,chs_sd_sigmaX_nPV[j1],chs_sd_sigmaY_nPV[j1],chs_sd_sigmaXerr_nPV[j1],chs_sd_sigmaYerr_nPV[j1]);
chs_graph_sdSigmavsNpv[j1]->SetLineColor(j1+1);
chs_graph_sdSigmavsNpv[j1]->SetLineWidth(linewidth);
chs_graph_sdSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
chs_graph_sdSigmavsNpv[j1]->SetMarkerColor(j1+1);
chs_graph_sdSigmavsNpv[j1]->SetMarkerSize(markersize);
puppi_graph_sdMeanvsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_sd_meanX_nPV[j1],puppi_sd_meanY_nPV[j1],puppi_sd_meanXerr_nPV[j1],puppi_sd_meanYerr_nPV[j1]);
puppi_graph_sdMeanvsNpv[j1]->SetLineColor(j1+1);
puppi_graph_sdMeanvsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_sdMeanvsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_sdMeanvsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_sdMeanvsNpv[j1]->SetMarkerSize(markersize);
puppi_graph_sdSigmavsNpv[j1]  = new TGraphErrors(nPVnumOfbins,puppi_sd_sigmaX_nPV[j1],puppi_sd_sigmaY_nPV[j1],puppi_sd_sigmaXerr_nPV[j1],puppi_sd_sigmaYerr_nPV[j1]);
puppi_graph_sdSigmavsNpv[j1]->SetLineColor(j1+1);
puppi_graph_sdSigmavsNpv[j1]->SetLineWidth(linewidth);
puppi_graph_sdSigmavsNpv[j1]->SetMarkerStyle(markerstyle);
puppi_graph_sdSigmavsNpv[j1]->SetMarkerColor(j1+1);
puppi_graph_sdSigmavsNpv[j1]->SetMarkerSize(markersize);


}









gROOT->ProcessLine(".L tdrstyle.C");
 setTDRStyle();


char Cname1[100];
 
TCanvas *c_1D[24]; 
for(int k0=0;k0<24;k0++){
sprintf(Cname1,"c_1D%i",k0);

c_1D[k0]=new TCanvas(Cname1,Cname1,1000,800);

}

char Legname1[100];


TLegend *leg_1D[24]; 

for(int k0=0;k0<24;k0++){
sprintf(Legname1,"leg_1D%i",k0);

leg_1D[k0]=new TLegend(0.6074297,0.6101036,0.8574297,0.9002591);//0.5,0.6,0.750,0.89);//0.6,0.2,0.95,0.85);
leg_1D[k0]->SetTextFont(42);
leg_1D[k0]->SetLineColor(1);
leg_1D[k0]->SetLineStyle(1);
leg_1D[k0]->SetLineWidth(3);
leg_1D[k0]->SetFillColor(0);
leg_1D[k0]->SetFillStyle(1001);
leg_1D[k0]->SetShadowColor(0);
leg_1D[k0]->SetDrawOption(0);
leg_1D[k0]->SetBorderSize(0);
leg_1D[k0]->SetTextSize(0.0388601);



}




char TMgname1[100];


TMultiGraph *TMg_1D[24]; 
 

for(int k0=0;k0<24;k0++){
sprintf(TMgname1,"TMg_1D%i",k0);

TMg_1D[k0]=new TMultiGraph();

}


for(int j0=0;j0<12;j0++){
TMg_1D[j0]->SetMinimum(-20);
TMg_1D[j0]->SetMaximum(60);

}

for(int j0=12;j0<24;j0++){
TMg_1D[j0]->SetMinimum(0);
TMg_1D[j0]->SetMaximum(60);

}


////////////////////////////////////////////////////////////////Mean vs nPV

/*
TLatex *tex1 = new TLatex(13.37756,60.34121,"CMS Simulation Preliminary, #sqrt{s}= 13 TeV");
   tex1->SetLineWidth(0.35);
   tex1->SetTextSize(0.045);
   
 TLatex *tex2 = new TLatex(14.77879,54.93871,"QCD, Anti-kT (R=0.8)");
   tex2->SetLineWidth(0.03);
   tex2->SetTextSize(0.04);

TLatex *tex3 = new TLatex(15.74516,50.09225,"300 <p_{T} <470 GeV");
   tex3->SetLineWidth(0.03);
   tex3->SetTextSize(0.04);

TLatex *tex4 = new TLatex(17.38798,45.4164," |#eta|< 2.5");
   tex4->SetLineWidth(0.03);
   tex4->SetTextSize(0.04);

*/


TLatex *   tex1 = new TLatex(13.79871,60.62555,"CMS");
   tex1->SetTextSize(0.0388601);
   tex1->SetLineWidth(2);
  
 TLatex *tex2 = new TLatex(17.23991,60.62555,"Simulation Preliminary");
   tex2->SetTextFont(52);
   tex2->SetTextSize(0.03238342);
   tex2->SetLineWidth(2);
   
 TLatex *tex3 = new TLatex(47.73935,60.62555,"13 TeV");
   tex3->SetTextFont(42);
   tex3->SetTextSize(0.03367876);
   tex3->SetLineWidth(2);
   
TLatex *tex4 = new TLatex(14.69436,55.22305,"QCD, Anti-kT (R=0.8)");
   tex4->SetTextFont(42);
   tex4->SetTextSize(0.04404145);
   tex4->SetLineWidth(0);
   
tex5 = new TLatex(14.69436,51.0527,"<n_{PU}>= 40");
   tex5->SetTextFont(42);
   tex5->SetTextSize(0.04404145);
   tex5->SetLineWidth(2);


 TLatex *tex6 = new TLatex(14.69436,46.598,"p_{T} > 300 GeV");
   tex6->SetTextFont(42);
   tex6->SetTextSize(0.04404145);
   tex6->SetLineWidth(0);
   
TLatex *tex7 = new TLatex(14.69436,41.47984,"|#eta|< 2.5");
   tex7->SetTextFont(42);
   tex7->SetTextSize(0.04404145);
   tex7->SetLineWidth(0);














TLatex *texc12 = new TLatex(18.35434,40.23506,"PF ");
   texc12->SetLineWidth(0.05);
   texc12->SetTextSize(0.04);

c_1D[0]->cd();
TMg_1D[0]->Add(pf_graph_trMeanvsNpv[0]);
TMg_1D[0]->Add(pf_graph_trMeanvsNpv[1]);
TMg_1D[0]->Add(pf_graph_trMeanvsNpv[2]);
TMg_1D[0]->Add(pf_graph_trMeanvsNpv[3]);

TMg_1D[0]->SetTitle(" ;n_{PV} ;Resolution (GeV) ");
TMg_1D[0]->Draw("AP");

leg_1D[0]->AddEntry(pf_graph_trMeanvsNpv[0],"Trimming (r_{sub}=0.2,pT_{frac}=0.05)","l");
leg_1D[0]->AddEntry(pf_graph_trMeanvsNpv[1],"Trimming (r_{sub}=0.1,pT_{frac}=0.03)","l");
leg_1D[0]->AddEntry(pf_graph_trMeanvsNpv[2],"Trimming (r_{sub}=0.2,pT_{frac}=0.03) ","l");
leg_1D[0]->AddEntry(pf_graph_trMeanvsNpv[3],"Trimming (r_{sub}=0.3,pT_{frac}=0.03)","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();



texc12->Draw();
leg_1D[0]->Draw();

////////////////////////////tr2

TLatex *texc13 = new TLatex(14.69436,40.23506,"PF+CHS with trimming");
   texc13->SetLineWidth(0.05);
   texc13->SetTextSize(0.04);

c_1D[1]->cd();
TMg_1D[1]->Add(chs_graph_trMeanvsNpv[0]);
TMg_1D[1]->Add(chs_graph_trMeanvsNpv[1]);
TMg_1D[1]->Add(chs_graph_trMeanvsNpv[2]);
TMg_1D[1]->Add(chs_graph_trMeanvsNpv[3]);
TMg_1D[1]->SetTitle(" ;n_{PV} ;Resolution (GeV) ");
TMg_1D[1]->Draw("AP");

leg_1D[1]->AddEntry(chs_graph_trMeanvsNpv[0],"Trimming (r_{sub}=0.2,pT_{frac}=0.05)","l");
leg_1D[1]->AddEntry(chs_graph_trMeanvsNpv[1],"Trimming (r_{sub}=0.1,pT_{frac}=0.03)","l");
leg_1D[1]->AddEntry(chs_graph_trMeanvsNpv[2],"Trimming (r_{sub}=0.2,pT_{frac}=0.03) ","l");
leg_1D[1]->AddEntry(chs_graph_trMeanvsNpv[3],"Trimming (r_{sub}=0.3,pT_{frac}=0.03)","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc13->Draw();
leg_1D[1]->Draw();


////////////////////////////////////////tr3


TLatex *texc14 = new TLatex(18.35434,40.23506,"PF+PUPPI ");
   texc14->SetLineWidth(0.05);
   texc14->SetTextSize(0.04);

c_1D[2]->cd();
TMg_1D[2]->Add(puppi_graph_trMeanvsNpv[0]);
TMg_1D[2]->Add(puppi_graph_trMeanvsNpv[1]);
TMg_1D[2]->Add(puppi_graph_trMeanvsNpv[2]);
TMg_1D[2]->Add(puppi_graph_trMeanvsNpv[3]);

TMg_1D[2]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[2]->Draw("AP");

leg_1D[2]->AddEntry(puppi_graph_trMeanvsNpv[0],"Trimming (r_{sub}=0.2,pT_{frac}=0.05)","l");
leg_1D[2]->AddEntry(puppi_graph_trMeanvsNpv[1],"Trimming (r_{sub}=0.1,pT_{frac}=0.03)","l");
leg_1D[2]->AddEntry(puppi_graph_trMeanvsNpv[2],"Trimming (r_{sub}=0.2,pT_{frac}=0.03) ","l");
leg_1D[2]->AddEntry(puppi_graph_trMeanvsNpv[3],"Trimming (r_{sub}=0.3,pT_{frac}=0.03)","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc14->Draw();
leg_1D[2]->Draw();

///////////////////////////////////////////////tr4
/*
TLatex *texc15 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.3,pT_{frac}=0.03)");
   texc15->SetLineWidth(0.05);
   texc15->SetTextSize(0.03);

c_1D[3]->cd();
TMg_1D[3]->Add(pf_graph_trMeanvsNpv[3]);
TMg_1D[3]->Add(chs_graph_trMeanvsNpv[3]);
TMg_1D[3]->Add(puppi_graph_trMeanvsNpv[3]);
TMg_1D[3]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[3]->Draw("AP");

leg_1D[3]->AddEntry(pf_graph_trMeanvsNpv[3],"PF  ","lp");
leg_1D[3]->AddEntry(chs_graph_trMeanvsNpv[3],"PF+CHS ","lp");
leg_1D[3]->AddEntry(puppi_graph_trMeanvsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc15->Draw();
leg_1D[3]->Draw();
*/
//////////////////////////////////////////////////pr1

TLatex *texc16 = new TLatex(18.35434,40.23506,"PF");
   texc16->SetLineWidth(0.05);
   texc16->SetTextSize(0.04);

c_1D[4]->cd();
TMg_1D[4]->Add(pf_graph_prMeanvsNpv[0]);
TMg_1D[4]->Add(pf_graph_prMeanvsNpv[1]);
TMg_1D[4]->Add(pf_graph_prMeanvsNpv[2]);
TMg_1D[4]->Add(pf_graph_prMeanvsNpv[3]);

TMg_1D[4]->SetTitle(" ;n_{PV} ;Resolution (GeV)");
TMg_1D[4]->Draw("AP");

leg_1D[4]->AddEntry(pf_graph_prMeanvsNpv[0],"Pruning (z_{cut}=0.1,r_{cut}=0.5) ","l");
leg_1D[4]->AddEntry(pf_graph_prMeanvsNpv[1],"Prun (z_{cut}=0.05,r_{cut}=0.5) ","l");
leg_1D[4]->AddEntry(pf_graph_prMeanvsNpv[2],"Pruned (z_{cut}=0.05,r_{cut}=0.75)","l");
leg_1D[4]->AddEntry(pf_graph_prMeanvsNpv[3],"Pruned (z_{cut}=0.1,r_{cut}=0.75)","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc16->Draw();
leg_1D[4]->Draw();

////////////////////////////////////////pr2
TLatex *texc17 = new TLatex(18.35434,40.23506,"PF+CHS");
   texc17->SetLineWidth(0.05);
   texc17->SetTextSize(0.04);

c_1D[5]->cd();
TMg_1D[5]->Add(chs_graph_prMeanvsNpv[0]);
TMg_1D[5]->Add(chs_graph_prMeanvsNpv[1]);
TMg_1D[5]->Add(chs_graph_prMeanvsNpv[2]);
TMg_1D[5]->Add(chs_graph_prMeanvsNpv[3]);

TMg_1D[5]->SetTitle(" ;n_{PV} ;Resolution (GeV)");
TMg_1D[5]->Draw("AP");

leg_1D[5]->AddEntry(chs_graph_prMeanvsNpv[0],"Pruned (z_{cut}=0.1,r_{cut}=0.5) ","l");
leg_1D[5]->AddEntry(chs_graph_prMeanvsNpv[1],"Pruned (z_{cut}=0.05,r_{cut}=0.5) ","l");
leg_1D[5]->AddEntry(chs_graph_prMeanvsNpv[2],"Pruned (z_{cut}=0.05,r_{cut}=0.75)","l");
leg_1D[5]->AddEntry(chs_graph_prMeanvsNpv[3],"Pruned (z_{cut}=0.1,r_{cut}=0.75)","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc17->Draw();
leg_1D[5]->Draw();



///////////////////////////////////////////pr3

TLatex *texc18 = new TLatex(18.35434,40.23506,"PF+PUPPI ");
   texc18->SetLineWidth(0.05);
   texc18->SetTextSize(0.04);

c_1D[6]->cd();
TMg_1D[6]->Add(puppi_graph_prMeanvsNpv[0]);
TMg_1D[6]->Add(puppi_graph_prMeanvsNpv[1]);
TMg_1D[6]->Add(puppi_graph_prMeanvsNpv[2]);
TMg_1D[6]->Add(puppi_graph_prMeanvsNpv[3]);

TMg_1D[6]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[6]->Draw("AP");

leg_1D[6]->AddEntry(puppi_graph_prMeanvsNpv[0],"Pruned (z_{cut}=0.1,r_{cut}=0.5) ","l");
leg_1D[6]->AddEntry(puppi_graph_prMeanvsNpv[1],"Pruned (z_{cut}=0.05,r_{cut}=0.5) ","l");
leg_1D[6]->AddEntry(puppi_graph_prMeanvsNpv[2],"Pruned (z_{cut}=0.05,r_{cut}=0.75)","l");
leg_1D[6]->AddEntry(puppi_graph_prMeanvsNpv[3],"Pruned (z_{cut}=0.1,r_{cut}=0.75)","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc18->Draw();
leg_1D[6]->Draw();







//////////////////////////////////////////pr4
/*
TLatex *texc19 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.05,r_{cut}=0.75)");
   texc19->SetLineWidth(0.05);
   texc19->SetTextSize(0.03);

c_1D[7]->cd();
TMg_1D[7]->Add(pf_graph_prMeanvsNpv[3]);
TMg_1D[7]->Add(chs_graph_prMeanvsNpv[3]);
TMg_1D[7]->Add(puppi_graph_prMeanvsNpv[3]);
TMg_1D[7]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[7]->Draw("AP");

leg_1D[7]->AddEntry(pf_graph_prMeanvsNpv[3],"PF  ","lp");
leg_1D[7]->AddEntry(chs_graph_prMeanvsNpv[3],"PF+CHS ","lp");
leg_1D[7]->AddEntry(puppi_graph_prMeanvsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc19->Draw();
leg_1D[7]->Draw();
*/
///////////////////////////////sd1

TLatex *texc20 = new TLatex(18.35434,40.23506,"PF ");
   texc20->SetLineWidth(0.05);
   texc20->SetTextSize(0.04);

c_1D[8]->cd();
TMg_1D[8]->Add(pf_graph_sdMeanvsNpv[0]);
TMg_1D[8]->Add(pf_graph_sdMeanvsNpv[1]);
TMg_1D[8]->Add(pf_graph_sdMeanvsNpv[2]);
//TMg_1D[8]->Add(pf_graph_sdMeanvsNpv[3]);

TMg_1D[8]->SetTitle(" ;n_{PV} ;Resolution (GeV) ");
TMg_1D[8]->Draw("AP");

leg_1D[8]->AddEntry(pf_graph_sdMeanvsNpv[0],"Softdrop #beta = 2 ","l");
leg_1D[8]->AddEntry(pf_graph_sdMeanvsNpv[1],"Softdrop #beta = 0 ","l");
leg_1D[8]->AddEntry(pf_graph_sdMeanvsNpv[2],"Softdrop #beta = 1 ","l");
//leg_1D[8]->AddEntry(pf_graph_sdMeanvsNpv[3],"Softdrop #beta = -1 ","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc20->Draw();
leg_1D[8]->Draw();

////////////////////////////////sd2

TLatex *texc21 = new TLatex(18.35434,40.23506,"PF+CHS ");
   texc21->SetLineWidth(0.05);
   texc21->SetTextSize(0.04);

c_1D[9]->cd();
TMg_1D[9]->Add(chs_graph_sdMeanvsNpv[0]);
TMg_1D[9]->Add(chs_graph_sdMeanvsNpv[1]);
TMg_1D[9]->Add(chs_graph_sdMeanvsNpv[2]);
//TMg_1D[9]->Add(chs_graph_sdMeanvsNpv[3]);

TMg_1D[9]->SetTitle(" ;n_{PV} ;Resolution (GeV) ");
TMg_1D[9]->Draw("AP");

leg_1D[9]->AddEntry(chs_graph_sdMeanvsNpv[0],"Softdrop #beta = 2 ","l");
leg_1D[9]->AddEntry(chs_graph_sdMeanvsNpv[1],"Softdrop #beta = 0 ","l");
leg_1D[9]->AddEntry(chs_graph_sdMeanvsNpv[2],"Softdrop #beta = 1 ","l");
//leg_1D[9]->AddEntry(chs_graph_sdMeanvsNpv[3],"Softdrop #beta = -1 ","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc21->Draw();
leg_1D[9]->Draw();

////////////////////////////////////sd3

TLatex *texc22 = new TLatex(18.35434,40.23506,"PF+PUPPI ");
   texc22->SetLineWidth(0.05);
   texc22->SetTextSize(0.04);

c_1D[10]->cd();
TMg_1D[10]->Add(puppi_graph_sdMeanvsNpv[0]);
TMg_1D[10]->Add(puppi_graph_sdMeanvsNpv[1]);
TMg_1D[10]->Add(puppi_graph_sdMeanvsNpv[2]);
//TMg_1D[10]->Add(puppi_graph_sdMeanvsNpv[3]);

TMg_1D[10]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[10]->Draw("AP");

leg_1D[10]->AddEntry(puppi_graph_sdMeanvsNpv[0],"Softdrop #beta = 2 ","l");
leg_1D[10]->AddEntry(puppi_graph_sdMeanvsNpv[1],"Softdrop #beta = 0 ","l");
leg_1D[10]->AddEntry(puppi_graph_sdMeanvsNpv[2],"Softdrop #beta = 1 ","l");
//leg_1D[10]->AddEntry(puppi_graph_sdMeanvsNpv[3],"Softdrop #beta = -1 ","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc22->Draw();
leg_1D[10]->Draw();



///////////////////////////////////////////////////////////////Mean vs nPV





/*
////////////////////////////tr1

TLatex *texc0 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.2,pT_{frac}=0.05)");
   texc0->SetLineWidth(0.05);
   texc0->SetTextSize(0.03);

c_1D[0]->cd();
TMg_1D[0]->Add(pf_graph_trMeanvsNpv[0]);
TMg_1D[0]->Add(chs_graph_trMeanvsNpv[0]);
TMg_1D[0]->Add(puppi_graph_trMeanvsNpv[0]);
TMg_1D[0]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[0]->Draw("AP");

leg_1D[0]->AddEntry(pf_graph_trMeanvsNpv[0],"PF  ","lp");
leg_1D[0]->AddEntry(chs_graph_trMeanvsNpv[0],"PF+CHS ","lp");
leg_1D[0]->AddEntry(puppi_graph_trMeanvsNpv[0],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
texc0->Draw();
leg_1D[0]->Draw();

////////////////////////////tr2

TLatex *texc1 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.1,pT_{frac}=0.03)");
   texc1->SetLineWidth(0.05);
   texc1->SetTextSize(0.03);

c_1D[1]->cd();
TMg_1D[1]->Add(pf_graph_trMeanvsNpv[1]);
TMg_1D[1]->Add(chs_graph_trMeanvsNpv[1]);
TMg_1D[1]->Add(puppi_graph_trMeanvsNpv[1]);
TMg_1D[1]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[1]->Draw("AP");

leg_1D[1]->AddEntry(pf_graph_trMeanvsNpv[1],"PF  ","lp");
leg_1D[1]->AddEntry(chs_graph_trMeanvsNpv[1],"PF+CHS ","lp");
leg_1D[1]->AddEntry(puppi_graph_trMeanvsNpv[1],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc1->Draw();
leg_1D[1]->Draw();


////////////////////////////////////////tr3


TLatex *texc2 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.2,pT_{frac}=0.03)");
   texc2->SetLineWidth(0.05);
   texc2->SetTextSize(0.03);

c_1D[2]->cd();
TMg_1D[2]->Add(pf_graph_trMeanvsNpv[2]);
TMg_1D[2]->Add(chs_graph_trMeanvsNpv[2]);
TMg_1D[2]->Add(puppi_graph_trMeanvsNpv[2]);
TMg_1D[2]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[2]->Draw("AP");

leg_1D[2]->AddEntry(pf_graph_trMeanvsNpv[2],"PF  ","lp");
leg_1D[2]->AddEntry(chs_graph_trMeanvsNpv[2],"PF+CHS ","lp");
leg_1D[2]->AddEntry(puppi_graph_trMeanvsNpv[2],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc2->Draw();
leg_1D[2]->Draw();

///////////////////////////////////////////////tr4

TLatex *texc3 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.3,pT_{frac}=0.03)");
   texc3->SetLineWidth(0.05);
   texc3->SetTextSize(0.03);

c_1D[3]->cd();
TMg_1D[3]->Add(pf_graph_trMeanvsNpv[3]);
TMg_1D[3]->Add(chs_graph_trMeanvsNpv[3]);
TMg_1D[3]->Add(puppi_graph_trMeanvsNpv[3]);
TMg_1D[3]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[3]->Draw("AP");

leg_1D[3]->AddEntry(pf_graph_trMeanvsNpv[3],"PF  ","lp");
leg_1D[3]->AddEntry(chs_graph_trMeanvsNpv[3],"PF+CHS ","lp");
leg_1D[3]->AddEntry(puppi_graph_trMeanvsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc3->Draw();
leg_1D[3]->Draw();
//////////////////////////////////////////////////pr1

TLatex *texc4 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.1,r_{cut}=0.5)");
   texc4->SetLineWidth(0.05);
   texc4->SetTextSize(0.03);

c_1D[4]->cd();
TMg_1D[4]->Add(pf_graph_prMeanvsNpv[0]);
TMg_1D[4]->Add(chs_graph_prMeanvsNpv[0]);
TMg_1D[4]->Add(puppi_graph_prMeanvsNpv[0]);
TMg_1D[4]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[4]->Draw("AP");

leg_1D[4]->AddEntry(pf_graph_prMeanvsNpv[0],"PF  ","lp");
leg_1D[4]->AddEntry(chs_graph_prMeanvsNpv[0],"PF+CHS ","lp");
leg_1D[4]->AddEntry(puppi_graph_prMeanvsNpv[0],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc4->Draw();
leg_1D[4]->Draw();

////////////////////////////////////////pr2
TLatex *texc5 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.05,r_{cut}=0.5)");
   texc5->SetLineWidth(0.05);
   texc5->SetTextSize(0.03);

c_1D[5]->cd();
TMg_1D[5]->Add(pf_graph_prMeanvsNpv[1]);
TMg_1D[5]->Add(chs_graph_prMeanvsNpv[1]);
TMg_1D[5]->Add(puppi_graph_prMeanvsNpv[1]);
TMg_1D[5]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[5]->Draw("AP");

leg_1D[5]->AddEntry(pf_graph_prMeanvsNpv[1],"PF  ","lp");
leg_1D[5]->AddEntry(chs_graph_prMeanvsNpv[1],"PF+CHS ","lp");
leg_1D[5]->AddEntry(puppi_graph_prMeanvsNpv[1],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc5->Draw();
leg_1D[5]->Draw();



///////////////////////////////////////////pr3

TLatex *texc6 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.05,r_{cut}=0.75)");
   texc6->SetLineWidth(0.05);
   texc6->SetTextSize(0.03);

c_1D[6]->cd();
TMg_1D[6]->Add(pf_graph_prMeanvsNpv[2]);
TMg_1D[6]->Add(chs_graph_prMeanvsNpv[2]);
TMg_1D[6]->Add(puppi_graph_prMeanvsNpv[2]);
TMg_1D[6]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[6]->Draw("AP");

leg_1D[6]->AddEntry(pf_graph_prMeanvsNpv[2],"PF  ","lp");
leg_1D[6]->AddEntry(chs_graph_prMeanvsNpv[2],"PF+CHS ","lp");
leg_1D[6]->AddEntry(puppi_graph_prMeanvsNpv[2],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc6->Draw();
leg_1D[6]->Draw();







//////////////////////////////////////////pr4

TLatex *texc7 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.05,r_{cut}=0.75)");
   texc7->SetLineWidth(0.05);
   texc7->SetTextSize(0.03);

c_1D[7]->cd();
TMg_1D[7]->Add(pf_graph_prMeanvsNpv[3]);
TMg_1D[7]->Add(chs_graph_prMeanvsNpv[3]);
TMg_1D[7]->Add(puppi_graph_prMeanvsNpv[3]);
TMg_1D[7]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[7]->Draw("AP");

leg_1D[7]->AddEntry(pf_graph_prMeanvsNpv[3],"PF  ","lp");
leg_1D[7]->AddEntry(chs_graph_prMeanvsNpv[3],"PF+CHS ","lp");
leg_1D[7]->AddEntry(puppi_graph_prMeanvsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc7->Draw();
leg_1D[7]->Draw();

///////////////////////////////sd1

TLatex *texc8 = new TLatex(20.45501,24.0062,"Softdrop #beta= 2");
   texc8->SetLineWidth(0.05);
   texc8->SetTextSize(0.03);

c_1D[8]->cd();
TMg_1D[8]->Add(pf_graph_sdMeanvsNpv[0]);
TMg_1D[8]->Add(chs_graph_sdMeanvsNpv[0]);
TMg_1D[8]->Add(puppi_graph_sdMeanvsNpv[0]);
TMg_1D[8]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[8]->Draw("AP");

leg_1D[8]->AddEntry(pf_graph_sdMeanvsNpv[0],"PF  ","lp");
leg_1D[8]->AddEntry(chs_graph_sdMeanvsNpv[0],"PF+CHS ","lp");
leg_1D[8]->AddEntry(puppi_graph_sdMeanvsNpv[0],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc8->Draw();
leg_1D[8]->Draw();

////////////////////////////////sd2

TLatex *texc9 = new TLatex(20.45501,24.0062,"Softdrop #beta= 0");
   texc9->SetLineWidth(0.05);
   texc9->SetTextSize(0.03);

c_1D[9]->cd();
TMg_1D[9]->Add(pf_graph_sdMeanvsNpv[1]);
TMg_1D[9]->Add(chs_graph_sdMeanvsNpv[1]);
TMg_1D[9]->Add(puppi_graph_sdMeanvsNpv[1]);
TMg_1D[9]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[9]->Draw("AP");

leg_1D[9]->AddEntry(pf_graph_sdMeanvsNpv[1],"PF  ","lp");
leg_1D[9]->AddEntry(chs_graph_sdMeanvsNpv[1],"PF+CHS ","lp");
leg_1D[9]->AddEntry(puppi_graph_sdMeanvsNpv[1],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc9->Draw();
leg_1D[9]->Draw();

////////////////////////////////////sd3

TLatex *texc10 = new TLatex(20.45501,24.0062,"Softdrop #beta= 1");
   texc10->SetLineWidth(0.05);
   texc10->SetTextSize(0.03);

c_1D[10]->cd();
TMg_1D[10]->Add(pf_graph_sdMeanvsNpv[2]);
TMg_1D[10]->Add(chs_graph_sdMeanvsNpv[2]);
TMg_1D[10]->Add(puppi_graph_sdMeanvsNpv[2]);
TMg_1D[10]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[10]->Draw("AP");

leg_1D[10]->AddEntry(pf_graph_sdMeanvsNpv[2],"PF  ","lp");
leg_1D[10]->AddEntry(chs_graph_sdMeanvsNpv[2],"PF+CHS ","lp");
leg_1D[10]->AddEntry(puppi_graph_sdMeanvsNpv[2],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc10->Draw();
leg_1D[10]->Draw();
////////////////////////////////////////////sd4
TLatex *texc11 = new TLatex(20.45501,24.0062,"Softdrop #beta= -1");
   texc11->SetLineWidth(0.05);
   texc11->SetTextSize(0.03);

c_1D[11]->cd();
TMg_1D[11]->Add(pf_graph_sdMeanvsNpv[3]);
TMg_1D[11]->Add(chs_graph_sdMeanvsNpv[3]);
TMg_1D[11]->Add(puppi_graph_sdMeanvsNpv[3]);
TMg_1D[11]->SetTitle(" ;n_{PV} ;Mean[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[11]->Draw("AP");

leg_1D[11]->AddEntry(pf_graph_sdMeanvsNpv[3],"PF  ","lp");
leg_1D[11]->AddEntry(chs_graph_sdMeanvsNpv[3],"PF+CHS ","lp");
leg_1D[11]->AddEntry(puppi_graph_sdMeanvsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc11->Draw();
leg_1D[11]->Draw();

*/

/////////////////////////////////////////////////RMS vs nPV plots
////////////////////////////tr1

/*
TLatex *tex1 = new TLatex(23.53793,60.84428,"CMS Simulation Preliminary, #sqrt{s}= 13 TeV");
   tex1->SetLineWidth(0.3);
   tex1->SetTextSize(0.04);
   
 TLatex *tex2 = new TLatex(25.07888,55.98329,"QCD, Anti-kT (R=0.8)");
   tex2->SetLineWidth(0.02);
   tex2->SetTextSize(0.03);

TLatex *tex3 = new TLatex(25.07888,52.78526,"   p_{T} > 300 GeV");
   tex3->SetLineWidth(0.02);
   tex3->SetTextSize(0.03);

TLatex *tex4 = new TLatex(27.33103,49.971," |#eta|< 2.5");
   tex4->SetLineWidth(0.02);
   tex4->SetTextSize(0.03);
*/



TLatex *texc12 = new TLatex(14.69436,36.36168,"PF with trimming");
   texc12->SetLineWidth(0.05);
   texc12->SetTextSize(0.04);
   texc12->SetTextFont(42);
c_1D[12]->cd();
TMg_1D[12]->Add(pf_graph_trSigmavsNpv[0]);
TMg_1D[12]->Add(pf_graph_trSigmavsNpv[1]);
TMg_1D[12]->Add(pf_graph_trSigmavsNpv[2]);
TMg_1D[12]->Add(pf_graph_trSigmavsNpv[3]);

TMg_1D[12]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[12]->Draw("AP");

leg_1D[12]->AddEntry(pf_graph_trSigmavsNpv[0],"r_{sub}=0.2,pT_{frac}=0.05","l");
leg_1D[12]->AddEntry(pf_graph_trSigmavsNpv[1],"r_{sub}=0.1,pT_{frac}=0.03","l");
leg_1D[12]->AddEntry(pf_graph_trSigmavsNpv[2],"r_{sub}=0.2,pT_{frac}=0.03","l");
leg_1D[12]->AddEntry(pf_graph_trSigmavsNpv[3],"r_{sub}=0.3,pT_{frac}=0.03","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();


texc12->Draw();
leg_1D[12]->Draw();

////////////////////////////tr2

TLatex *texc13 = new TLatex(14.69436,36.36168,"PF+CHS with trimming");
   texc13->SetLineWidth(0.05);
   texc13->SetTextSize(0.04);
   texc13->SetTextFont(42);

c_1D[13]->cd();
TMg_1D[13]->Add(chs_graph_trSigmavsNpv[0]);
TMg_1D[13]->Add(chs_graph_trSigmavsNpv[1]);
TMg_1D[13]->Add(chs_graph_trSigmavsNpv[2]);
TMg_1D[13]->Add(chs_graph_trSigmavsNpv[3]);
TMg_1D[13]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[13]->Draw("AP");

leg_1D[13]->AddEntry(chs_graph_trSigmavsNpv[0],"r_{sub}=0.2,pT_{frac}=0.05","l");
leg_1D[13]->AddEntry(chs_graph_trSigmavsNpv[1],"r_{sub}=0.1,pT_{frac}=0.03","l");
leg_1D[13]->AddEntry(chs_graph_trSigmavsNpv[2],"r_{sub}=0.2,pT_{frac}=0.03","l");
leg_1D[13]->AddEntry(chs_graph_trSigmavsNpv[3],"r_{sub}=0.3,pT_{frac}=0.03","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();
texc13->Draw();
leg_1D[13]->Draw();


////////////////////////////////////////tr3


TLatex *texc14 = new TLatex(14.69436,36.36168,"PF+PUPPI with trimming");
   texc14->SetLineWidth(0.05);
   texc14->SetTextSize(0.04);
   texc14->SetTextFont(42);
c_1D[14]->cd();
TMg_1D[14]->Add(puppi_graph_trSigmavsNpv[0]);
TMg_1D[14]->Add(puppi_graph_trSigmavsNpv[1]);
TMg_1D[14]->Add(puppi_graph_trSigmavsNpv[2]);
TMg_1D[14]->Add(puppi_graph_trSigmavsNpv[3]);

TMg_1D[14]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[14]->Draw("AP");

leg_1D[14]->AddEntry(puppi_graph_trSigmavsNpv[0],"r_{sub}=0.2,pT_{frac}=0.05","l");
leg_1D[14]->AddEntry(puppi_graph_trSigmavsNpv[1],"r_{sub}=0.1,pT_{frac}=0.03","l");
leg_1D[14]->AddEntry(puppi_graph_trSigmavsNpv[2],"r_{sub}=0.2,pT_{frac}=0.03","l");
leg_1D[14]->AddEntry(puppi_graph_trSigmavsNpv[3],"r_{sub}=0.3,pT_{frac}=0.03","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texc14->Draw();
leg_1D[14]->Draw();

///////////////////////////////////////////////tr4
/*
TLatex *texc15 = new TLatex(20.45501,24.0062,"Trimmed (r_{sub}=0.3,pT_{frac}=0.03)");
   texc15->SetLineWidth(0.05);
   texc15->SetTextSize(0.03);

c_1D[15]->cd();
TMg_1D[15]->Add(pf_graph_trSigmavsNpv[3]);
TMg_1D[15]->Add(chs_graph_trSigmavsNpv[3]);
TMg_1D[15]->Add(puppi_graph_trSigmavsNpv[3]);
TMg_1D[15]->SetTitle(" ;n_{PV} ;RMS[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[15]->Draw("AP");

leg_1D[15]->AddEntry(pf_graph_trSigmavsNpv[3],"PF  ","lp");
leg_1D[15]->AddEntry(chs_graph_trSigmavsNpv[3],"PF+CHS ","lp");
leg_1D[15]->AddEntry(puppi_graph_trSigmavsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc15->Draw();
leg_1D[15]->Draw();
*/
//////////////////////////////////////////////////pr1

TLatex *texc16 = new TLatex(14.69436,36.36168,"PF with pruning");
   texc16->SetLineWidth(0.05);
   texc16->SetTextSize(0.04);
   texc16->SetTextFont(42);

c_1D[16]->cd();
TMg_1D[16]->Add(pf_graph_prSigmavsNpv[0]);
TMg_1D[16]->Add(pf_graph_prSigmavsNpv[1]);
TMg_1D[16]->Add(pf_graph_prSigmavsNpv[2]);
TMg_1D[16]->Add(pf_graph_prSigmavsNpv[3]);

TMg_1D[16]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[16]->Draw("AP");

leg_1D[16]->AddEntry(pf_graph_prSigmavsNpv[0],"z_{cut}=0.1,r_{cut}=0.5","l");
leg_1D[16]->AddEntry(pf_graph_prSigmavsNpv[1],"z_{cut}=0.05,r_{cut}=0.5","l");
leg_1D[16]->AddEntry(pf_graph_prSigmavsNpv[2],"z_{cut}=0.05,r_{cut}=0.75","l");
leg_1D[16]->AddEntry(pf_graph_prSigmavsNpv[3],"z_{cut}=0.1,r_{cut}=0.75","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();
texc16->Draw();
leg_1D[16]->Draw();

////////////////////////////////////////pr2
TLatex *texc17 = new TLatex(14.69436,36.36168,"PF+CHS with pruning");
   texc17->SetLineWidth(0.05);
   texc17->SetTextSize(0.04);
   texc17->SetTextFont(42);
c_1D[17]->cd();
TMg_1D[17]->Add(chs_graph_prSigmavsNpv[0]);
TMg_1D[17]->Add(chs_graph_prSigmavsNpv[1]);
TMg_1D[17]->Add(chs_graph_prSigmavsNpv[2]);
TMg_1D[17]->Add(chs_graph_prSigmavsNpv[3]);

TMg_1D[17]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[17]->Draw("AP");

leg_1D[17]->AddEntry(chs_graph_prSigmavsNpv[0],"z_{cut}=0.1,r_{cut}=0.5","l");
leg_1D[17]->AddEntry(chs_graph_prSigmavsNpv[1],"z_{cut}=0.05,r_{cut}=0.5","l");
leg_1D[17]->AddEntry(chs_graph_prSigmavsNpv[2],"z_{cut}=0.05,r_{cut}=0.75","l");
leg_1D[17]->AddEntry(chs_graph_prSigmavsNpv[3],"z_{cut}=0.1,r_{cut}=0.75","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texc17->Draw();
leg_1D[17]->Draw();



///////////////////////////////////////////pr3

TLatex *texc18 = new TLatex(14.69436,36.36168,"PF+PUPPI with pruning");
   texc18->SetLineWidth(0.05);
   texc18->SetTextSize(0.04);
   texc18->SetTextFont(42);
c_1D[18]->cd();
TMg_1D[18]->Add(puppi_graph_prSigmavsNpv[0]);
TMg_1D[18]->Add(puppi_graph_prSigmavsNpv[1]);
TMg_1D[18]->Add(puppi_graph_prSigmavsNpv[2]);
TMg_1D[18]->Add(puppi_graph_prSigmavsNpv[3]);

TMg_1D[18]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[18]->Draw("AP");

leg_1D[18]->AddEntry(puppi_graph_prSigmavsNpv[0],"z_{cut}=0.1,r_{cut}=0.5","l");
leg_1D[18]->AddEntry(puppi_graph_prSigmavsNpv[1],"z_{cut}=0.05,r_{cut}=0.5","l");
leg_1D[18]->AddEntry(puppi_graph_prSigmavsNpv[2],"z_{cut}=0.05,r_{cut}=0.75","l");
leg_1D[18]->AddEntry(puppi_graph_prSigmavsNpv[3],"z_{cut}=0.1,r_{cut}=0.75","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();
texc18->Draw();
leg_1D[18]->Draw();







//////////////////////////////////////////pr4
/*
TLatex *texc19 = new TLatex(20.45501,24.0062,"Pruned (z_{cut}=0.05,r_{cut}=0.75)");
   texc19->SetLineWidth(0.05);
   texc19->SetTextSize(0.03);

c_1D[19]->cd();
TMg_1D[19]->Add(pf_graph_prSigmavsNpv[3]);
TMg_1D[19]->Add(chs_graph_prSigmavsNpv[3]);
TMg_1D[19]->Add(puppi_graph_prSigmavsNpv[3]);
TMg_1D[19]->SetTitle(" ;n_{PV} ;RMS[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[19]->Draw("AP");

leg_1D[19]->AddEntry(pf_graph_prSigmavsNpv[3],"PF  ","lp");
leg_1D[19]->AddEntry(chs_graph_prSigmavsNpv[3],"PF+CHS ","lp");
leg_1D[19]->AddEntry(puppi_graph_prSigmavsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc19->Draw();
leg_1D[19]->Draw();
*/
///////////////////////////////sd1

TLatex *texc20 = new TLatex(14.69436,36.36168,"PF with softdrop");
   texc20->SetLineWidth(0.05);
   texc20->SetTextSize(0.04);
   texc20->SetTextFont(42);
c_1D[20]->cd();
TMg_1D[20]->Add(pf_graph_sdSigmavsNpv[0]);
TMg_1D[20]->Add(pf_graph_sdSigmavsNpv[1]);
TMg_1D[20]->Add(pf_graph_sdSigmavsNpv[2]);
//TMg_1D[20]->Add(pf_graph_sdSigmavsNpv[3]);

TMg_1D[20]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[20]->Draw("AP");

leg_1D[20]->AddEntry(pf_graph_sdSigmavsNpv[0],"#beta = 2 ","l");
leg_1D[20]->AddEntry(pf_graph_sdSigmavsNpv[1],"#beta = 0 ","l");
leg_1D[20]->AddEntry(pf_graph_sdSigmavsNpv[2],"#beta = 1 ","l");
//leg_1D[20]->AddEntry(pf_graph_sdSigmavsNpv[3],"Softdrop #beta = -1 ","l");

tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texc20->Draw();
leg_1D[20]->Draw();

////////////////////////////////sd2

TLatex *texc21 = new TLatex(14.69436,36.36168,"PF+CHS with softdrop");
   texc21->SetLineWidth(0.05);
   texc21->SetTextSize(0.04);
   texc21->SetTextFont(42);
c_1D[21]->cd();
TMg_1D[21]->Add(chs_graph_sdSigmavsNpv[0]);
TMg_1D[21]->Add(chs_graph_sdSigmavsNpv[1]);
TMg_1D[21]->Add(chs_graph_sdSigmavsNpv[2]);
//TMg_1D[21]->Add(chs_graph_sdSigmavsNpv[3]);

TMg_1D[21]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[21]->Draw("AP");

leg_1D[21]->AddEntry(chs_graph_sdSigmavsNpv[0],"#beta = 2 ","l");
leg_1D[21]->AddEntry(chs_graph_sdSigmavsNpv[1],"#beta = 0 ","l");
leg_1D[21]->AddEntry(chs_graph_sdSigmavsNpv[2],"#beta = 1 ","l");
//leg_1D[21]->AddEntry(chs_graph_sdSigmavsNpv[3],"Softdrop #beta = -1 ","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texc21->Draw();
leg_1D[21]->Draw();

////////////////////////////////////sd3

TLatex *texc22 = new TLatex(14.69436,36.36168,"PF+PUPPI with softdrop");
   texc22->SetLineWidth(0.05);
   texc22->SetTextSize(0.04);
   texc22->SetTextFont(42);
c_1D[22]->cd();
TMg_1D[22]->Add(puppi_graph_sdSigmavsNpv[0]);
TMg_1D[22]->Add(puppi_graph_sdSigmavsNpv[1]);
TMg_1D[22]->Add(puppi_graph_sdSigmavsNpv[2]);
//TMg_1D[22]->Add(puppi_graph_sdSigmavsNpv[3]);

TMg_1D[22]->SetTitle(" ;n_{PV} ;Resolution(m_{reco}-m_{gen}) [GeV]");
TMg_1D[22]->Draw("AP");

leg_1D[22]->AddEntry(puppi_graph_sdSigmavsNpv[0],"#beta = 2 ","l");
leg_1D[22]->AddEntry(puppi_graph_sdSigmavsNpv[1],"#beta = 0 ","l");
leg_1D[22]->AddEntry(puppi_graph_sdSigmavsNpv[2],"#beta = 1 ","l");
//leg_1D[22]->AddEntry(puppi_graph_sdSigmavsNpv[3],"Softdrop #beta = -1 ","l");


tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();
texc22->Draw();
leg_1D[22]->Draw();

char name[100];
char namepng[100];

for(int j1=0;j1<12;j1++){
sprintf(name,"Mean_Vs_nPV_15to50_%ig.pdf",j1);
cout<<"name = "<<name<<endl;
//c_1D[j1]->SaveAs(name);
c_1D[j1]->Close();

}

for(int j1=12;j1<24;j1++){
sprintf(name,"Resolution_Vs_nPV_15to50_%ig.pdf",j1);
sprintf(namepng,"Resolution_Vs_nPV_15to50_%ig.png",j1);

cout<<"name = "<<name<<endl;
c_1D[j1]->SaveAs(name);
c_1D[j1]->SaveAs(namepng);
//c_1D[j1]->Close();

}







////////////////////////////////////////////sd4
/*
TLatex *texc23 = new TLatex(20.45501,24.0062,"Softdrop #beta= -1");
   texc23->SetLineWidth(0.05);
   texc23->SetTextSize(0.03);

c_1D[23]->cd();
TMg_1D[23]->Add(pf_graph_sdSigmavsNpv[3]);
TMg_1D[23]->Add(chs_graph_sdSigmavsNpv[3]);
TMg_1D[23]->Add(puppi_graph_sdSigmavsNpv[3]);
TMg_1D[23]->SetTitle(" ;n_{PV} ;RMS[m_{RecoJet}-m_{GenJet}] (GeV) ");
TMg_1D[23]->Draw("AP");

leg_1D[23]->AddEntry(pf_graph_sdSigmavsNpv[3],"PF  ","lp");
leg_1D[23]->AddEntry(chs_graph_sdSigmavsNpv[3],"PF+CHS ","lp");
leg_1D[23]->AddEntry(puppi_graph_sdSigmavsNpv[3],"PF+PUPPI ","lp");

tex1->Draw();
tex2->Draw();
tex3->Draw();
texc23->Draw();
leg_1D[23]->Draw();

*/














































}//main

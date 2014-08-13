#include<iostream>
#include <TH1.h>
#include "TStyle.h"
#include <TCanvas.h>
#include<Riostream.h>
#include "TROOT.h"
#include<TGraph.h>
#include "TGraphErrors.h"
using namespace std;

void 1DmassCode(){

gStyle->SetOptStat(0000);

int maxEvents=40000;

double pTmax=10000;
double pTmin=300;

const unsigned int nPVmin=20;
const unsigned int nPVmax=60;
const unsigned double labelsize=0.05;
const unsigned double titlesize=0.05;
const unsigned double titleOffset=1.5;
const unsigned int nPVnumOfbins=20;//should be a divisible number of (nPVmax-nPVmin)

const unsigned int NumtrPar=4;
const unsigned int NumprPar=4;
const unsigned int NumsdPar=4;

 char histnamePF[100];
 char histnameCHS[100];
 char histnamePUPPI[100];

  

  double rangeMin=0;
  double rangeMax=150;
 
  int NoOfbins=50;

  char *nstTr[4];
  char *nstPr[4];
  char *nstSd[4];

  nstTr[0]="tr1";
  nstTr[1]="tr2";
  nstTr[2]="tr3";
  nstTr[3]="tr4";

  nstPr[0]="pr1";
  nstPr[1]="pr2";
  nstPr[2]="pr3";
  nstPr[3]="pr4";

  nstSd[0]="sd1";
  nstSd[1]="sd2";
  nstSd[2]="sd3";
  nstSd[3]="sd4";


TH1F *hpf=new TH1F("hpf","hpf",50,0,150);
        hpf->SetLineColor(6);
        hpf->SetLineWidth(2);
        hpf->GetXaxis()->SetTitle("m_{jet}");
        hpf->GetXaxis()->SetTitleSize(titlesize);
        hpf->GetXaxis()->SetLabelSize(labelsize);
        hpf->GetYaxis()->SetTitle("Normalized to unity");
        hpf->GetYaxis()->SetTitleSize(titlesize);
        hpf->GetYaxis()->SetLabelSize(labelsize);
        hpf->GetYaxis()->SetTitleOffset(titleOffset);


TH1F *hchs=new TH1F("hchs","hchs",50,0,150);

        hchs->SetLineColor(6);
        hchs->SetLineWidth(2);
        hchs->GetXaxis()->SetTitle("m_{jet}");
        hchs->GetXaxis()->SetTitleSize(titlesize);
        hchs->GetXaxis()->SetLabelSize(labelsize);
        hchs->GetYaxis()->SetTitle("Normalized to unity");
        hchs->GetYaxis()->SetTitleSize(titlesize);
        hchs->GetYaxis()->SetLabelSize(labelsize);
        hchs->GetYaxis()->SetTitleOffset(titleOffset);




  TH1F *JetRsAKTrPF_h[NumtrPar];
  TH1F *JetRsAKTrCHS_h[NumtrPar];
  TH1F *JetRsAKTrPUPPI_h[NumtrPar];

  TH1F *JetRsAKPrPF_h[NumprPar];
  TH1F *JetRsAKPrCHS_h[NumprPar];
  TH1F *JetRsAKPrPUPPI_h[NumprPar];

  TH1F *JetRsAKSdPF_h[NumsdPar];
  TH1F *JetRsAKSdCHS_h[NumsdPar];
  TH1F *JetRsAKSdPUPPI_h[NumsdPar];


  for (int jj=0;jj<NumtrPar;jj++) {
  
        sprintf(histnamePF,"ResponseAK8PF_Tr_%i",jj);
        JetRsAKTrPF_h[jj] = new TH1F(histnamePF,histnamePF,NoOfbins,rangeMin,rangeMax);
        JetRsAKTrPF_h[jj]->SetLineColor(jj+1);
        JetRsAKTrPF_h[jj]->SetLineWidth(4);
        JetRsAKTrPF_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKTrPF_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKTrPF_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKTrPF_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKTrPF_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKTrPF_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKTrPF_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


        sprintf(histnameCHS,"ResponseAK8CHS_Tr_%i",jj);
        JetRsAKTrCHS_h[jj] = new TH1F(histnameCHS,histnameCHS,NoOfbins,rangeMin,rangeMax);
        JetRsAKTrCHS_h[jj]->SetLineColor(jj+2);
        JetRsAKTrCHS_h[jj]->SetLineWidth(4);
        JetRsAKTrCHS_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKTrCHS_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKTrCHS_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKTrCHS_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKTrCHS_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKTrCHS_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKTrCHS_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


        sprintf(histnamePUPPI,"ResponseAK8PUPPI_Tr_%i",jj);
        JetRsAKTrPUPPI_h[jj] = new TH1F(histnamePUPPI,histnamePUPPI,NoOfbins,rangeMin,rangeMax);
        JetRsAKTrPUPPI_h[jj]->SetLineColor(jj+3);
        JetRsAKTrPUPPI_h[jj]->SetLineWidth(4);
        JetRsAKTrPUPPI_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKTrPUPPI_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKTrPUPPI_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKTrPUPPI_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKTrPUPPI_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKTrPUPPI_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKTrPUPPI_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


     }


 for (int jj=0;jj<NumprPar;jj++) {
  
        sprintf(histnamePF,"ResponseAK8PF_Pr_%i",jj);
        JetRsAKPrPF_h[jj] = new TH1F(histnamePF,histnamePF,NoOfbins,rangeMin,rangeMax);
        JetRsAKPrPF_h[jj]->SetLineColor(jj+1);
        JetRsAKPrPF_h[jj]->SetLineWidth(4);
        JetRsAKPrPF_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKPrPF_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKPrPF_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKPrPF_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKPrPF_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKPrPF_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKPrPF_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


        sprintf(histnameCHS,"ResponseAK8CHS_Pr_%i",jj);
        JetRsAKPrCHS_h[jj] = new TH1F(histnameCHS,histnameCHS,NoOfbins,rangeMin,rangeMax);
        JetRsAKPrCHS_h[jj]->SetLineColor(jj+2);
        JetRsAKPrCHS_h[jj]->SetLineWidth(4);
        JetRsAKPrCHS_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKPrCHS_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKPrCHS_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKPrCHS_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKPrCHS_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKPrCHS_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKPrCHS_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


        sprintf(histnamePUPPI,"ResponseAK8PUPPI_Pr_%i",jj);
        JetRsAKPrPUPPI_h[jj] = new TH1F(histnamePUPPI,histnamePUPPI,NoOfbins,rangeMin,rangeMax);
        JetRsAKPrPUPPI_h[jj]->SetLineColor(jj+3);
        JetRsAKPrPUPPI_h[jj]->SetLineWidth(4);
        JetRsAKPrPUPPI_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKPrPUPPI_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKPrPUPPI_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKPrPUPPI_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKPrPUPPI_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKPrPUPPI_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKPrPUPPI_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);
     }

      for (int jj=0;jj<NumsdPar;jj++) {
  
        sprintf(histnamePF,"ResponseAK8PF_Sd_%i",jj);
        JetRsAKSdPF_h[jj] = new TH1F(histnamePF,histnamePF,NoOfbins,rangeMin,rangeMax);
        JetRsAKSdPF_h[jj]->SetLineColor(jj+1);
        JetRsAKSdPF_h[jj]->SetLineWidth(4);
        JetRsAKSdPF_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKSdPF_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKSdPF_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKSdPF_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKSdPF_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKSdPF_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKSdPF_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);

        sprintf(histnameCHS,"ResponseAK8CHS_Sd_%i",jj);
        JetRsAKSdCHS_h[jj] = new TH1F(histnameCHS,histnameCHS,NoOfbins,rangeMin,rangeMax);
        JetRsAKSdCHS_h[jj]->SetLineColor(jj+2);
        JetRsAKSdCHS_h[jj]->SetLineWidth(4);
        JetRsAKSdCHS_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKSdCHS_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKSdCHS_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKSdCHS_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKSdCHS_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKSdCHS_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKSdCHS_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);

        sprintf(histnamePUPPI,"ResponseAK8PUPPI_Sd_%i",jj);
        JetRsAKSdPUPPI_h[jj] = new TH1F(histnamePUPPI,histnamePUPPI,NoOfbins,rangeMin,rangeMax);
        JetRsAKSdPUPPI_h[jj]->SetLineColor(jj+3);
        JetRsAKSdPUPPI_h[jj]->SetLineWidth(4);
        JetRsAKSdPUPPI_h[jj]->GetXaxis()->SetTitle("m_{jet} (GeV)");
        JetRsAKSdPUPPI_h[jj]->GetXaxis()->SetTitleSize(titlesize);
        JetRsAKSdPUPPI_h[jj]->GetXaxis()->SetLabelSize(labelsize);
        JetRsAKSdPUPPI_h[jj]->GetYaxis()->SetTitle("Normalized to unity");
        JetRsAKSdPUPPI_h[jj]->GetYaxis()->SetTitleSize(titlesize);
        JetRsAKSdPUPPI_h[jj]->GetYaxis()->SetLabelSize(labelsize);
        JetRsAKSdPUPPI_h[jj]->GetYaxis()->SetTitleOffset(titleOffset);


       }



////////////////////////////////set histogram color
for(int i1=0;i1<4;i1++){
JetRsAKTrPF_h[i1]->SetLineColor(1);
JetRsAKPrPF_h[i1]->SetLineColor(1);
JetRsAKSdPF_h[i1]->SetLineColor(1);


JetRsAKTrCHS_h[i1]->SetLineColor(2);
JetRsAKPrCHS_h[i1]->SetLineColor(2);
JetRsAKSdCHS_h[i1]->SetLineColor(2);


JetRsAKTrPUPPI_h[i1]->SetLineColor(3);
JetRsAKPrPUPPI_h[i1]->SetLineColor(3);
JetRsAKSdPUPPI_h[i1]->SetLineColor(3);


}








TFile *f1 = new TFile("QCD300to470csa14Puppiv3.root","READ","My root file1");
    TTree *trpf = (TTree*)f1->Get("chs");
    TTree *trchs = (TTree*)f1->Get("pf");
    TTree *trpuppi = (TTree*)f1->Get("puppi");
    TTree *trsk = (TTree*)f1->Get("softkiller");
    TTree *trgen = (TTree*)f1->Get("gen");






   int   npu;
  //////////////////////////////////////PF Tree
   trpf->SetBranchAddress("npv", &npu);
   vector<float>   *pf_m;
   vector<float>   *pf_ptcorrphil;

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
trpf->SetBranchAddress("m",&pf_m);
trpf->SetBranchAddress("ptcorrphil",&pf_ptcorrphil);
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
   vector<float>   *chs_m;
   vector<float>   *chs_ptcorrphil;
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

trchs->SetBranchAddress("m",&chs_m);
trchs->SetBranchAddress("ptcorrphil",&chs_ptcorrphil);
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


///////////////////////////////gen

int n = (int) trpf->GetEntries();
    

    cout<<"Total Entries = "<<n<<endl;
    cout<<"Code will run on "<<maxEvents<<" Entries"<<endl;
    for (int i = 0; i < n; ++i) {//loop over entries 

     if(i==maxEvents){
     break;
                     }

    trpf->GetEntry(i);
    trchs->GetEntry(i);
    trpuppi->GetEntry(i);
    trsk->GetEntry(i);
    trgen->GetEntry(i);



if(pf_ptcorrphil->at(0) > pTmin && pf_ptcorrphil->at(0) < pTmax ){
hpf->Fill(pf_m->at(0));

}

if(chs_ptcorrphil->at(0) > pTmin && chs_ptcorrphil->at(0) < pTmax ){
hchs->Fill(chs_m->at(0));

}






// PF
   
      if(pf_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin  && pf_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax ){
      JetRsAKTrPF_h[0]->Fill(pf_mtrimsafe_Rtrim_020_Ptfrac_005->at(0));
          }

if(pf_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) > pTmin && pf_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrPF_h[1]->Fill(pf_mtrimsafe_Rtrim_010_Ptfrac_003->at(0));
 }

if(pf_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) > pTmin && pf_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrPF_h[2]->Fill(pf_mtrimsafe_Rtrim_020_Ptfrac_003->at(0));
 }

if(pf_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) > pTmin &&  pf_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrPF_h[3]->Fill(pf_mtrimsafe_Rtrim_030_Ptfrac_003->at(0));
 }
//

if(pf_ptprunedsafe_zcut_010_R_cut_050->at(0) > pTmin &&  pf_ptprunedsafe_zcut_010_R_cut_050->at(0) < pTmax ){
       JetRsAKPrPF_h[0]->Fill(pf_mprunedsafe_zcut_010_R_cut_050->at(0));
 }

if(pf_ptprunedsafe_zcut_005_R_cut_050->at(0) > pTmin &&  pf_ptprunedsafe_zcut_005_R_cut_050->at(0) < pTmax ){
       JetRsAKPrPF_h[1]->Fill(pf_mprunedsafe_zcut_005_R_cut_050->at(0));
 }
if(pf_ptprunedsafe_zcut_005_R_cut_075->at(0) > pTmin  && pf_ptprunedsafe_zcut_005_R_cut_075->at(0) < pTmax ){
       JetRsAKPrPF_h[2]->Fill(pf_mprunedsafe_zcut_005_R_cut_075->at(0));
 }
if(pf_ptprunedsafe_zcut_010_R_cut_075->at(0) > pTmin && pf_ptprunedsafe_zcut_010_R_cut_075->at(0) < pTmax ){
       JetRsAKPrPF_h[3]->Fill(pf_mprunedsafe_zcut_010_R_cut_075->at(0));
 }

//

if(pf_ptsoftdropsafe_beta20->at(0) > pTmin && pf_ptsoftdropsafe_beta20->at(0) < pTmax ){
       JetRsAKSdPF_h[0]->Fill(pf_msoftdropsafe_beta20->at(0));
 }
if(pf_ptsoftdropsafe_beta00->at(0) > pTmin && pf_ptsoftdropsafe_beta00->at(0) < pTmax ){
       JetRsAKSdPF_h[1]->Fill(pf_msoftdropsafe_beta00->at(0));
 }
if(pf_ptsoftdropsafe_beta10->at(0) > pTmin  && pf_ptsoftdropsafe_beta10->at(0) < pTmax ){
       JetRsAKSdPF_h[2]->Fill(pf_msoftdropsafe_beta10->at(0));
 }
if(pf_ptsoftdropsafe_betam1->at(0) > pTmin && pf_ptsoftdropsafe_betam1->at(0) < pTmax ){
       JetRsAKSdPF_h[3]->Fill(pf_msoftdropsafe_betam1->at(0));
 }

//PF+CHS
if(chs_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin &&  chs_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax ){
      JetRsAKTrCHS_h[0]->Fill(chs_mtrimsafe_Rtrim_020_Ptfrac_005->at(0));
 }

if(chs_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) > pTmin &&  chs_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrCHS_h[1]->Fill(chs_mtrimsafe_Rtrim_010_Ptfrac_003->at(0));
 }

if(chs_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) > pTmin && chs_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrCHS_h[2]->Fill(chs_mtrimsafe_Rtrim_020_Ptfrac_003->at(0));
 }

if(chs_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) > pTmin &&  chs_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) < pTmax ){
      JetRsAKTrCHS_h[3]->Fill(chs_mtrimsafe_Rtrim_030_Ptfrac_003->at(0));
 }
//

if(chs_ptprunedsafe_zcut_010_R_cut_050->at(0) > pTmin &&  chs_ptprunedsafe_zcut_010_R_cut_050->at(0) < pTmax ){
       JetRsAKPrCHS_h[0]->Fill(chs_mprunedsafe_zcut_010_R_cut_050->at(0));
 }
if(chs_ptprunedsafe_zcut_005_R_cut_050->at(0) > pTmin  && chs_ptprunedsafe_zcut_005_R_cut_050->at(0) < pTmax ){
       JetRsAKPrCHS_h[1]->Fill(chs_mprunedsafe_zcut_005_R_cut_050->at(0));
 }
if(chs_ptprunedsafe_zcut_005_R_cut_075->at(0) > pTmin  && chs_ptprunedsafe_zcut_005_R_cut_075->at(0) < pTmax ){
       JetRsAKPrCHS_h[2]->Fill(chs_mprunedsafe_zcut_005_R_cut_075->at(0));
 }
if(chs_ptprunedsafe_zcut_010_R_cut_075->at(0) > pTmin  && chs_ptprunedsafe_zcut_010_R_cut_075->at(0) < pTmax ){
       JetRsAKPrCHS_h[3]->Fill(chs_mprunedsafe_zcut_010_R_cut_075->at(0));
 }

// SD

if(chs_ptsoftdropsafe_beta20->at(0) > pTmin  && chs_ptsoftdropsafe_beta20->at(0) < pTmax ){
     JetRsAKSdCHS_h[0]->Fill(chs_msoftdropsafe_beta20->at(0));
 }
if(chs_ptsoftdropsafe_beta00->at(0) > pTmin  && chs_ptsoftdropsafe_beta00->at(0) < pTmax ){
       JetRsAKSdCHS_h[1]->Fill(chs_msoftdropsafe_beta00->at(0));
}
if(chs_ptsoftdropsafe_beta10->at(0) > pTmin  && chs_ptsoftdropsafe_beta10->at(0) ){
       JetRsAKSdCHS_h[2]->Fill(chs_msoftdropsafe_beta10->at(0));
}
if(chs_ptsoftdropsafe_betam1->at(0) > pTmin  && chs_ptsoftdropsafe_betam1->at(0) < pTmax ){
       JetRsAKSdCHS_h[3]->Fill(chs_msoftdropsafe_betam1->at(0));
}


// Puppi

/*

if(puppi_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) > pTmin &&  puppi_pttrimsafe_Rtrim_020_Ptfrac_005->at(0) < pTmax ){
      JetRsAKTrPUPPI_h[0]->Fill(puppi_mtrim_Rtrim_020_Ptfrac_005->at(0)-gen_mtrim_Rtrim_020_Ptfrac_005->at(0));
 }

if(puppi_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) > pTmin && puppi_pttrimsafe_Rtrim_010_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_010_Ptfrac_003->at(0) < pTmax){
      JetRsAKTrPUPPI_h[1]->Fill(puppi_mtrim_Rtrim_010_Ptfrac_003->at(0)-gen_mtrim_Rtrim_010_Ptfrac_003->at(0));
 }


if(puppi_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) > pTmin && puppi_pttrimsafe_Rtrim_020_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_020_Ptfrac_003->at(0) < pTmax){
      JetRsAKTrPUPPI_h[2]->Fill(puppi_mtrim_Rtrim_020_Ptfrac_003->at(0)-gen_mtrim_Rtrim_020_Ptfrac_003->at(0));
 }
if(puppi_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) > pTmin && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) > pTmin && puppi_pttrimsafe_Rtrim_030_Ptfrac_003->at(0) < pTmax && gen_pttrim_Rtrim_030_Ptfrac_003->at(0) < pTmax){
      JetRsAKTrPUPPI_h[3]->Fill(puppi_mtrim_Rtrim_030_Ptfrac_003->at(0)-gen_mtrim_Rtrim_030_Ptfrac_003->at(0));
 }


if(puppi_ptprunedsafe_zcut_010_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_050->at(0) > pTmin && puppi_ptprunedsafe_zcut_010_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_050->at(0) < pTmax){
       JetRsAKPrPUPPI_h[0]->Fill(puppi_mpruned_zcut_010_R_cut_050->at(0)-gen_mpruned_zcut_010_R_cut_050->at(0));
 }
if(puppi_ptprunedsafe_zcut_005_R_cut_050->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_050->at(0) > pTmin && puppi_ptprunedsafe_zcut_005_R_cut_050->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_050->at(0) < pTmax){
       JetRsAKPrPUPPI_h[1]->Fill(puppi_mpruned_zcut_005_R_cut_050->at(0)-gen_mpruned_zcut_005_R_cut_050->at(0));
 }
if(puppi_ptprunedsafe_zcut_005_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_005_R_cut_075->at(0) > pTmin && puppi_ptprunedsafe_zcut_005_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_005_R_cut_075->at(0) < pTmax){
       JetRsAKPrPUPPI_h[2]->Fill(puppi_mpruned_zcut_005_R_cut_075->at(0)-gen_mpruned_zcut_005_R_cut_075->at(0));
 }
if(puppi_ptprunedsafe_zcut_010_R_cut_075->at(0) > pTmin && gen_ptpruned_zcut_010_R_cut_075->at(0) > pTmin && puppi_ptprunedsafe_zcut_010_R_cut_075->at(0) < pTmax && gen_ptpruned_zcut_010_R_cut_075->at(0) < pTmax){
       JetRsAKPrPUPPI_h[3]->Fill(puppi_mpruned_zcut_010_R_cut_075->at(0)-gen_mpruned_zcut_010_R_cut_075->at(0));


if(puppi_ptsoftdropsafe_beta20->at(0) > pTmin && puppi_ptsoftdropsafe_beta20->at(0) < pTmax ){
       JetRsAKSdPUPPI_h[0]->Fill(puppi_msoftdrop_beta20->at(0));
}
if(puppi_ptsoftdropsafe_beta00->at(0) > pTmin &&  puppi_ptsoftdropsafe_beta00->at(0) < pTmax && gen_ptsoftdrop_beta00->at(0) < pTmax){
       JetRsAKSdPUPPI_h[1]->Fill(puppi_msoftdrop_beta00->at(0)-gen_msoftdrop_beta00->at(0));
}
if(puppi_ptsoftdropsafe_beta10->at(0) > pTmin && gen_ptsoftdrop_beta10->at(0) > pTmin && puppi_ptsoftdropsafe_beta10->at(0) < pTmax && gen_ptsoftdrop_beta10->at(0) < pTmax){
       JetRsAKSdPUPPI_h[2]->Fill(puppi_msoftdrop_beta10->at(0)-gen_msoftdrop_beta10->at(0));
}
if(puppi_ptsoftdropsafe_betam1->at(0) > pTmin && gen_ptsoftdrop_betam1->at(0) > pTmin && puppi_ptsoftdropsafe_betam1->at(0) < pTmax && gen_ptsoftdrop_betam1->at(0) < pTmax){
       JetRsAKSdPUPPI_h[3]->Fill(puppi_msoftdrop_betam1->at(0)-gen_msoftdrop_betam1->at(0));



*/



}//event loop entries



//Normalize histograms


hpf->Scale(1/hpf->Integral());

hchs->Scale(1/hchs->Integral());

JetRsAKTrPF_h[0]->Scale(1/JetRsAKTrPF_h[0]->Integral());
JetRsAKTrPF_h[1]->Scale(1/JetRsAKTrPF_h[1]->Integral());
JetRsAKTrPF_h[2]->Scale(1/JetRsAKTrPF_h[2]->Integral());
JetRsAKTrPF_h[3]->Scale(1/JetRsAKTrPF_h[3]->Integral());

JetRsAKTrCHS_h[0]->Scale(1/JetRsAKTrCHS_h[0]->Integral());
JetRsAKTrCHS_h[1]->Scale(1/JetRsAKTrCHS_h[1]->Integral());
JetRsAKTrCHS_h[2]->Scale(1/JetRsAKTrCHS_h[2]->Integral());
JetRsAKTrCHS_h[3]->Scale(1/JetRsAKTrCHS_h[3]->Integral());

JetRsAKTrPUPPI_h[0]->Scale(1/JetRsAKTrPUPPI_h[0]->Integral());
JetRsAKTrPUPPI_h[1]->Scale(1/JetRsAKTrPUPPI_h[1]->Integral());
JetRsAKTrPUPPI_h[2]->Scale(1/JetRsAKTrPUPPI_h[2]->Integral());
JetRsAKTrPUPPI_h[3]->Scale(1/JetRsAKTrPUPPI_h[3]->Integral());

JetRsAKPrPF_h[0]->Scale(1/JetRsAKPrPF_h[0]->Integral());
JetRsAKPrPF_h[1]->Scale(1/JetRsAKPrPF_h[1]->Integral());
JetRsAKPrPF_h[2]->Scale(1/JetRsAKPrPF_h[2]->Integral());
JetRsAKPrPF_h[3]->Scale(1/JetRsAKPrPF_h[3]->Integral());

JetRsAKPrCHS_h[0]->Scale(1/JetRsAKPrCHS_h[0]->Integral());
JetRsAKPrCHS_h[1]->Scale(1/JetRsAKPrCHS_h[1]->Integral());
JetRsAKPrCHS_h[2]->Scale(1/JetRsAKPrCHS_h[2]->Integral());
JetRsAKPrCHS_h[3]->Scale(1/JetRsAKPrCHS_h[3]->Integral());

JetRsAKPrPUPPI_h[0]->Scale(1/JetRsAKPrPUPPI_h[0]->Integral());
JetRsAKPrPUPPI_h[1]->Scale(1/JetRsAKPrPUPPI_h[1]->Integral());
JetRsAKPrPUPPI_h[2]->Scale(1/JetRsAKPrPUPPI_h[2]->Integral());
JetRsAKPrPUPPI_h[3]->Scale(1/JetRsAKPrPUPPI_h[3]->Integral());

JetRsAKSdPF_h[0]->Scale(1/JetRsAKSdPF_h[0]->Integral());
JetRsAKSdPF_h[1]->Scale(1/JetRsAKSdPF_h[1]->Integral());
JetRsAKSdPF_h[2]->Scale(1/JetRsAKSdPF_h[2]->Integral());
JetRsAKSdPF_h[3]->Scale(1/JetRsAKSdPF_h[3]->Integral());

JetRsAKSdCHS_h[0]->Scale(1/JetRsAKSdCHS_h[0]->Integral());
JetRsAKSdCHS_h[1]->Scale(1/JetRsAKSdCHS_h[1]->Integral());
JetRsAKSdCHS_h[2]->Scale(1/JetRsAKSdCHS_h[2]->Integral());
JetRsAKSdCHS_h[3]->Scale(1/JetRsAKSdCHS_h[3]->Integral());

JetRsAKSdPUPPI_h[0]->Scale(1/JetRsAKSdPUPPI_h[0]->Integral());
JetRsAKSdPUPPI_h[1]->Scale(1/JetRsAKSdPUPPI_h[1]->Integral());
JetRsAKSdPUPPI_h[2]->Scale(1/JetRsAKSdPUPPI_h[2]->Integral());
JetRsAKSdPUPPI_h[3]->Scale(1/JetRsAKSdPUPPI_h[3]->Integral());





gROOT->ProcessLine(".L tdrstyle.C");
setTDRStyle();

 TLatex *tex1 = new TLatex(3.78634,0.172621,"CMS");
   tex1->SetLineWidth(2);

     TLatex *tex2 = new TLatex(23.86,0.172621,"Simulation Preliminary");
   tex2->SetTextFont(52);
   tex2->SetTextSize(0.04323763);
   tex2->SetLineWidth(2);
   
     TLatex *tex3 = new TLatex(123.3919,0.172621,"13 TeV");
   tex3->SetTextFont(42);
   tex3->SetTextSize(0.04323763);
   tex3->SetLineWidth(2);
  ///////////////////////////////////////////////////////////////////
      TLatex *tex4 = new TLatex(6.71375,0.1569342,"QCD, Anti-kT (R=0.8)");
   tex4->SetTextFont(42);
   tex4->SetLineWidth(0);
   tex4->SetTextSize(0.0475614);

      TLatex *tex5 = new TLatex(8.804756,0.1224232,"p_{T} >300 GeV");
   tex5->SetTextFont(42);
   tex5->SetTextSize(0.0475614);
   tex5->SetLineWidth(0);
  
      TLatex *tex6 = new TLatex(6.295548,0.1354209," |#eta|< 2.5");
   tex6->SetTextFont(42);
   tex6->SetLineWidth(0);
   tex6->SetTextSize(0.0475614);
   
     
   
      TLatex *tex7 = new TLatex(7.131951,0.1466257,"<n_{PU}> = 40");
   tex7->SetTextFont(42);
   tex7->SetLineWidth(2);
  

   TLatex *texTr = new TLatex(63.44973,0.1439366,"PF with trimming");
   texTr->SetTextFont(42);
   texTr->SetTextSize(0.04539);
   texTr->SetLineWidth(2);


   TLatex *texSd = new TLatex(63.44973,0.1439366,"PF with softdrop");
   texSd->SetTextFont(42);
   texSd->SetTextSize(0.04539);
   texSd->SetLineWidth(2);

   TLatex *texPr = new TLatex(63.44973,0.1439366,"PF with pruning");
   texPr->SetTextFont(42);
   texPr->SetTextSize(0.04539);
   texPr->SetLineWidth(2);


//////////////////////////////////////////////////////////////////////////////////////

char Cname1[100];
 
TCanvas *c_1D[15]; 
for(int k0=0;k0<15;k0++){
sprintf(Cname1,"c_1D%i",k0);

c_1D[k0]=new TCanvas(Cname1,Cname1,408,129,1400,500);


}

char Legname1[100];


TLegend *leg_1D[15]; 

for(int k0=0;k0<15;k0++){
sprintf(Legname1,"leg_1D%i",k0);

//leg_1D[k0]=new TLegend(0.6,0.7,0.89,0.89);
//leg_1D[k0]= new TLegend(0.35083826,0.6766091,0.7,0.9085201);
leg_1D[k0]=new TLegend(0.5182893,0.4718955,0.9572334,0.8026634);//0.5,0.5864753,0.93899,0.9172432);
leg_1D[k0]->SetTextFont(42);
leg_1D[k0]->SetTextSize(0.04);
leg_1D[k0]->SetLineColor(1);
leg_1D[k0]->SetLineStyle(1);
leg_1D[k0]->SetLineWidth(1);
leg_1D[k0]->SetFillColor(0);
leg_1D[k0]->SetFillStyle(1001);
leg_1D[k0]->SetShadowColor(0);
leg_1D[k0]->SetDrawOption(0);
leg_1D[k0]->SetBorderSize(0);


}

//trim 1

//TLatex *texc0 = new TLatex(11.58537,0.09060423,"Trimmed (r_{filt}=0.2,pT_{frac}=0.05)");
//  texc0->SetLineWidth(0.01);
//  texc0->SetTextSize(0.02);
c_1D[0]->Divide(3,1);
c_1D[0]->Range(-140,-0.02378049,110,0.1591463);
/*
leg_1D[0]->AddEntry(JetRsAKTrCHS_h[0],"Trimmed (r_{filt}=0.2,pT_{frac}=0.05)","lp");
leg_1D[0]->AddEntry(JetRsAKTrCHS_h[1],"Trimmed (r_{filt}=0.1,pT_{frac}=0.03)","lp");
leg_1D[0]->AddEntry(JetRsAKTrCHS_h[2],"Trimmed (r_{filt}=0.2,pT_{frac}=0.03)","lp");
leg_1D[0]->AddEntry(JetRsAKTrCHS_h[3],"Trimmed (r_{filt}=0.3,pT_{frac}=0.03)","lp");
leg_1D[0]->AddEntry(JetRsAKPrCHS_h[0],"Pruned (z_{cut}=0.1,r_{cut}=0.5)","lp");
leg_1D[0]->AddEntry(JetRsAKPrCHS_h[1],"Pruned (z_{cut}=0.05,r_{cut}=0.5)","lp");
leg_1D[0]->AddEntry(JetRsAKPrCHS_h[2],"Pruned (z_{cut}=0.05,r_{cut}=0.75)","lp");
leg_1D[0]->AddEntry(JetRsAKPrCHS_h[3],"Pruned (z_{cut}=0.1,r_{cut}=0.75)","lp");
leg_1D[0]->AddEntry(JetRsAKSdCHS_h[0],"Softdrop #beta = 2  ","lp");
leg_1D[0]->AddEntry(JetRsAKSdCHS_h[1],"Softdrop #beta = 0  ","lp");
leg_1D[0]->AddEntry(JetRsAKSdCHS_h[2],"Softdrop #beta = 1 ","lp");
*/

//leg_1D[0]->AddEntry(JetRsAKSdCHS_h[3],"Softdrop #beta = -1  ","lp");



JetRsAKTrCHS_h[0]->SetLineColor(1);
JetRsAKTrCHS_h[1]->SetLineColor(2);
JetRsAKTrCHS_h[2]->SetLineColor(3);
JetRsAKTrCHS_h[3]->SetLineColor(4);

JetRsAKTrCHS_h[0]->SetLineWidth(2);
JetRsAKTrCHS_h[1]->SetLineWidth(2);
JetRsAKTrCHS_h[2]->SetLineWidth(2);
JetRsAKTrCHS_h[3]->SetLineWidth(2);

//JetRsAKTrCHS_h[0]->SetLineStyle(1);
//JetRsAKTrCHS_h[1]->SetLineStyle(2);
//JetRsAKTrCHS_h[2]->SetLineStyle(3);
//JetRsAKTrCHS_h[3]->SetLineStyle(5);

//JetRsAKTrCHS_h[0]->SetMarkerStyle(2);
//JetRsAKTrCHS_h[1]->SetMarkerStyle(3);
//JetRsAKTrCHS_h[2]->SetMarkerStyle(4);
//JetRsAKTrCHS_h[3]->SetMarkerStyle(5);


JetRsAKPrCHS_h[0]->SetLineColor(1);
JetRsAKPrCHS_h[1]->SetLineColor(2);
JetRsAKPrCHS_h[2]->SetLineColor(3);
JetRsAKPrCHS_h[3]->SetLineColor(4);

JetRsAKPrCHS_h[0]->SetLineWidth(2);
JetRsAKPrCHS_h[1]->SetLineWidth(2);
JetRsAKPrCHS_h[2]->SetLineWidth(2);
JetRsAKPrCHS_h[3]->SetLineWidth(2);

//JetRsAKPrCHS_h[0]->SetLineStyle(1);
//JetRsAKPrCHS_h[1]->SetLineStyle(2);
//JetRsAKPrCHS_h[2]->SetLineStyle(3);
//JetRsAKPrCHS_h[3]->SetLineStyle(5);

//JetRsAKPrCHS_h[0]->SetMarkerStyle(2);
//JetRsAKPrCHS_h[1]->SetMarkerStyle(3);
//JetRsAKPrCHS_h[2]->SetMarkerStyle(4);
//JetRsAKPrCHS_h[3]->SetMarkerStyle(5);


JetRsAKSdCHS_h[0]->SetLineColor(1);
JetRsAKSdCHS_h[1]->SetLineColor(2);
JetRsAKSdCHS_h[2]->SetLineColor(3);

JetRsAKSdCHS_h[0]->SetLineWidth(2);
JetRsAKSdCHS_h[1]->SetLineWidth(2);
JetRsAKSdCHS_h[2]->SetLineWidth(2);


//JetRsAKSdCHS_h[0]->SetLineStyle(1);
//JetRsAKSdCHS_h[1]->SetLineStyle(2);
//JetRsAKSdCHS_h[2]->SetLineStyle(3);

//JetRsAKSdCHS_h[0]->SetMarkerStyle(2);
//JetRsAKSdCHS_h[1]->SetMarkerStyle(3);
//JetRsAKSdCHS_h[2]->SetMarkerStyle(4);

c_1D[0]->cd(1);
leg_1D[1]->AddEntry(JetRsAKTrCHS_h[0],"r_{sub}=0.2,pT_{frac}=0.05","l");
leg_1D[1]->AddEntry(JetRsAKTrCHS_h[1],"r_{sub}=0.1,pT_{frac}=0.03","l");
leg_1D[1]->AddEntry(JetRsAKTrCHS_h[2],"r_{sub}=0.2,pT_{frac}=0.03","l");
leg_1D[1]->AddEntry(JetRsAKTrCHS_h[3],"r_{sub}=0.3,pT_{frac}=0.03","l");
leg_1D[1]->AddEntry(hchs,"ungroomed","l");

JetRsAKTrCHS_h[0]->Draw();
JetRsAKTrCHS_h[0]->GetYaxis()->SetRangeUser(0,0.17);

//JetRsAKTrCHS_h[0]->GetXaxis()->SetTitleSize(0.5);
//JetRsAKTrCHS_h[0]->GetXaxis()->SetLabelSize(0.5);
//JetRsAKTrCHS_h[0]->GetXaxis()->SetTitleSize(1.5);


//JetRsAKTrPF_h[0]->GetYaxis()->SetLabelSize(0.5);
//JetRsAKTrPF_h[0]->GetYaxis()->SetTitleSize(0.5);

JetRsAKTrCHS_h[1]->Draw("SAME");
JetRsAKTrCHS_h[2]->Draw("SAME");
JetRsAKTrCHS_h[3]->Draw("SAME");
hchs->Draw("SAME");

leg_1D[1]->Draw();
tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texTr->Draw();


c_1D[0]->cd(2);
leg_1D[2]->AddEntry(JetRsAKSdCHS_h[0],"#beta = 2  ","lp");
leg_1D[2]->AddEntry(JetRsAKSdCHS_h[1],"#beta = 0  ","lp");
leg_1D[2]->AddEntry(JetRsAKSdCHS_h[2],"#beta = 1 ","lp");
leg_1D[2]->AddEntry(hchs,"ungroomed","l");

JetRsAKSdCHS_h[0]->Draw("");
JetRsAKSdCHS_h[0]->GetYaxis()->SetRangeUser(0,0.17);
JetRsAKSdCHS_h[0]->GetYaxis()->SetLabelOffset(0);
JetRsAKSdCHS_h[1]->Draw("SAME");
JetRsAKSdCHS_h[2]->Draw("SAME");
hchs->Draw("SAME");

//JetRsAKSdCHS_h[3]->Draw("SAME");
leg_1D[2]->Draw();
tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texSd->Draw();
c_1D[0]->cd(3);
leg_1D[3]->AddEntry(JetRsAKPrCHS_h[0],"z_{cut}=0.1,r_{cut}=0.5","lp");
leg_1D[3]->AddEntry(JetRsAKPrCHS_h[1],"z_{cut}=0.05,r_{cut}=0.5","lp");
leg_1D[3]->AddEntry(JetRsAKPrCHS_h[2],"z_{cut}=0.05,r_{cut}=0.75","lp");
leg_1D[3]->AddEntry(JetRsAKPrCHS_h[3],"z_{cut}=0.1,r_{cut}=0.75","lp");
leg_1D[3]->AddEntry(hchs,"ungroomed","l");

JetRsAKPrCHS_h[0]->Draw("");
JetRsAKPrCHS_h[0]->GetYaxis()->SetRangeUser(0,0.17);

JetRsAKPrCHS_h[1]->Draw("SAME");
JetRsAKPrCHS_h[2]->Draw("SAME");
JetRsAKPrCHS_h[3]->Draw("SAME");
hchs->Draw("SAME");

leg_1D[3]->Draw();
tex1->Draw();
tex2->Draw();
tex3->Draw();
tex4->Draw();
tex5->Draw();
tex6->Draw();
tex7->Draw();

texPr->Draw();

/*

char name[100];



for(int j1=1;j1<14;j1++){
sprintf(name,"plot%ig.pdf",j1);
cout<<"name = "<<name<<endl;
c_1D[j1]->SaveAs(name);
c_1D[j1]->Close();

}
*/

for(int j1=1;j1<14;j1++){
c_1D[j1]->Close();
}


c_1D[0]->SaveAs("1DPF_QCD.pdf");
c_1D[0]->SaveAs("1DPF_QCD.png");


}//main




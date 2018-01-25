#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <iostream>
#include <bitset>
#include <TCanvas.h>
#include <algorithm>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <TColor.h>
#include <TLine.h>

int noHlts = 13;

std::string hltsName[13] = {"HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
                            "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi",
                            "HLT_Mu20_TkMu0_Phi",
                            "HLT_Dimuon14_Phi_Barrel_Seagulls",
                            "HLT_Mu25_TkMu0_Phi",
                            "HLT_Dimuon24_Phi_noCorrL1",
                            "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                            "HLT_DoubleMu4_JpsiTrk_Displaced",
                            "HLT_DoubleMu4_Jpsi_Displaced",
                            "HLT_DoubleMu4_3_Jpsi_Displaced",
                            "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
                            "HLT_Dimuon25_Jpsi",
                            "HLT_Dimuon0_Jpsi"};

int skimXTree(std::string path, std::string filename, std::string treename = "xTree", std::string dirname = "rootuple")
{

   TFile *oldfile = TFile::Open((path+filename).data());
   TDirectory *directory = (TDirectory*)oldfile->Get(dirname.data());
   TTree *oldtree = (TTree*)directory->Get(treename.data());
   Long64_t nentries = oldtree->GetEntries();
   ULong64_t event   = 0;
   oldtree->SetBranchAddress("event",&event);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile((treename + "_skim_" + filename).data(),"RECREATE");
   TTree *newtree = oldtree->CloneTree();

   newtree->Print();
   newtree->Write();

   return 0;

}



int selectXTree()
{

   TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmedNP.root");

   TTree *oldtree = (TTree*)oldfile->Get("jTree");
   Long64_t nentries = oldtree->GetEntries();
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0;

   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("l_xy",&xyl);
   oldtree->SetBranchAddress("lErr_xy",&xylErr);
   oldtree->SetBranchAddress("cosAlpha",&cosA);
   oldtree->SetBranchAddress("phi_M",&phiM);
   oldtree->SetBranchAddress("jpsi_M",&jPsiM);

   oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
   oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

   oldtree->SetBranchAddress("ctauPV",&ctau);
   oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

   oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
   oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("skimmedNPCos.root","RECREATE");
   TTree *newtree = oldtree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      std::bitset<16> pM(phiMType);
      std::bitset<16> pP(phiPType);
      // std::cout << phiMType << "-" << binary << std::endl;
      // std::cout<<vProb <<std::endl;
      //if (vProb > 0.0) newtree->Fill();
      if (jPsiM > 2.8 && phiM < 1.1 && cosA > 0.995 && phiM > 0.95 && pP.test(1) && pM.test(1) && vProb > 0.1  ) newtree->Fill();
   }

   newtree->Print();
   newtree->Write();

   return 0;

}

int drawPTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

   TFile *oldfile = TFile::Open(path.data());
   TTree *oldtree = (TTree*)oldfile->Get(treename.data());

   Long64_t nentries = oldtree->GetEntries();
   Double_t pM = 0.0;
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

   oldtree->SetBranchAddress("pM",&pM);
   oldtree->SetBranchAddress("p_vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);

   //Create a new file + a clone of old tree in new file
   TCanvas c("c","c",1200,1600);

   TFile *newfile = new TFile("drawSkim.root","RECREATE");
   // for(int j = 0; j < 1; j++)
   // {


     Double_t xmin = 5.0, xmax = 6.0;
     Int_t xBin = ((xmax - xmin)/0.01);

     TTree *newtree = oldtree->CloneTree(0);
     // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
     TH1F* phiHist = new TH1F("phiHist","phiHist",250,0.0,1.25);

     std::vector<TH1F*> phiHists;


     for (size_t i = 0; i < 13; i++)
          phiHists.push_back(new TH1F((hltsName[i] + "_phi").data(),(hltsName[i] + "_phi").data(),200,0.25,1.25));


     for (Long64_t i=0;i<nentries; i++) {
        oldtree->GetEntry(i);
        std::bitset<16> tB(trigger);
        // std::bitset<16> pM(phiMType);
        // std::bitset<16> pP(phiPType);
       bool test = false;
       // bool jpsimass = jPsiM < 3.15 && jPsiM > 3.0;
       // bool phimass = phiM < 1.06 && phiM > 0.98;
       for (int j = 0; j < 13; j++){
          // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
          if (tB.test(j))
          {
            test = true;
            phiHists[j]->Fill(pM);

          }
        }
	      // if(cosA > 0.995 && phimass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
        if (test)
        phiHist->Fill(pM);

     }


     //newtree->Draw("phi_M","","same");
   // }
   phiHist->SetMinimum(1.0);
   phiHist->SetMaximum(phiHist->GetMaximum()*5.0);
   //oldtree->Draw("phi_M");
   // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

   phiHist->SetLineColor(kBlue);
   phiHist->Write();
   phiHist->Draw();

   TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
   leg->AddEntry(phiHist,(phiHist->GetName()),"l");
   for (size_t i = 0; i < 13; i++)
    {
      phiHists[i]->SetLineColor(colors[i]);
      phiHists[i]->SetLineWidth(2);
      if(i>5) phiHists[i]->SetLineStyle(kDashed);
      phiHists[i]->Draw("same");
      leg->AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
      phiHists[i]->Write();
    }

    leg->Draw();
    c.SetLogy(1);
    c.SaveAs("phitriggerCheck.png");
    c.SaveAs("phitriggerCheck.eps");
    c.SaveAs("phitriggerCheck.root");

   return 0;

}


int drawJTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

   TFile *oldfile = TFile::Open(path.data());
   TTree *oldtree = (TTree*)oldfile->Get(treename.data());

   Long64_t nentries = oldtree->GetEntries();
   Double_t pM = 0.0;
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

   oldtree->SetBranchAddress("jM",&jM);
   oldtree->SetBranchAddress("p_vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);

   //Create a new file + a clone of old tree in new file
   TCanvas c("c","c",1200,1600);

   TFile *newfile = new TFile("drawSkim.root","RECREATE");
   // for(int j = 0; j < 1; j++)
   // {


     Double_t xmin = 5.0, xmax = 6.0;
     Int_t xBin = ((xmax - xmin)/0.01);

     TTree *newtree = oldtree->CloneTree(0);
     // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
     TH1F* jpsiHist = new TH1F("jpsiHist","jpsiHist",140,2.6,3.3);

     std::vector<TH1F*> phiHists;


     for (size_t i = 0; i < 13; i++)
          jpsiHists.push_back(new TH1F((hltsName[i] + "_jpsi").data(),(hltsName[i] + "_jpsi").data(),140,2.6,3.3));


     for (Long64_t i=0;i<nentries; i++) {
        oldtree->GetEntry(i);
        std::bitset<16> tB(trigger);
        // std::bitset<16> pM(phiMType);
        // std::bitset<16> pP(phiPType);
       bool test = false;
       // bool jpsimass = jPsiM < 3.15 && jPsiM > 3.0;
       // bool phimass = phiM < 1.06 && phiM > 0.98;
       for (int j = 0; j < 13; j++){
          // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
          if (tB.test(j))
          {
            test = true;
            jpsiHists[j]->Fill(jM);

          }
        }
	      // if(cosA > 0.995 && phimass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
        if (test)
        jpsiHist->Fill(jM);

     }


     //newtree->Draw("phi_M","","same");
   // }
   jpsiHist->SetMinimum(1.0);
   jpsiHist->SetMaximum(jpsiHist->GetMaximum()*5.0);
   //oldtree->Draw("phi_M");
   // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

   jpsiHist->SetLineColor(kBlue);
   jpsiHist->Write();
   jpsiHist->Draw();

   TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
   leg->AddEntry(jpsiHist,(jpsiHist->GetName()),"l");
   for (size_t i = 0; i < 13; i++)
    {
      jpsiHists[i]->SetLineColor(colors[i]);
      jpsiHists[i]->SetLineWidth(2);
      if(i>5) phiHists[i]->SetLineStyle(kDashed);
      jpsiHists[i]->Draw("same");
      leg->AddEntry(jpsiHists[i],(jpsiHists[i]->GetName()),"l");
      jpsiHists[i]->Write();
    }

    leg->Draw();
    c.SetLogy(1);
    c.SaveAs("phitriggerCheck.png");
    c.SaveAs("phitriggerCheck.eps");
    c.SaveAs("phitriggerCheck.root");

   return 0;

}


int drawXTree(std::string path = "/Users/adrianodiflorio/Documents/Git/X4140/iPythons/xTree.root",std::string treename = "xTree")
{

  UInt_t colors[13] = {1,2,3,6,7,8,30,40,46,38,29,34,9};

   TFile *oldfile = TFile::Open(path.data());
   TTree *oldtree = (TTree*)oldfile->Get(treename.data());

   Long64_t nentries = oldtree->GetEntries();
   Double_t xM = 0.0;
   Double_t xyl   = 0.0;
   Double_t xylErr   = 0.0;
   Double_t cosA  = 0.0;

   Double_t phiM  = 0.0;
   Double_t jPsiM  = 0.0;

   Double_t ctau  = 0.0;
   Double_t ctauErr  = 0.0;

   Double_t vProb  = 0.0;

   Int_t phiMType = 0, phiPType = 0;
   UInt_t phi_trigger = 0, jpsi_trigger = 0, trigger = 0;

   oldtree->SetBranchAddress("xM",&xM);
   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);
   oldtree->SetBranchAddress("l_xy",&xyl);
   oldtree->SetBranchAddress("lErr_xy",&xylErr);
   oldtree->SetBranchAddress("cosAlpha",&cosA);
   oldtree->SetBranchAddress("phi_M",&phiM);
   oldtree->SetBranchAddress("jpsi_M",&jPsiM);

   oldtree->SetBranchAddress("phi_muonM_type",&phiMType);
   oldtree->SetBranchAddress("phi_muonP_type",&phiPType);

   oldtree->SetBranchAddress("ctauPV",&ctau);
   oldtree->SetBranchAddress("ctauErrPV",&ctauErr);

   oldtree->SetBranchAddress("phi_trigger",&phi_trigger);
   oldtree->SetBranchAddress("jpsi_trigger",&jpsi_trigger);
   //Create a new file + a clone of old tree in new file
   TCanvas c("c","c",1200,1600);

   TFile *newfile = new TFile("drawSkim.root","RECREATE");
   // for(int j = 0; j < 1; j++)
   // {


     Double_t xmin = 5.0, xmax = 6.0;
     Int_t xBin = ((xmax - xmin)/0.01);

     TTree *newtree = oldtree->CloneTree(0);
     // TH1F* phi_triggrHist = new TH1F("phi_triggrHist","phi_triggrHist",600,0.6,1.2);
     TH1F* phiHist = new TH1F("phiHist","phiHist",250,0.0,1.25);
     TH1F* jpsiHist = new TH1F("jpsiHist","jpsiHist",140,2.6,3.3);
     TH1F* xHist = new TH1F("xHist","xHist",xBin,xmin,xmax);

     std::vector<TH1F*> phiHists;
     std::vector<TH1F*> jpsiHists;
     std::vector<TH1F*> xHists;


     for (size_t i = 0; i < noHlts; i++)
     {
          phiHists.push_back(new TH1F((hltsName[i] + "_phi").data(),(hltsName[i] + "_phi").data(),200,0.25,1.25));
          jpsiHists.push_back(new TH1F((hltsName[i] + "_jpsi").data(),(hltsName[i] + "_jpsi").data(),140,2.6,3.3));
          xHists.push_back(new TH1F((hltsName[i] + "_x").data(),(hltsName[i] + "_x").data(),xBin,xmin,xmax));
     }


     for (Long64_t i=0;i<nentries; i++) {
        oldtree->GetEntry(i);
        std::bitset<16> tB(trigger);
        // std::bitset<16> pM(phiMType);
        // std::bitset<16> pP(phiPType);
        bool test = false;
       bool jpsimass = jPsiM < 3.15 && jPsiM > 3.0;
       bool phimass = phiM < 1.06 && phiM > 0.98;
       for (int j = 0; j < 13; j++){
          // if (tB.test(j) && cosA > 0.995 && vProb > 0.01 && xyl/xylErr > 2.0 && trigger > 0)
          if (tB.test(j))
          {
            test = true;
            phiHists[j]->Fill(phiM);
            jpsiHists[j]->Fill(jPsiM);
   	    //  if(phimass && jpsimass && cosA > 0.995 && vProb>0.02 && xyl/xylErr > 2.0)

            xHists[j]->Fill(xM);
          }
        }
	      // if(cosA > 0.995 && phimass && jpsimass && vProb>0.02 && xyl/xylErr > 2.0)
        if (test){
        xHist->Fill(xM);
        phiHist->Fill(phiM);
        jpsiHist->Fill(jPsiM);
      }
     }


     //newtree->Draw("phi_M","","same");
   // }
   phiHist->SetMinimum(1.0);
   phiHist->SetMaximum(phiHist->GetMaximum()*5.0);
   //oldtree->Draw("phi_M");
   // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

   phiHist->SetLineColor(kBlue);
   phiHist->Write();
   phiHist->Draw();

   TLegend* leg = new TLegend(0.1,0.5,0.45,0.9);
   leg->AddEntry(phiHist,(phiHist->GetName()),"l");
   for (size_t i = 0; i < 13; i++)
    {
      phiHists[i]->SetLineColor(colors[i]);
      phiHists[i]->SetLineWidth(2);
      if(i>5) phiHists[i]->SetLineStyle(kDashed);
      phiHists[i]->Draw("same");
      leg->AddEntry(phiHists[i],(phiHists[i]->GetName()),"l");
      phiHists[i]->Write();
    }

    leg->Draw();
    c.SetLogy(1);
    c.SaveAs("phitriggerCheck.png");
    c.SaveAs("phitriggerCheck.eps");
    c.SaveAs("phitriggerCheck.root");

    jpsiHist->SetMinimum(1.0);
    jpsiHist->SetMaximum(jpsiHist->GetMaximum()*5.0);
    //oldtree->Draw("phi_M");
    // TH1F* phi_triggrHist = (TH1F*)gDirectory->Get("phi_triggrHist");

    jpsiHist->SetLineColor(kBlue);
    jpsiHist->Write();
    jpsiHist->Draw();

    leg = new TLegend(0.1,0.5,0.4,0.9);
    leg->AddEntry(jpsiHist,(jpsiHist->GetName()),"l");
    for (size_t i = 0; i < 13; i++)
     {
       jpsiHists[i]->SetLineColor(colors[i]);
       jpsiHists[i]->SetLineWidth(2);
       if(i>5) jpsiHists[i]->SetLineStyle(kDashed);
       jpsiHists[i]->Draw("same");
       leg->AddEntry(jpsiHists[i],(jpsiHists[i]->GetName()),"l");
       jpsiHists[i]->Write();
     }


   // phi_triggrHist->Draw("same");
   leg->Draw();
   c.SetLogy(1);
   c.SaveAs("jpsitriggerCheck.png");
   c.SaveAs("jpsitriggerCheck.eps");
   c.SaveAs("jpsitriggerCheck.root");

   xHist->SetMinimum(0.1);
   xHist->SetMaximum(xHist->GetMaximum());//*5.0);

   TLine line (5.35,0.1,5.35,xHist->GetMaximum());
   line.SetLineColor(kRed);
   line.SetLineWidth(2);

   xHist->SetLineColor(kBlue);
   xHist->Write();
   xHist->Draw();

   leg = new TLegend(0.1,0.55,0.35,0.9);
   leg->AddEntry(xHist,(xHist->GetName()),"l");
   for (size_t i = 0; i < 13; i++)
    {
      xHists[i]->SetLineColor(colors[i]);
      xHists[i]->SetLineWidth(2);
      if(i>5) xHists[i]->SetLineStyle(kDashed);
      xHists[i]->Draw("same");
      leg->AddEntry(xHists[i],(xHists[i]->GetName()),"l");
      xHists[i]->Write();
    }

  // line.Draw();
  // phi_triggrHist->Draw("same");
  leg->Draw();
  c.SetLogy(0);
  c.SaveAs("xtriggerCheck.png");
  c.SaveAs("xtriggerCheck.eps");
  c.SaveAs("xtriggerCheck.root");

   return 0;

}

int selectXTreeHLT()
{

   TFile *oldfile = TFile::Open("/Users/adrianodiflorio/Documents/Git/X4140/iPythons/skimmed_vProb.root");

   TTree *oldtree = (TTree*)oldfile->Get("xTree");
   Long64_t nentries = oldtree->GetEntries();
   Double_t vProb   = 0.0;
   UInt_t trigger = 0;
   oldtree->SetBranchAddress("vProb",&vProb);
   oldtree->SetBranchAddress("trigger",&trigger);
   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("skimmed.root","RECREATE");
   TTree *newtree = oldtree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      // std::bitset<16> binary(trigger);
      // std::cout << trigger << "-" << binary << std::endl;
      if(trigger>0)
      newtree->Fill();
   }

   newtree->Print();
   newtree->Write();

   return 0;

}

/*
# coding: utf-8

# In[1]:


import ROOT
from ROOT import TFile,TH1,TH1F,TCanvas,TNtuple,TTreeReader,TTreeReaderValue
from ROOT import RooFit
from ROOT.RooFit import Layout
from ROOT import RooStats
from ROOT import RooAbsData
RooAbsData.setDefaultStorageType ( RooAbsData.Tree )
from array import array
import sys


# In[3]:


from ROOT import RooRealVar,RooAbsPdf,RooChebychev,RooExponential,RooGaussian,RooAbsPdf,RooPlot,RooAddPdf,RooDataHist,RooArgSet,RooArgList
from ROOT import kGreen,kRed,kBlack,kBlue,kDashed,kDotted,kMagenta,RooVoigtian
from ROOT.RooFit import Components,LineColor,LineStyle,Name,Normalization,Range,SelectVars
from ROOT import RooDataSet,RooFormulaVar,RooLinkedList


# In[7]:


no_hlts = 13


# In[8]:


#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/09Jan2017.root"
#rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/allphi_DataF.root"
rootfile = "/Users/adrianodiflorio/Desktop/mmkk2017/phiJpsiTriggersBCDEF.root"
inputfile = TFile(rootfile,"READ")
inputfile.ls()
xTupleDir = (inputfile.Get("rootuple"))
xTupleDir.ls()
#pTuple = (xTupleDir.Get("pTree"))
xTuple = (xTupleDir.Get("xTree"))
#jTuple = (xTupleDir.Get("jTree"))

event = 2

xTuple.SetBranchAddress("event",event);
newfile = TFile("small.root","recreate");
newtree = xTuple.CloneTree();
newtree.CopyEntries(xTuple);

newtree.Write()

sys.exit()

# In[16]:

file = TFile("newFile.root","RECREATE")
canvas = TCanvas("canvas","canvas",1200,1000)
mass = RooRealVar("xM","M(#mu#mu#mu#mu)[GeV]",5.15,5.55)
trigger = RooRealVar("trigger","trigger",0.0,10000)
vProb = RooRealVar("vProb","vProb",-1.0,1.0)
alldata = RooDataSet("alldata","alldata",xTuple,RooArgSet(mass), RooFormulaVar("vProb","vProb","vProb>0.01",RooArgList(vProb)))#,cutFormula)
frame = mass.frame(Range(5.15,5.55))
alldata.plotOn(frame,RooLinkedList())
alldata.Write()
frame.Draw()


# In[ ]:


canvas.SaveAs("testCanvas.eps")
*/

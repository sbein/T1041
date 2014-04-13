#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TBEvent.h"
#include "shashlik.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

const int MAXADC=4095;

// inputs data file and event in file to display (default is to integrate all)
void calDisplay(TString fdat, int ndisplay=0){
  TFile *f = new TFile(fdat);
  if (f->IsZombie()){
    cout << "Cannot open file: " << fdat << endl;
    return;
  }
  Bool_t singleEvent=ndisplay>0;
  Mapper *mapper=Mapper::Instance();
  TH2F *hModU=new TH2F("hModU","Modules UpSteam RO",4,0.5,4.5,4,0.5,4.5);
  TH2F *hModD=new TH2F("hModD","Modules DownStream RO",4,0.5,4.5,4,0.5,4.5);
  TH2F *hChanU=new TH2F("hChanU","Channels UpStream RO",8,0.5,4.5,8,0.5,4.5);
  TH2F *hChanD=new TH2F("hChanD","Channels DownStream RO",8,0.5,4.5,8,0.5,4.5);
  hModU->SetMinimum(0);
  hModD->SetMinimum(0);
  hChanU->SetMinimum(0);
  hChanU->SetMinimum(0);
  if (singleEvent){
    hModU->SetMaximum(MAXADC*4);
    hModD->SetMaximum(MAXADC*4);
    hChanU->SetMaximum(MAXADC);
    hChanU->SetMaximum(MAXADC);
  }
  gStyle->SetOptStat(0);

  TBEvent *event = new TBEvent();
  TTree *t1041 = (TTree*)f->Get("t1041"); 
  TBranch *bevent = t1041->GetBranch("tbevent");
  bevent->SetAddress(&event);

  Int_t start=0; Int_t end=t1041->GetEntries();
  if (singleEvent) {
    start=ndisplay-1;
    end=ndisplay;
  }
  for (Int_t i=start; i<end; i++) {
    t1041->GetEntry(i);
    for (Int_t j=0; j<event->NPadeChan(); j++){
      PadeChannel pch=event->GetPadeChan(j);
      // loop over ADC samples
      UShort_t* wform=pch.GetWform();
      UShort_t max=0;
      // find the value to plot (just using the peak sample value now)
       for (Int_t k=0; k<event->GetPadeChan(j).__SAMPLES(); k++){
	if (wform[k]>200 && wform[k]>max) max=wform[k];
      }
      ///////////////////////////////////////////////////////////////
      int module,fiber;
      mapper->Pade2Fiber(pch.GetBoardID(), pch.GetChannelID(), module, fiber);
      int xm,ym;
      mapper->ModuleXY(module,xm,ym);
      if (module<0) hModU->Fill(xm,ym,max);
      else hModD->Fill(xm,ym,max);
      
      float xf,yf;
      int fiberID=module*100+fiber;
      mapper->FiberXY(fiberID, xf, yf);
      if (module<0) hChanU->Fill(xf,yf,max);
      else hChanD->Fill(xf,yf,max);

    }
  }
  TCanvas *c1=new TCanvas("c","c",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hModD->Draw("colz");
  c1->cd(2);
  hModU->Draw("colz");
  c1->cd(3)->SetGrid();
  hChanD->Draw("colz");
  c1->cd(4)->SetGrid();
  hChanU->Draw("colz");
}



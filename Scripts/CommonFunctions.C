#ifndef COMMONFUNCTIONS__C
#define COMMONFUNCTIONS__C
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include "TObject.h"

TH1 *getHistRebinned(TH1 *inh, Int_t nNewBins) {
  TH1 *reth = (TH1*)inh->Clone(Form("%s_Rebinned",inh->GetName()));
  reth->SetDirectory(0);
  reth->Reset();
  reth->RebinX(nNewBins);
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc = inh->GetBinContent(i);
    Double_t be = inh->GetBinError(i);
    if(be==0) continue;
    Int_t bind = reth->FindBin(inh->GetBinCenter(i));
    Double_t sumw = 1./(be*be);
    reth->SetBinContent(bind,reth->GetBinContent(bind)+bc*sumw);
    reth->SetBinError(bind,reth->GetBinError(bind)+sumw);
  }
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    Double_t bc = reth->GetBinContent(i);
    Double_t be = reth->GetBinError(i);
    if(be==0) continue;
    reth->SetBinContent(i,bc/be);
    reth->SetBinError(i,TMath::Sqrt(1./be));
  };
  return reth;
}
TH1 *getHistRebinned(TH1 *inh, TH1 *hTarget) {
  TH1 *reth = (TH1*)hTarget->Clone(Form("%s_Rebinned",inh->GetName()));
  reth->SetDirectory(0);
  reth->Reset();
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc = inh->GetBinContent(i);
    Double_t be = inh->GetBinError(i);
    if(be==0) continue;
    Int_t bind = reth->FindBin(inh->GetBinCenter(i));
    Double_t sumw = 1./(be*be);
    reth->SetBinContent(bind,reth->GetBinContent(bind)+bc*sumw);
    reth->SetBinError(bind,reth->GetBinError(bind)+sumw);
  }
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    Double_t bc = reth->GetBinContent(i);
    Double_t be = reth->GetBinError(i);
    if(be==0) continue;
    reth->SetBinContent(i,bc/be);
    reth->SetBinError(i,TMath::Sqrt(1./be));
  };
  return reth;
}

TLegend *Legend(Double_t x0=0.5, Double_t y0=0.5, Double_t x1=0.8, Double_t y1=0.8) {
  TLegend *tleg = new TLegend(x0,y0,x1,y1);
  tleg->SetTextFont(43);
  tleg->SetTextSize(30);
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  return tleg;
}
TObject *getObj(TObject *indf, TString objName, Bool_t getClone=kFALSE) {
  TObject *obj = (TObject*)((TDirectory*)indf)->Get(objName.Data());
  if(!indf) { printf("Directory provided does not exist!\n"); return 0; };
  if(!obj) { printf("Object %s not found in %s!\n", objName.Data(), indf->GetName()); return 0; };
  if(getClone) return (TObject*)obj->Clone(objName.Data());
  return obj;
}
TFile *getFile(TString fiName, TString fiMode="READ") {
  TFile *retFile = new TFile(fiName.Data(),fiMode.Data());
  if(!retFile) { printf("File %s not found!\n", fiName.Data()); return 0; };
  if(retFile->IsZombie()) { printf("File %s is a zombie!\n",fiName.Data()); delete retFile; return 0; };
  return retFile;
}
#endif

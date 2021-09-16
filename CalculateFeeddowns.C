#include "Scripts/FD.h"
#include "TList.h"
#include "Scripts/CommonFunctions.C"

TH1 *CalculateOneFeeddown(TList *dataList, TList *mcList) {
  TH3D *hmc = (TH3D*)mcList->FindObject("DCAxy_noChi2");
  TH2D *hwdc = (TH2D*)mcList->FindObject("WithinDCA_withChi2");
  TH3D *hdt = (TH3D*)dataList->FindObject("DCAxy_noChi2");
  FD *fd = new FD(hdt,0,hmc,hwdc);
  TH1 *reth = fd->PerformFit(0.2,3);
  delete fd;
  return reth;
};
TH1 *CalculateEfficiency(TList *mcList) {
  TH2D *hrec = (TH2D*)mcList->FindObject("Spectra_ch");
  TH2D *hgen = (TH2D*)mcList->FindObject("Spectra_ch_Gen");
  TH1 *reth = hrec->ProjectionX("effres");
  TH1 *den  = hgen->ProjectionX("hgen");
  reth->Divide(den);
  delete den;
  return reth;
};
void CalculateAllFeeddowns(TString dataFileName, TString mcFileName, TString outputFileName) {
  TFile *dtFile = getFile(dataFileName);
  TFile *mcFile = getFile(mcFileName);
  TList *lok = dtFile->GetListOfKeys();
  TList *outList = new TList();
  outList->SetName("FDList");
  outList->SetOwner(kTRUE);
  TList *effList = new TList();
  effList->SetName("EffAndFD");
  for(Int_t i=0;i<lok->GetEntries();i++) {
    printf("Processing %i\n",i);
    TString lnm(lok->At(i)->GetName());
    if(!lnm.Contains("SpectraList")) continue;
    if(lnm.EqualTo("SpectraList_1")) continue;
    TList *dtList = (TList*)getObj(dtFile,lnm);
    TList *mcList = (TList*)getObj(mcFile,lnm);
    if(!dtList || !mcList) continue;
    TH1 *htoadd = CalculateOneFeeddown(dtList,mcList);
    TString hName("FD");
    if(lnm.Contains("_") && !lnm.EqualTo("SpectraList_0")) hName += lnm(lnm.Index("_"),lnm.Length());
    htoadd->SetName(hName.Data());
    outList->Add(htoadd);
    TH1 *heff = CalculateEfficiency(mcList);
    heff->Multiply(htoadd);
    hName="EffRescaled_Cent0";
    if(lnm.Contains("_") && !lnm.EqualTo("SpectraList_0")) hName += "_SystFlag"+lnm(lnm.Index("_")+1,lnm.Length()) + "_";
    heff->SetName(hName.Data());
    effList->Add(heff);
  };
  TFile *outFile = getFile(outputFileName,"RECREATE");
  outList->Write("FDList",TObject::kSingleKey);
  effList->Write("EffAndFD",TObject::kSingleKey);
  outFile->Close();
}

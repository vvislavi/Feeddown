/* Note 2021-09-16:
  -- FD is estimated using 2015 dataset only. Fitting 2018 is rather unstable and fluctuating quite a bit.
     Also, in some cases, the fit fails at higher pt. Possibly more statistics would help, but one should then
     likely run on train.
     Checking the available fd's, it looks like the difference between 2015 and 2018 is < 1% anyways. Maybe
     I can just add an extra % to the systematics.
  */
#include "Scripts/FD.h"
#include "TList.h"
#include "Scripts/CommonFunctions.C"
#include "TF1.h"
Double_t g_xmax=3;
TH1 *CalculateOneFeeddown(TList *dataList, TList *mcList, Bool_t quitPrematurely=kFALSE) {
  TH3D *hmc = (TH3D*)mcList->FindObject("DCAxy_noChi2");
  TH2D *hwdc = (TH2D*)mcList->FindObject("WithinDCA_withChi2");
  TH3D *hdt = (TH3D*)dataList->FindObject("DCAxy_noChi2");
  FD *fd = new FD(hdt,0,hmc,hwdc);
  if(quitPrematurely) return 0;
  TH1 *reth = fd->PerformFit(0.2,g_xmax);
  delete fd;
  // delete hmc;
  // delete hwdc;
  // delete hdt;
  Double_t lval = reth->GetBinContent(reth->FindBin(g_xmax-(1e-6)));
  for(Int_t i=reth->FindBin(g_xmax);i<=reth->FindBin(3.0 - 1e-6);i++) reth->SetBinContent(i,lval);
  // if(g_xmax<3) {
  //   Double_t yval=0;
  //   Int_t i=reth->FindBin(3- (1e-6));
  //   do {i--; yval=reth->GetBinContent(i); } while((i>1) && (yval=0));
  //   for(Int_t j=i;j<=reth->FindBin(3- (1e-6)); j++) reth->SetBinContent(j,yval);
  // }
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
void ModifyList(TList *inlist) {
  TH3 *hwchi = (TH3*)inlist->FindObject("DCAxy_withChi2");
  TH2 *hwithinDCA = (TH2*)hwchi->Project3D("yx");
  inlist->Remove(inlist->FindObject("WithinDCA_withChi2"));
  hwithinDCA->SetName("WithinDCA_withChi2");
  inlist->Add(hwithinDCA);
  inlist->Remove(inlist->FindObject("Spectra_ch"));
  hwchi->GetYaxis()->SetRange(1,1);
  TH2 *hch = (TH2*)hwchi->Project3D("zx");
  hch->SetName("Spectra_ch");
  inlist->Add(hch);
}
TF1 *getIncreasingFunction(Double_t increment) {
  Double_t a = increment/2.8;
  Double_t b = 1-0.2*a;
  TF1 *retf = new TF1("incrFunc","[0]*x+[1]",0,10);
  retf->SetParameters(a,b);
  return retf;
};
// TList *dtList;
// TList *mcList;
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
    if(i==6 || i==14 || i==15 || i==16 || i==17) g_xmax=2.7;
    else if(i==18) g_xmax = 2.4;
    else g_xmax=3.0;
    printf("Processing %i\n",i);
    TString lnm(lok->At(i)->GetName());
    if(!lnm.Contains("SpectraList")) continue;
    TList *dtList = (TList*)getObj(dtFile,lnm);
    TList *mcList = (TList*)getObj(mcFile,lnm);
    if(!dtList || !mcList) continue;
    if(lnm.EqualTo("SpectraList_1")) ModifyList(mcList);
    TH1 *htoadd = CalculateOneFeeddown(dtList,mcList);
    TString hName("FD");
    if(lnm.Contains("_") && !lnm.EqualTo("SpectraList_0")) hName += lnm(lnm.Index("_"),lnm.Length());
    htoadd->SetName(hName.Data());
    outList->Add(htoadd);
    delete dtList;
    delete mcList;
  };
  TFile *outFile = getFile(outputFileName,"RECREATE");
  outList->Write("FDList",TObject::kSingleKey);
  outFile->Close();
}
void ReplaceAxis(TH1 **target, TH1 *source) {
  TH1 *trg = *target;
  TH1 *reth = (TH1*)source->Clone("temph");
  reth->Reset();
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    Double_t bcent = reth->GetBinCenter(i);
    Int_t tbin = trg->FindBin(bcent);
    Double_t bcon = trg->GetBinContent(tbin);
    if(!bcon) continue;
    reth->SetBinContent(i,bcon);
  }
  TString nameBU(trg->GetName());
  delete trg;
  reth->SetName(nameBU.Data());
  (*target) = reth;
}
void CalculateAllEfficiencies(TString fdFileName, TString mcFileName, TString outputFileName) {
  TFile *fdFile = getFile(fdFileName);
  TFile *mcFile = getFile(mcFileName);
  TList *fdList = (TList*)getObj(fdFile,"FDList");
  if(!fdList) return;
  TList *effList = new TList();
  effList->SetName("EffAndFD");
  for(Int_t i=0;i<fdList->GetEntries();i++) {
    printf("Processing %i\n",i);
    TString fdnm(fdList->At(i)->GetName());
    printf("FD name is %s\n",fdnm.Data());
    if(!fdnm.Contains("FD")) continue;
    TString mcnm("SpectraList");
    if(fdnm.Contains("_")) mcnm.Append(fdnm(fdnm.Index("_"),fdnm.Length())); else mcnm.Append("_0");
    printf("mcnm is %s\n",mcnm.Data());
    TList *mcList = (TList*)getObj(mcFile,mcnm);
    if(!fdList || !mcList) continue;
    if(mcnm.EqualTo("SpectraList_1")) ModifyList(mcList);
    TH1 *htoadd = (TH1*)fdList->At(i);
    TH1 *heff = CalculateEfficiency(mcList);
    ReplaceAxis(&heff,htoadd);
    heff->Multiply(htoadd);
    TString hName="EffRescaled_Cent0";
    if(fdnm.Contains("_")) hName += "_SystFlag"+fdnm(fdnm.Index("_")+1,fdnm.Length()) + "_";
    heff->SetName(hName.Data());
    effList->Add(heff);
    if(mcnm.EqualTo("SpectraList_0")) { //Efficiency + 4% at pt 3 GeV/c
      TF1 *ftemp = getIncreasingFunction(0.04);
      TH1 *heffMod1 = (TH1*)heff->Clone("EffRescaled_Cent0_SystFlag25_");
      heffMod1->Multiply(ftemp);
      delete ftemp;
      TH1 *heffMod2 = (TH1*)heff->Clone("EffRescaled_Cent0_SystFlag26_");
      ftemp = getIncreasingFunction(-0.04);
      heffMod2->Multiply(ftemp);
      effList->Add(heffMod1);
      effList->Add(heffMod2);
    }
  };
  TFile *outFile = getFile(outputFileName,"RECREATE");
  effList->Write("EffAndFD",TObject::kSingleKey);
  outFile->Close();
}

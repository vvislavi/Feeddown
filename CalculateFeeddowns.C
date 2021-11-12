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
#include "AliEffFDContainer.h"
TList *getSomeList(TString infi) {
  TFile *tf = getFile(infi);
  AliEffFDContainer *fd = (AliEffFDContainer*)getObj(tf,"Eff_YS_FB32_EffAndFDEtaM0806");
  return fd->fOutList;
}
AliEffFDContainer *getSomeEF(TString objname, TString infi) {
  TFile *tf = getFile(infi);
  AliEffFDContainer *fd = (AliEffFDContainer*)getObj(tf,objname);
  tf->Close();
  return fd;
}
TString g_mcFile   = "~/AliceData/pPb_Efficiencies_All/LHC20f11c.root";//DataYukoMerged.root";
TString g_dataFile = "~/AliceData/pPb_Efficiencies_All/LHC16_pass1/merged.root";//MCYukoMerged.root";
TString g_objName;//  = "Eff_YS_FB32_EffAndFDEtaM0806";
//Name collections
TString etaPFs[]  = {"0806","0604","0402","0200","0002","0204","0406","0608"};
TString namePFs[] = {"VV_FB96","YS_FB32","ZZ_FB96"};
void ConfigObjName(TString lNamePF, TString lEtaPF) {
  g_objName = Form("Eff_%s_EffAndFDEtaM%s",lNamePF.Data(),lEtaPF.Data());
};

AliEffFDContainer *getDataList(TString objName=g_objName, TString infi=g_dataFile) { return getSomeEF(objName,infi); };
AliEffFDContainer *getMCList(TString objName=g_objName, TString infi=g_mcFile) { return getSomeEF(objName,infi); };
AliEffFDContainer *efdt, *efmc;
void Init() {
  efdt = getDataList(g_objName,g_dataFile);
  efmc = getMCList(g_objName,g_mcFile);
};
Double_t g_xmax=2.5;
FD *fd;
TH1 *CalculateOneFeeddown(AliEffFDContainer *dtef, AliEffFDContainer *mcef, Bool_t quitPrematurely=kFALSE, Bool_t wChi=kFALSE, Int_t centBin=0) {
  TString wChiFlag = wChi?"with":"no";
  TH3D *hmc = dynamic_cast<TH3D*>(mcef->fetchObj("DCAxy_"+wChiFlag+"Chi2"));
  TH3D *hdt = dynamic_cast<TH3D*>(dtef->fetchObj("DCAxy_"+wChiFlag+"Chi2"));
  TH2D *hwdc = dynamic_cast<TH2D*>(mcef->fetchObj("WithinDCA_withChi2"));
  fd = new FD(hdt,0,hmc,hwdc);
  if(quitPrematurely) return 0;
  TH1 *reth = fd->PerformFit(0.2,g_xmax,"",0,centBin);
  delete fd;
  Double_t lval = reth->GetBinContent(reth->FindBin(g_xmax-(1e-6)));
  for(Int_t i=reth->FindBin(g_xmax);i<=reth->FindBin(3.0 - 1e-6);i++) reth->SetBinContent(i,lval);
  return reth;
};
TH2 *CalculateOneFeeddown2D(AliEffFDContainer *dtef, AliEffFDContainer *mcef) {
  TH1 *htemp = CalculateOneFeeddown(dtef,mcef,kTRUE,kFALSE,0);
  TH3D *hdt = dynamic_cast<TH3D*>(dtef->fetchObj("DCAxy_noChi2"));
  TH2D *rh = (TH2D*)hdt->Project3D("yx");
  rh->SetName("FDvsCent");
  rh->Reset();
  for(Int_t i=1;i<=rh->GetNbinsY();i++) {
    TH1 *htemp = CalculateOneFeeddown(efdt,efmc,kFALSE,kFALSE,i);
    for(Int_t j=1;j<=htemp->GetNbinsX();j++) {
      if(!htemp->GetBinContent(j)) continue;
      rh->SetBinContent(j,i,htemp->GetBinContent(j));
    };
  };
  return rh;
}
TH1 *CalculateEfficiency(AliEffFDContainer *mcef, Int_t centBin=0) {
  TH2D *hrec = dynamic_cast<TH2D*>(mcef->fetchObj("nChRec_Weighted"));
  TH2D *hgen = dynamic_cast<TH2D*>(mcef->fetchObj("nChGen_Weighted"));
  Int_t cLow = centBin?centBin:1;
  Int_t cHigh = centBin?centBin:hrec->GetNbinsY();
  TH1 *reth = hrec->ProjectionX("effres",cLow,cHigh);
  TH1 *den  = hgen->ProjectionX("hgen",cLow,cHigh);
  reth->Divide(den);
  delete den;
  return reth;
};
TH2 *CalculateEfficiency2D(AliEffFDContainer *mcef) {
  TH2D *hrec = dynamic_cast<TH2D*>(mcef->fetchObj("nChRec_Weighted"));
  TH2D *hgen = dynamic_cast<TH2D*>(mcef->fetchObj("nChGen_Weighted"));
  hrec = (TH2D*)hrec->Clone("eff2D");
  hrec->Divide(hgen);
  return hrec;
};
void WriteOneSet(TString outFile, TString outName) {
  Init();
  TH2 *fd2D = CalculateOneFeeddown2D(efdt,efmc);
  TH1 *fdMB = CalculateOneFeeddown(efdt,efmc,kFALSE,kFALSE,0);
  TH2 *ef2D = CalculateEfficiency2D(efmc);
  TH1 *efMB = CalculateEfficiency(efmc,0);
  fd2D->SetName(Form("FDvsCent_%s",outName.Data()));
  fdMB->SetName(Form("FD_MB_%s",outName.Data()));
  ef2D->SetName(Form("EFvsCent_%s",outName.Data()));
  efMB->SetName(Form("EF_MB_%s",outName.Data()));
  delete efdt;
  delete efmc;
  // delete fd;
  TFile *tf = new TFile(outFile.Data(),"UPDATE");
  fd2D->Write(0,TObject::kOverwrite);
  fdMB->Write(0,TObject::kOverwrite);
  ef2D->Write(0,TObject::kOverwrite);
  efMB->Write(0,TObject::kOverwrite);
  tf->Close();
};
void WriteOneConfig(TString namePF) {
  for(Int_t i=0;i<8;i++) {
    ConfigObjName(namePF,etaPFs[i]);
    WriteOneSet(Form("CombinedOutput/%s.root",namePF.Data()),Form("%s%s",i<4?"EtaNeg":"EtaPos",etaPFs[i].Data()));
  }
}

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
/*void CalculateAllFeeddowns(TString dataFileName, TString mcFileName, TString outputFileName) {
  TFile *dtFile = getFile(dataFileName);
  TFile *mcFile = getFile(mcFileName);
  TList *lok = dtFile->GetListOfKeys();
  TList *outList = new TList();
  outList->SetName("FDList");
  outList->SetOwner(kTRUE);
  TList *effList = new TList();
  effList->SetName("EffAndFD");
  for(Int_t i=0;i<lok->GetEntries();i++) {
    if(i>0) continue;
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
}*/
/*void ReplaceAxis(TH1 **target, TH1 *source) {
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
}*/
/*void CalculateAllEfficiencies(TString fdFileName, TString mcFileName, TString outputFileName) {
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
}*/

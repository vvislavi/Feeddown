#include "FD.h"
using namespace RooFit;
FD::FD():
  TNamed("",""),
  fDCAData(0),
  fDCAMC(0),
  fWithinDCAData(0),
  fWithinDCAMC(0)
{
};
FD::~FD() {
  delete fDCAData;
  delete fDCAMC;
  delete fWithinDCAData;
  delete fWithinDCAMC;
};
FD::FD(TH3D *DCAData, TH2D *withinDCAData, TH3D *DCAMC, TH2D *withinDCAMC):
  TNamed("",""),
  fDCAData(0),
  fDCAMC(0),
  fWithinDCAData(0),
  fWithinDCAMC(0)
{
  Initialize(DCAData,withinDCAData,DCAMC,withinDCAMC);
};
void FD::Initialize(TH3D *DCAData, TH2D *withinDCAData, TH3D *DCAMC, TH2D *withinDCAMC) {
  if(DCAData) { delete fDCAData; fDCAData = (TH3D*)DCAData->Clone(Form("%s_DCAData",this->GetName())); fDCAData->SetDirectory(0); };
  if(DCAMC) { delete fDCAMC; fDCAMC = (TH3D*)DCAMC->Clone(Form("%s_DCAMC",this->GetName())); fDCAMC->SetDirectory(0); };
  if(withinDCAData) { delete fWithinDCAData; fWithinDCAData = (TH2D*)withinDCAData->Clone(Form("%s_withinDCAData",this->GetName())); fWithinDCAData->SetDirectory(0); };
  if(withinDCAMC) { delete fWithinDCAMC; fWithinDCAMC = (TH2D*)withinDCAMC->Clone(Form("%s_withinDCAMC",this->GetName())); fWithinDCAMC->SetDirectory(0); };
};
TH1 **FD::getMCProjections(Int_t PtBinNo) {
  if(!fDCAMC) {printf("MC template not set!\n"); return 0; };
  TH1 **hRet = new TH1*[fDCAMC->GetNbinsY()];
  fDCAMC->GetXaxis()->SetRange(PtBinNo,PtBinNo);
  for(Int_t i=0;i<fDCAMC->GetNbinsY();i++) {
    fDCAMC->GetYaxis()->SetRange(i+1,i+1);
    hRet[i] = (TH1*)fDCAMC->Project3D("z")->Clone(Form("Proj_%i",i));
    hRet[i]->SetDirectory(0);
    Integerize(hRet[i]);
  }
  fDCAMC->GetXaxis()->UnZoom();
  fDCAMC->GetYaxis()->UnZoom();
  return hRet;
};
TH1 *FD::getDataProjection(Int_t PtBinNo) {
  if(!fDCAData) {printf("MC template not set!\n"); return 0; };
  TH1 *hRet=0;
  fDCAData->GetXaxis()->SetRange(PtBinNo,PtBinNo);
  fDCAData->GetYaxis()->SetRange(1,1);
  hRet = (TH1*)fDCAData->Project3D("z")->Clone("ProjData");
  hRet->SetDirectory(0);
  fDCAData->GetXaxis()->UnZoom();
  fDCAData->GetYaxis()->UnZoom();
  return hRet;
};
void FD::Integerize(TH1 *inh) {
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bconp = inh->GetBinContent(i);
    Double_t bcon = TMath::Floor(bconp);
    inh->SetBinContent(i,bcon);
    inh->SetBinError(i,TMath::Sqrt(bcon));
  }
}
Double_t FD::FitOneBin(Int_t ptBin, TString outFile, Double_t *outFracs) {
  printf("Fitting bin %i\n\n\n\n\n\n\n\n\n",ptBin);
  RooRealVar x("x", "x", -3, 3);
  TH1 *l_signal = getDataProjection(ptBin);
  TH1 **l_temps = getMCProjections(ptBin);

  RooDataHist dsig("dsig", "dsig", x, Import(*l_signal));
  RooDataHist mcprim("mcprim", "mcprim", x, Import(*l_temps[0]));
  RooDataHist mcsec("mcsec", "mcsec", x, Import(*l_temps[1]));
  RooDataHist mcmat("mcmat", "mcmat", x, Import(*l_temps[2]));

  RooHistFunc f_prim("f_prim","f_prim",x,mcprim);
  RooHistFunc f_sec("f_sec","f_sec",x,mcsec);
  RooHistFunc f_mat("f_mat","f_mat",x,mcmat);

  RooRealVar Aprim("Aprim","Aprim",0.7,0.0,100);
  RooRealVar Asec("Asec","Asec",0.2,0.0,100);
  RooRealVar Amat("Amat","Amat",0.1,0,100);

  RooRealSumPdf model0("model0","model0",
     RooArgList(f_prim,f_sec,f_mat),
     RooArgList(Aprim, Asec, Amat),
     true);

  auto result0 = model0.fitTo(dsig, PrintLevel(0), Save());

  if(!outFile.IsNull()) {
    TCanvas* can = new TCanvas("can", "", 600, 600);

    TPaveText pt(-19.5, 1, -2, 25);
    pt.SetFillStyle(0);
    pt.SetBorderSize(0);

    auto frame = x.frame(Title("No template uncertainties"));
    dsig.plotOn(frame);
    model0.plotOn(frame, LineColor(kGreen+2));
    model0.plotOn(frame, Components(f_prim), LineColor(kAzure));
    model0.plotOn(frame, Components(f_sec), LineColor(kRed));
    frame->Draw();
    can->Print(outFile.Data());
  };
  if(!fWithinDCAMC) {printf("fWithinDCAMC not defined, will not return a value!\n"); return 0; };
  Double_t primV = fWithinDCAMC->GetBinContent(ptBin,1);//l_temps[0]->Integral();//
  Double_t secwV = fWithinDCAMC->GetBinContent(ptBin,2);//l_temps[1]->Integral();//
  Double_t secmV = fWithinDCAMC->GetBinContent(ptBin,3);//l_temps[2]->Integral();//
  delete [] l_temps;
  delete l_signal;
  primV*=Aprim.getVal();
  secwV*=Asec.getVal();
  secmV*=Amat.getVal();
  if(outFracs) { outFracs[0] = Aprim.getVal(); outFracs[1] = Asec.getVal(); outFracs[2] = Amat.getVal(); };
  return primV/(primV+secwV+secmV);
}
TH1 *FD::PerformFit(Double_t ptMin, Double_t ptMax, TString outFile, TH1 **outHists) {
  TH1 *retH = fDCAData->Project3D("x");
  retH->SetName("PrimOverAll");
  retH->Reset();
  if(outHists) {
    outHists[0] = (TH1*)retH->Clone("primaries");
    outHists[1] = (TH1*)retH->Clone("secondaries");
    outHists[2] = (TH1*)retH->Clone("material");
  }
  Double_t *outFracs = outHists?(new Double_t[3]):0;
  for(Int_t i=retH->FindBin(ptMin+(1e-6)); i<=retH->FindBin(ptMax-1e-6); i++) {
    Double_t nval = FitOneBin(i,outFile,outFracs);
    retH->SetBinContent(i,nval);
    if(outHists) for(Int_t j=0;j<3;j++) outHists[j]->SetBinContent(i,outFracs[j]);
  }
  return retH;
}

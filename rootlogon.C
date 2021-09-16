{
  gROOT->ProcessLine(".x Scripts/LoadLibs.C");
  TFile *tf1 = new TFile("../AliceData/Efficiency_LHC20j6a_v3/LHC20j6a.root","READ");
  TList *tl1 = (TList*)tf1->Get("SpectraList_1");
  TH3D *hmc = (TH3D*)tl1->FindObject("DCAxy_noChi2");
  TH2D *hwdc = (TH2D*)tl1->FindObject("WithinDCA_withChi2");
  // TH2D *thmc[2];
  // thmc[0] = (TH2D*)tl1->FindObject("WithinDCA_withChi2");
  // thmc[1] = (TH2D*)tl1->FindObject("WithinDCA_noChi2");
  TFile *tf2 = new TFile("../AliceData/LHC15o_pass2_DCAxy_6Runs_v3/merged.root");
  TList *tl2 = (TList*)tf2->Get("SpectraList_1");
  TH3D *hdt = (TH3D*)tl2->FindObject("DCAxy_noChi2");
  // TH2D *thdt[2];
  // thdt[0] = (TH2D*)tl2->FindObject("WithinDCA_withChi2");
  // thdt[1] = (TH2D*)tl2->FindObject("WithinDCA_noChi2");
  // TH1D *hall[3];
  // TFile *tfin = new TFile("fracs.root","READ");
  // hall[0] = (TH1D*)tfin->Get("prim");
  // hall[1] = (TH1D*)tfin->Get("sec");
  // hall[2] = (TH1D*)tfin->Get("mat");
  // hall[3] = (TH1D*)tfin->Get("Signal");
  // TObjArray *mcar = new TObjArray();
  // for(Int_t i=0;i<3;i++) mcar->Add(hall[i]);
  // gROOT->ProcessLine(".L Scripts/FeedDown.cxx+");
  // FeedDown *fd = new FeedDown(hdt,hmc,hwdc);
}

#ifndef FD__H
#define FD__H
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TNamed.h"
#include "TMath.h"


class FD:public TNamed {
  public:
    FD();
    ~FD();
    FD(TH3D *DCAdata, TH2D *withinDCAdata, TH3D *fDCAMC, TH2D *withinDCAMC);
    TH1 **getMCProjections(Int_t PtBinNo);
    TH1 *getDataProjection(Int_t PtBinNo, Int_t centBin);
    void Initialize(TH3D *DCAData, TH2D *withinDCAData, TH3D *fDCAMC, TH2D *withinDCAMC);
    Double_t FitOneBin(Int_t ptBin, TString outFile="", Double_t *outFracs=0, Int_t centBin=0);
    TH1 *PerformFit(Double_t ptMin=0.2, Double_t ptMax=5, TString outFile="", TH1 **outFracs=0, Int_t centBin=0);
    void NormalizeByBinSize(TH1 *inh);
    TH1 *rebinDCA(TH1*);
    ClassDef(FD, 1);
  private:
    TH3D *fDCAData;
    TH3D *fDCAMC;
    TH2D *fWithinDCAData;
    TH2D *fWithinDCAMC;
    void Integerize(TH1 *inh);
};
#endif

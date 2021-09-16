#include "TAxis.h"
#include "TGaxis.h"
#include "TPad.h"

void ProcessAxisPx(TAxis *ax, Int_t ls, Double_t lof, Int_t ts, Double_t tof, Int_t nDiv, Double_t ltick) {
  ax->SetNdivisions(nDiv);
  ax->SetTickLength(ltick);
  ax->SetTitleFont(43);
  ax->SetLabelFont(43);
  ax->SetLabelSize(ls);
  ax->SetLabelOffset(lof);
  ax->SetTitleSize(ts);
  ax->SetTitleOffset(tof);
};


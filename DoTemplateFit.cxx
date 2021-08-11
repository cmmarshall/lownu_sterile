#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"

int main()
{
  TFile *f = new TFile("output.root","READ");
  TH2D* h = (TH2D*)f->Get("hElepVsEv");
  TH1D * templates[10];
  for(int i=0; i<10; i++) { 
    templates[i] = (TH1D*)h->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  }

  TH1D * intrinsic = (TH1D*)f->Get("hLepE_sm");
  TH1D * target = (TH1D*)f->Get("hElep_w");

  TemplateFitter tf( templates, intrinsic, target );
  double energy_bins[11];
  for( int b = 0; b <= 10; ++b ) {
    energy_bins[b] = h->GetXaxis()->GetBinLowEdge(b+1);
  }
  tf.setEnergyBins( energy_bins );

  double bf_dm2, bf_theta;
  bool isOK = tf.doFit( bf_theta, bf_dm2 );

  printf( "Best-fit theta = %f, dm2 = %f\n", bf_theta, bf_dm2 );

}







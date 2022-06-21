#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

int main()
{
  char var[20] = "ElepReco";
  int nuCut, EvCut;
  nuCut = 3;
  EvCut = 0;


  TFile *CC_f  = new TFile("/dune/app/users/qvuong/lownu_sterile/gen_data/new_tgt/output.root","READ");
  TFile *nue_f = new TFile("/dune/app/users/qvuong/lownu_sterile/gen_data/new_tgt/nue_output.root","READ");

  TH2D* CC_hm    = (TH2D*)CC_f->Get(Form("m_h%sVsEv%d",var,nuCut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_m_h%sVsEv%d",var,nuCut));
  TH2D* CC_he    = (TH2D*)CC_f->Get(Form("e_h%sVsEv%d",var,nuCut));
  TH2D* nue_hm   = (TH2D*)nue_f->Get(Form("m_h%sVsEv%d",var,EvCut));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("m_h%sVsEv%d_w",var,EvCut));
  TH2D* nue_he   = (TH2D*)nue_f->Get(Form("e_h%sVsEv%d",var,EvCut));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("e_h%sVsEv%d_w",var,EvCut));

  TH1D * CC_templates_m[480];
  TH1D * CC_templates_m_nc[480];
  TH1D * CC_templates_e[480];
  TH1D * nue_templates_m[480];
  TH1D * nue_templates_m_w[480];
  TH1D * nue_templates_e[480];
  TH1D * nue_templates_e_w[480];

  for(int i=0; i<480; i++) {
    CC_templates_m[i]    = (TH1D*)CC_hm->ProjectionY(Form("CC_m_bin%d",i+1),i+1,i+1);
    CC_templates_m_nc[i] = (TH1D*)CC_hm_nc->ProjectionY(Form("CC_nc_m_bin%d",i+1),i+1,i+1);
    CC_templates_e[i]    = (TH1D*)CC_he->ProjectionY(Form("CC_e_bin%d",i+1),i+1,i+1);
    nue_templates_m[i]    = (TH1D*)nue_hm->ProjectionY(Form("nue_m_bin%d",i+1),i+1,i+1);
    nue_templates_m_w[i]  = (TH1D*)nue_hm_w->ProjectionY(Form("nue_w_m_bin%d",i+1),i+1,i+1);
    nue_templates_e[i]    = (TH1D*)nue_he->ProjectionY(Form("nue_e_bin%d",i+1),i+1,i+1);
    nue_templates_e_w[i]  = (TH1D*)nue_he_w->ProjectionY(Form("nue_w_e_bin%d",i+1),i+1,i+1);
  }

  TemplateFitter tf( CC_templates_m, CC_templates_m_nc, CC_templates_e, nue_templates_m, nue_templates_m_w, nue_templates_e, nue_templates_e_w );

  double energy_bins[481];
  for( int b = 0; b <= 480; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  } 
  tf.setEnergyBins( energy_bins );


  double oscpar[3], oscpar_max[3], oscpar_min[3], stepsize[3], seed[3];
  int N = 1;

  oscpar_max[0] = 0.04; oscpar_min[0] = 0.01;	//Uee2 range
  oscpar_max[1] = 0.01; oscpar_min[1] = 0.01;	//Umm2 range
  oscpar_max[2] = 10.0; oscpar_min[2] = 1.0;	//dm2 range

  double s2mue2_max, s2mue2_min;  		//sin(2theta_mue)^2 range
  s2mue2_min = (oscpar_min[0]*oscpar_min[0]) * (oscpar_min[1]*oscpar_min[1]);
  s2mue2_max = (oscpar_max[0]*oscpar_max[0])  * (oscpar_max[1]*oscpar_max[1]);

  for(int j = 0; j < 3; j++) {
    stepsize[j] = (oscpar_max[j] - oscpar_min[j])/N;
  }

  TH2D *h = new TH2D("h","",N+1,s2mue2_min,s2mue2_max,N+1,oscpar_min[2],oscpar_max[2]);

/*
  for(int i2 = 0; i2 < N+1; i2++) {
    oscpar[2] = oscpar_min[2] + i2*stepsize[2];
*/

  {
    oscpar[2] = oscpar_min[2];				//choose dm2 = 1.0
    oscpar[1] = oscpar_min[1];				//choose Umm2 = 0.01
    for(int i0 = 0; i0 < N+1; i0++) {
      oscpar[0] = oscpar_min[0] + i0*stepsize[0];	//vary Uee2 from 0.01 to 0.04


/* 
// Putting in the oscillation parameter separately!!
  {
    {
      oscpar[0] = 0.01;
      oscpar[1] = 0.01;
      oscpar[2] = 1.0;
*/
    
      std::cout << oscpar[0] << "\t" << oscpar[1] << "\t" << oscpar[2] << "\n";

      for(int j = 0; j < 3; j++) {
        seed[j] = oscpar[j]; 
      }
      std::cout << seed[0] << "\t" << seed[1] << "\t" << seed[2] << "\n";

      double s2mue2 = (oscpar[0]*oscpar[0]) * (oscpar[1]*oscpar[1]);
  
      tf.setPara( var, oscpar, nuCut, EvCut, seed );

      tf.getTarget( oscpar );

      double bf_dm2, bf_Uee2, bf_Umm2;
      double par[3];
      bool isOK = tf.doFit( bf_Uee2, bf_Umm2 , bf_dm2);
      par[0] = bf_Uee2;
      par[1] = bf_Umm2;
      par[2] = bf_dm2;
      double chi2 = tf.bfChi2( par  );
      printf( "nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2, chi2 );

      h->Fill(s2mue2,oscpar[2]);
    }
  }
}




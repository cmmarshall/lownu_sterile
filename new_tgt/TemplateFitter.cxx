#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TRandom3.h>

TemplateFitter::TemplateFitter(TH1D * CC_templates_m[480], TH1D * CC_templates_m_nc[480], TH1D * CC_templates_e[480], TH1D * nue_templates_m[480], TH1D * nue_templates_m_w[480], TH1D * nue_templates_e[480], TH1D * nue_templates_e_w[480] )
{
  for( int i = 0; i < 480; ++i ) {
  CC_m_templates[i] = CC_templates_m[i];
  CC_nc_m_templates[i] = CC_templates_m_nc[i];
  CC_e_templates[i] = CC_templates_e[i];
  nue_m_templates[i] = nue_templates_m[i];
  nue_w_m_templates[i] = nue_templates_m_w[i];
  nue_e_templates[i] = nue_templates_e[i];
  nue_w_e_templates[i] = nue_templates_e_w[i];
  }

  //CC_m_target = CC_target_m;
  //CC_e_target = CC_target_e;
  //target_nue = nue_target;

}

void TemplateFitter::setEnergyBins( double bins[481] )
{
  for( int i = 0; i < 481; ++ i ) {m_energy_bins[i] = bins[i];}
}


double TemplateFitter::getPmue( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mue2 = 4 * Uee2 * Umm2;
  double prob = s2mue2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPee( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2ee2 = 4 * Uee2 * (1 - Uee2);
  double prob = 1.0 - s2ee2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPmm( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mm2 = 4 * Umm2 * (1 - Umm2);
  double prob = 1.0 - s2mm2  * pow(sin(del),2);
  return prob;
}

double b0,b1,b2;
double s0,s1,s2;
void TemplateFitter::setPara( char var[20], double oscpar[3], int nuCut, int EvCut, double seed[3] )
{
  name  = var;
  cutNu = nuCut;
  cutEv = EvCut;
  for(int i = 0; i < 3; i++){
    ospar[i] = oscpar[i];
    std::cout << ospar[i];
  }
  b0 = oscpar[0];
  b1 = oscpar[1];
  b2 = oscpar[2];

  s0 = seed[0];
  s1 = seed[1];
  s2 = seed[2];
}

const int N = 40;
TH2D *h0  = new TH2D("h0","",N,0,0.1, N,0,12.0);
TH2D *h1  = new TH2D("h1","",N,0,12.0,N,0,0.1);
TH2D *h2  = new TH2D("h2","",N,0,0.1, N,0,0.1);
TH2D *h0z = new TH2D("h0z","",N,0.0*b1,2.0*b1,N,0.5*b2,1.5*b2);
TH2D *h1z = new TH2D("h1z","",N,0.5*b2,1.5*b2,N,0.0*b0,2.0*b0);
TH2D *h2z = new TH2D("h2z","",N,0.0*b0,2.0*b0,N,0.0*b1,2.0*b1);

TH1D * CC_e_target = new TH1D("CC_e_target","",100,0,16);
TH1D * CC_m_target = new TH1D("CC_m_target","",100,0,16);
TH1D * target_nue = new TH1D("target_nue","",100,0,16);

void TemplateFitter::getTarget( double *oscpar )
{
  TRandom3 *rando = new TRandom3(8888); 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();
/*
  TH1D * CC_e_target = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_m_target = (TH1D*) CC_tp_e->Clone();
  TH1D * target_nue = (TH1D*) CC_tp_e->Clone();
*/
  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, ospar[0], ospar[1], ospar[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, ospar[0], ospar[1], ospar[2]);
      mm  = mm  + getPmm(e, ospar[0], ospar[1], ospar[2]);
    }
    double Pmue = mue/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);


  for(int bx = 1; bx <= CC_tp_e->GetNbinsX(); bx++){
    double mean = CC_tp_e->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CC_e_target->AddBinContent(bx, fl_bc);
    //std::cout << CC_tp_e->GetBinContent(bx) << "\t" << CC_e_target->GetBinContent(bx) << "\n";
  }

  for(int bx = 1; bx <= CC_tp_m->GetNbinsX(); bx++){
    double mean = CC_tp_m->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CC_m_target->AddBinContent(bx, fl_bc);
    //std::cout << CC_tp_m->GetBinContent(bx) << "\t" << CC_m_target->GetBinContent(bx) << "\n";
  }

  for(int bx =  1; bx <= nue_tp->GetNbinsX(); bx++){
    double mean = nue_tp->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    target_nue->AddBinContent(bx, fl_bc);
    //std::cout << nue_tp->GetBinContent(bx) << "\t" << target_nue->GetBinContent(bx) << "\n";
  }

}

// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( const double * par )
{
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    }
    double Pmue = mue/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);
/*
  for(int bx = 1; bx <= CC_tp_e->GetNbinsX(); bx++){
    double mean = CC_tp_e->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CC_e_target->AddBinContent(bx, fl_bc);
    std::cout << CC_tp_e->GetBinContent(bx) << "\t" << CC_e_target->GetBinContent(bx) << "\n";

  }

  //std::cout << CC_tp_m->GetNbinsX() << "\t" << CC_m_target->GetNbinsX() << "\n";
  //std::cout << nue_tp->GetNbinsX() << "\t" << target_nue->GetNbinsX() << "\n";

  for(int bx = 1; bx <= CC_tp_m->GetNbinsX(); bx++){
    double mean = CC_tp_m->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CC_m_target->AddBinContent(bx, fl_bc);
    std::cout << CC_tp_m->GetBinContent(bx) << "\t" << CC_m_target->GetBinContent(bx) << "\n";
  }

  for(int bx =  1; bx <= nue_tp->GetNbinsX(); bx++){
    double mean = nue_tp->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    target_nue->AddBinContent(bx, fl_bc);
    std::cout << nue_tp->GetBinContent(bx) << "\t" << target_nue->GetBinContent(bx) << "\n";
  }
*/
  // calculate the chi2 with the "data" target
  const int nbins_E = 100;
  double chi2 = 0.0;

  TMatrixD target(3*nbins_E, 1);
  TMatrixD temp(3*nbins_E, 1);
  TMatrixD diff(3*nbins_E, 1);
  TMatrixD diff_T(1, 3*nbins_E);
  TMatrixD unc(3*nbins_E, 3*nbins_E);
  TMatrixD cov(3*nbins_E, 3*nbins_E);
  TMatrixD covmtr_tot(3*nbins_E, 3*nbins_E);

  for( int bx = 0; bx < nbins_E; bx++ ) {
    temp[bx][0] = CC_tp_m->GetBinContent(bx+1); 
    temp[bx+nbins_E][0] = CC_tp_e->GetBinContent(bx+1);
    temp[bx+2*nbins_E][0] = nue_tp->GetBinContent(bx+1);

    target[bx][0] = CC_m_target->GetBinContent(bx+1); 
    target[bx+nbins_E][0] = CC_e_target->GetBinContent(bx+1);
    target[bx+2*nbins_E][0] = target_nue->GetBinContent(bx+1);

    diff[bx][0] = temp[bx][0] - target[bx][0];
    diff[bx+nbins_E][0] = temp[bx+nbins_E][0] - target[bx+nbins_E][0];
    diff[bx+2*nbins_E][0] = temp[bx+2*nbins_E][0] - target[bx+2*nbins_E][0];
 
    diff_T[0][bx] = diff[bx][0]; 
    diff_T[0][bx+nbins_E] = diff[bx+nbins_E][0];
    diff_T[0][bx+2*nbins_E] = diff[bx+2*nbins_E][0];

    for( int by = 0; by < nbins_E; by++ ) {
      unc[bx][by] 	       = 0;
      unc[bx][by+nbins_E]   = 0;
      unc[bx][by+2*nbins_E] = 0;

      unc[bx+nbins_E][by]              = 0;
      unc[bx+nbins_E][by+nbins_E]   = 0;
      unc[bx+nbins_E][by+2*nbins_E] = 0;

      unc[bx+2*nbins_E][by]              = 0;
      unc[bx+2*nbins_E][by+nbins_E]   = 0;
      unc[bx+2*nbins_E][by+2*nbins_E] = 0;

      if(bx == by) {
        unc[bx][by] 			      = CC_m_target->GetBinContent(bx+1);
        unc[bx+nbins_E][by+nbins_E]     = CC_e_target->GetBinContent(bx+1);
        unc[bx+2*nbins_E][by+2*nbins_E] = target_nue->GetBinContent(bx+1);
      }
    }
  }

  covmtr_tot = unc;
  TDecompSVD svd(covmtr_tot);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, 3*nbins_E);
  for( int bx = 0; bx < 3*nbins_E; bx++ ) {
    for( int by = 0; by < 3*nbins_E; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";

  h0->Fill(par[1], par[2]);
  h1->Fill(par[2], par[0]);
  h2->Fill(par[0], par[1]);

  h0z->Fill(par[1], par[2]);
  h1z->Fill(par[2], par[0]);
  h2z->Fill(par[0], par[1]);

  return chi2;
}

double bf0[1], bf1[1], bf2[1];

bool TemplateFitter::doFit( double &Uee2, double &Umm2, double &dm2 )
{
/*
  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed1[1][0] = 0.001;
  seed1[1][1] = 0.001;
  seed1[1][2] = 0.1;
  seed1[2][0] = 0.1;
  seed1[2][1] = 0.1;
  seed1[2][2] = 10.0;

  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  seed2[1][0] = 0.001;
  seed2[1][1] = 0.001;
  seed2[1][2] = 0.1;
  seed2[2][0] = 0.1;
  seed2[2][1] = 0.1;
  seed2[2][2] = 10.0;

  double p0, p1, p2;
  if(para == 1) {
    p0 = seed1[s][0];
    p1 = seed1[s][1];
    p2 = seed1[s][2];
  }
  else if(para == 2) {
    p0 = seed2[s][0];
    p1 = seed2[s][1];
    p2 = seed2[s][2];
  }
  std::cout << "seed = " << s << "\t" << p0 << "\t" << p1 << "\t" << p2 << "\n";
*/
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(1000000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(1000000);
  fitter->SetTolerance(0.1); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  fitter->SetVariable( 0, "Uee2", s0, 0.00001 );
  fitter->SetVariable( 1, "Umm2", s1, 0.00001 );
  fitter->SetVariable( 2, "dm2",  s2, 0.01 );
  fitter->SetVariableLowerLimit(0, 0.0);
  fitter->SetVariableLowerLimit(1, 0.0);
  fitter->SetVariableLowerLimit(2, 0.0);

  // 3 free parameters = theta, dm2
  ROOT::Math::Functor lf( this, &TemplateFitter::getChi2, 3 );
  ROOT::Math::Functor functor( lf, 3 );
  fitter->SetFunction( functor );

  // Go!
  fitter->Minimize();

/*
  if( fitter->Status() != 0 ) {
    std::cout << "Something bad happened" << std::endl;
    return false;
  }
*/

  const double *bestfit = fitter->X();
  Uee2 = bestfit[0];
  Umm2 = bestfit[1];
  dm2 = bestfit[2];
  bf0[0] = bestfit[0];
  bf1[0] = bestfit[1];
  bf2[0] = bestfit[2];
  
  double chi2 = fitter->MinValue();
  //h->Fill(Uee2, Umm2, dm2);
  return true;

}
/*
void TemplateFitter::Draw()
{

  double par[3];
  par[0] = bf0[0];
  par[1] = bf1[0];
  par[2] = bf2[0];

  std::cout << "best_fit: " << par[0] << "\t" << par[1] << "\t" << par[2] << "\n";

  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    }
    double Pmue = mue/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CC_e_target->SetName("e_target");
  CC_e_target->SetLineColor(kBlack);
  CC_m_target->SetName("m_target");
  CC_m_target->SetLineColor(kBlack);
  target_nue->SetName("nue_target");
  target_nue->SetLineColor(kBlack);

  CC_tp_e->SetName("tp_e");
  CC_tp_e->SetMarkerStyle(21);
  CC_tp_e->SetMarkerColor(kRed);
  CC_tp_e->SetMarkerSize(0.5);
  CC_tp_me->SetName("tp_me");
  CC_tp_me->SetMarkerStyle(21);
  CC_tp_me->SetMarkerColor(kBlue);
  CC_tp_me->SetMarkerSize(0.5);

  CC_tp_m->SetName("tp_m");
  CC_tp_m->SetMarkerStyle(21);
  CC_tp_m->SetMarkerColor(kRed);
  CC_tp_m->SetMarkerSize(0.5);
  CC_tp_em->SetName("tp_em");
  CC_tp_em->SetMarkerStyle(21);
  CC_tp_em->SetMarkerColor(kBlue);
  CC_tp_em->SetMarkerSize(0.5);

  nue_tp->SetName("tp_nue");
  nue_tp->SetMarkerStyle(21);
  nue_tp->SetMarkerColor(kRed);
  nue_tp->SetMarkerSize(0.5);
  nue_tp_os->SetName("tp_nue_os");
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kBlue);
  nue_tp_os->SetMarkerSize(0.5);

  TCanvas *c_e = new TCanvas("c_e","",900,700);
  CC_e_target->Draw();
  CC_tp_me->Draw("same");
  CC_tp_e->Draw("same");
  TLegend *legend_e = new TLegend(0.65,0.70,0.9,0.9);
  legend_e->AddEntry(CC_e_target,"fluctuated target");
  legend_e->AddEntry(CC_tp_me,"oscillated mu->e");
  legend_e->AddEntry(CC_tp_e,"templates at bestfit");
  legend_e->Draw();
  c_e->SaveAs(Form("%s_fit_CC_e_nCov_%d%d%d_%d_ft.png",name,para,cutNu,cutEv,s));

  TCanvas *c_m = new TCanvas("c_m","",900,700);
  CC_m_target->Draw();
  CC_tp_em->Draw("same");
  CC_tp_m->Draw("same");
  TLegend *legend_m = new TLegend(0.65,0.70,0.9,0.9);
  legend_m->AddEntry(CC_m_target,"fluctuated target");
  legend_m->AddEntry(CC_tp_em,"oscillated e->mu");
  legend_m->AddEntry(CC_tp_m,"templates at bestfit");
  legend_m->Draw();
  c_m->SaveAs(Form("%s_fit_CC_m_nCov_%d%d%d_%d_ft.png",name,para,cutNu,cutEv,s));

  TCanvas *c = new TCanvas("c","",900,700);
  gPad->SetLogy();
  target_nue->Draw();
  nue_tp_os->Draw("same");
  nue_tp->Draw("same");
  TLegend *legend_nue = new TLegend(0.65,0.78,0.9,0.9);
  legend_nue->AddEntry(target_nue,"fluctuated target");
  legend_nue->AddEntry(nue_tp_os,"oscillated mu->e & e->mu");
  legend_nue->AddEntry(nue_tp,"templates at bestfit");
  legend_nue->Draw();
  c->SaveAs(Form("%s_fit_nue_nCov_%d%d%d_%d_ft_Log.png",name,para,cutNu,cutEv,s));

  TFile *f = new TFile(Form("/dune/app/users/qvuong/lownu/Elep_combine/chi2/%s_chi2_nCov_%d%d%d.root",name,para,cutNu,cutEv));
  std::cout << 1 << "\n";
  TH2D *hc0 = (TH2D*)f->Get("h0");
  TH2D *hc1 = (TH2D*)f->Get("h1");
  TH2D *hc2 = (TH2D*)f->Get("h2");
  TH2D *hc0z = (TH2D*)f->Get("h0z");
  TH2D *hc1z = (TH2D*)f->Get("h1z");
  TH2D *hc2z = (TH2D*)f->Get("h2z");
  TH2D *hc0d = (TH2D*)f->Get("h0d");
  TH2D *hc1d = (TH2D*)f->Get("h1d");
  TH2D *hc2d = (TH2D*)f->Get("h2d");
  TH2D *hc0dz = (TH2D*)f->Get("h0dz");
  TH2D *hc1dz = (TH2D*)f->Get("h1dz");
  TH2D *hc2dz = (TH2D*)f->Get("h2dz");

  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed1[1][0] = 0.001;
  seed1[1][1] = 0.001;
  seed1[1][2] = 0.1;
  seed1[2][0] = 0.1;
  seed1[2][1] = 0.1;
  seed1[2][2] = 10.0;

  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  seed2[1][0] = 0.001;
  seed2[1][1] = 0.001;
  seed2[1][2] = 0.1;
  seed2[2][0] = 0.1;
  seed2[2][1] = 0.1;
  seed2[2][2] = 10.0;

  double p0[1], p1[1], p2[1]; //true values
  double s0[1], s1[1], s2[1];
  if(para == 1) {
    p0[0] = seed1[0][0];
    p1[0] = seed1[0][1];
    p2[0] = seed1[0][2];

    s0[0] = seed1[s][0];
    s1[0] = seed1[s][1];
    s2[0] = seed1[s][2];
  }
  else if(para == 2) {
    p0[0] = seed2[0][0];
    p1[0] = seed2[0][1];
    p2[0] = seed2[0][2];

    s0[0] = seed2[s][0];
    s1[0] = seed2[s][1];
    s2[0] = seed2[s][2];
  }

  int nG = 1;
  TGraph *g0   = new TGraph(nG,p1,p2);
  TGraph *g1   = new TGraph(nG,p2,p0);
  TGraph *g2   = new TGraph(nG,p0,p1);
  TGraph *g0s  = new TGraph(nG,s1,s2);
  TGraph *g1s  = new TGraph(nG,s2,s0);
  TGraph *g2s  = new TGraph(nG,s0,s1);
  TGraph *g0BF = new TGraph(nG,bf1,bf2);
  TGraph *g1BF = new TGraph(nG,bf2,bf0);
  TGraph *g2BF = new TGraph(nG,bf0,bf1);

  g0->SetName("g0");
  g0->SetMarkerStyle(29);
  g0->SetMarkerColor(kRed);
  g0->SetMarkerSize(3);
  g1->SetName("g1");
  g1->SetMarkerStyle(29);
  g1->SetMarkerColor(kRed);
  g1->SetMarkerSize(3);
  g2->SetName("g2");
  g2->SetMarkerStyle(29);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerSize(3);

  g0s->SetName("g0s");
  g0s->SetMarkerStyle(3);
  g0s->SetMarkerColor(kGreen);
  g0s->SetMarkerSize(3);
  g1s->SetName("g1s");
  g1s->SetMarkerStyle(3);
  g1s->SetMarkerColor(kGreen);
  g1s->SetMarkerSize(3);
  g2s->SetName("g2s");
  g2s->SetMarkerStyle(3);
  g2s->SetMarkerColor(kGreen);
  g2s->SetMarkerSize(3);

  g0BF->SetName("g0BF");
  g0BF->SetMarkerStyle(8);
  g0BF->SetMarkerColor(kBlue);
  g0BF->SetMarkerSize(1.1);
  g1BF->SetName("g1BF");
  g1BF->SetMarkerStyle(8);
  g1BF->SetMarkerColor(kBlue);
  g1BF->SetMarkerSize(1.1);
  g2BF->SetName("g2BF");
  g2BF->SetMarkerStyle(8);
  g2BF->SetMarkerColor(kBlue);
  g2BF->SetMarkerSize(1.1);

  h0->GetXaxis()->SetTitle("Umm2");
  h0->GetYaxis()->SetTitle("dm2");
  h0->SetStats(0);
  h1->GetXaxis()->SetTitle("dm2");
  h1->GetYaxis()->SetTitle("Uee2");
  h1->SetStats(0);
  h2->GetXaxis()->SetTitle("Uee2");
  h2->GetYaxis()->SetTitle("Umm2");
  h2->SetStats(0);

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetLogz(0);
  hc0->Draw("colz");
  h0->Draw("same");
  g0->Draw("same P");
  g0s->Draw("same P");
  g0BF->Draw("same P");
  TLegend *legend0 = new TLegend(0.65,0.70,0.9,0.9);
  legend0->AddEntry("g0","  True values");
  legend0->AddEntry("g0s","  Initial Seeds");
  legend0->AddEntry("g0BF","  Bestfit values");
  legend0->Draw();
  c0->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d0.png",name,para,cutNu,cutEv,s));
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogz(0);
  hc1->Draw("colz");
  h1->Draw("same");
  g1->Draw("same P");
  g1s->Draw("same P");
  g1BF->Draw("same P");
  TLegend *legend1 = new TLegend(0.65,0.70,0.9,0.9);
  legend1->AddEntry("g1","  True values");
  legend1->AddEntry("g1s","  Initial Seeds");
  legend1->AddEntry("g1BF","  Bestfit values");
  legend1->Draw();
  c1->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d1.png",name,para,cutNu,cutEv,s));
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogz(0);
  hc2->Draw("colz");
  h2->Draw("same");
  g2->Draw("same P");
  g2s->Draw("same P");
  g2BF->Draw("same P");
  TLegend *legend2 = new TLegend(0.65,0.70,0.9,0.9);
  legend2->AddEntry("g2","  True values");
  legend2->AddEntry("g2s","  Initial Seeds");
  legend2->AddEntry("g2BF","  Bestfit values");
  legend2->Draw();
  c2->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d2.png",name,para,cutNu,cutEv,s));
  TCanvas *c0z = new TCanvas("c0z","",800,600);
  c0z->SetLogz(0);
  hc0z->Draw("colz");
  h0z->Draw("same");
  g0->Draw("same P");
  g0s->Draw("same P");
  g0BF->Draw("same P");
  TLegend *legend0z = new TLegend(0.65,0.70,0.9,0.9);
  legend0z->AddEntry("g0","  True values");
  legend0z->AddEntry("g0s","  Initial Seeds");
  legend0z->AddEntry("g0BF","  Bestfit values");
  legend0z->Draw();
  c0z->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d0_zoom.png",name,para,cutNu,cutEv,s));
  TCanvas *c1z = new TCanvas("c1z","",800,600);
  c1z->SetLogz(0);
  hc1z->Draw("colz");
  h1z->Draw("same");
  g1->Draw("same P");
  g1s->Draw("same P");
  g1BF->Draw("same P");
  TLegend *legend1z = new TLegend(0.65,0.70,0.9,0.9);
  legend1z->AddEntry("g1","  True values");
  legend1z->AddEntry("g1s","  Initial Seeds");
  legend1z->AddEntry("g1BF","  Bestfit values");
  legend1z->Draw();
  c1z->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d1_zoom.png",name,para,cutNu,cutEv,s));
  TCanvas *c2z = new TCanvas("c2z","",800,600);
  c2z->SetLogz(0);
  hc2z->Draw("colz");
  h2z->Draw("same");
  g2->Draw("same P");
  g2s->Draw("same P");
  g2BF->Draw("same P");
  TLegend *legend2z = new TLegend(0.65,0.70,0.9,0.9);
  legend2z->AddEntry("g2","  True values");
  legend2z->AddEntry("g2s","  Initial Seeds");
  legend2z->AddEntry("g2BF","  Bestfit values");
  legend2z->Draw();
  c2z->SaveAs(Form("%s_fit_parDraw_nCov_%d%d%d_%d2_zoom.png",name,para,cutNu,cutEv,s));

  TCanvas *c0d = new TCanvas("c0d","",800,600);
  c0d->SetLogz(0);
  hc0d->Draw("colz");
  h0->Draw("same");
  g0->Draw("same P");
  g0s->Draw("same P");
  g0BF->Draw("same P");
  TLegend *legend0d = new TLegend(0.65,0.70,0.9,0.9);
  legend0d->AddEntry("g0","  True values");
  legend0d->AddEntry("g0s","  Initial Seeds");
  legend0d->AddEntry("g0BF","  Bestfit values");
  legend0d->Draw();
  c0d->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d0.png",name,para,cutNu,cutEv,s));
  TCanvas *c1d = new TCanvas("c1d","",800,600);
  c1d->SetLogz(0);
  hc1d->Draw("colz");
  h1->Draw("same");
  g1->Draw("same P");
  g1s->Draw("same P");
  g1BF->Draw("same P");
  TLegend *legend1d = new TLegend(0.65,0.70,0.9,0.9);
  legend1d->AddEntry("g1","  True values");
  legend1d->AddEntry("g1s","  Initial Seeds");
  legend1d->AddEntry("g1BF","  Bestfit values");
  legend1d->Draw();
  c1d->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d1.png",name,para,cutNu,cutEv,s));
  TCanvas *c2d = new TCanvas("c2d","",800,600);
  c2d->SetLogz(0);
  hc2d->Draw("colz");
  h2->Draw("same");
  g2->Draw("same P");
  g2s->Draw("same P");
  g2BF->Draw("same P");
  TLegend *legend2d = new TLegend(0.65,0.70,0.9,0.9);
  legend2d->AddEntry("g2","  True values");
  legend2d->AddEntry("g2s","  Initial Seeds");
  legend2d->AddEntry("g2BF","  Bestfit values");
  legend2d->Draw();
  c2d->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d2.png",name,para,cutNu,cutEv,s));
  TCanvas *c0dz = new TCanvas("c0dz","",800,600);
  c0dz->SetLogz(0);
  hc0dz->Draw("colz");
  h0z->Draw("same");
  g0->Draw("same P");
  g0s->Draw("same P");
  g0BF->Draw("same P");
  TLegend *legend0dz = new TLegend(0.65,0.70,0.9,0.9);
  legend0dz->AddEntry("g0","  True values");
  legend0dz->AddEntry("g0s","  Initial Seeds");
  legend0dz->AddEntry("g0BF","  Bestfit values");
  legend0dz->Draw();
  c0dz->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d0_zoom.png",name,para,cutNu,cutEv,s));
  TCanvas *c1dz = new TCanvas("c1dz","",800,600);
  c1dz->SetLogz(0);
  hc1dz->Draw("colz");
  h1z->Draw("same");
  g1->Draw("same P");
  g1s->Draw("same P");
  g1BF->Draw("same P");
  TLegend *legend1dz = new TLegend(0.65,0.70,0.9,0.9);
  legend1dz->AddEntry("g1","  True values");
  legend1dz->AddEntry("g1s","  Initial Seeds");
  legend1dz->AddEntry("g1BF","  Bestfit values");
  legend1dz->Draw();
  c1dz->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d1_zoom.png",name,para,cutNu,cutEv,s));
  TCanvas *c2dz = new TCanvas("c2dz","",800,600);
  c2dz->SetLogz(0);
  hc2dz->Draw("colz");
  h2z->Draw("same");
  g2->Draw("same P");
  g2s->Draw("same P");
  g2BF->Draw("same P");
  TLegend *legend2dz = new TLegend(0.65,0.70,0.9,0.9);
  legend2dz->AddEntry("g2","  True values");
  legend2dz->AddEntry("g2s","  Initial Seeds");
  legend2dz->AddEntry("g2BF","  Bestfit values");
  legend2dz->Draw();
  c2dz->SaveAs(Form("%s_fit_parDraw_diff_nCov_%d%d%d_%d2_zoom.png",name,para,cutNu,cutEv,s));

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TFile *out = new TFile(Form("%s_fitResults_nCov_%d%d%d_%d.root",name,para,cutNu,cutEv,s),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  g0BF->Write();
  g1BF->Write();
  g2BF->Write();
  out->Close();
}

void TemplateFitter::TrueDraw()
{
  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;
  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;
  double par[3];
  if(para == 1) {
    par[0] = seed1[0][0];
    par[1] = seed1[0][1];
    par[2] = seed1[0][2];
  }
  else if(para == 2) {
    par[0] = seed2[0][0];
    par[1] = seed2[0][1];
    par[2] = seed2[0][2];
  }

  std::cout << "true values: " << par[0] << "\t" << par[1] << "\t" << par[2] << "\n";

  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    }
    double Pmue = mue/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CC_e_target->SetName("e_target");
  CC_e_target->SetLineColor(kBlack);
  CC_m_target->SetName("m_target");
  CC_m_target->SetLineColor(kBlack);
  target_nue->SetName("nue_target");
  target_nue->SetLineColor(kBlack);

  CC_tp_e->SetName("tp_e");
  CC_tp_e->SetMarkerStyle(21);
  CC_tp_e->SetMarkerColor(kRed);
  CC_tp_e->SetMarkerSize(0.5);
  CC_tp_me->SetName("tp_me");
  CC_tp_me->SetMarkerStyle(21);
  CC_tp_me->SetMarkerColor(kBlue);
  CC_tp_me->SetMarkerSize(0.5);

  CC_tp_m->SetName("tp_m");
  CC_tp_m->SetMarkerStyle(21);
  CC_tp_m->SetMarkerColor(kRed);
  CC_tp_m->SetMarkerSize(0.5);
  CC_tp_em->SetName("tp_em");
  CC_tp_em->SetMarkerStyle(21);
  CC_tp_em->SetMarkerColor(kBlue);
  CC_tp_em->SetMarkerSize(0.5);

  nue_tp->SetName("tp_nue");
  nue_tp->SetMarkerStyle(21);
  nue_tp->SetMarkerColor(kRed);
  nue_tp->SetMarkerSize(0.5);
  nue_tp_os->SetName("tp_nue_os");
  nue_tp_os->SetMarkerStyle(21);
  nue_tp_os->SetMarkerColor(kBlue);
  nue_tp_os->SetMarkerSize(0.5);

  TCanvas *c_e = new TCanvas("c_e","",900,700);
  CC_e_target->Draw();
  CC_tp_me->Draw("same");
  CC_tp_e->Draw("same");
  TLegend *legend_e = new TLegend(0.65,0.70,0.9,0.9);
  legend_e->AddEntry(CC_e_target,"fluctuated target");
  legend_e->AddEntry(CC_tp_me,"true oscillated mu->e");
  legend_e->AddEntry(CC_tp_e,"templates at true values");
  legend_e->Draw();
  c_e->SaveAs(Form("%s_TrueDraw_CC_e_%d%d%d_ft.png",name,para,cutNu,cutEv));

  TCanvas *c_m = new TCanvas("c_m","",900,700);
  CC_m_target->Draw();
  CC_tp_em->Draw("same");
  CC_tp_m->Draw("same");
  TLegend *legend_m = new TLegend(0.65,0.70,0.9,0.9);
  legend_m->AddEntry(CC_m_target,"fluctuated target");
  legend_m->AddEntry(CC_tp_em,"true oscillated e->mu");
  legend_m->AddEntry(CC_tp_m,"templates at true values");
  legend_m->Draw();
  c_m->SaveAs(Form("%s_TrueDraw_CC_m_%d%d%d_ft.png",name,para,cutNu,cutEv));

  TCanvas *c = new TCanvas("c","",900,700);
  gPad->SetLogy();
  target_nue->Draw();
  nue_tp_os->Draw("same");
  nue_tp->Draw("same");
  TLegend *legend_nue = new TLegend(0.65,0.78,0.9,0.9);
  legend_nue->AddEntry(target_nue,"fluctuated target");
  legend_nue->AddEntry(nue_tp_os,"oscillated mu->e & e->mu");
  legend_nue->AddEntry(nue_tp,"templates at true values");
  legend_nue->Draw();
  c->SaveAs(Form("%s_TrueDraw_nue_%d%d%d_ft_Log.png",name,para,cutNu,cutEv));
}
*/

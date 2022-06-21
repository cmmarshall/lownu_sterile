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


double TemplateFitter::bfChi2 ( double *par )
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

  std::cout << "bestfit values: " << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";

  return chi2;
}

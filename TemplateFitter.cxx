#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"

TemplateFitter::TemplateFitter(TH1D * templates[10], TH1D * intrinsic, TH1D * target)
{
  for( int i = 0; i < 10; ++i ) m_templates[i] = templates[i];
  m_intrinsic = intrinsic;
  m_target = target;
}

void TemplateFitter::setEnergyBins( double bins[11] )
{
  for( int i = 0; i < 11; ++ i ) m_energy_bins[i] = bins[i];
}

double TemplateFitter::getOscProb( double energy, double theta, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double prob = pow(sin(2*theta),2)  * pow(sin(del),2);
  return prob;
}

// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( const double * par )
{

  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * temp = (TH1D*) m_intrinsic->Clone("tmp");

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 0; i < 10; ++i ) {
    double e = (m_energy_bins[i+1] - m_energy_bins[i])/2.;
    double oscProb = getOscProb(e, par[0], par[1]); // par[0] = theta, par[1] = dm2
    temp->Add( m_templates[i], oscProb );
  }
  // now we have temp = intrinsic + oscillated nu_e CC

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  for( int bx = 1; bx <= m_target->GetNbinsX(); ++bx ) {
    double tgt = m_target->GetBinContent(bx);
    double diff = temp->GetBinContent(bx) - tgt;
    if( tgt > 0. ) chi2 += (diff*diff) / tgt;
  }

  return chi2;

}



bool TemplateFitter::doFit( double &theta, double &dm2 )
{
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(10000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(10000);
  fitter->SetTolerance(0.001); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  fitter->SetVariable( 0, "theta", 0.001, 0.00001 );
  fitter->SetVariable( 1, "dm2", 5.0, 0.01 );

  // 2 free parameters = theta, dm2
  ROOT::Math::Functor lf( this, &TemplateFitter::getChi2, 2 );
  ROOT::Math::Functor functor( lf, 2 );
  fitter->SetFunction( functor );

  // Go!
  fitter->Minimize();

  if( fitter->Status() != 0 ) {
    std::cout << "Something bad happened" << std::endl;
    return false;
  }

  const double *bestfit = fitter->X();
  theta = bestfit[0];
  dm2 = bestfit[1];
  return true;

}





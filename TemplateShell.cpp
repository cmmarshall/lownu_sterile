#include <string>

// The templates are some histogram for each bin of neutrino energy
TH2D * templates[10];

// this is the thing you are trying to fit to, i.e. the data distribution
TH2D * target;

bool doFit( double bestfit_norms[10] )
{
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(10000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(10000);
  fitter->SetTolerance(0.001); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  for( int i = 0; i < 10; ++i ) {
    std::string varname = "temp" + std::to_string(i);
    fitter->SetVariable( i, varname, 1.0, 0.001 );
  }

  // 10 free paraemters = 10 normalizations
  ROOT::Math::Functor lf( this, getChi2, 10 );
  ROOT::Math::Functor functor( lf, 10 );
  fitter->SetFunction( functor );

  // Go!
  fitter->Minimize();

  if( fitter->Status() != 0 ) {
    std::cout << "Something bad happened" << std::endl;
    return false;
  }

  const double *bestfit = fitter->X();
  for( int i = 0; i < 10; ++i ) {
    bestfit_norms[i] = bestfit[0];
  }
  return true;

}

// function whose return Minuit mimizes, must take const double* and return double
double getChi2( const double * par ) const
{

  // Create a histogram "temp" from the templates, scaled based on the normalizations in par
  TH2D * temp = templates[0]->Clone("tmp");
  temp->Scale( par[0] );
  for( int i = 1; i < 10; +=i ) temp->Add( templates[i], par[i] );

  double chi2 = 0.0;
  for( int bx = 1; bx <= target->GetNbinsX(); ++bx ) {
    for( int by = 1; by <= target->GetNbinsY(); ++by ) {
      double tgt = target->GetBinContent(bx,by);
      double diff = temp->GetBinContent(bx,by) - tgt;
      if( tgt > 0. ) chi2 += (diff*diff) / tgt;
    }
  }

  return chi2;

}



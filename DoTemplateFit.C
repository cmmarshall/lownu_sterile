#include <string>

// The templates are reconstructed lepton energy, in a slice of true neutrino energy
TH1D * templates[10];
// The template for the intrinsic nue
TH1D * intrinsic; 
// define the energy bins used in the template

// this is the thing you are trying to fit to, i.e. the data distribution
TH1D * target;




double getOscProb( double energy, double theta, double dm2 )
{
    double L = 0.5;
    double del = 1.27*L*dm2/energy;
    double prob = pow(sin(2*theta),2)  * pow(sin(del),2);
    return prob;
}

// function whose return Minuit mimizes, must take const double* and return double
double getChi2( const double * par )
{

  double energy_max = 10;
  double energy_min = 0;
  int N = 10;
  double size = (energy_max-energy_min)/N;
  double energy_range[N];
  for(int i = 0; i < N; i++){
    energy_range[i] = i*size;
  }

  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * temp = (TH1D*)intrinsic->Clone("tmp");

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 0; i < N; ++i ) {
    double e = (energy_range[i+1] - energy_range[i])/2.;
    double oscProb = getOscProb(e, par[0], par[1]); // par[0] = theta, par[1] = dm2
    temp->Add( templates[i], oscProb );
  }
  // now we have temp = intrinsic + oscillated nu_e CC

  // calculate the chi2 with the "data" target
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



bool doFit( double &theta, double &dm2 )
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
  ROOT::Math::Functor functor( &getChi2, 2 );
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


// main function
void DoTemplateFit()
{
  TFile *f = new TFile("output.root","READ");
  TH2D* h = (TH2D*)f->Get("hElepVsEv");
  for(int i=0; i<10; i++) { 
    templates[i] = (TH1D*)h->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  }

  intrinsic = (TH1D*)f->Get("hLepE_sm");
  target = (TH1D*)f->Get("hElep_w");

  /*
  double par[2];
  par[0] = 0.001;
  par[1] = 0.5;

  cout << getChi2(par) << endl;

  templates[2]->Draw();
  */

  double bf_dm2, bf_theta;
  bool isOK = doFit( bf_theta, bf_dm2 );

  printf( "Best-fit theta = %f, dm2 = %f\n", bf_theta, bf_dm2 );

}







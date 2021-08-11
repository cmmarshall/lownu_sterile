#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 1

#include "Math/Minimizer.h"
#include "TH1.h"

class TemplateFitter {

  public:

    TemplateFitter(TH1D * templates[10], TH1D * intrinsic, TH1D * target);
    ~TemplateFitter(){};
    void setEnergyBins(double bins[11]);
    bool doFit(double &theta, double &dm2);

  private:

    double getOscProb(double energy, double theta, double dm2);
    double getChi2(const double * par);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * m_templates[10];
    // The template for the intrinsic nue
    TH1D * m_intrinsic; 
    // define the energy bins used in the template
    double m_energy_bins[11];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * m_target;
};

#endif

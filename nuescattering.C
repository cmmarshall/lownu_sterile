    #include <string>

void nuescattering()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "tree", "tree" );
  for(int i = 10; i<12; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_0%d.root",i) );
  std::cout << "Nue File number:" << i << "\n";
  }
  
  //tree->Add( "nueFHC_000.root" );

  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,40);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,16);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,40);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,16);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,40);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,16);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,40);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,40);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,40);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,16);
  TH1D *hElep_w0   = new TH1D("hElep_w0","",100,0,40);
  TH1D *hElep_w3   = new TH1D("hElep_w3","",100,0,16);
  TH1D *hElep_w0_test = new TH1D("hElep_w0_test","",100,0,40);
  TH1D *os_hElep_w0   = new TH1D("os_hElep_w0","",100,0,40);
  TH1D *os_hElep_w3   = new TH1D("os_hElep_w3","",100,0,16);
  TH1D *unos_hElep_w0   = new TH1D("unos_hElep_w0","",100,0,40);
  TH1D *unos_hElep_w3   = new TH1D("unos_hElep_w3","",100,0,16);

  TH2D *BthetaVsEe = new TH2D("BthetaVsEe","",400,0.5,60,20,0,150);  

  const Int_t nbinsX0 = 400; const Int_t nbinsY0 = 100;
  Double_t xEdges0[nbinsX0+1], yEdges0[nbinsY0+1];
  xEdges0[0]=yEdges0[0]=0;
  for(int i=0; i<nbinsX0+1; i++)
  {
    if(i<200)                xEdges0[i+1] = xEdges0[i] + 0.02;
    else if(i>=200 && i<240) xEdges0[i+1] = xEdges0[i] + 0.1;
    //else if(i>=240 && i<280) xEdges0[i+1] = xEdges0[i] + 0.2;
    else                     xEdges0[i+1] = xEdges0[i] + 0.2;
  }
  for(int i=0; i<nbinsY0+1; i++)
  {
    yEdges0[i+1]=yEdges0[i]+0.4;
  }

  const Int_t nbinsX = 280; const Int_t nbinsY = 100;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else                     xEdges[i+1] = xEdges[i] + 0.2;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.16;
  }
   
  TH2D *m_hElepVsEv0 = new TH2D("m_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *m_hElepVsEv3 = new TH2D("m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0 = new TH2D("e_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *e_hElepVsEv3 = new TH2D("e_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);

  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  double Enu, E[2], p[2]; // the energy of the neutrino, in GeV
  double px[2], py[2], pz[2];
  double best_px[2], best_py[2], best_pz[2];
  int pdg[2];
  double BthetaX, BthetaY, Btheta;
  double PthetaX, PthetaY, Ptheta;
  double BthetaX_sm, BthetaY_sm, Btheta_sm;

  tree->SetBranchAddress( "Enu", &Enu );
  tree->SetBranchAddress( "E", &E );
  tree->SetBranchAddress( "pdg", &pdg );
  tree->SetBranchAddress( "px", &px );
  tree->SetBranchAddress( "py", &py );
  tree->SetBranchAddress( "pz", &pz );
  tree->SetBranchAddress( "best_px", &best_px );
  tree->SetBranchAddress( "best_py", &best_py );
  tree->SetBranchAddress( "best_pz", &best_pz );

  TF1 *tsmear1 = new TF1( "tsmear1", "3.29 + 3.485*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmear2 = new TF1( "tsmear2", "10.287 + 4.889*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmearRatio = new TF1( "tsmearRatio", "0.039 + 0.551*pow(x,-1.) - 0.268*pow(x,-0.5)", 0., 999.9 );
  TF1 *doubleGaus = new TF1( "dg", "[0]*TMath::Exp(-0.5*pow(x/[1],2)) + [2]*TMath::Exp(-0.5*pow(x/[3],2))", -1000., 1000. );

  double s2tW = 0.23;
  double m_C_LL = -1./2. + s2tW; 
  double m_C_LR = s2tW;
  double e_C_LL = 1./2. + s2tW;
  double e_C_LR = s2tW;
  double y, sigma_m, sigma_e;

  double Uee2 = 0.01;
  double Umm2 = 0.0016;
  double s2ee2 = 4*Uee2*(1 - Uee2); 
  double s2mm2 = 4*Umm2*(1 - Umm2); 
  double s2me2 = 4*Uee2*Umm2;
  double dm2 = 1.3; //eV2
  double L = 500; //m
  double del, pmumu, pmue, pee;

  const int N = tree->GetEntries();
  //const int N = 1000;
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);


    //for(int i=0; i<2; i++){
    if( pdg[1]==14 || pdg[1]==12 ){
    //if(E[i]<.5) continue;

    BthetaX = atan2(best_px[0],best_pz[0])*1E3;
    BthetaY = atan2(best_py[0],best_pz[0])*1E3;
    Btheta = atan2(sqrt(best_py[0]*best_py[0]+best_px[0]*best_px[0]),best_pz[0])*1E3;
 /*
    PthetaX = atan2(perf_px[i],perf_pz[i])*1E3;
    PthetaY = atan2(perf_py[i],perf_pz[i])*1E3;
    Ptheta = atan2(sqrt(perf_py[i]*perf_py[i]+perf_px[i]*perf_px[i]),perf_pz[i])*1E3;
*/

    double p1 = tsmear1->Eval(E[0]);
    double p3 = tsmear2->Eval(E[0]);
    double ratio = tsmearRatio->Eval(E[0]);

    doubleGaus->SetParameter( 1, p1 );
    doubleGaus->SetParameter( 3, p3 );
    doubleGaus->SetParameter( 0, 1. );
    doubleGaus->SetParameter( 2, ratio );

    BthetaX_sm = BthetaX + doubleGaus->GetRandom();
    BthetaY_sm = BthetaY + doubleGaus->GetRandom();
    Btheta_sm = atan(sqrt(tan(BthetaX_sm/1E3)*tan(BthetaX_sm/1E3) + tan(BthetaY_sm/1E3)*tan(BthetaY_sm/1E3)))*1E3;    
      
    //std::cout << E[0]*Btheta_sm*Btheta_sm << "\n";
    if(E[0]*Btheta_sm*Btheta_sm/1E3<3){
      BthetaVsEe->Fill(E[0],Btheta_sm);
    
      del = 1.27/1E3*dm2*L/Enu;
      pmue = s2me2 * pow(sin(del),2);
      pmumu = 1 - s2mm2 * pow(sin(del),2);   
      pee = 1 - s2ee2 * pow(sin(del),2);
     
      for(int i=0; i<2; i++){
        p[i]=sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
      }
      y = E[0]/Enu;
      sigma_m = (m_C_LL*m_C_LL + m_C_LR*m_C_LR*(1.-y)*(1.-y));
      sigma_e = (e_C_LL*e_C_LL + e_C_LR*e_C_LR*(1.-y)*(1.-y));

      if(pdg[1] == 14){
        //muon to muon ("survival")
        mm_hElep_w0->Fill(E[0],pmumu);
        mm_hElep_w3->Fill(E[0],pmumu);
      
        //muon to electron
        me_hElep_w0->Fill(E[0],pmue*sigma_e/sigma_m);
        me_hElep_w3->Fill(E[0],pmue*sigma_e/sigma_m);

        //muon target weighted     
        m_hElepVsEv0->Fill(Enu,E[0],sigma_e/sigma_m);
        m_hElepVsEv3->Fill(Enu,E[0],sigma_e/sigma_m);

        //muon target     
        m_hElepVsEv0->Fill(Enu,E[0]);
        m_hElepVsEv3->Fill(Enu,E[0]);
      }

      if(pdg[1] == 12){
        //electron to electron ("survival")
        ee_hElep_w0->Fill(E[0],pee);
        ee_hElep_w3->Fill(E[0],pee);
      
        //electron to muon 
        em_hElep_w0->Fill(E[0],pmue*sigma_m/sigma_e);
        em_hElep_w3->Fill(E[0],pmue*sigma_m/sigma_e);
      
        //electron target weighted
        e_hElepVsEv0->Fill(Enu,E[0],sigma_m/sigma_e);
        e_hElepVsEv3->Fill(Enu,E[0],sigma_m/sigma_e);

        //electron target
        e_hElepVsEv0->Fill(Enu,E[0]);
        e_hElepVsEv3->Fill(Enu,E[0]);
      }
    }
  }

  }
  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);
  hElep_w0->Add(e_hElep_w0); hElep_w0->Add(m_hElep_w0);  
  hElep_w3->Add(e_hElep_w3); hElep_w3->Add(m_hElep_w3);  
  os_hElep_w0->Add(me_hElep_w0); os_hElep_w0->Add(em_hElep_w0);
  os_hElep_w3->Add(me_hElep_w3); os_hElep_w3->Add(em_hElep_w3);
  unos_hElep_w0->Add(ee_hElep_w0); unos_hElep_w0->Add(mm_hElep_w0);
  unos_hElep_w3->Add(ee_hElep_w3); unos_hElep_w3->Add(mm_hElep_w3);

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  THStack *e_hsElep_w0 = new THStack("e_hsElep_w0","nu+e oscillated + unoscillated nue weighted Elep_reco");
  me_hElep_w0->SetMarkerStyle(21);
  me_hElep_w0->SetMarkerSize(0.5);
  me_hElep_w0->SetMarkerColor(kRed);
  me_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w0->Add(me_hElep_w0);
  ee_hElep_w0->SetMarkerStyle(21);
  ee_hElep_w0->SetMarkerSize(0.5);
  ee_hElep_w0->SetMarkerColor(kBlue);
  ee_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w0->Add(ee_hElep_w0);
  THStack *e_hsElep_w3 = new THStack("e_hsElep_w3","nu+e oscillated + unoscillated nue weighted Elep_reco (zoom)");
  me_hElep_w3->SetMarkerStyle(21);
  me_hElep_w3->SetMarkerSize(0.5);
  me_hElep_w3->SetMarkerColor(kRed);
  me_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w3->Add(me_hElep_w3);
  ee_hElep_w3->SetMarkerStyle(21);
  ee_hElep_w3->SetMarkerSize(0.5);
  ee_hElep_w3->SetMarkerColor(kBlue);
  ee_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w3->Add(ee_hElep_w3);
  THStack *m_hsElep_w0 = new THStack("m_hsElep_w0","nu+e oscillated + unoscillated numu weighted Elep_reco");
  em_hElep_w0->SetMarkerStyle(21);
  em_hElep_w0->SetMarkerSize(0.5);
  em_hElep_w0->SetMarkerColor(kRed);
  em_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w0->Add(em_hElep_w0);
  mm_hElep_w0->SetMarkerStyle(21);
  mm_hElep_w0->SetMarkerSize(0.5);
  mm_hElep_w0->SetMarkerColor(kBlue);
  mm_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w0->Add(mm_hElep_w0);
  THStack *m_hsElep_w3 = new THStack("m_hsElep_w3","nu+e oscillated + unoscillated numu weighted Elep_reco (zoom)");
  em_hElep_w3->SetMarkerStyle(21);
  em_hElep_w3->SetMarkerSize(0.5);
  em_hElep_w3->SetMarkerColor(kRed);
  em_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w3->Add(em_hElep_w3);
  mm_hElep_w3->SetMarkerStyle(21);
  mm_hElep_w3->SetMarkerSize(0.5);
  mm_hElep_w3->SetMarkerColor(kBlue);
  mm_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w3->Add(mm_hElep_w3);
  TCanvas *cElep_w = new TCanvas("cElep_w","",1600,1200);
  cElep_w->Divide(2,2);
  cElep_w->cd(1);
  e_hsElep_w0->Draw();
  cElep_w->cd(2);
  m_hsElep_w0->Draw();
  cElep_w->cd(3);
  e_hsElep_w3->Draw();
  cElep_w->cd(4);
  m_hsElep_w3->Draw();
  cElep_w->SaveAs("true_nue_Elep_reco_w_1.pdf");
  
  THStack *hsElep_w0 = new THStack("hsElep_w0","oscillated + unoscillated weighted Elep_reco");
  os_hElep_w0->SetMarkerStyle(21);
  os_hElep_w0->SetMarkerSize(0.5);
  os_hElep_w0->SetMarkerColor(kRed);
  os_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  hsElep_w0->Add(os_hElep_w0);
  unos_hElep_w0->SetMarkerStyle(21);
  unos_hElep_w0->SetMarkerSize(0.5);
  unos_hElep_w0->SetMarkerColor(kBlue);
  unos_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  hsElep_w0->Add(unos_hElep_w0);
  THStack *hsElep_w3 = new THStack("hsElep_w3","oscillated + unoscillated weighted Elep_reco");
  os_hElep_w3->SetMarkerStyle(21);
  os_hElep_w3->SetMarkerSize(0.5);
  os_hElep_w3->SetMarkerColor(kRed);
  os_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  hsElep_w3->Add(os_hElep_w3);
  unos_hElep_w3->SetMarkerStyle(21);
  unos_hElep_w3->SetMarkerSize(0.5);
  unos_hElep_w3->SetMarkerColor(kBlue);
  unos_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  hsElep_w3->Add(unos_hElep_w3);
  TCanvas *cElep_w_sum = new TCanvas("cElep_w_sum","",1600,600);
  cElep_w_sum->Divide(2,1);
  cElep_w_sum->cd(1);
  hsElep_w0->Draw();
  cElep_w_sum->cd(2);
  hsElep_w3->Draw();
  cElep_w_sum->SaveAs("true_nue_Elep_reco_w_1_sum.pdf");
  
  e_hElep_w0->SetTitle("nu+e nue weighted Elep_reco");
  e_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w0->SetTitle("nu+e numu weighted Elep_reco");
  m_hElep_w0->GetXaxis()->SetTitle("Elep_reco");

  m_hElepVsEv0->SetTitle("nu+e numu Elep_reco Vs true Ev");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  m_hElepVsEv3->SetTitle("nu+e numu Elep_reco Vs true Ev (zoom)");
  m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv0->SetTitle("nu+e nue Elep_reco Vs true Ev");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv3->SetTitle("nu+e nue Elep_reco Vs true Ev (zoom)");
  e_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",1600,1200);
  cElepVsEv->Divide(2,2);
  cElepVsEv->cd(1);
  m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(3);
  m_hElepVsEv3->Draw("colz");
  cElepVsEv->cd(4);
  e_hElepVsEv3->Draw("colz");
  cElepVsEv->SaveAs("true_nue_ElepVsEv_1.pdf");

  TF1 *cut = new TF1("cut","sqrt(3*1E3/x)",0.5,60);
  TCanvas *cBthetaVsEe = new TCanvas("cBthetaVsEe","",800,600);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut->Draw("same");
  cBthetaVsEe->SaveAs("nue_BthetaVsEe.pdf");

  TFile *out = new TFile("nue_output_1.root","RECREATE");
  e_hElep_w0->Write();
  e_hElep_w3->Write();
  me_hElep_w0->Write();
  me_hElep_w3->Write();
  ee_hElep_w0->Write();
  ee_hElep_w3->Write();
  m_hElep_w0->Write();
  m_hElep_w3->Write();
  em_hElep_w0->Write();
  em_hElep_w3->Write();
  mm_hElep_w0->Write();
  mm_hElep_w3->Write();
  //os_hElep_w0->Write();
  //unos_hElep_w0->Write();
  hElep_w0->Write();
  hElep_w3->Write();
  m_hElepVsEv0->Write();
  m_hElepVsEv3->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv3->Write();
  out->Close();

}


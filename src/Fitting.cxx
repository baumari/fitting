#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TApplication.h>
#include <TVirtualPad.h>
#include <TString.h>
#include <TF1.h>

#include <KTF1Spline.hh>
#include <KUtil.hh>
#include <KExpdata.hh>
#include <KTheodata.hh>

#define ROOT
//#define ADDERR 

int main(int argc, char* argv[]){
#ifdef ADDERR
  #ifdef ROOT
  std::cout << "You cannot use root fit in adderr option!!" << std::endl;
  std::exit(EXIT_FAILURE);
  #endif
#endif

  /** Check environment **/
#ifndef ROOT
  if(argc != 3){
    fprintf(stderr, "Usage: %s expfile theofile\n",argv[0]);
    fprintf(stderr, "Multiple theofiles are only allowed in ROOT-mode.\n");
    std::exit(EXIT_FAILURE);
  }  
#endif  
  if(argc < 3){
    fprintf(stderr, "Usage: %s expfile theofile1 theofile2 ...\n",argv[0]);
    std::exit(EXIT_FAILURE);
  }
  /** Read experimental data **/
  KExpdataCS Exp(argv[1]);
  /** Read Theoretical data **/  
  std::vector<KTheodata> vTheo;
  for(int itheo = 2; itheo != argc; ++itheo){
    vTheo.push_back(KTheodata(argv[itheo]));
  }

#ifdef ROOT
  TApplication app("app", &argc, argv);
  auto c = new TCanvas("c","c");
  c->Divide(1,2);  
  auto gr = new TGraphErrors(Exp.GetN(), &Exp.fx[0], &Exp.fy[0], &Exp.fx_err[0], &Exp.fy_err[0]);
  gr->SetTitle("Experimental Data");
  std::vector<KTF1Spline*> vSp;
  for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
    vSp.push_back(new KTF1Spline(Form("sp%d",itheo), vTheo[itheo], false));
    vSp.back()->Getf()->SetName(Form("f%d", itheo));
  }
  std::string lambda_func("[&](double *x, double *p){return p[0]*f0->Eval(x[0])");
  for(int itheo = 1; itheo != vTheo.size(); ++itheo){
    std::string s(Form(" + p[%d]*f%d->Eval(x[0])", itheo, itheo));
    lambda_func += s;
  }
  lambda_func += ";}";
  double xmin, xmax;
  vSp[0]->Getf()->GetRange(xmin, xmax);
  auto fitf = new TF1("fitf", lambda_func.c_str(), xmin, xmax, vTheo.size());
  for(std::size_t ipar = 0; ipar != vTheo.size(); ++ipar) fitf->SetParLimits(ipar, 0., 1.);
  gr->Fit(fitf, "Q0", "", vTheo[0].GetfxMin(), vTheo[0].GetfxMax());
  gr->SetMarkerColor(kBlack);
  gr->SetMarkerStyle(3);
  gr->SetMarkerSize(1);
  gr->Draw("AP");
  fitf->SetLineWidth(1);
  fitf->SetLineColor(kBlack);
  fitf->Draw("lsame");
  std::vector<TF1*> vfResult;
  for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
    vfResult.push_back(new TF1(std::string(Form("fRes%d", itheo)).c_str(),
   			       std::string(Form("[&](double *x, double *p){return p[0]*f%d->Eval(x[0]);}",itheo)).c_str(),
   			       xmin, xmax, 1));
    vfResult.back()->FixParameter(0, fitf->GetParameter(itheo));
  }
  for(std::size_t icolor = 2; auto&& x : vfResult){
    x->SetLineWidth(1);
    x->SetLineColor(icolor++);
    x->Draw("lsame");
  }
  c->Update();
  std::cout << "Scale Factor(Theo_to_Exp): ";
  for(int ipar = 0; ipar != vTheo.size(); ++ipar)
    std::cout << Form("p[%d]: %lf , error: %lf ", ipar, fitf->GetParameter(ipar), fitf->GetParError(ipar));
  std::cout << std::endl;
  std::cout << "Reduced-chisq: " << fitf->GetChisquare()/fitf->GetNDF() << std::endl;

  std::cout << "Additional error-estimation (mizumashi)" << std::endl;
  for(auto &&err : Exp.fy_err) err*=sqrt(fitf->GetChisquare()/fitf->GetNDF());
  auto gr_additional = new TGraphErrors(Exp.GetN(), &Exp.fx[0], &Exp.fy[0], &Exp.fx_err[0], &Exp.fy_err[0]);
  gr_additional->Fit(fitf, "Q0", "", vTheo[0].GetfxMin(), vTheo[0].GetfxMax());
  if(vTheo.size() == 1){
    for(int ipar = 0; ipar != vTheo.size(); ++ipar)
    std::cout << Form("p[%d]: %lf , error: %lf ", ipar, fitf->GetParameter(ipar), fitf->GetParError(ipar));
    std::cout << std::endl;
    std::cout << "Reduced-chisq: " << fitf->GetChisquare()/fitf->GetNDF() << std::endl;
  }else{
    const double delta = 0.01; // varying ratio of coefficient
    const double limit = 0.2; // varying in 10% is acceptable <= you can change
    double MinChisq = fitf->GetChisquare();
    std::vector<double> vMinCoeff;
    for(std::size_t idx = 0; idx != vTheo.size(); ++idx) vMinCoeff.push_back(fitf->GetParameter(idx));      
    std::vector<std::vector<double> > vcoeff, vchisq;
    vcoeff.resize(vTheo.size()); vchisq.resize(vTheo.size());
    for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
      for(double ratio = 1.-limit; ratio <= 1.+limit; ratio+=delta){
	fitf->FixParameter(itheo, vMinCoeff[itheo]*ratio);
	gr_additional->Fit(fitf, "Q0", "", vTheo[0].GetfxMin(), vTheo[0].GetfxMax());
	vcoeff[itheo].push_back(fitf->GetParameter(itheo)); vchisq[itheo].push_back(fitf->GetChisquare());
	fitf->ReleaseParameter(itheo); // release parameter
	fitf->SetParLimits(itheo, 0.,1.);
      }
    }
    auto cc = new TCanvas("cc","Chisquare Contour");
    cc->DivideSquare(vTheo.size());
    std::vector<TGraph*> vGrerr;
    std::vector<TLine*> vLine;
    for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
      vGrerr.push_back(new TGraph(vcoeff[itheo].size(), &vcoeff[itheo][0], &vchisq[itheo][0]));
      vGrerr.back()->SetMarkerColor(kBlack);
      vGrerr.back()->SetMarkerStyle(3);
      vGrerr.back()->SetMarkerSize(1);
      vGrerr.back()->SetTitle(Form("p[%d]", itheo));
      cc->cd(itheo+1);
      vGrerr.back()->Draw("AP");
      cc->Update();      
      vLine.push_back(new TLine(gPad->GetUxmin(), MinChisq+1, gPad->GetUxmax(), MinChisq+1));
      vLine.back()->SetLineColor(kBlack);
      vLine.back()->SetLineWidth(1);
      vLine.back()->Draw("lsame");
      cc->Update();            
    }
  }

  // double MinChisq = f->GetChisquare();
  // double MinCoeff = f->GetParameter(2*Theo.GetN()+1);
  // double chisq = MinChisq;
  // double coeff;
  // const double delta = 0.001;
  // std::vector<double> vCoeff, vChisq;
  // double limit = MinCoeff*0.10; // varying in 20% is acceptable <= you can change
  // for(coeff = MinCoeff-limit; coeff <= MinCoeff+limit;){
  //   f->FixParameter(2*Theo.GetN()+1, coeff);
  //   gr->Fit(f, "Q0", "", Theo.GetfxMin(), Theo.GetfxMax());
  //   vCoeff.push_back(f->GetParameter(2*Theo.GetN()+1));
  //   vChisq.push_back(f->GetChisquare());
  //   coeff += delta*MinCoeff;
  // }
  // c->cd(2);
  // auto grcont = new TGraph((int)vCoeff.size(), &vCoeff[0], &vChisq[0]);
  // grcont->SetTitle("Chisquare contour;Coefficient;Chisquare");  
  // grcont->SetMarkerColor(kBlack);
  // grcont->SetMarkerStyle(3);
  // grcont->SetMarkerSize(1);
  // grcont->Draw("AP");
  // c->Update(); 
  
  // if(vChisq[0] < MinChisq+1 || vChisq.back() < MinChisq+1){
  //   std::cout << "Calcualtion spaces are small!!" << std::endl;
  //   app.Run();    
  //   return 0;
  // }
  // double err_coeff[2];
  // int flag = 0;
  // for(std::size_t idx = 0; idx != vChisq.size(); ++idx){
  //   if(flag == 0 && (vChisq[idx] < MinChisq+1)){
  //     flag = 1;
  //     err_coeff[0] = vCoeff[idx];
  //   }
  //   if(flag == 1 && (vChisq[idx] > MinChisq+1)){
  //     err_coeff[1] = vCoeff[idx];
  //     break;
  //   }
  // }
  // std::cout << "Error(chisq-cont): " << fabs(err_coeff[0]-MinCoeff) << std::endl;
  app.Run();
#endif
  return 0;
}





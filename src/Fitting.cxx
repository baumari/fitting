#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <algorithm>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TApplication.h>
#include <TVirtualPad.h>
#include <TString.h>
#include <TF1.h>
#include <TAxis.h>

#include <KTF1Spline.hh>
#include <KUtil.hh>
#include <KExpdata.hh>
#include <KTheodata.hh>
#include <KOptions.hh>

void usage(){
  std::cout << "fitting expfile theofile1 theofile2 ..." << std::endl;
}


int main(int argc, char* argv[]){
  KOptions opt;
  opt.Add("h","help","Print this help.");
  opt.Add("e","err","Estimation of error bands.");
  opt.Add("o","out","result.out","Output file. "
	  "If this option is selected, graphs are not displayed.");
  if(!opt.Check(argc, argv)){
    usage();
    std::exit(EXIT_FAILURE);
  }
  if(opt.Exist("h")){
    usage();
    opt.Description();
    std::exit(EXIT_SUCCESS);
  }
  if(opt.nArg() < 2){
    usage();
    std::exit(EXIT_FAILURE);
  }
  
  /** Read experimental data **/
  KExpdataCS Exp(argv[opt.LeadArg()]);
  /** Read Theoretical data **/  
  std::vector<KTheodata> vTheo;
  for(int itheo = 1; itheo != opt.nArg(); ++itheo){
    vTheo.push_back(KTheodata(argv[opt.LeadArg() + itheo]));
  }
  /*** Fitting Routine ****/
  TApplication app("app", &argc, argv);
  auto gr = new TGraphErrors(Exp.GetN(), &Exp.fx[0], &Exp.fy[0], &Exp.fx_err[0], &Exp.fy_err[0]);
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
  gr->Fit(fitf, "Q0", "", xmin, xmax);
  std::vector<TF1*> vfResult;
  for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
    vfResult.push_back(new TF1(std::string(Form("fRes%d", itheo)).c_str(),
   			       std::string(Form("[&](double *x, double *p){return p[0]*f%d->Eval(x[0]);}",itheo)).c_str(),
   			       xmin, xmax, 1));
    vfResult.back()->FixParameter(0, fitf->GetParameter(itheo));
  }
  std::cout << "Scale Factor(Theo_to_Exp): ";
  for(int ipar = 0; ipar != vTheo.size(); ++ipar)
    std::cout << Form("p[%d]: %lf , error: %lf ", ipar, fitf->GetParameter(ipar), fitf->GetParError(ipar));
  std::cout << std::endl;
  std::cout << "Reduced-chisq: " << fitf->GetChisquare()/fitf->GetNDF() << std::endl;
  if(!opt.Exist("err")){
    if(opt.Exist("o")){
      std::ofstream ofs(opt.Get<std::string>("o").c_str());
      ofs << fitf->GetChisquare()/fitf->GetNDF() << " ";
      for(int ipar = 0; ipar != vTheo.size(); ++ipar)
	ofs << fitf->GetParameter(ipar) << " " << fitf->GetParError(ipar) << " ";
      ofs << std::endl;
      delete gr; delete fitf;
      for(auto && x : vSp) delete x; vSp.clear();
      for(auto && x : vfResult) delete x; vfResult.clear();
      std::exit(EXIT_SUCCESS);
    }else{
      auto c = new TCanvas("c","c");
      gr->SetTitle("Experimental Data");
      gr->SetMarkerColor(kBlack);
      gr->Draw("AP");
      gr->SetMarkerStyle(3);
      gr->SetMarkerSize(1);
      gr->GetXaxis()->SetLimits(xmin, xmax);
      gr->GetYaxis()->SetRangeUser(fitf->GetMinimum(xmin, xmax)*0.8,
				   fitf->GetMaximum(xmin, xmax)*1.2);
      fitf->SetLineWidth(1);
      fitf->SetLineColor(kBlack);
      fitf->Draw("lsame");
      for(std::size_t icolor = 2; auto&& x : vfResult){
	x->SetLineWidth(1);
	x->SetLineColor(icolor++);
	x->Draw("lsame");
      }
      c->Update();      
      app.Run();
    }
  }

  /*** Estimation of error bands ***/
  std::cout << "Additional error-estimation (mizumashi)" << std::endl;
  double reduced_chi = fitf->GetChisquare()/fitf->GetNDF();
  for(auto &&err : Exp.fy_err) err*=sqrt(fitf->GetChisquare()/fitf->GetNDF());
  auto gr_additional
    = new TGraphErrors(Exp.GetN(), &Exp.fx[0], &Exp.fy[0], &Exp.fx_err[0], &Exp.fy_err[0]);
  gr_additional->Fit(fitf, "Q0", "", xmin, xmax);
  if(vTheo.size() == 1){
    std::cout << Form("p[%d]: %lf , error: %lf ",
		      0, fitf->GetParameter(0), fitf->GetParError(0));
    std::cout << std::endl;
    std::cout << "Reduced-chisq: " << fitf->GetChisquare()/fitf->GetNDF() << std::endl;
    if(opt.Exist("o")){
      std::ofstream ofs(opt.Get<std::string>("o").c_str());
      ofs << fitf->GetChisquare()/fitf->GetNDF() << " "
	  << fitf->GetParameter(0) << " " << fitf->GetParError(0) << std::endl;
    }
    delete gr_additional;
    std::exit(EXIT_SUCCESS);
  }else{ // vTheo.size() >= 2
    double MinChisq = fitf->GetChisquare();
    std::vector<double> vMinCoeff;
    for(std::size_t idx = 0; idx != vTheo.size(); ++idx)
      vMinCoeff.push_back(fitf->GetParameter(idx));
    std::vector<double> vErrorLow(vTheo.size());
    std::vector<double> vErrorHigh(vTheo.size());
    std::vector<std::vector<double>> vChisq(vTheo.size()), vCoeff(vTheo.size());
    // only for graphics
    double count = 0; // counter for binary search
    double bin_min, bin_max; // binary fit range
    double dchisq;
    int zero = 0;
    double delta = 0.2;
    for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){ //binary search
      do{ // low
	if(vMinCoeff[itheo]*(1.-count*delta) < 0){
	  zero = 1;
	  break;
	}
 	// initialization of parameters
	for(int ipar = 0; ipar != vTheo.size(); ++ipar){
	  fitf->ReleaseParameter(ipar);
	  fitf->SetParLimits(ipar, 0, 1);
	}	
	fitf->FixParameter(itheo, vMinCoeff[itheo]*(1.-count*delta));
	gr_additional->Fit(fitf, "Q0", "", xmin, xmax);
	count += 1;
	vCoeff[itheo].push_back(vMinCoeff[itheo]*(1.-count*delta));
	vChisq[itheo].push_back(fitf->GetChisquare());
      }while(fitf->GetChisquare() < MinChisq + 1); // fit range determined here
      if(zero) vErrorLow[itheo] = 0;
      else{
	bin_min = vMinCoeff[itheo]*(1.-count*delta);
	bin_max = vMinCoeff[itheo];
	while(1){
	  // initialization of parameters
	  for(int ipar = 0; ipar != vTheo.size(); ++ipar){
	    fitf->ReleaseParameter(ipar);
	    fitf->SetParLimits(ipar, 0, 1);
	  }	
	  fitf->FixParameter(itheo, (bin_min+bin_max)/2.);
	  gr_additional->Fit(fitf, "Q0", "", xmin, xmax);
	  dchisq = fitf->GetChisquare();
	  vCoeff[itheo].push_back((bin_min+bin_max)/2.);
	  vChisq[itheo].push_back(dchisq);	
	  if(fabs(dchisq-(MinChisq+1)) < 0.0001*MinChisq) break;
	  if(dchisq-(MinChisq+1) > 0) bin_min = (bin_min+bin_max)/2.;
	  else bin_max = (bin_min+bin_max)/2.;
	}
	vErrorLow[itheo] = (bin_min+bin_max)/2.;
      }
      fitf->ReleaseParameter(itheo);	
      count = 0; zero = 0;
      do{ // high
	if(vMinCoeff[itheo]*(1.+count*delta) > 1){
	  zero = 1;
	  break;
	}
	// initialization of parameters
	for(int ipar = 0; ipar != vTheo.size(); ++ipar){
	  fitf->ReleaseParameter(ipar);
	  fitf->SetParLimits(ipar, 0, 1);
	}	
	fitf->FixParameter(itheo, vMinCoeff[itheo]*(1.+count*delta));
	gr_additional->Fit(fitf, "Q0", "", xmin, xmax);
	count += 1;
	vCoeff[itheo].push_back(vMinCoeff[itheo]*(1.+count*delta));
	vChisq[itheo].push_back(fitf->GetChisquare());		
      }while(fitf->GetChisquare() < MinChisq + 1); // fit range determined here
      if(zero) vErrorHigh[itheo] = 1;
      else{
	bin_min = vMinCoeff[itheo];
	bin_max = vMinCoeff[itheo]*(1.+count*delta);
	while(1){
	  // initialization of parameters
	  for(int ipar = 0; ipar != vTheo.size(); ++ipar){
	    fitf->ReleaseParameter(ipar);
	    fitf->SetParLimits(ipar, 0, 1);
	  }	
	  fitf->FixParameter(itheo, (bin_min+bin_max)/2.);
	  gr_additional->Fit(fitf, "Q0", "", xmin, xmax);
	  dchisq = fitf->GetChisquare();
	  vCoeff[itheo].push_back((bin_min+bin_max)/2.);
	  vChisq[itheo].push_back(dchisq);		
	  if(fabs(dchisq-(MinChisq+1)) < 0.0001*MinChisq) break;
	  if(dchisq-(MinChisq+1) > 0) bin_max = (bin_min+bin_max)/2.;
	  else bin_min = (bin_min+bin_max)/2.;
	}
	vErrorHigh[itheo] = (bin_min+bin_max)/2.;
      }
      count = 0; zero = 0;
      fitf->ReleaseParameter(itheo);	
    }
    for(std::size_t idx = 0; idx != vTheo.size(); ++idx)
      std::cout << Form("p[%d]: %lf error_low: %lf error_high: %lf ", idx, vMinCoeff[idx],
			vErrorLow[idx], vErrorHigh[idx]);
    std::cout << std::endl;
    if(!opt.Exist("o")){
      auto cc = new TCanvas("cc","Chisquare Contour");
      cc->DivideSquare(vTheo.size());
      std::vector<TGraph*> vGrerr;
      std::vector<TLine*> vLine;
      for(std::size_t itheo = 0; itheo != vTheo.size(); ++itheo){
	vGrerr.push_back(new TGraph(vCoeff[itheo].size(),
				    &vCoeff[itheo][0], &vChisq[itheo][0]));
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
      app.Run();
    }else{
      std::ofstream ofs(opt.Get<std::string>("o").c_str());
      ofs << "# minchisq reduced_chi coeff errlow errhigh" << std::endl;
      ofs << MinChisq << " " << reduced_chi << " ";
      for(std::size_t idx = 0; idx != vTheo.size(); ++idx)
	ofs << vMinCoeff[idx] << " " <<  vErrorLow[idx] << " " <<  vErrorHigh[idx] << " ";
      ofs << std::endl;
      delete gr_additional;
    }
  }
  return 0;
}





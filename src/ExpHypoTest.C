/*************************************************************
 * @author   Triston Ruiseco
 * @file     ExpHypoTest.cpp
 * @date     03/22/2021
 * @brief    Analyzes two exponentially distributed data samples and simulates
 *           highly configurable tests for various measurements
 *           per experiment then plots the signifiance of the tests vs the
 *           measurements per experiment used to generate them.
 *************************************************************/

// Standard library includes
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>


// ROOT includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TFile.h"
#include "TFormula.h"
#include "TF1.h"

// Local includes
#include "LT.h"

// Directives and declarations for namespaces and namespace members
using std::cout, std::cin, std::string, std::vector, std::stoi, std::ifstream,
			std::sort, std::erf, std::sqrt, std::stringstream, std::exp, std::abs,
			std::min, std::max, std::to_string, std::stod, std::pow;


// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
bool strsame(string a, string b);

/**
 * Exponential distribution probability density function
 * param x: data
 * param rate: rate parameter of exp distribution
 * return probability of observing x in exp distribution
 */
double ExpPDF(double x, double rate);

/**
 * Equal type I and type II error probability finding function
 * param arr0: vector of some distribution for some hypothesis 0
 * param arr1: vector of some distribution for some hypothesis 1
 * return significance level of test such that alpha = beta
 */
double FindABSame(const vector<double> &arr0, const vector<double> &arr1);

/**
 * return GaussianPDF(x,0,1)
 */
double sigma(int x);

/**
 * return LomaxPDF(x,a,b)
 */
double pLom(double x, double a, double b);

/**
 * return ExponentialPDF(x,l)
 */
double pExp(double x, double l);

/**
 * return first index of arr such that arr[index] is of lesser value than y
 */
int FirstIndexLess(const vector<double>& arr, double y);

/**
 * Confidence vs measurements plot generation and exportation function
 */
void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title);

int CanvasNumber = 0;
TCanvas* PlotHypotheses(vector<double>& array0, vector<double>& array1,
			string& title, string& xlabel, int order, double alpha = 0.);


// Begin primary program function
int main(int argc, char** argv){

  // Command line option parsing variables
  bool argexists = 0;
  bool printhelp = 0;

  // Command line option storage variables
  bool haveInput[2] = {0,0};
  string InputFile[2] = {"\0","\0"};
  int Nexp = 0;
  int mpe = 0;
	double rate = 0.05131;
<<<<<<< HEAD
	double alpha = 50.0;
	double beta = 996.0;
=======
	double alpha = 10.0;
	double beta = 199.2;
>>>>>>> parent of 91c6502... first commit

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    argexists = 0;
    if(strsame(argv[i],"--help") || strsame(argv[i],"-h")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"--mpe")){
      argexists = 1;
      mpe = stoi(argv[++i]);
      if(mpe <= 0){
        printhelp = 1;
        cout << "--mpe must be a positive integer\n";
      }
    }
    if(strsame(argv[i],"--Nexp")){
      argexists = 1;
      Nexp = stoi(argv[++i]);
      if(Nexp <= 0){
        printhelp = 1;
        cout << "--Nexp must be a positive integer\n";
      }
    }
    if(strsame(argv[i],"--H0")){
      argexists = 1;
      InputFile[0] = string(argv[++i]);
      haveInput[0] = 1;
    }
    if(strsame(argv[i],"--H1")){
      argexists = 1;
      InputFile[1] = string(argv[++i]);
      haveInput[1] = 1;
    }
		if(strsame(argv[i],"--rate")){
      argexists = 1;
      double r = stod(argv[++i]);
      if(r > 0.0){
        rate = r;
      }
    }
		if(strsame(argv[i],"--gamma")){
      argexists = 1;
      double a = stod(argv[++i]);
      if(a > 0.0){
         alpha = a;
      }
			double b = stod(argv[++i]);
			if(b > 0.0){
				 beta = b;
			}
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }

  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide required input, or provides an undefined option */
  if(printhelp || !haveInput[0] || !haveInput[1]){
    cout << "Usage: " << argv[0] << " [options] --Nexp [integer] --mpe [integer] --H0 [filename] --H1 [filename]\n"
         << "  descriptions:\n"
         << "   --Nexp [integer]      number of experiments per test\n"
         << "   --mpe [integer]       max measurements/experiment\n"
         << "   --H0 [filename]       input data for hypothesis 0\n"
         << "   --H1 [filename]       input data for hypothesis 1\n"
         << "  options:\n"
				 << "		--rate [num]				  rate paramater for H0 exp dist (0.05131)\n"
				 << "		--gamma [num] [num]		shape, rate for H1 gamma dist (50,996)\n"
         << "   --help(-h)            print options\n";

    return 0;
  }

  // Experimental data storage objects
  vector<double> times[2];
  string testr = "\0";
  double temptime = 0.0;

  // Read in experimental data
  for(int h = 0; h < 2; ++h){
    ifstream inFile(InputFile[h]);

    // Check if file can be opened
    if(!inFile.is_open()){
      cout << "Failed to open " << InputFile[h] << "\n";
      return 0;
    }

    // Store all experimental data
    LT inputLT(Nexp*mpe, "Reading data set: ");
    for(int i = 0; i < Nexp*mpe && inFile >> temptime; ++i){
      times[h].push_back(temptime);
      inputLT.track();
    }
    inFile.close();

    // Check that enough data was found to perform requested calculations
    if(int(times[h].size()) < (Nexp*mpe)){
      cout << InputFile[h] << " contains too few measurements to complete indicated analysis\n";
      return 0;
    }
  }

  // First, we use our simulated data to build the *probability distribution* for the
  // for a single measurement of the time between decays (integrating over all nuisances implicitly)
  // we represent these probability distributions with histograms
  TH1D* pdf0 = new TH1D("pdf0", "pdf0",
         1000, 0.0, 1000);
  TH1D* pdf1 = new TH1D("pdf1", "pdf1",
         1000, 0.0, 1000);

	double maxtime = 0.0;
  int num_vals = times[1].size();
  for(int i = 0; i < num_vals; ++i){
      pdf0->Fill(times[0][i]);
      pdf1->Fill(times[1][i]);
			if(times[1][i] > maxtime){
				maxtime = times[1][i];
			}
  }
	cout << "\nLargest time value sampled from H1: " << maxtime << " with probability " << pLom(maxtime, alpha, beta) << "\n";

  // Normalize histograms to make PDFs
  pdf0->Scale(1./pdf0->Integral());
  pdf1->Scale(1./pdf1->Integral());

  // Data analysis storage objects
  vector<double> LLR[2];    // log-likelihood of experiment for each h
  vector<double> CR;        // significance level of each generated test
  vector<double> MPE;       // keeps index-wise track of  measures/exp for CR
  double tLLR = 0.0;
	double ab = 0.0;
	string title = "\0";
	string xlabel = "\0";
	string file = "\0";

  /* Construct vector of log likelihood ratios for each experiment,
     for each hypothesis, for each value of measurements/experiment, */
  for(int M = 1; M <= mpe; M = M*1.1 + 1){    // loop exponentially over many measurements/experiment
    for(int h = 0; h < 2; ++h){        				// loop over all hypotheses
      LLR[h].clear();
      for(int e = 0; e < Nexp; ++e){   				// loop over all experiments
        tLLR = 0;
        for(int m = 0; m < M; ++m){    				// loop over all measurements
          tLLR += log(pLom(times[h][(M*e)+m], alpha, beta)); // M*e used to prevent reusing data
          tLLR -= log(pExp(times[h][(M*e)+m], rate));        // across experiments within the same test
        }
        LLR[h].push_back(tLLR);
      }
    }
    // Sort test distributions
    sort(LLR[0].begin(),LLR[0].end());
    sort(LLR[1].begin(),LLR[1].end());

		// Find and record alpha
		ab = FindABSame(LLR[0], LLR[1]);
		CR.push_back(ab);
		MPE.push_back(double(M)/1000.0);

    // Draw distribution of LLRs
    title = to_string(M)+" measurements / experiment";
    xlabel = "#lambda = log [ #it{L}(H1) / #it{L}(H0) ]";
    TCanvas* can1 = PlotHypotheses(LLR[0], LLR[1], title, xlabel, 1, ab);
		file = "Hypotheses" + to_string(M) + ".png";
    can1->SaveAs(file.c_str());
  }

	// Sort time intervals
	sort(times[0].begin(), times[0].end());
	sort(times[1].begin(), times[1].end());

	// Draw distribution of intervals
	title = "Distribution of Time Between Decays";
	xlabel = "Inter-decay Interval (s)";
	TCanvas* can0 = PlotHypotheses(times[0], times[1], title, xlabel, 0);
	can0->SaveAs("Times.png");

  // Plot and save results
  stringstream title2;
  title2 << Nexp << " experiments with #lambda_{0} = " << rate << "s^{-1}, #lambda_{1} ~ Gamma(" << 50 << "," << 996 << ") s^{-1}";
  PlotConfidence(MPE, CR, title2.str());

  return 0;
}

// Program-specific helper function definitions
double sigma(int x){
  return erf(double(x)/sqrt(2));
}

double pLom(double x, double a, double b){
	double power = -(a+1.0);

	double base = 1.0;
	base += x/b;

	double coeff = a/b;

	return (coeff*pow(base, power));
}

double pExp(double x, double l){
	return (l*exp(-x*l));
}

int FirstIndexLess(const vector<double>& arr, double y){
  int n = arr.size();
  for(int i = 0; i < n; ++i){
    if(arr[i] < y){
      return i;
    }
  }
  return n;
}

bool strsame(string a, string b){
  if(a.length()==b.length()){
    int n = a.length();
    for(int i = 0; i < n; i++){
      if(a.at(i)!=b.at(i)){
        return 0;
      }
    }
    return 1;
  }
  return 0;
}

double ExpPDF(double x, double rate){
  return (rate*exp(-rate*x));
}

double FindABSame(const vector<double>& arr0, const vector<double>& arr1){
  int n0 = arr0.size();
  int n1 = arr1.size();
  vector<int> x[2];

  int b = 0;
  for(int a = n0 - 1; a > -1; --a){
    b = 0;
    for(int i = 0; i < n1; ++i){
      if(arr1[i+1] < arr0[a]){
        ++b;
      }
    }
    x[0].push_back(a);
    x[1].push_back(abs(b + a - n0));
  }
  int k = x[1].size();
  int mindex = 0;
  for(int i = 0; i < k; ++i){
    if(x[1][i] < x[1][mindex]){
      mindex = i;
    }
  }
  return (1.0 - (double(x[0][mindex])/double(n0)));
}

void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title){
  TCanvas* canvas = (TCanvas*) new TCanvas("canvas", "Canvas_Title",200,10,500,400);
  canvas->SetGrid();

	// Configure Canvas
  double lm = 0.15;
  double rm = 0.04;
  double bm = 0.15;
  double tm = 0.07;
  canvas->SetLeftMargin(lm);
  canvas->SetRightMargin(rm);
  canvas->SetBottomMargin(bm);
  canvas->SetTopMargin(tm);
  canvas->SetLogy();
  canvas->Draw();
  canvas->Update();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// Configure Graph
  int N = arr0.size();
  TGraph* graph = new TGraph(N, &(arr0[0]), &(arr1[0]));
  graph->SetLineColor(kAzure+1);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(0);

	// Configure X
	graph->GetXaxis()->SetTitleFont(42);
	graph->GetXaxis()->SetTitleSize(0.05);
	graph->GetXaxis()->SetTitleOffset(1.1);
	graph->GetXaxis()->SetLabelFont(42);
	graph->GetXaxis()->SetLabelSize(0.04);
	graph->GetXaxis()->SetTickSize(0.);
  graph->GetXaxis()->SetTitle("1000s of Measurements/Experiment");
	graph->GetXaxis()->CenterTitle();

	// Configure Y
	graph->GetYaxis()->SetTitleFont(42);
	graph->GetYaxis()->SetTitleSize(0.05);
	graph->GetYaxis()->SetTitleOffset(1.1);
	graph->GetYaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetTitle("Test Significance #alpha (= #beta)");
	graph->GetYaxis()->CenterTitle();

	// Draw configured graph on canvas
  graph->Draw();
  canvas->Update();

	// Format default latex text drawing settings
  TLatex text;
	text.SetTextFont(42);
  text.SetTextAlign(21);
  text.SetTextSize(0.04);
  text.SetNDC();

	// Draw title
  text.DrawLatex((1.-rm+lm)/2., 1.-tm+0.012, title.c_str());

	// Draw lines at sigma levels of significance and label them
  TLine* line = new TLine();
  line->SetLineWidth(2);

  int n = arr1.size();
  int i = 1;
  int k = FirstIndexLess(arr1, 1 - sigma(i));

	text.SetTextSize(0.035);
	text.SetTextAlign(33);
	text.SetTextAngle(90);

  while(k < n && i < 8){
    double xndc = (1.-rm-lm)*((arr0[k]-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+lm;
    cout << i << "sigma at " << arr0[k] << "\n";
    line->SetLineColor(kRed+3-i);
    line->DrawLineNDC(xndc,bm,xndc,1.-tm);
    text.DrawLatex(xndc+0.005, 1-tm-0.01, Form("x = %.3f, #alpha = %d #sigma", arr0[k] , i ));
    ++i;
    k = FirstIndexLess(arr1, 1 - sigma(i));
  }

	// Export in log scale as well as normal
  canvas->SaveAs("confidenceLogY.png");
  canvas->SetLogy(0);
  canvas->SaveAs("confidence.png");
}

// Plot an arrays of doubles as a histograms (and return canvas)
TCanvas* PlotHypotheses(vector<double>& array0, vector<double>& array1,
			string& title, string& xlabel, int order, double alpha){
  int N0 = array0.size();
  int N1 = array1.size();
  double hmin = std::min(array0[0], array1[0]);
  double hmax = std::max(array0[N0-1], array1[N1-1]);
  // create histograms
  TH1D* hist = new TH1D(Form("hist0_%d", CanvasNumber),
			Form("hist0_%d", CanvasNumber),
			100, hmin, hmax);
  TH1D* hist1 = new TH1D(Form("hist1_%d", CanvasNumber),
			 Form("hist1_%d", CanvasNumber),
			 100, hmin, hmax);

  for(int i = 0; i < N0; i++){
    hist->Fill(array0[i]);
  }
  for(int i = 0; i < N1; i++){
    hist1->Fill(array1[i]);
  }

  // Normalize to unit area to make a probability distributions
  hist->Scale(1./hist->Integral());
  hist1->Scale(1./hist1->Integral());

  // some formating settings
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas* canvas = (TCanvas*) new TCanvas(Form("canvas_%d", CanvasNumber),
					   Form("canvas_%d", CanvasNumber),500,400);
  double hlo = 0.15;
  double hhi = 0.04;
  double hbo = 0.15;
  double hto = 0.07;
  canvas->SetLeftMargin(hlo);
  canvas->SetRightMargin(hhi);
  canvas->SetBottomMargin(hbo);
  canvas->SetTopMargin(hto);
  canvas->SetGridy();
  canvas->SetLogy();
  canvas->Draw();
  canvas->cd();

	TH1D* front = nullptr;
	TH1D* back = nullptr;
	if(order == 0){
		back = hist1;
		front = hist;
	} else if(order == 1){
		back = hist;
		front = hist1;
	}

  back->GetXaxis()->CenterTitle();
  back->GetXaxis()->SetTitleFont(42);
  back->GetXaxis()->SetTitleSize(0.05);
  back->GetXaxis()->SetTitleOffset(1.1);
  back->GetXaxis()->SetLabelFont(42);
  back->GetXaxis()->SetLabelSize(0.04);
  back->GetXaxis()->SetTitle(xlabel.c_str());
  back->GetXaxis()->SetTickSize(0.);
  back->GetYaxis()->CenterTitle();
  back->GetYaxis()->SetTitleFont(42);
  back->GetYaxis()->SetTitleSize(0.05);
  back->GetYaxis()->SetTitleOffset(1.1);
  back->GetYaxis()->SetLabelFont(42);
  back->GetYaxis()->SetLabelSize(0.035);
  back->GetYaxis()->SetTitle("Probability");
  back->GetYaxis()->SetRangeUser(0.5*std::min(1./double(N0),1./double(N1)),
				 2*std::max(hist->GetMaximum(),hist1->GetMaximum()));

	hist->SetLineColor(kBlue+2);
	hist->SetFillColor(kBlue+1);
	hist->SetFillStyle(3004);

	hist1->SetLineColor(kGreen+2);
	hist1->SetFillColor(kGreen+1);
	hist1->SetFillStyle(3004);

	back->Draw("hist");
	front->Draw("hist same");

  TLegend* leg = new TLegend(0.76,0.68,0.93,0.8947);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
	if(order == 0){
		leg->AddEntry(hist, "H0");
	  leg->AddEntry(hist1, "H1");
	} else if(order == 1){
		leg->AddEntry(hist, "P(#lambda | H0)");
	  leg->AddEntry(hist1, "P(#lambda | H1)");
	}
  leg->Draw();

  // Draw some text on canvas
  TLatex l;
  l.SetTextFont(42);
  l.SetTextAlign(21);
  l.SetTextSize(0.04);
  l.SetNDC();

  l.DrawLatex((1.-hhi+hlo)/2., 1.-hto+0.012, title.c_str());

  if(alpha > 1./double(N0)){
    double lambda_crit = array0[std::min(int((1.-alpha)*N0),N0-1)];
    double beta = 0.;

    TH1D* histp = new TH1D(Form("histp0_%d", CanvasNumber),
			   Form("histp0_%d", CanvasNumber),
			   100, hmin, hmax);
    TH1D* histp1 = new TH1D(Form("histp1_%d", CanvasNumber),
			    Form("histp1_%d", CanvasNumber),
			    100, hmin, hmax);

    for(int i = 0; i < N0; i++){
      if(array0[i] > lambda_crit){
	      histp->Fill(array0[i]);
      }
    }
    for(int i = 0; i < N1; i++){
      if(array1[i] < lambda_crit){
      	histp1->Fill(array1[i]);
	      beta++;
      }
    }

    beta /= double(N1);

    histp->Scale(1./double(N0));
    histp1->Scale(1./double(N1));

    histp->SetLineColor(kBlue+2);
    histp->SetFillColor(kBlue+2);
    histp->SetFillStyle(3104);
    histp1->SetLineColor(kGreen+2);
    histp1->SetFillColor(kGreen+2);
    histp1->SetFillStyle(3104);
    histp1->Draw("hist same");
    histp->Draw("hist same");

    l.SetTextAlign(11);
    l.SetTextSize(0.045);
    l.SetTextColor(kBlue+3);
    if(alpha > 0.01){
      l.DrawLatex(0.8,0.63, Form("#alpha = %.2f", alpha));
    } else {
      l.DrawLatex(0.8,0.63, Form("#alpha = %.3f", alpha));
    }
    l.SetTextColor(kGreen+3);
    if(beta > 0.001){
      l.DrawLatex(0.8,0.55, Form("#beta = %.3f", beta));
    } else {
      l.DrawLatex(0.8,0.55, Form("#beta = %.1e", beta));
    }
    TLine* line = new TLine();
    line->SetLineWidth(2);

    double crit = hlo + (1.-hhi-hlo)*(lambda_crit-hmin)/(hmax-hmin);
    line->DrawLineNDC(crit, hbo, crit, 1-hto);
    l.SetTextSize(0.035);
    l.SetTextAlign(33);
    l.SetTextAngle(90);
    l.SetTextColor(kBlack);
    l.DrawLatex(crit+0.005, 1-hto-0.01, Form("#scale[1.4]{#lambda_{#alpha}} = %.3f", lambda_crit));

  }


  CanvasNumber++;

  return canvas;
}

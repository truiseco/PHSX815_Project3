/*************************************************************
* @author   Triston Ruiseco
* @file     ExpHypoSim.cpp
* @date     02/22/2021
* @brief    Generates an exponentially distributed data sample and exports it
            to a file according to users specifications.
*************************************************************/

// Std Includes
#include <iostream>
#include <fstream>

// Local Includes
#include "Random.h"
#include "LT.h"

// Directives and declarations for namespaces and namespace members
using std::string, std::stoi, std::stod, std::stol, std::cout, std::ofstream;

// Global variables and objects
Random rng(314159);
int mode;
int alpha;
double beta;
double rate;

// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
bool strsame(string a, string b);

/**
 * Model-dependent rate sampling function
 * return a an exponential rate parameter according to model in use
          (fixed if model 0, gamma distributed if model 1)
 */
double getRate();

// Begin primary program function
int main(int argc, char** argv){

  // Command line option parsing variables
  bool argexists = 0;
  bool printhelp = 0;

  // Command line option storage variables
  string filename = "data.txt";
  bool output = 0;
  int num_meas = 1;
  mode = 0;
  rate = 0.05131;
<<<<<<< HEAD
  alpha = 50;
  beta = 996;
=======
  alpha = 10;
  beta = 199.2;
>>>>>>> parent of 91c6502... first commit

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    if(strsame(argv[i],"--help")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"-f")){
      argexists = 1;
      mode = 0;
    }
    if(strsame(argv[i],"--rate")){
      double r = stod(argv[++i]);
      if(r > 0.0){
        rate = r;
      }
    }
    if(strsame(argv[i],"-d")){
      argexists = 1;
      mode = 1;
    }
    if(strsame(argv[i],"--alpha")){
      int a = stoi(argv[++i]);
      if(a > 0){
         alpha = a;
      }
    }
    if(strsame(argv[i],"--beta")){
      double b = stod(argv[++i]);
      if(b > 0.0){
         beta = b;
      }
    }
    if(strsame(argv[i],"--measures")){
      argexists = 1;
      int arg_num = stoi(argv[++i]);
      if(arg_num > 0){
	       num_meas = arg_num;
      }
    }
    if(strsame(argv[i],"--output")){
      argexists = 1;
      filename = string(argv[++i]);
      output = 1;
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }

  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide required input, or provides an undefined option */
  if(printhelp || !output){
  cout << "\nUsage: " << argv[0] << " --output [file] [options]\n"
       << "  options and descriptions:\n"
       << "   --help(-h)            print options\n"
       << "   --measures [int]      number of time measurements (1)\n"
       << "   --output [file]       name of output file\n"
       << "== FIXED MODE ========================================================\n"
       << "   -f                    run in fixed mode: simulate exponential data\n"
       << "                            with fixed rate parameter (default mode)\n"
       << "   --rate [number]       rate parameter of exponential dist (0.05131)\n"
       << "== DISTRIBUTED MODE ==================================================\n"
       << "   -d                    run in distributed mode: sample exponential\n"
       << "                            rate froma gamma distribution with shape\n"
       << "                            and rate parameters alpha and beta\n"
       << "   --alpha [int]         shape parameter of gamma distribution (50)\n"
       << "   --beta [number]       rate parameter of gamma distribution (996)\n\n";

  return 0;
  }
  cout << "\n";

  // Generate and write to file an exponentially distributed sample
  LT myLT(double(num_meas), "File output progress: ");  // Loop progress tracker

  ofstream outFile;
  outFile.open(filename);
  for(int m = 0; m < num_meas; ++m){
    outFile << rng.Exponential(getRate()) << " ";
    myLT.track();
  }
  outFile.close();
  cout << "\n\n";

  return 0;
}

// Program-specific helper function definitions
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

double getRate(){
  double lambda = rate; // Fixed rate (model 0)
  if(mode == 1){ // Distributed rate (model 1)
    lambda = rng.Gamma(alpha, beta);
  }
  return lambda;
}

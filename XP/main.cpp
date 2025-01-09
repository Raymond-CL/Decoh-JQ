#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>


const  int pTnum=100;
double pTmin=40.;
double pTmax=1040;
double pT,pTbin,pTbinmin,pTbinmax;
const  int Pnum=10;
const  int num=100;
std::vector<std::vector<double>> P (pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<double>> Su(pTnum, std::vector<double>(Pnum, 0.0));
std::vector<std::vector<std::vector<double>>> PTable (pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
std::vector<std::vector<std::vector<double>>> SuTable(pTnum, std::vector<std::vector<double>>(Pnum, std::vector<double>(num, 0.0)));
double Pval=0.0;
int    ival,jval;
double Rsize=0.4;
double Qmed=39.99;  //medium scale
double LambdaQCD=0.245748;  //keep consistent with the _lambdaQCD value in the Main/src/QCD.cpp
double lambda=std::log(Qmed/LambdaQCD);
unsigned int Nc=3;
unsigned int nf=3;
double b=11./3.*Nc-2./3.*nf;
double ymin,ymax,ybin;
double PSum,nMul;

//------------------------------------------------------------
int Factorial(const int n){
  if(n<=1) return 1;        //0!=1, 1!=1
  return n*Factorial(n-1);  //n!
}

//------------------------------------------------------------
double LinearInterpolation(const double yp){
  int ylow=(int)(std::floor(yp/ybin));
  int yhigh=ylow+1;
  if(ylow<0 || yhigh>=num) return 0.0;
  return PTable[ival][jval][ylow]+(PTable[ival][jval][yhigh]-PTable[ival][jval][ylow])*(yp-ylow*ybin)/ybin;
}

//------------------------------------------------------------
double Function(double yp, void *params){
  // (void)(params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  double Pval=LinearInterpolation(yp);
  int    nFactorial=Factorial(jval+1);
  return (ymax-yp)*gamma0square*Pval*nFactorial;
}

//------------------------------------------------------------
double FunctionBase(double yp, void *params){
  // (void)(params);
  double gamma0square=4.*Nc/(b*(yp+lambda));
  return (ymax-yp)*gamma0square;
}


//------------------------------------------------------------
void CalProbability(){

  std::stringstream ss;
  ss << "./out-" << Qmed << ".dat";
  std::string OutputString=ss.str();
  std::ofstream OutputFile;
  OutputFile.open(OutputString);
  OutputFile << "# pT P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 <Mul>" << std::endl;

  pTbin=(pTmax-pTmin)/pTnum;
  for(int i=0; i<pTnum; i++){
    ival=i;
    pT=pTmin+pTbin/2.+pTbin*i;
    //pTbinmin=pT-pTbin/2.;
    pTbinmax=pT+pTbin/2.;
    pTbinmax=pTbinmax*Rsize;  //attention Rsize
    OutputFile << pT+pTbin/2. << " ";

    PSum=0.0;
    nMul=0.0;
    for(int j=0; j<Pnum; j++){
      jval=j;
      ybin=std::log(pTbinmax/Qmed)/num;

      double res,err;

      //Index
      if(j==0){
        for(int k=0; k<num; k++){
          ymin=0.0;
          ymax=ybin*(k+1);  //update current y value
          //QAG adaptive integration
          gsl_function FunBase;
          FunBase.function = &FunctionBase;
          FunBase.params = nullptr;
          gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
          res=0.0; err=0.0;
          gsl_integration_qag(&FunBase, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
          gsl_integration_workspace_free(workspace);
          PTable[i][j][k]=exp(-res);
        }
      }
      //Index
      for(int k=0; k<num; k++){
        if(j>0){
          for(int jj=0; jj<j; jj++){
            PTable[i][j][k]+=(double)(jj+1)/j*PTable[i][j-1-jj][k]*SuTable[i][jj][k]/Factorial(jj+1);
          }
        }

        ymin=0.0;
        ymax=ybin*(k+1);  //update current y value
        //QAG adaptive integration
        gsl_function Fun;
        Fun.function = &Function;
        Fun.params = nullptr;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&Fun, ymin, ymax, 0, 1e-3, 10000, 6, workspace, &res, &err);
        gsl_integration_workspace_free(workspace);
        SuTable[i][j][k]=res;
      }

      
      //Probility value
      if(j==0){
        ymin=0.0;
        ymax=std::log(pTbinmax/Qmed);
        //QAG adaptive integration
        gsl_function FBase;
        FBase.function = &FunctionBase;
        FBase.params = nullptr;
        gsl_integration_workspace *ws = gsl_integration_workspace_alloc(10000);
        res=0.0; err=0.0;
        gsl_integration_qag(&FBase, ymin, ymax, 0, 1e-3, 10000, 6, ws, &res, &err);
        gsl_integration_workspace_free(ws);
        P[i][j]=exp(-res);
      }
      if(j>0){
        for(int jj=0; jj<j; jj++){
          P[i][j]+=(double)(jj+1)/j*P[i][j-1-jj]*Su[i][jj]/Factorial(jj+1);
        }
      }
      
      ymin=0.0;
      ymax=std::log(pTbinmax/Qmed);
      //QAG adaptive integration
      gsl_function F;
      F.function = &Function;
      F.params = nullptr;
      gsl_integration_workspace *wsp = gsl_integration_workspace_alloc(10000);
      res=0.0; err=0.0;
      gsl_integration_qag(&F, ymin, ymax, 0, 1e-3, 10000, 6, wsp, &res, &err);
      gsl_integration_workspace_free(wsp);
      Su[i][j]=res;

      PSum+=P[i][j];
      nMul+=(j+1)*P[i][j];
      OutputFile << P[i][j] << " ";
      // std::cout << "pT=" << pT+pTbin/2. << ",  n=" << j+1 << ",  Pn=" << P[i][j] << ",  PSum" << PSum << std::endl;
    }
    OutputFile << nMul << std::endl;
    // std::cout << "pT=" << pT+pTbin/2. << ",  PSum=" << PSum << ",  nMul=" << nMul << std::endl;
  }

  OutputFile.close();
}



//####################################################################################################
int main(){

  CalProbability();

  return 0;
}

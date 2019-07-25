#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector Filter_Rule(NumericVector P, double x, double y) {
  int N = P.length();
  IntegerVector Signals(N);
  double hi=P[0], lo=P[0];
  
  int S=0;
  for(int k=1; k<N; k++) {
    switch(S) {
    default: // case 0
      if( (P[k]-lo)/lo > x ) {
        S=1;
        hi=P[k];
      } else if( (P[k]-hi)/hi < -x ) {
        S=-1;
        lo=P[k];
      } else {
        S=0;
        lo=std::min(lo,P[k]);
        hi=std::max(hi,P[k]);
      } break;
    case 1:
      if( (P[k]-hi)/hi < -y ) {
        S=0;  // long position is liquidated
        lo=P[k];
        hi=P[k];
        if( (P[k]-hi)/hi < -x ) {
          S=-1; // long position is reversed to short
          lo=P[k];  // reset the low
        } 
      } else {
        S=1;  // maintain long position
        hi=std::max(hi,P[k]);
      } break;
    case -1:
      if( (P[k]-lo)/lo > y ) {
        S=0; 
        lo=P[k];
        hi=P[k];
        if( (P[k]-lo)/lo > x ) {
          S=1;
          hi=P[k];
        }
      } else {
        S=-1;
        lo=std::min(lo,P[k]);
      } break;
    }
    Signals[k] = S;
  }
  return Signals;
}

// [[Rcpp::export]]
IntegerVector MA_Crossover(NumericVector P, int fast_n, int slow_n, double b) {
  int N = P.length();
  IntegerVector Signals(N);
  
  double P_fast = mean( P[seq(slow_n-fast_n,slow_n-1)] );
  double P_slow = mean( P[seq(0,slow_n-1)] );
  
  int S=0;
  for(int i=slow_n-1; i<N; i++) {
    switch(S){
    default: // case 0
      if(P_fast > (1+b)*P_slow) {
        S=1;
      } else if (P_fast < (1-b)*P_slow) {
        S=-1;
      } else {
        S=0;
      } break;
    case 1:
      if(P_fast < P_slow) {
        S=0;
        if(P_fast < (1-b)*P_slow) {
          S=-1;
        }
      } else { S=1; } break;
    case -1:
      if(P_fast > P_slow) {
        S=0;
        if(P_fast > (1+b)*P_slow) {
          S=1;
        }
      } else { S=-1; } break;
    }
    P_fast = P_fast + (P[i+1] - P[i-fast_n+1])/fast_n;
    P_slow = P_slow + (P[i+1] - P[i-slow_n+1])/slow_n;
    Signals[i] = S;
  }
  return Signals;
}


// [[Rcpp::export]]
IntegerVector Channel_Breakout(NumericVector P, int win, double x, double b, int c) {
  int N = P.length();
  IntegerVector Signals(N);
  int count=1;
  
  double Pmax=max( P[seq(0,win-1)] );
  double Pmin=min( P[seq(0,win-1)] );
  
  int S=0;
  for(int k=win; k<N; k++) {
    if(count==1) {
      if( Pmax/Pmin-1 <= x ) {
        if( P[k]/Pmax-1 > b ) {
          S=1;
          count = c;
        } else if( P[k]/Pmin-1 < -b ) {
          S=-1;
          count = c;
        } else {
          S=0;
        }
      } else {
        S=0;
      }
    } else {
      count--;
    }
    
    if( P[k-win] < Pmax ) {
      Pmax = std::max(Pmax, P[k]);
    } else {
      Pmax = max( P[seq(k-win+1,k)] );
    }
    if( P[k-win] > Pmin ) {
      Pmin = std::min(Pmin, P[k]);
    } else {
      Pmin = min( P[seq(k-win+1,k)] ); 
    }
    Signals[k] = S;
  }
  return Signals;
}


// [[Rcpp::export]]
IntegerVector RSI_Rule(NumericVector P, int h, int v, int d) {
  int N = P.length();
  IntegerVector Signals(N);

  double RSI=0, U=0, D=0;
  for(int i=1; i<=h; i++) {
    U = (P[i]>P[i-1]) ? (U + P[i]-P[i-1]) : U;
    D = (P[i]<P[i-1]) ? (D + P[i-1]-P[i]) : D;
    RSI=100*U/(U+D);
  }
  
  int overbought=0;
  int oversold=0;
  int S=0;
  for(int k=h+1; k<N; k++) {
    if(RSI>50+v) overbought=overbought+1;
    if(RSI<50-v) oversold=oversold+1;
    if( (overbought>=d) & (RSI<50+v) ) {
      S=-1;
      overbought=0;
    }
    if( (oversold>=d) & (RSI>50-v) ) {
      S=1;
      oversold=0;
    }
    Signals[k]=S;
    U = (P[k-h]>P[k-h-1]) ? (U - P[k-h]+P[k-h-1]) : U;
    U = (P[k]>P[k-1]) ? (U + P[k]-P[k-1]) : U;
    D = (P[k-h]<P[k-h-1]) ? (D - P[k-h-1]+P[k-h]) : D;
    D = (P[k]<P[k-1]) ? (D + P[k-1]-P[k]) : D;
    RSI=100*U/(U+D);
  }
  return Signals;
}




#include <iostream>
#include <algorithm>
//#include <Rcpp.h>
#include <vector>
#include <submat.h>

/**
 *Copyright (C) 2011  Ruben Dezeure
 *Contact: dezeurer@student.ethz.ch
 *
 *This program is free software; you can redistribute it and/or
 *modify it under the terms of the GNU General Public License
 *as published by the Free Software Foundation; either version 2
 *of the License, or (at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 **/

using namespace std;
using namespace Rcpp;

RcppExport SEXP main(SEXP a)
{
  // now this main is just used for testing the other functions
  
  std::vector<int> k(4);
  k[0]=1;
  k[1]=2;
  k[2]=3;
  k[3]=4;
  std::vector<int> nextset;
  
  for (int i = 0; i < 10; ++i)
    {
      nextset = getNextSet(100,4,k);
      cout << "{" << nextset[0] <<","<<nextset[1]<<","<<nextset[2]<<","<<nextset[3]<<"}" << endl;
    }
  
  return a;
}

/**
 *This function calculates partial correlation of i and j given the set k
 * C is the correlation matrix among nodes
 */
double pcorOrder(int i,int j,std::vector<int> k,NumericMatrix Corr)
{
  double r;
  double cutat = 0.99999;
  
  if (k.size() == 0)
    {
      r = Corr(i,j);
      
    } 
  /**
   //optimization for the case k.size() ==1 
    else if(k.size() ==1)
    {
      //no match for call to Rcpp::NumericMatrix with (int,vector<int>)
      
      //use vector product calls and power taking and so on ...
      //maybe use armadillo for this?
      
      r = (Corr(i,j)-Corr(i,k)*Corr(j,k))/sqrt((1-pow(Corr(j,k),2))*(1-Corr(i,k)^2));
      }*/ 
  else
    {
      // push_front only works on integervector of rcpp library, 
      // maybe better to use stl vectors
      //k.push_front(j);
      //k.push_front(i);
      //NumericMatrix sub = C(k,k);
      
      int m=Corr.nrow(),n=Corr.ncol();
      arma::mat C(Corr.begin(),m,n,false);//reuses memory and avoids extra copy
      //need an efficient way to get the submatrix off of this. Problem is that i j k not 
      //really represent a range of consecutive rows and columns :/
      //use Rinside? look at presentation a LOT of overhead
      
      //      arma::mat sub = C[c(i,j,k),c(i,j,k)]
      std::vector<int> rows(k.size()+2);
      int l =0; // index in rows vector
      rows[l] = i; 
      l++;
      rows[l] = j; 
      l++;
      
      std::vector<int>::iterator row;
 
      for (row = k.begin(); row !=k.end();++row)
	{
	  rows[l] = *row;
	  l++;
	}
      
      std::vector<int> cols = rows;
      
      //might be the big performance stumble block right here, 
      //possibly needs to be optimized
      arma::mat sub = submat(C,rows.begin(),rows.end(),cols.begin(),cols.end());
      arma::mat PM;
      
      try
	{
	  PM = arma::inv(sub);
	}
      catch(runtime_error re)
	{
	  cout << "Caught error yes mam!" << endl;
	  //the matrix appears to be singular
	  cout << sub << endl;
	}
      
      //the correlation matrix is always a positive semi definite matrix
      //inverse can be done faster if the matrix is a positive definite symmetric matrix
      //we specify this so: inv( sympd(sub) )
		    
      //PM <- pseudoinverse(C(c(i,j,k),c(i,j,k)))
      //return -PM[1,2]/sqrt(PM[1,1]*PM[2,2]);
      r = -PM(1,2)/sqrt(PM(1,1)*PM(2,2));      
      //r <- -PM(1,2)/sqrt(PM(1,1)*PM(2,2))
      //invert matrix better done by rcpparmadillo instead of calling the R function
      //not much improvement expected calling the R function :p. maybe good for comparison
    }
  //if(is.na(r)) r<-0
  if(R_IsNA(r))
    r = 0;
  
  //min(cut.at,max(-cut.at,r))
  return min(cutat,max(-cutat,r));
}

  
/**
 * Generate the next set in a list of all possible sets of size k out of 1:n
 * 
 */
std::vector<int> getNextSet(int n, int k,std::vector<int> previous)
{
 /** initial implementation purely based on the R code, might be a faster way to do this*/
  int sum = 0;
  std::vector<int>::iterator row;
  int iter = n-k+1;
  
  for (row = previous.begin(); row !=previous.end();++iter,++row)
    {
      sum += (iter - *row == 0);
    }

  int chInd = k-sum;
  chInd = chInd -1; //working with c++ indexing, not R here
  cout << "chInd = "<< chInd << endl;
  cout << "k= " << k << " sum = " << sum << endl;
  
  if(chInd == 0)
    {
      //was last set to check
      previous[0]=-1; //marks finished
    }else
    {
      //there is still a set to go
      cout << "there is still a set to go"<<endl;
      previous[chInd] =  previous[chInd] +1;
      //do we need this really? Yes to cover all subsets!
      if (chInd < k)
      {
	for (int i = chInd+1; i < k; ++i)
	  {
	    previous[i]=previous[i-1] + 1;
	  }
      }
    }
  return previous;
}



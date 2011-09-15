library(inline)
bod <-
  '
 // now this main is just used for testing the other functions
  
  std::vector<int> k(4);
  k[0]=1;
  k[1]=2;
  k[2]=99;
  k[3]=100;
  std::vector<int> nextset;
  nextset = k;
  for (int i = 0; i < 10; ++i)
    {
      nextset = getNextSet(100,4,nextset);
      cout << "{" << nextset[0] <<","<<nextset[1]<<","<<nextset[2]<<","<<nextset[3]<<"}" << endl;
    }
  
  return a;
'
inc <- '
#include <iostream>
#include <algorithm>
//#include <Rcpp.h>
#include <vector>

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

'
fun <- cxxfunction(signature(a="numeric"),body = bod,includes=inc,plugin="RcppArmadillo")

fun(0)

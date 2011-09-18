#include <iostream>
#include <algorithm>
//#include <Rcpp.h>
#include <vector>
#include <submat.h>
#include <graphfuncts.h>


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

/**
 * function for computing the skeleton in our DAG
 * arguments
 * n number of samples
 * p number of variables
 * m_max max order to keep looking for
 */
RcppExport SEXP main(SEXP p,SEXP alpha,SEXP m_max, SEXP C)
{
  // now this main will be implemented as the skeleton function! :p
  
  int p = Rcpp::as<int>(p);
  double alpha = Rcpp::as<double>(alpha);
  int m_max = Rcpp::as<int>(m_max);
  NumericMatrix Corr(C);
  
  //using c++ datatypes and trying to writing all functions myself, 
  //i.e. not calling R from c++
  //how to store graph? Use matrix?
  bool G[p*p];// only store upper triangle? -> store as vector?
  
  //all connections exist
  initialiseGraph(G,p);

  int row = 0;// the current row we are studying the connections of
  int connections[p-1]; // stores the connections the current row has (maximum p-1)
  //static allocation ==> efficient
  int sizeothers,x,y;
  double pval;
  
  for (int ord = 0; ord <= m_max; ++ord)
    {
      //look for next row with connections then iterate over the remaining connections
      //alternative: save all connections explicitly in double array and loop over those.
      
      row = getNextRowWithConnections(row,G,p); // row == x
      x = row;
      
      while(row != -1)
	{
	  connections = getRowConnections(row,G,p); // getting the connections belonging to the row
	  for (int i = 0; i < p-1; ++i)
	    {
	      //one of the remaining edge tests
	      y=connections[i];
	      
	      //so now we are gonna check the correlation between row and connections[i]
	      //in respect to every other subset of lengths ord of the remaining connections
	      if(y == -1) break;//reached end of connections
	      //y == connections[i]
	      std::vector others = getOtherConnections(i,connections);
	      sizeothers = others.size();
	      if (sizeothers < ord)
		{
		  continue;//goto next loop iteration
		}
	      
	      //initial subset, TODO is there a more efficient way? Builtin way for this?
	      std::vector subset = getSeqVector(ord);
	      
	     
	      while(subset[0] != -1)
		{
		  std::vector k = others[subset] //does this work? check! otherwise write own function that does this.
		  //pval = pcorOrder(x,y,k,C)
		  pval = pcorOrder(x,y,k,Corr);
		  
		  if (pval >= alpha)
		    {
		      //independent
		      G[x,y]=false;G[y,x]=false;
		      break; //no more checking to be done
		    }
		  
		  subset = getNextSet(sizeothers,ord,subset);
		}
	      
	    }
	  row = getNextRowWithConnections(row+1,G,p);
	}
    }
  
  return wrap(G); // return graph matrix, in later stage return a more complete object
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



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
 */
using namespace std;
using namespace Rcpp;

/**
 * This file contains functions doing operations on the graph
 *
 */

/**
 * initialise our graph with all connections on
 *
 *
 */
void initialiseGraph(bool G[],int p)
{
  //cout << "initialisegraph" << endl;
  
  for (int i = 0; i < p*p; ++i)
    {
      G[i] = true;
    }

  for (int i = 0; i < p; ++i)
    {
      //diagonal connections obviously don\'t exist
      G[i*p+i] = false;
    }
}

/**
 * any method looks if there is still a connection in the graph matrix
 */
bool any(bool G[],int p)
{
  //cout << "any" << endl;
  
  for (int i = 0; i < p*p; ++i)
    {
      if (G[i] == true)
	{
	  return true;
	}   
    }   
}

/**
 * row = the row from which we want the connections in the graph
 * G = the graph
 * p = the number of rows
 * returns a integer vector containing the connections of the next 
 * row that still has connections.
 * -1 signals end of connections
 */
void getRowConnections(int row,bool G[],int p,int* connections)
{
  //cout << "getrowconnections" << endl;
  
  int index = 0;//keeping track of the index in the connections vector
  
  //there is still a row >= startrow with connections
  for (int j = row; j < p; ++j)
    {
      if (G[row*p + j] == true)
	{
	  //row*p + row + j, we are only looking from the diagonal and on, otherwise we have doubles in our connections, Graph is symmetric remember!
	  connections[index] = j;
	  index++;  
	}
    }
  
  if (index != p-1)
    {
      //we didn\'t have p-1 connections
      connections[index] = -1; // signal end
    }

  return;
}

/**
 * startrow = the row from which to start looking for a row with connections in the graph
 * G = the graph
 * p = the number of rows
 * returns an integer value that says which row still has connections
 */
int getNextRowWithConnections(int startrow,bool G[],int p)
{
  //cout << "getNextRowWithConnections" << endl;
  
  int returnrow = -1;
  //cout << "startrow = " << startrow << endl;
  //cout << " p = " << p << endl;
  
  
  for (int row = startrow; row < p; row++)
    {
      
      for (int i = row; i < p; i++)
	{
	  
	  //only search through the upper triangular
	  if(G[row*p+i] == 1)
	    {
	      returnrow = row;
	      //break; // this wouldn't exit the whole loop only the inner one!
	      return returnrow;
	    }
	}
    }
  
  //seems like there are no more rows >= startrows containing connections
  return returnrow;
}

/**
 * get the other connection when excluding y
 */
void getOtherConnections(std::vector<int>* others,int j,int* connections,int p)
{
  //cout << "getOtherConnections" << endl;
  (*others).resize(0);
  for (int i = 0; i < p-1; ++i)
    {
      if (connections[i] == -1) break;

      //we don\'t want y==connections[j] in here
      if (i != j)
	{
	  (*others).push_back(connections[i]);
	}
    }
  
  return;
}

/**
 * Creates a vector of size ord with elements 0:ord-1
 *
 */
void getSeqVector(std::vector<int>* subset,int ord)
{
  //cout << "getSeqVector" << endl;
  
  (*subset).resize(ord);

  for (int i = 0; i < ord; ++i)
    {
      (*subset)[i] = i;
    }

  return;
}

/**
 *Get a subset of a vector with the subset signalled by the indices in the vector subset
 *
 */
std::vector<int> getSubset(std::vector<int> set,std::vector<int> subsetind)
{
  //cout << "getsubset" << endl;//DEBUG

  std::vector<int> subset(subsetind.size(),0);
  
  for (int i = 0; i < subsetind.size(); ++i)
    {
      //cout << "i = " << i << " subsetind[i] = " << subsetind[i] << " set[subsetind[i]] = " << set[subsetind[i]] << endl;
      
      subset[i] = set[subsetind[i]];
    }
  return subset;
}



/**
 *Copyright (C) 2011 Alain Hauser
 * code from Alain Hauser
 *to get submatrices
 *not very fast, try to speed up/ optimise yourself!
 *ALTERNATIVE: use the method eleme of the armadillo matrix class
 *from mail {"elem" returns a vector, so you should then reshape its elements into a matrix 
 *using one of the member functions "reshape" or "set_size". 
 *I have, however, never used that approach, but I think it should work as described.}
 **/
namespace arma
{
  
  template <typename T, typename InputIterator> Mat<T> submat(const Mat<T>& input,
							      InputIterator firstRow, InputIterator lastRow,
							      InputIterator firstCol, InputIterator lastCol)
  {
    //cout << "submat" << endl;//DEBUG
    Mat<T> result(std::distance(firstRow, lastRow), std::distance(firstCol, lastCol));
    InputIterator row, col;
    unsigned int i = 0;
    unsigned int j = 0;

    for (row = firstRow; row != lastRow; ++i, ++row) {
      j = 0;
      for (col = firstCol; col != lastCol; ++j, ++col)
	result(i, j) = input(*row, *col);
    }

    return result;
  }

  /**
   * Help function to extract arbitrary subvectors
   */
  template <typename T, typename InputIterator> Col<T> subvec(const Mat<T>& input,
							      InputIterator firstRow, InputIterator lastRow,
							      const unsigned int colind)
  {
    //cout << "subvec" << endl;//DEBUG

    Col<T> result(std::distance(firstRow, lastRow));
    unsigned int i = 0;

    for (; firstRow != lastRow; ++i, ++firstRow)
      result(i) = input(*firstRow, colind);

    return result;
  }
}


/**
 *converts the boolian matrix to a logical matrix so we can easily return it to R
 */
LogicalMatrix convertToLogical(bool G[],int p)
{
  cout << "converttological" << endl;

  LogicalMatrix log(p,p);
  
  for (int i = 0; i < p; ++i)
    {
      for (int j = 0; j<p; ++j)
	{
	  log(i,j) = G[i*p+j];
	}
    }



  return log;
  
}

/**
 *This function calculates partial correlation of i and j given the set k
 * C is the correlation matrix among nodes
 */
double pcorOrder(int i,int j,std::vector<int> k,NumericMatrix Corr)
{  
  //cout << "pcororder" << endl;//DEBUG

  double r;
  double cutat = 0.99999;
  
  if (k.size() == 0)
    {
      r = Corr(i,j);
    }
   //optimization for the case k.size() ==1 
   else if(k.size() ==1)
   {
     r = (Corr(i,j)-Corr(i,k[0])*Corr(j,k[0]))/sqrt((1-pow(Corr(j,k[0]),2))*(1-pow(Corr(i,k[0]),2)));
   } 
  else
    { 
      int m=Corr.nrow(),n=Corr.ncol();
      
      arma::mat C(Corr.begin(),m,n,false);//reuses memory and avoids extra copy
      //need an efficient way to get the submatrix off of this. Problem is that i j k not 
      //really represent a range of consecutive rows and columns :/
      //use Rinside? look at presentation a LOT of overhead
      
      //      arma::mat sub = C[c(i,j,k),c(i,j,k)]
      //cout << "C["<<i<<","<<j<<"] = " << C[i,j] << endl;
      
      std::vector<int> rows(k.size()+2);
      int l =0; // index in rows vector
      rows[l] = i; 
      l++;
      rows[l] = j; 
      l++;
      //cout << "i and j " << i << " and " << j << endl;
      
      std::vector<int>::iterator row;
 
      for (row = k.begin(); row !=k.end();++row)
	{
	  rows[l] = *row;
	  //cout << "k elem " << *row << endl;
	  
	  l++;
	}
      
      std::vector<int> cols = rows;
      
      //might be the big performance stumble block right here,
      //possibly needs to be optimized
      //cout << "C matrix " << m << "," << n<< endl;
      //cout << "row " << rows.front() << "," << rows.back() << endl;
      //cout << "end of k " << k.back() << endl;
      
      
      arma::mat sub = arma::submat(C,rows.begin(),rows.end(),cols.begin(),cols.end());

      arma::mat PM;
      
      try
	{
	  PM = arma::inv(sub);
	}
      catch(runtime_error re)
	{
	  //the matrix appears to be singular
	  cout << "Caught error yes mam!" << endl;
	  cout << "The matrix appears to be singular :s" << endl;

	  cout << "Some DEBUG info: "<<endl;
	  cout << "i = " << i << " j = " << j << endl;

	  for(int l = 0; l< k.size();l++)
	    {
	      cout << k[l] << " , ";
	    }
	  cout << endl;

	  cout << Corr << endl;

	  cout << sub << endl;
	  //otherwise we would just get a matrix index out of bounds error right below here.
	}
      //cout << "PM size " << PM.size() << endl;
      
      //the correlation matrix is always a positive semi definite matrix
      //inverse can be done faster if the matrix is a positive definite symmetric matrix
      //we specify this so: inv( sympd(sub) )
		    
      //PM <- pseudoinverse(C(c(i,j,k),c(i,j,k)))
      //return -PM[1,2]/sqrt(PM[1,1]*PM[2,2]);
      cout << "PM " << endl << PM << endl;
      
      r = -PM(0,1)/sqrt(PM(0,0)*PM(1,1));
      
    }
  //if(is.na(r)) r<-0
  if(R_IsNA(r))
    r = 0;
  
  //min(cut.at,max(-cut.at,r))
  return min(cutat,max(-cutat,r));
}

double zStat(int i,int j,std::vector<int> k,NumericMatrix Corr,long n)
{
  double r = pcorOrder(i,j,k,Corr);
  r = sqrt(n-k.size()-3.0)*(0.5*log((1.0+r)/(1.0-r)));
  //check if na? see R code
  return r;
}

double gaussCItest(int i,int j,std::vector<int> k,NumericMatrix Corr,long n)
{
  double z = zStat(i,j,k,Corr,n);
  return 2*stats::pnorm_0(fabs(z),0,0);
  //x=abs(z), lower.tail=FALSE, log.p=FALSE
}
  
/**
 * Generate the next set in a list of all possible sets of size k out of 1:n
 * 
 */
std::vector<int> getNextSet(int n, int k,std::vector<int> previous)
{
  //cout << "getnextset" << endl;//DEBUG

  /** initial implementation purely based on the R code, might be a faster way to do this*/
  int sum = 0;
  std::vector<int>::iterator row;
  int iter = n-k;
  //cout << "iter = " << iter << endl;
  
  for (row = previous.begin(); row !=previous.end();++iter,++row)
    {
      sum += ((iter - *row) == 0);
    }

  int chInd = k-sum;
  chInd = chInd -1; //working with c++ indexing, not R here
  //cout << "chInd = "<< chInd << endl;
  //cout << "k= " << k << " sum = " << sum << endl;
  
  if(chInd == -1 || k == 0)
    {
      //was last set to check, k == 0 there was no set actually ^^
      previous.resize(1); 
      previous[0] = -1; //marks finished
      //      cout << "=-=FINISHED WITH SETS" << endl;
      
    }else
    {
      //there is still a set to go
      //cout << "there is still a set to go"<<endl;
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


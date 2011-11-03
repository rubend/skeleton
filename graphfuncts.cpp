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
  
  for (int j = 0; j < p; ++j)
    {
      if (G[row*p + j] == true)
	{
	  connections[index] = j;
	  index++;  
	}
    }
  
  if (index != p-1)
    {
      //we didn\'t have p-1 connections
      connections[index] = -1; // signal end of connections
    }

  return;
}

/**
 * find a row with enough neighbours as we need (neighboursneeded),Return the row otherwise return -1 ==> couldn't find a row that satisfies this rule
 */
int getRowWithEnoughConnections(bool G[],int p,int neighboursneeded)
{ 
  int connectionsInRow = 0;
  //cout << "neighboursneeded = " << neighboursneeded << endl;
  for (int i = 0; i < p; ++i)
    {
      for (int j = 0; j < p; ++j)
	{
	  if(G[i*p+j] == true)
	    {
	      connectionsInRow++;
	    }
	}      
      if (connectionsInRow >= neighboursneeded)
	{
	  return i;
	}
      connectionsInRow = 0;
    }
  return -1; //means there isn't any!
}

/**
 * startrow = the row from which to start looking for a row with connections in the graph
 * G = the graph
 * p = the number of rows
 * returns an integer value that says which row still has connections
 */
int getNextRowWithConnections(int startrow,bool G[],int p)
{ 
  int returnrow = -1;
  for (int row = startrow; row < p; row++)
    {
      
      for (int i = 0; i < p; i++)
	{
	  
	  if(G[row*p+i] == 1)
	    {
	      returnrow = row;
	      return returnrow;
	    }
	}
    }
  
  //seems like there are no more rows >= startrows containing connections
  return returnrow;
}

/**
 * get the maximum number of neighbours a node in the graph has.
 * Might be faster to  keep track of how many neighbours each node has and updating and looking for new max when one decreases and so on ...
 */
int getMaxConnectionsLeft(bool G[],int p)
{
  int currentMax = 0;
  int connectionsInRow = 0;
  for (int i = 0; i < p; ++i)
    {
      for (int j = 0; j < p; ++j)
	{
	  if(G[i*p+j] == true)
	    {
	      connectionsInRow++;
	    }
	}
      if (connectionsInRow > currentMax)
	{
	  currentMax = connectionsInRow;
	}
      connectionsInRow =0;
    }
  return currentMax;  
}

/**
 * find a row with more neighbours then neighboursneeded, if you do find one return true. Else false.
 */
bool graphContainsNodeWithMoreNeighbours(bool G[],int p,int neighboursneeded)
{
  int connectionsInRow = 0;
  for (int i = 0; i < p; ++i)
    {
      for (int j = 0; j < p; ++j)
	{
	  if(G[i*p+j] == true)
	    {
	      connectionsInRow++;
	    }
	}
      if (connectionsInRow >= neighboursneeded)
	{
	  return true;
	}
      connectionsInRow =0;
    }
  return false;  
}


/**
 * get the other connection when excluding y
 */
void getOtherConnections(std::vector<int>* others,int j,int* connections,int p)
{
  (*others).resize(0);
  for (int i = 0; i < p-1; ++i)
    {
      if (connections[i] == -1) break;

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
 
  std::vector<int> subset(subsetind.size(),0);
  
  for (int i = 0; i < subsetind.size(); ++i)
    {     
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
	  //cout << "k elem " << *row << endl;
	  
	  l++;
	}
      
      std::vector<int> cols = rows;

      //might be the big performance stumble block right here,TODO check with callgrind
      arma::mat sub = arma::submat(C,rows.begin(),rows.end(),cols.begin(),cols.end());

      arma::mat PM;
      //arma::mat rhs;
      //rhs.eye(sub.n_cols,2);//we only need the first 2 columns of the inverse
      
      try
	{

	  //the correlation matrix is always a positive semi definite matrix
	  //inverse can be done faster if the matrix is a positive definite symmetric matrix
	  //we specify this so: inv( sympd(sub) )
	  PM = arma::inv(arma::sympd(sub));
	  //PM = arma::solve(arma::sympd(sub),rhs);
	}
      catch(runtime_error re)
	{
	  //the matrix appears to be singular
	  cout << "Caught error yes m'am!" << endl;
	  cout << "The matrix appears to be singular :s" << endl;

	  cout << "Some DEBUG info: "<<endl;
	  cout << "i = " << i << " j = " << j << endl;
	  cout << "k = ";
	  for(int l = 0; l< k.size();l++)
	    {
	      cout << k[l] << " , ";
	    }
	  cout << endl;

	  cout << Corr << endl;

	  cout << sub << endl;
	  //otherwise we would just get a matrix index out of bounds error right below here.
	}

      // cout << "PM " << endl << PM << endl;      
      r = -PM(0,1)/sqrt(PM(0,0)*PM(1,1));
      
    }
  if(R_IsNA(r))
    r = 0;
  
  return min(cutat,max(-cutat,r));
}

double zStat(int i,int j,std::vector<int> k,NumericMatrix Corr,long n)
{
  double r = pcorOrder(i,j,k,Corr);
  r = sqrt(n-k.size()-3.0)*(0.5*log((1.0+r)/(1.0-r)));
  if(r != r)//r is NaN
    r = 0.0;
  
  return r;
}

double gaussCItest(int i,int j,std::vector<int> k,NumericMatrix Corr,long n)
{
  double z = zStat(i,j,k,Corr,n);
  return 2*stats::pnorm_0(fabs(z),0,0);
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
  int iter = n-k;
  
  for (row = previous.begin(); row !=previous.end();++iter,++row)
    {
      sum += ((iter - *row) == 0);
    }

  int chInd = k-sum;
  chInd = chInd -1; //working with c++ indexing, not R here
  
  if(chInd == -1 || k == 0)
    {
      //was last set to check, k == 0 there was no set actually ^^
      previous.resize(1); 
      previous[0] = -1; //marks finished     
    }else
    {
      //there is still a set to go
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

void printOutG(bool G[],int p)
{
  cout << "Printout G:" << endl;
  
  for (int i = 0; i < p; ++i)
    {
      for (int j = 0; j < p; ++j)
	{
	  cout << G[i*p+j] << " " << endl;
	}
      cout << endl;
    }
  cout << "end of G" << endl;
}


#include <iostream>
#include <algorithm>
#include <Rcpp.h>
#include <Rcpparmadillo.h>
#include <vector>
//#include <submat.h>
//#include <graphfuncts.h>

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

double pcorOrder(int,int,std::vector<int>,NumericMatrix);
std::vector<int> getNextSet(int,int,std::vector<int>);


/**
 * function for computing the skeleton in our DAG
 * arguments
 * n number of samples
 * p number of variables
 * m_max max order to keep looking for
 */
RcppExport SEXP skeleton(SEXP pt,SEXP alphat,SEXP m_maxt, SEXP C)
{
  // now this main will be implemented as the skeleton function! :p
  
  int p = as<int>(pt);
  double alpha = as<double>(alphat);
  int m_max = as<int>(m_maxt);
  NumericMatrix Corr(C);
  
  //using c++ datatypes and trying to write all functions myself, 
  //i.e. not calling R from c++
  //how to store graph? Use matrix?
  bool G[p*p];// only store upper triangle? -> store as vector?
  
  //all connections exist
  cout << "Got till the initialisation" << endl;
  
  initialiseGraph(G,p);

  int row = 0;// the current row we are studying the connections of
  int connections[p-1]; // stores the connections the current row has (maximum p-1)
  //static allocation ==> efficient
  int sizeothers,x,y;
  double pval;
  std::vector<int> subset(0);//declare in beginning and increase size every loop run
  std::vector<int> others(0);
  for (int ord = 0; ord <= m_max; ++ord)
    {
      //look for next row with connections then iterate over the remaining connections
      //alternative: save all connections explicitly in double array and loop over those.
      
      row = getNextRowWithConnections(row,G,p); // row == x
      x = row;
      
      while(row != -1)
	{
	  getRowConnections(row,G,p,connections); // getting the connections belonging to the row
	  for (int i = 0; i < p-1; ++i)
	    {
	      //one of the remaining edge tests
	      y=connections[i];
	      
	      //so now we are gonna check the correlation between row and connections[i]
	      //in respect to every other subset of lengths ord of the remaining connections
	      if(y == -1) break;//reached end of connections
	      //y == connections[i]
	      getOtherConnections(&others,i,connections,p);
	      sizeothers = others.size();
	      if (sizeothers < ord)
		{
		  continue;//goto next loop iteration
		}
	      
	      //initial subset, TODO is there a more efficient way? Builtin way for this?
	      getSeqVector(&subset,ord);//get sequence vector of size ord in subset 
	      
	     
	      while(subset[0] != -1)
		{
		  std::vector<int> k = getSubset(others,subset); //does this work? check! otherwise write own function that does this.
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
  
  //convert G to a logicalMatrix for returning to R?
  
  
  return wrap(convertToLogical(G,p)); // return graph matrix, in later stage return a more complete object
}
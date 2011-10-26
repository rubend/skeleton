// now this main will be implemented as the skeleton function! :p
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
int p = as<int>(pt);
double alpha = as<double>(alphat);
int m_max = as<int>(m_maxt);
NumericMatrix Corr(C);
long n = as<long>(nt);
  
bool G[p*p];// only store upper triangle? -> store as vector?
  
initialiseGraph(G,p);

int row;// the current row we are studying the connections of
int connections[p-1]; // stores the connections the current row has (maximum p-1)
//static allocation ==> efficient
int sizeothers,x,y;
double pval;
std::vector<int> subset(0,0);//declare in beginning and increase size every loop run
std::vector<int> others(0);
int iter =0;

//SPEEDUP: parallelise with openmp!
for (int ord = 0; ord <= m_max; ++ord)
  //additional check for if any connections are left might be useful, especially if mmax is huge
  {
    //cout << "ord = " << ord << endl;
    
    //look for next row with connections then iterate over the remaining connections
    //alternative: save all connections explicitly in double array and loop over those.
    row = 0; // reset to zero
    
    row = getNextRowWithConnections(row,G,p); // row == x
      
    while(row != -1)
      {
	x = row; // just temporary variable for clarity when calling pcororder
	getRowConnections(row,G,p,connections); // getting the connections belonging to the row
	for (int i = 0; i < p-1; ++i)
	  {
	    //one of the remaining edge tests
	    y=connections[i];
	      
	    //so now we are gonna check the correlation between row and connections[i]
	    //in respect to every other subset of lengths ord of the remaining connections
	    if(y == -1) break;//reached end of connections
	    getOtherConnections(&others,i,connections,p);
	    sizeothers = others.size();

	    if (sizeothers < ord)
	      {
		continue;//goto next loop iteration
	      }
	      
	      
	    getSeqVector(&subset,ord);//get sequence vector of size ord in subset

	    while(subset.size() == 0 || subset[0] != -1)
	      {
		//subset.size() == 0 is possible if ord == 0, subset[0] == -1 is the stop condition
		
		std::vector<int> k = getSubset(others,subset); 

		pval = gaussCItest(x,y,k,Corr,n);
		
		iter++;
		if (pval >=  alpha)
		  {
		    //independent
		    G[p*x+y]=false;G[p*y+x]=false;
	
		    break; //no more checking to be done
		  }
		subset = getNextSet(sizeothers,ord,subset);  
	      }
	  }
	
	row = getNextRowWithConnections(row+1,G,p);
      }
  }

return wrap(convertToLogical(G,p)); // return graph matrix, in later stage return a more complete object

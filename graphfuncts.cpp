
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

/**
 * This file contains functions doing operations on the graph
 *
 */

/**
 * initialise our graph with all connections on
 *
 *
 */
void initialiseGraph(bool* G,int p)
{
  for (int i = 0; i < p*p; ++i)
    {
      G[i] = true;
    }

  for (int i = 0; i < p; ++i)
    {
      //diagonal connections obviously don't exist
      G[i*p+i] = false;
    }
}

/**
 * any method looks if there is still a connection in the graph matrix
 */
bool any(bool *G,int p)
{
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
void getRowConnections(int row,bool* G,int p,int* connections)
{
  int index = 0;//keeping track of the index in the connections vector
  
  //there is still a row >= startrow with connections
  for (int j = 0; j < p; ++j)
    {
      if (G[row*p+j] == true)
	{
	  connections[index] = j;
	  index++;  
	}
    }
  
  if (index != p-1)
    {
      //we didn't have p-1 connections
      connections[index+1] = -1; // signal end
    }

  return;
}

/**
 * startrow = the row from which to start looking for a row with connections in the graph
 * G = the graph
 * p = the number of rows
 * returns an integer value that says which row still has connections
 */
int getNextRowWithConnections(int startrow,bool* G,int p)
{
  for (int row = startrow; row < p; ++row)
    {
      for (int i = row; i < p; ++i)
	{
	  //only search through the upper triangular
	  if(G[row*p+i] == true)
	    {
	      //we have a hit!
	      return row;
	    }
	}
    }
  
  //seems like there are no more rows >= startrows containing connections
  return -1;
}

/**
 * get the other connection when excluding y
 */
std::vector<int> getOtherConnections(int j,int* connections,int p)
{
  std::vector<int> others(0);

  for (int i = 0; i < p-1; ++i)
    {
      if (connections[i] == -1) break;

      //we don't want y==connections[j] in here
      if (i != j)
	{
	  others.push_back(connections[i]);
	}
    }
  
  return others;
}

/**
 * Creates a vector of size ord with elements 1:ord
 *
 */
std::vector<int> getSeqVector(int ord)
{
  std::vector<int> seq(ord);
  for (int i = 0; i < ord; ++i)
    {
      seq[i] = i+1;
    }
  return seq;
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
 *converts the boolian matrix to a logical matrix so we can easily return it to R
 */
LogicalMatrix convertToLogical(bool* G,int p)
{
  LogicalMatrix log(p,p);
  
  for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j<p; ++j)
	{
	  log(i,j) = G[i*p+j];
	}
    }
  return log;
  
}

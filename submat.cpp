#include <Rcpp.h>

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

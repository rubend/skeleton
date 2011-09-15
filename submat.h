template <typename T, typename InputIterator> Mat<T> submat(const Mat<T>& input,InputIterator firstRow, InputIterator lastRow,InputIterator firstCol, InputIterator lastCol);

  template <typename T, typename InputIterator> Col<T> subvec(const Mat<T>& input,InputIterator firstRow, InputIterator lastRow,const unsigned int colind);




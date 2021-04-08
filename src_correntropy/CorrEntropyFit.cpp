#include "correntropy_regression.h"



int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CorrEntropyFit [DATAI] [DATAO]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      return -1;
    }
    using T = double;
    std::string FileI = argv[1];
    std::string FileO = argv[2];
    std::ifstream is(FileI);
    T sigma;
    is >> sigma;
    MyVector<T> B = ReadVector<T>(is);
    MyMatrix<T> M = ReadMatrix<T>(is);
    MyVector<T> Beta = Compute_CorrEntropy_Regression(B, M, sigma);
    //
    std::ofstream os(FileO);
    int len = Beta.size();
    for (int i=0; i<len; i++)
      os << Beta(i) << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }

}

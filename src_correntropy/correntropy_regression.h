#ifndef INCLUDE_CORRENTROPY_REGRESSION_H
#define INCLUDE_CORRENTROPY_REGRESSION_H


#include "MAT_Matrix.h"




/*
  We have a number of function f_j for 1 <= j <= m
  A number of values x_i for 1 <= i <= n and a fitting value b_i
  ---
  So, we want to find the coefficients beta_j such that
  b_i is near sum_j beta_j f_j(x_i)
  The matrix encoding is M(i,j) = f_j(x_i)
  ---
  There are many ways to do this fitting, here we consider the correntropy.
  That is we want to maximize
  sum_i g(b_i, sum_j beta_j f_j(x_i))
  with g(x, y) = exp(-(x-y)^2 / 2 \sigma^2)
  ---
  The classical algorithm is the half-quadratic optimization.
  However, it may not be relevant here as we do not have a second term.
  Half-Quadratic is a generalized case.
  ---
  F(beta) = sum_i g(b_i - sum_j beta_j f_j(x_i))
          = sum_i g(b_i - sum_j beta_j M(i,j) )
  dF/dbeta_k = sum_i - M(i,k) g'(b_i - sum_j beta_j M(i,j) )
  d^2F(dbeta_k dbeta_l) = sum_i  M(i,k) M(i,l) g''(b_i - sum_j beta_j M(i,j) )
  ---
  The function is expressed as
  f(h) = a + b h + Q[h] / 2
  f(h_0 + h) = a + b h_0 + bh + Q[h_0] + 2h_0 Q h + Q[h]
             = C + (b + h_0 Q) h + Q[h]
  The direction is computed as - b Q^{-1}
 */
template<typename T>
MyVector<T> Compute_CorrEntropy_Regression(MyVector<T> const& B, MyMatrix<T> const& M, double const& sigma)
{
  int n = M.rows();
  int m = M.cols();
  int n_B = B.size();
  std::cerr << "n=" << n << " m=" << m << " n_B=" << n_B << "\n";
  if (n != n_B) {
    std::cerr << "We should have n=" << n << " n_B=" << n_B << " equal\n";
    throw TerminalException{1};
  }
  //
  double alpha = 1 / (2 * sigma * sigma);
  auto g=[&](T const& x) -> T {
    return exp(- alpha  * x * x);
  };
  auto gD1=[&](T const& x) -> T {
    return -2 * alpha * x * exp(- alpha  * x * x);
  };
  auto gD2=[&](T const& x) -> T {
    return (-2 * alpha  + 4 * alpha * alpha * x * x) * exp(- alpha  * x * x);
  };
  //
  auto f=[&](MyVector<T> const& beta) -> T {
    T sum=0;
    for (int i=0; i<n; i++) {
      T delta = B(i);
      for (int j=0; j<m; j++)
        delta -= beta(j) * M(i,j);
      sum += g(delta);
    }
    return sum;
  };
  auto ComputeGradient=[&](MyVector<T> const& beta) -> MyVector<T> {
    MyVector<T> grad(m);
    for (int k=0; k<m; k++) {
      T sum = 0;
      for (int i=0; i<n; i++) {
        T delta = B(i);
        for (int j=0; j<m; j++)
          delta -= beta(j) * M(i,j);
        sum += - M(i,k) * gD1(delta);
      }
      grad(k) = sum;
    }
    return grad;
  };
  auto ComputeHessian=[&](MyVector<T> const& beta) -> MyMatrix<T> {
    MyMatrix<T> Hess(m,m);
    for (int k=0; k<m; k++) {
      for (int l=0; l<m; l++) {
        T sum = 0;
        for (int i=0; i<n; i++) {
          T delta = B(i);
          for (int j=0; j<m; j++)
            delta -= beta(j) * M(i,j);
          sum += M(i,k) * M(i,l) * gD2(delta);
        }
        Hess(k,l) = sum;
      }
    }
    return Hess;
  };
  //
  // The iteration algorithm
  //
  MyVector<T> beta = ZeroVector<T>(m);
  int n_iter = 100;
  int i_iter=0;
  while(true) {
    T eVal = f(beta);
    std::cerr << "i_iter=" << i_iter << " eVal=" << eVal << "\n";
    T eFact = 0.8;
    auto DirUpgrade=[&](MyVector<T> const& edir) -> bool {
      T eScal = 1;
      int itermax = 1000;
      int iter=0;
      while(true) {
        MyVector<T> betaN = beta + eScal * edir;
        T eValN = f(betaN);
        //        std::cerr << "   eScal=" << eScal << " eVal=" << eVal << " eValN=" << eValN << "\n";
        if (eValN > eVal) {
          beta = betaN;
          return true;
        }
        eScal *= eFact;
        iter++;
        if (iter > itermax)
          return false;
      }
    };
    //
    MyVector<T> grad = ComputeGradient(beta);
    MyMatrix<T> hess = ComputeHessian(beta);
    MyVector<T> dir = - hess.colPivHouseholderQr().solve(grad);
    //    std::cerr << "We have dir\n";
    T eScal = dir.dot(grad);
    //    std::cerr << "We have eScal\n";
    bool reply;
    if (eScal > 0)
      reply = DirUpgrade(dir);
    else
      reply = DirUpgrade(grad);
    if (!reply)
      break;
    //
    i_iter++;
    if (i_iter == n_iter)
      break;
  }
  return beta;
}





#endif

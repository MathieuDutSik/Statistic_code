#ifndef INCLUDE_CORRENTROPY_REGRESSION_H
#define INCLUDE_CORRENTROPY_REGRESSION_H







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
 */
template<typename T>
std::vector<T> Computre_Regression(MyVector<T> const& B, MyMatrix<T> const& M, double const& sigma)
{
  double alpha = -1 / (2 * sigma * sigma);
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
  int n = M.rows();
  int m = M.cols();
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

}





#endif

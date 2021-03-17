#ifndef INCLUDE_CORRENTROPY_REGRESSION_H
#define INCLUDE_CORRENTROPY_REGRESSION_H







/*
  We have a number of function f_j for 1 <= j <= m
  A number of values x_i for 1 <= i <= n and a fitting value b_i
  ---
  So, we want to find the coefficients beta_j such that
  b_i is near sum_j beta_j f_j(x_i)
  ---
  There are many ways to do this fitting, here we consider the correntropy.
  That is we want to maximize
  sum_i g(b_i, sum_j beta_j f_j(x_i))
  with g(x, y) = exp(-(x-y)^2 / 2 \sigma^2)
  ---
  The classical algorithm is the half-quadratic optimization.
  However, it may not be relevant here as we do not have a second term.
  Half-Quadratic is only for 
 */
template<typename T>
std::vector<T> Computre_Regression(MyVector<T> const& ListValue, MyMatrix<T> const& ListFct, double const& sigma)
{
  auto ComputeGradient=[&](MyVector<T> const& beta) -> MyVector<T> {
                                                                    
  };

}





#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            HOMEWORK #3                  %%%
%%%                       Donggyu Kim       %%%
%%%-----------------------------------------%%%
%%%                      Language: Octave   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read data
if (opt == 1)  
  %------------------%
  %   Problem 1 (e)  %
  %------------------%
  load classification.test;
else
  %------------------------%
  %   Problem 1 (a) - (d)  %
  %------------------------%
  load classification.dat;
end

%------------------%
%   Problem 1 (a)  %
%------------------%


[N, m] = size(classification);

x1 = classification(:,1);
x2 = classification(:,2);
y = classification(:,3);

A = zeros(0,0); % plotted as "o"
B = zeros(0,0); % plotted as "x"
A_row = 1;
B_row = 1;

for i = 1:N
  if (y(i) == 0)
    A(A_row, 1) = x1(i);
    A(A_row, 2) = x2(i);
    A_row = A_row + 1;
  else
    B(B_row, 1) = x1(i);
    B(B_row, 2) = x2(i);
    B_row = B_row + 1;
  endif
endfor

hold on;
plot(A(:,1), A(:,2), 'o');
plot(B(:,1), B(:,2), 'x');

%------------------%
%   Problem 1 (b)  %
%------------------%

% Maximum likelyhood estimates
pi_ML = mean(y);
mu11_ML = sum(y .* x1) / sum(y);
mu21_ML = sum(y .* x2) / sum(y);
mu10_ML = sum((ones(N, 1) - y).*x1) / sum(ones(N, 1) - y);
mu20_ML = sum((ones(N, 1) - y).*x2) / sum(ones(N, 1) - y);
sigma1_ML = (sum(y .* ((x1 - mu11_ML).^2)) + sum((ones(N, 1) - y) .* (x1 - mu10_ML).^2)) / N;
sigma2_ML = (sum(y .* ((x2 - mu21_ML).^2)) + sum((ones(N, 1) - y) .* (x2 - mu20_ML).^2)) / N;
sigma12_ML = (sum(y .* (x1 - mu11_ML) .* (x2 - mu21_ML)) + sum((ones(N, 1) - y) .* (x1 - mu10_ML) .* mu20_ML)) / N;

%{
problem 1 (b)
-------------
pi_ML =  0.50000
mu11_ML = -4.4740
mu21_ML = -4.4585
mu10_ML =  0.070360
mu20_ML = -0.045253
sigma1_ML =  7.4447
sigma2_ML =  19.172
sigma12_ML =  9.5255

problem 1 (e)
-------------
pi_ML =  0.52000
mu11_ML = -4.4462
mu21_ML = -4.5071
mu10_ML = -0.19159
mu20_ML = -0.15322
sigma1_ML =  7.2495
sigma2_ML =  17.959
sigma12_ML =  9.2556
%}

plot(mu11_ML, mu21_ML, 'r');
plot(mu10_ML, mu20_ML, 'g');

global Sigma = [sigma1_ML, sigma12_ML ; sigma12_ML, sigma2_ML];
global mu1 = [mu11_ML ; mu21_ML];
global mu0 = [mu10_ML ; mu20_ML];
global beta_ = (Sigma^-1) * (mu1 - mu0);
global gamma_ = -0.5 * (mu1 - mu0)' * (Sigma^-1) * (mu1 + mu0) - log(pi_ML / (1 - pi_ML));

%{
Problem 1 (b)
-------------
beta_ =

  -0.61042
  -0.23020

gamma_ = -1.8624

Problem 1 (e)
-------------
beta_ =

  -0.58689
  -0.24243

gamma_ = -2.0059
%}

% Posterior prob
function retval = gen_p(x)
  global Sigma;
  global mu1;
  global mu0;
  global beta_;
  global gamma_;
  retval = 1 / (1 + exp(-beta_' * x - gamma_));
endfunction

% Find points whose posterior prob is 0.5
function points = gen_p_eq_0_5(X1, X2)
  points = zeros(0, 0);
  points_idx = 1;
  [m, n] = size(X1);
  for i = 1:m
    for j = 1:n
      x = [X1(i, j) ; X2(i, j)];
      p = gen_p(x);
      if (abs(p - 0.5) < 0.1^3)
        points(points_idx, 1) = x(1);
        points(points_idx, 2) = x(2);
        points_idx = points_idx + 1; 
      endif
    endfor
  endfor
endfunction

% Draw a red line
x1_range = linspace(-15,5);
x2_range = linspace(-25,5);
[X1, X2] = meshgrid(x1_range, x2_range);
line_b = gen_p_eq_0_5(X1, X2);
plot(line_b(:,1), line_b(:,2), 'r');

%------------------%
%   Problem 1 (c)  %
%------------------%

function z = eta(theta, x)
  z = x * theta;
endfunction

function z = mu(theta, x)
  z = 1 / (1 + exp(-eta(theta, x)));
endfunction

% IRLS algorithm
function theta = logistic_regression(X, y)
  theta = zeros(3,1);
  theta_0 = zeros(3,1);
  lim = 10^-6;
  N = size(y);
  t = 0;
  while t == 0 || norm(theta - theta_0) > lim
    % compute logistic functions
    mu_ = zeros(0);
    for n = 1:N
      x_n = X(n,:);
      mu_(n, 1) = mu(theta, x_n);
    endfor
    % compute W(t)
    W = zeros(0, 0);
    for n = 1:N
      W(n,n) = mu_(n, 1) * (1 - mu_(n, 1));
    endfor
    % compute theta(t+1)
    t = t + 1;
    theta_0 = theta;
    theta = ((X'*W*X)^-1)*(X'*W*X*theta + X'*(y - mu_));
  end
endfunction

X = [x1, x2, ones(size(y))];
global theta = logistic_regression(X, y);

%{
Problem 1 (c)
-------------
theta =

  -3.65885
  -0.22892
  -4.76961

Problem 1 (e)
-------------
theta =

  -3.36062
  -0.28787
  -4.36076
%}

function points = logistic_p_eq_0_5(X1, X2)
  global theta;
  points = zeros(0, 0);
  points_idx = 1;
  [m, n] = size(X1);
  for i = 1:m
    for j = 1:n
      x = [X1(i, j), X2(i, j), 1];
      p = mu(theta, x);
      if (abs(p - 0.5) < 0.1^2)
        points(points_idx, 1) = x(1);
        points(points_idx, 2) = x(2);
        points_idx = points_idx + 1; 
      endif
    endfor
  endfor
endfunction

% Draw a green line
line_c = logistic_p_eq_0_5(X1, X2);
plot(line_c(:,1), line_c(:,2), 'g');

%------------------%
%   Problem 1 (d)  %
%------------------%

% solving the normal equation
function theta = linear_regression(X, y)
  theta = ((X'*X)^-1)*X'*y;
endfunction

global theta_ = linear_regression(X, y);

%{
Problem 1 (d)
-------------
theta_ =

  -0.119579
   0.025781
   0.294764

Problem 1 (e)
-------------
theta_ =

  -0.121700
   0.026638
   0.291827
%}

function points = linear_p_eq_0_5(X1, X2)
  global theta_;
  points = zeros(0, 0);
  points_idx = 1;
  [m, n] = size(X1);
  for i = 1:m
    for j = 1:n
      x = [X1(i, j) ; X2(i, j) ; 1];
      p = theta_'*x;
      if (abs(p - 0.5) < 0.1^3)
        points(points_idx, 1) = x(1);
        points(points_idx, 2) = x(2);
        points_idx = points_idx + 1; 
      endif
    endfor
  endfor
endfunction

% Draw a black line
line_d = linear_p_eq_0_5(X1, X2);
plot(line_d(:,1), line_d(:,2), 'k');

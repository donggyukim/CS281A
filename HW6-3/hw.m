%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            HOMEWORK #6-3                %%%
%%%                       Donggyu Kim       %%%
%%%-----------------------------------------%%%
%%%                      Language: Octave   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read data
load hmm_gauss.dat;
load hmm_test.dat;
global M = 4; % state #
global n = 0;
[T, n] = size(hmm_gauss);

% gaussian emission
function retval = p(y_t, mu_i, sigma2_i)
  global n;
  retval = exp(-((y_t - mu_i)' * (y_t - mu_i))/(2*sigma2_i)) / sqrt(((2*pi) ^ n) * sigma2_i);
endfunction

%%%%%%%%%%%%%%%%%%%%%%
%% Prob 3 (b)       %%
%% HMM EM algorithm %%
%%%%%%%%%%%%%%%%%%%%%%
% M x M transition matrix
A = ones(M,M) / M;
% initial transition
pi_ = ones(1,M) / M;
% mean parameters
mu = zeros(n,M);
prev_mu = zeros(n,M);
% variance parameters
sigma2 = ones(1,M);
prev_sigma2 = zeros(1,M);

% gaussian emission
function retval = p(y_t, mu_i, sigma2_i)
  global n;
  retval = exp(-((y_t - mu_i)' * (y_t - mu_i))/(2*sigma2_i)) / sqrt(((2*pi) ^ n) * sigma2_i);
endfunction

% EM algorithm
delta = 1;
lim = 0.1^3;
while (delta > lim)
  %%%%%%%%%%%%%%
  %%% E step %%%
  %%%%%%%%%%%%%%

  % compute alphas
  alpha = zeros(0);
  y_0 = hmm_gauss(1,:)';
  for i = 1 : M
    alpha(1, i) = pi_(i) * p(y_0, mu(:,i), sigma2(i));
  endfor
  % normalization
  alpha(1, :) = alpha(1, :) / sum(alpha(1, :));

  for t = 1 : T-1
    y_t_1 = hmm_gauss(t+1,:)';
    for j = 1 : M
      alpha(t+1, j) = 0;
      for i = 1 : M
        alpha(t+1, j) = alpha(t+1, j) + alpha(t, i) * A(i, j) * p(y_t_1, mu(:,j), sigma2(j)); 
      endfor
    endfor
    alpha(t+1,:) = alpha(t+1,:) / sum(alpha(t+1,:));
  endfor

  % compute gammas
  gamma = zeros(0);
  gamma(T,:) = alpha(T,:);
  for k = 1 : T-1
    t = T - k;
    for i = 1 : M
      gamma(t, i) = 0;
      for j = 1 : M
        sum_alpha_t_a = 0;
        for l = 1 : M
          sum_alpha_t_a = sum_alpha_t_a + alpha(t, l) * A(l, j);
        endfor
        gamma(t, i) = gamma(t, i) + (alpha(t, i) * A(i, j) * gamma(t+1, j)) / sum_alpha_t_a;
      endfor
    endfor
  endfor

  % compute xi
  xi = zeros(0);
  for t = 1 : T - 1
    for i = 1 : M
      for j = 1 : M
        y_t_1 = hmm_gauss(t+1,:)';
        xi(i, j, t) = alpha(t, i) * gamma(t+1, j) * p(y_t_1, mu(:,j), sigma2(j)) * A(i, j);
        xi(i, j, t) = xi(i, j, t) / alpha(t+1, j);
      endfor
    endfor
  endfor

  %%%%%%%%%%%%%%
  %%% M step %%%
  %%%%%%%%%%%%%%
  delta = 0;

  % compute pi_
  for i = 1 : M
    pi_(i) = gamma(1, i);
  endfor

  % compute A
  for i = 1 : M
    for j = 1 : M
      sum_xi = 0;
      sum_gamma = 0;
      for t = 1 : T - 1
        sum_xi = sum_xi + xi(i, j, t);
        sum_gamma = sum_gamma + gamma(t, i);
      endfor
      A(i, j) = sum_xi / sum_gamma;
    endfor
  endfor

  % compute mu
  for i = 1 : M
    prev_mu(:,i) = mu(:,i);
    sum_gamma_y_t = [0 ; 0];
    sum_gamma = 0;
    for t = 1 : T
      y_t = hmm_gauss(t,:)';
      sum_gamma_y_t = sum_gamma_y_t + gamma(t, i) * y_t;
      sum_gamma = sum_gamma + gamma(t, i);
    endfor
    mu(:,i) = sum_gamma_y_t / sum_gamma;
    delta += norm(mu(:,i) - prev_mu(:,i));
  endfor

  % compute sigma2
  for i = 1 : M
    prev_sigma2(i) = sigma2(i);
    sum_gamma_y_t_mu_i = 0;
    sum_gamma = 0;
    for t = 1 : T
      y_t = hmm_gauss(t, :)';
      mu_i = mu(:,i);
      sum_gamma_y_t_mu_i = sum_gamma_y_t_mu_i + gamma(t, i) * (y_t - mu_i)' * (y_t - mu_i);
      sum_gamma = sum_gamma + gamma(t, i);
    endfor
    sigma2(i) = sum_gamma_y_t_mu_i / sum_gamma;
    delta += norm(sigma2(i) - prev_sigma2(i));
  endfor
end

mu
sigma2

%{
mu =

  -0.252091  -0.252091  -0.252091  -0.252091
   0.029432   0.029432   0.029432   0.029432

sigma2 =

   3.1308   3.1308   3.1308   3.1308
%}
   

hold on;
plot(hmm_test(:,1), hmm_test(:,2), 'x');
plot(mu(:,1), mu(:,2), 'o');


%%%%%%%%%%%%%%%%%%%%%%
%% Prob 3 (c)       %%
%% GMM EM algorithm %%
%%%%%%%%%%%%%%%%%%%%%%

% initial transition
pi_ = rand(1,M);
pi_ = pi_ / sum(pi_);
% mean parameters
mu = rand(n,M)*10;
% variance parameters
sigma2 = rand(1,M)*10;
% expectation
tau = zeros(T,M);
tau_prev = zeros(T,M);
delta = 1;
lim = 0.1^6;

while (delta > lim)
  %%%%%%%%%%%%%%
  %%% E step %%%
  %%%%%%%%%%%%%%
  delta = 0;

  % compute tau
  for t = 1 : T
    tau_prev(t,:) = tau(t,:);
    for i = 1 : M
      y_t = hmm_gauss(t,:)';
      tau(t,i) = pi_(i) * p(y_t, mu(:,i), sigma2(i));
    endfor
    delta = delta + norm(tau(t,:) - tau_prev(t,:));
  endfor

  %%%%%%%%%%%%%%
  %%% M step %%%
  %%%%%%%%%%%%%%
  %delta = 0;

  % compute pi_
  for i = 1 : M
    pi_(i) = sum(tau(:,i)) / T;
  endfor

  % compute mu
  for i = 1 : M
    prev_mu(:,i) = mu(:,i);
    sum_tau_y = [0;0]; 
    sum_tau = 0;
    for t = 1 : T
      y_t = hmm_gauss(t,:)';
      sum_tau_y = sum_tau_y + tau(t,i) * y_t;
      sum_tau = sum_tau + tau(t,i);
    endfor
    mu(:,i) = sum_tau_y / sum_tau;
    delta += norm(mu(:,i) - prev_mu(:,i));
  endfor

  % compute sigma2
  for i = 1 : M
    prev_sigma2(i) = sigma2(i);
    sum_tau_y_t_mu_i = 0;
    sum_tau = 0;
    for t = 1 : T
      y_t = hmm_gauss(t, :)';
      mu_i = mu(:,i);
      sum_tau_y_t_mu_i = sum_tau_y_t_mu_i + tau(t, i) * (y_t - mu_i)' * (y_t - mu_i);
      sum_tau = sum_tau + tau(t, i);
    endfor
    sigma2(i) = sum_tau_y_t_mu_i / sum_tau;
    delta += norm(sigma2(i) - prev_sigma2(i));
  endfor
end

mu
sigma2

%{
mu =

  -0.603154  -0.604435  -0.624944  -0.559772
  -0.060845  -0.050238  -0.075793  -0.046821

sigma2 =

   1.8837   1.8818   1.8381   1.9670
%}

plot(mu(:,1), mu(:,2), 'o');
hold off;

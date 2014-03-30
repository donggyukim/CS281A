%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            HOMEWORK #3                  %%%
%%%                       Donggyu Kim       %%%
%%%-----------------------------------------%%%
%%%                      Language: Octave   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------%
%   Problem 1 (a)  %
%------------------%

% read data
load lms.dat;

X = lms(:,1:2);
y = lms(:,3);

% solve normal equation
theta_gold = inv(transpose(X) * X) * transpose(X) * y;

%{
theta_gold =

   5.0420
  -1.4006
%}

%------------------%
%   Problem 1 (b)  %
%------------------%

% covariant matrix = X^T * X
cov_mat = transpose(X) * X;

% lambdas: eigenvalues
lambdas = eig(cov_mat)

%{
lambdas =

    8.0934
   22.9518
%}

% V: eigenvectros
[V, D] = eig(cov_mat)

%{
V =

  -0.84687  -0.53181
  -0.53181   0.84687

D =

Diagonal Matrix

    8.0934         0
         0   22.9518
%}

% Cost function J
function z = J(X, y, theta)
	Y = y - X * theta;
	z = transpose(Y) * Y;
endfunction

% Contours of J in the parameter space
l1 = linspace(-5,15);
l2 = linspace(-10,10);
[A, B] = meshgrid(l1, l2);

figure('Name', 'Parameter Space Contour of J');
contour(A, B, J(X, y, [l1; l2]))

%------------------%
%   Problem 1 (c)  %
%------------------%

% LMS algorithm
function theta = LMS(X, y, step_size);
	figure;
	lim = 10^-9;
	theta = [0; 0];
	delta = [1; 1];
	while norm(delta) > lim
		delta = [0; 0];
		for n = 1 : size(y);
			x_n = transpose(X(n,:));
			y_n = y(n);
			delta = delta + (y_n - transpose(theta) * x_n) * x_n;
		endfor
		delta = step_size * delta;
		theta = theta + delta;
		plot(theta)
	end	
endfunction

LMS(X, y, 1/max(lambdas))
%{
ans =

   5.0420
  -1.4006
%}

LMS(X, y, 1/(2*max(lambdas)))
%{
ans =

   5.0420
  -1.4006
%}

LMS(X, y, 1/(4*max(lambdas)))
%{
ans =

   5.0420
  -1.4006
%}
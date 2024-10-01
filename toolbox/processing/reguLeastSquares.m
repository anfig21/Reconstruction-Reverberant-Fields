function [x,lambda] = reguLeastSquares(H,p,method,plotFlag)
%[x,lambda] = reguLeastSquares(H,p,method,plotFlag) Obtains the estimated
%coefficients to the problem
%                   min ||p - Hx||  s.t. ||x||_2
%via Tikhonov regularisation.
%   Input:
%       - H         : dictionary of orthonormal columns. M x N
%       - p         : measurement vector. M x 1
%       - method    : specifies the regularisation method:
%                       'lcurve' (Default)
%                       'gcv'
%                       'quasiopt'
%       - plotFlag  : 'true' to L-curve
%                     'false' to avoid plotting. Default value
%   Output:
%       - x         : column vector of coefficients. N x 1
%       - lambda    : optimal regularisation parameter. Scalar
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 2, error('reguLeastSquares Error: Not enough input parameters.'), end
if nargin < 3 || isempty(method), method = 'lcurve'; end
if nargin < 4, plotFlag = false; end

%% MAIN CODE
[U,s,V] = csvd(H);

switch method
    case 'lcurve', [lambda,~,~,~] = l_curve(U,s,p,'Tikh',[],[],plotFlag);
    case 'gcv', [lambda,~,~] = gcv(U,s,p,'Tikh');
    case 'quasiopt', [lambda,~,~] = quasiopt(U,s,p,'Tikh');
end

[x,~,~] = tikhonov(U,s,V,p,lambda);

% Picard condition
% figure, picard(U,s,p);

end


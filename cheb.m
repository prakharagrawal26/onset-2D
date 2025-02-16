%It gives differentiation matrix in physical space. Will be used only if
%using collocation scheme.
% CHEB compute D = differentiation matrix, x = Chebyshev grid
function [D,x] = cheb(N)
% N=5;
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D - diag(sum(D'));
% S=[1,0,-0.5,0,0,0;
%    0,0.5,0,-0.5,0,0;
%    0,0,0.5,0,-0.5,0;
%    0,0,0,0.5,0,-0.5;
%    0,0,0,0,0.5,0;
%    0,0,0,0,0,0.5];
function [am] = FourierApprox (ak, m, n,N) 
%FourierApprox gives an approximation to the inverse transformation of 
%Fourier coefficients by taking M organs.
%
%INPUT:
%  ak = Fourier coefficient vector.
%  m = Number of organs
%  n = vector of time.
%  N = period time.
%OUTPUT:
%  am = A signal approximation to the original signal
    k = -m:m;
    matrix = exp(1i*(2*pi/N)*k'*n) ;
    new_ak = ak(((N-1)/2)-m:((N-1)/2)+m);
    am = new_ak*matrix;


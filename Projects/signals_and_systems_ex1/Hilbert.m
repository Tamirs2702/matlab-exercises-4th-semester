function [hk] = Hilbert (gk ,N)
%Hilbert Returns a vector that is a multiple of gk with -j when k is positive
%and a multiple with j when k is negative
%
%INPUT:
%  gk= The original vector
%  N = period time.
%
%OUTPUT:
%  hk = A new vector used as a Hilbert filter
    h1 = 1i.*ones(1,(N-1)/2);
    h2 = -1i.*ones(1,(N-1)/2);
    h3 = [h1 , 1i , h2];
    hk = (h3.*gk);
% calculates the likelihood of a given non T-to-C (A-to-C or G-to-C) data
%
% Written by Pejman Mohammadi and Monica Golumbeanu
% pejman.mohammadi@bsse.ethz.ch
% monica.golumbeanu@bsse.ethz.ch
%
function [nLL, GnLL] = sGet_likeN(X,Y, W, epsilon, nu)
nLL =  Calc_nLL(X,Y, W, epsilon, nu);
EPS=1E-5;
GnLL = [
    Calc_nLL(X,Y, W, epsilon+EPS, nu)
    nLL
    Calc_nLL(X,Y, W, epsilon, nu+EPS)
    nLL
    ];
GnLL = (GnLL-nLL)/EPS;
end

function nLL = Calc_nLL(X,Y, W, epsilon, nu)
epsl = log(realmin('double'));
% Calculate the likelihood
XC = Y-X;
L = ...
       nu  *            c_binopdf(X, XC, 1-3*epsilon) + ...
    (1-nu) *            c_binopdf(X, XC,   epsilon) ;
nLL = -sum( W .* max(epsl, log(L)));
end

% scaled binomial
function p = c_binopdf(x, xc, p)
p = p.^x .* (1-p).^xc;
end
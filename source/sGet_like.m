% calculates the likelihood for a given T-to-C data
% 
% written by Pejman Mohammadi and Monica Golumbeanu
% pejman.mohammadi@bsse.ethz.ch
% monica.golumbeanu@bsse.ethz.ch
%
function [nLL, GnLL] = sGet_like(X,Y, W, epsilon, gamma, nu, psi)
nLL =  Calc_nLL(X,Y, W, epsilon, gamma, nu, psi);
EPS=1E-5;
GnLL = [
    Calc_nLL(X,Y, W, epsilon+EPS, gamma, nu, psi)
    Calc_nLL(X,Y, W, epsilon, gamma+EPS, nu, psi)
    Calc_nLL(X,Y, W, epsilon, gamma, nu+EPS, psi)
    Calc_nLL(X,Y, W, epsilon, gamma, nu, psi+EPS)
    ];
GnLL = (GnLL-nLL)/EPS;
end

function nLL = Calc_nLL(X,Y, W, epsilon, gamma, nu, psi)
%% Calculate the likelihood
epsl = log(realmin('double'));
theta = (1-gamma)*epsilon + (1-3*epsilon)*gamma;
XC = Y-X;
L = ...
    nu  *            c_binopdf(X, XC, 1-3*epsilon) + ...
    (1-nu) * (1-psi) * c_binopdf(X, XC, epsilon) + ...
    (1-nu) *    psi  * c_binopdf(X, XC, theta);
nLL = -sum( W .* max(epsl, log(L)));
end

% scaled binomial
function p = c_binopdf(x, xc, p)
p = p.^x .* (1-p).^xc;
end
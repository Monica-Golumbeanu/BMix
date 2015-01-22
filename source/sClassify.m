% Function sClassify computes the posterior probability for the three
% classes and chooses the class with the maximum posterior.
%
% Written by Monica Golumbeanu and Pejman Mohammadi
% pejman.mohammadi@bsse.ethz.ch
% monica.golumbeanu@bsse.ethz.ch
%
function [MLE_ClassI2, MLE_pC2] = sClassify(X,Y, epsilon, gamma, nu, psi)
theta = (1-gamma)*epsilon + (1-3*epsilon)*gamma;
Xc = Y-X; clear Y

lpC2(:,3) = log(1-nu) +    log(psi) + c_Lbino(X, Xc, theta);
lpC2(:,1) = log(nu)   +                c_Lbino(X, Xc, 1-3*epsilon);
lpC2(:,2) = log(1-nu) + log(1-psi)  + c_Lbino(X, Xc, epsilon);

lpC2 = lpC2 - repmat(max(lpC2,[],2),1,3);
pC2 = exp(lpC2); clear lpC2
pC2 = pC2 ./ repmat(sum(pC2,2),1, 3); % Normalize the evidences
[MLE_pC2, MLE_ClassI2] = max(pC2, [], 2);
end

function LL = c_Lbino(x, xc, p)
LL = log(p).*x + log(1-p).*xc;
end
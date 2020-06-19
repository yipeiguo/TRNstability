% Here we write another function to obtain the Jacobian matrix using
% simplified, more intuitive expressions for the elements.

function [JacobianMat,meandlogfdcj,dlogfijdcj_mat] = CreateJacMat_method2(c_ss,...
    IntParamsMat,phivec,mu)

N = length(phivec);

gamma_ij = IntParamsMat(:,1);
K_ij = IntParamsMat(:,2);
n_ij = IntParamsMat(:,3);
%cj = repelem(c_ss,N);
cj = repmat(c_ss,N,1);
dlogfijdcj = (n_ij.*gamma_ij.*cj.^(n_ij-1).*K_ij.^n_ij)./...
    ((K_ij.^n_ij+cj.^n_ij).*(K_ij.^n_ij+(1+gamma_ij).*cj.^n_ij));

dlogfijdcj_mat = reshape(dlogfijdcj,[N,N])';
meandlogfdcj = sum(repmat(c_ss,1,N).*dlogfijdcj_mat,1);

JacobianMat = -(mu*c_ss(end)).*(eye(N) + ...
    repmat(c_ss,1,N).*(repmat(meandlogfdcj,N,1)-dlogfijdcj_mat));

end

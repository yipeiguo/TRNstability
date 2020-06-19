% This is a function that defines the system of nonlinear equations for
% dc/dt. This set of equations will be solved to find steady-state
% solutions. Here, the x's are concentrations

function [F,fivec] = SetofEqns_v2(x,IntParamsMat,phivec,fijintfunc)

N = length(phivec);
% Define f_ij for every pair of (directed) interactions, where f_ij
% represents out j affects i.
fijvec = fijintfunc(IntParamsMat(:,1),IntParamsMat(:,2),IntParamsMat(:,3),...
    repmat(x,N,1));

% Calculat f_i = prod_j f_ij, i.e. how phi_i is being regulated by all
% other genes.
fivec = prod(reshape(fijvec,[N,N]),1)';

% Partition function
Z = sum(phivec.*fivec);

% Set of dynamical equations (RHS)
%F = mu.*x(end).*(phivec.*fivec./Z - x);
F = phivec.*fivec./Z - x;


end


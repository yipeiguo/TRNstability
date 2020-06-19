% This is a function that defines the system of nonlinear equations for
% dc/dt. This set of equations will be solved using odesolvers to calculate 
% the trajectory to steady-state solutions. Here, the y's are concentrations

function dydt = SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,mu)

N = length(phivec);
% Define f_ij for every pair of possible (directed) interactions, where f_ij
% represents out j affects i.
fijvec = fijintfunc(IntParamsMat(:,1),IntParamsMat(:,2),IntParamsMat(:,3),...
    repmat(y,N,1));

% Calculat f_i = prod_j f_ij, i.e. how phi_i is being regulated by all
% other genes.
fivec = prod(reshape(fijvec,[N,N]),1)';

% Partition function
Z = sum(phivec.*fivec);

% Set of dynamical equations (RHS)
cr = y(N,1);
dydt = (mu*cr).*(phivec.*fivec./Z - y);
%F = phivec.*fivec./Z - y;


end


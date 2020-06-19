% This is an event function that allows ode solver to stop when the event
% is reached. Here the event is when steady state is reached

% function [x,isterm,dir] = eventfun(t,y,IntPairs,IntParamsMat,phivec,fijintfunc,eps)
function [x,isterm,dir] = eventfun(t,y,IntPairs,IntParamsMat,phivec,fijintfunc,eps,tstart,tmax)

N = length(phivec);
% Define f_ij for every pair of (directed) interactions, where f_ij
% represents out j affects i.
ymat = sparse(IntPairs(:,1),IntPairs(:,2),y(IntPairs(:,1)),N,N);
fijvec = fijintfunc(IntParamsMat(:,1),IntParamsMat(:,2),IntParamsMat(:,3),...
    ymat(:));

% Calculat f_i = prod_j f_ij, i.e. how phi_i is being regulated by all
% other genes.
fivec = prod(reshape(fijvec,[N,N]),1)';

% Partition function
Z = sum(phivec.*fivec);

% At steady state:
F = phivec.*fivec./Z - y;

% Set of dynamical equations (RHS)
%cr = y(N,1);
%dydt = (mu*cr).*(phivec.*fivec./Z - y);

if (~exist('tstart', 'var')) || (~exist('tmax', 'var'))
    x = max(abs(F))-eps;
    %x = max(abs(dydt))>1e-6;
    %x = max(abs(dydt))-1e-7;
    %x = norm(dydt)-1e-8;
else
    x(1) = max(abs(F))-eps;
    x(2) = toc(tstart)<tmax;

end
isterm = true(size(x));
dir = zeros(size(x));  %or 0, doesn't matter

% isterm = 1;
% dir = 0;  %or 0, doesn't matter

end
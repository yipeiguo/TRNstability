% Here we look at how the maximum eigenvalue in our system scale with N if
% we assume there are no loops (except possibly self-loops)

close all; clear varriables;

%% Define parameters
Nscan = 100:100:1000;
%Nscan = [100,200:200:1000]; % number of genes (including ribosomes)
%Nscan = [1500,2000]; % number of genes (including ribosomes)
numtrials = 10; % number of trials (random networks) for each N and gamma
density = 0.01;
mu = 1;
phidist = 1; % 1= uniform, 2=uniform in log-space
%phir = 0.1;
maxfcscan = 1.5;
n = 1;
pdown = 0.5; % probability of being upregulating (pup = 1-pdown)

% fold change distribution
fcdist = 'uniform'; % either 'uniform','1overx' or '1overx2'

% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

eps = 1e-5;
tspan = [0 1e11];
%% Specify filename (where data will be stored)
runindx = 1;
namestr = sprintf('EigScaling_randomDAG_phidist%d_rho%.2f_n%d_Nmin%d_Nmax%d_fcmin%.1f_fcmax%.1f_pdown%.1f_numtrials%d_run%d',...
    [phidist,density,n,min(Nscan),max(Nscan),min(maxfcscan),max(maxfcscan),pdown,numtrials,runindx]);    
namestr = strrep(namestr,'.','pt');

%% Storage matrices
% for investigating how stability varies:
meanlambdaMat_Jtilde = zeros(length(Nscan),length(maxfcscan),numtrials);
minlambdaMat_Jtilde = zeros(length(Nscan),length(maxfcscan),numtrials);
Jtilde_varMat = zeros(length(Nscan),length(maxfcscan),numtrials);
Jtilde_meanMat = zeros(length(Nscan),length(maxfcscan),numtrials);

% For checking convergence:
maxfvalmat = zeros(length(Nscan),length(maxfcscan),numtrials);
exitflagmat = zeros(length(Nscan),length(maxfcscan),numtrials);

% Properties of solutions and jacobian
var_css = zeros(length(Nscan),length(maxfcscan),numtrials);
var_meandlogfdcj = zeros(length(Nscan),length(maxfcscan),numtrials);
var_cmeandlogfdcj = zeros(length(Nscan),length(maxfcscan),numtrials);

%% Main simulation
tic
for Nindx = 1:length(Nscan)    
    N = Nscan(Nindx);
    fprintf('N = %d \n',N);
    maxiter = N*N;
    for trialindx = 1:numtrials
        fprintf('trialindx = %d \n',trialindx);
        % If connections are drawn randomly:
        g = CreateRandomDAG(N,density,maxiter);
        disp('DAG created!');
        IntTypeMat = rand(N).*g; 
        IntTypeMat(IntTypeMat~=0) = IntTypeMat(IntTypeMat~=0) - pdown;
        IntTypeMat = sign(IntTypeMat);
        IntTypeVec = IntTypeMat(:);
        % IntType = 1 for upregulation; -1  for downregulation

        if phidist == 1
            phivec = rand(N,1);
        elseif phidist == 2
            logphi = rand(N,1).*10;
            phivec = exp(logphi);
        else % specify custome values
            %phivec = [1;2;3];
            phivec = ones(N,1)./(N);
        end
        %phivec = [sort(phivec,'descend')./sum(phivec).*(1-phir);phir];
        phivec = phivec./sum(phivec);
        
        % Interaction parameters 
        IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                                     % 2nd col: Kd, 3rd col: Hill's coefficient n
        IntParamsMat(:,2) = repmat(phivec,N,1);
        IntParamsMat(:,3) = 1;
        
        c_guess = phivec;
        for fcindx = 1:length(maxfcscan)
            maxfc = maxfcscan(fcindx);
            fprintf('maxfc = %d \n',maxfc);
            
            if strcmp(fcdist,'uniform')
                foldchangeVec = (rand(N*N,1).*(maxfc-1)).*(IntTypeVec~=0)+1;
            elseif strcmp(fcdist,'1overx')
                logfcVec = (rand(N*N,1).*log(maxfc)).*(IntTypeVec~=0);
                foldchangeVec = exp(logfcVec);
            elseif strcmp(fcdist,'1overx2')
                inversefcVec = 1-(rand(N*N,1).*(1-1/maxfc)).*(IntTypeVec~=0);
                foldchangeVec = 1./inversefcVec;
            end
            gammaVec = zeros(N*N,1);
            gammaVec(IntTypeVec==1) = foldchangeVec(IntTypeVec==1)-1;
            gammaVec(IntTypeVec==-1) = 1./foldchangeVec(IntTypeVec==-1)-1;
            
            IntParamsMat(:,1) = gammaVec;
            
            % Calculate steady state concentrations
            [c_ss,fval,exitflag,~] = ...
                fsolve(@(x) SetofEqns_v2(x,IntParamsMat,phivec,fijintfunc),c_guess);

            exitflagmat(Nindx,fcindx,trialindx) = exitflag;
            maxfvalmat(Nindx,fcindx,trialindx) = max(abs(fval));
        
            % Calculate stability of steady state:
            [JacobianMat,meandlogfdcj,dlogfijdcj_mat] = CreateJacMat_method2(c_ss,...
                IntParamsMat,phivec,mu);
            lambdavec = eig(JacobianMat);
            JtildeMat = JacobianMat./(-mu*c_ss(end))-eye(N);
            lambdavec_Jtilde = lambdavec./(-mu.*c_ss(end))-1;
            meanlambdaMat_Jtilde(Nindx,fcindx,trialindx) = mean(lambdavec_Jtilde);
            minlambdaMat_Jtilde(Nindx,fcindx,trialindx) = min(real(lambdavec_Jtilde));
            Jtilde_varMat(Nindx,fcindx,trialindx) = var(JtildeMat(:));
            Jtilde_meanMat(Nindx,fcindx,trialindx) = mean(JtildeMat(:));
            
            if exitflag>0
                c_guess = c_ss;
            end
        end
    end
    save(namestr);
end
toc


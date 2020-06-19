% In this script, we investigate how the mean and maximum eigenvalue of the
% Jacobian vary with the number of genes N for a fixed sparsity of
% connections (e.g. fully connected network etc.)

close all; 
clear variables;
rng('shuffle');

%% Define parameters
Nscan = 100:100:1000; % number of genes (including ribosomes)
numtrials = 10; % number of trials (random networks) for each N
%Nscan = ones(1,6).*200;
sparsity = 0.01;
mu = 1;
phidist = 1; % 1= uniform, 2=uniform in log-space
%phir = 0.1;
maxfcscan = 1.5;
n = 1;
%gammascan = gammamax;
pdown = 0.5; % probability of being upregulating (pup = 1-pdown)

% fold change distribution
fcdist = '1overx'; % either 'uniform','1overx' or '1overx2'

% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

% specify method of finding css
iffsolve = true;
ifodesolver = true;

% Optimization parameters
if iffsolve == true
    options_fsolve = optimoptions('fsolve','MaxFunctionEvaluations',1e5,'MaxIterations',1000);
end
if ifodesolver == true
    eps = 1e-5;
    tspan = [0 1e11];
end


%% Storage matrices
%phimat = zeros(max(Nscan),length(Nscan));
%cssmat = zeros(max(Nscan),length(Nscan)); 
% for investigating how stability varies with N:
meanlambdaMat_Jtilde = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
minlambdaMat_Jtilde = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
Jtilde_varMat = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);

% For checking convergence:
maxfvalmat = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
exitflagmat = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
if ifodesolver == true
    temat = zeros(length(Nscan),length(maxfcscan),numtrials);
end

% Properties of solutions and jacobian
var_css = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
var_meandlogfdcj = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);
var_cmeandlogfdcj = zeros(length(Nscan),length(maxfcscan),numtrials,iffsolve+ifodesolver);

%% Specify name of file
runindx = 1;
namestr = sprintf('EigScaling_fullrandom_phidist%d_rho%.2f_n%d_Nmin%d_Nmax%d_fcmin%.2f_fcmax%.2f_pdown%.2f_numtrials%d_fsolve%d_ode%d_run%d',...
    [phidist,sparsity,n,min(Nscan),max(Nscan),min(maxfcscan),max(maxfcscan),pdown,numtrials,iffsolve,ifodesolver,runindx]);    
namestr = strrep(namestr,'.','pt');

%% Main Calculation
tic
for Nindx = 1:length(Nscan)
    N = Nscan(Nindx);
    fprintf('N = %d \n',N);
    
    for trialindx = 1:numtrials
        fprintf('trialindx = %d \n',trialindx);
    
        %IntTypeMat = randn(N,N);
        IntTypeMat = sprand(N,N,sparsity); 
        %IntTypeMat(1:N+1:end) = 0; % no self-regulation
        IntTypeMat(end,:) = 0; % ribosomes don't regulate other genes
        %IntTypeMat = IntTypeMat~=0;
        IntTypeMat(IntTypeMat~=0) = IntTypeMat(IntTypeMat~=0) - pdown;
        IntTypeMat = sign(IntTypeMat);
        [IntPairs_i,IntPairs_j] = find(IntTypeMat);
        IntPairs = [IntPairs_i,IntPairs_j];
        IntTypeVec = full(IntTypeMat(:));
        % IntType = 1 for upregulation; -1  for downregulation
        
        if phidist == 1
            %phivec = rand(N-1,1);
            phivec = rand(N,1);
        elseif phidist == 2
            logphi = rand(N,1).*8;
            %logphi = rand(N-1,1).*10;
            phivec = exp(logphi);
        else % specify custome values
            %phivec = [1;2;3];
            %phivec = ones(N-1,1)./(N-1);
            phivec = ones(N,1)./(N);
        end
        %phivec = [sort(phivec,'descend')./sum(phivec).*(1-phir);phir];
        phivec = phivec./sum(phivec);
        
        % Interaction parameters 
        IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                                     % 2nd col: Kd, 3rd col: Hill's coefficient n
        IntParamsMat(:,2) = repmat(phivec,N,1);
        IntParamsMat(:,3) = n;

        cguess = phivec;
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
            if iffsolve == true
                [c_ss,fval,exitflag,~] = ...
                   fsolve(@(x) SetofEqns_v2(x,IntParamsMat,phivec,fijintfunc),cguess,options_fsolve);
                var_css(Nindx,fcindx,trialindx,1) = var(c_ss);
                exitflagmat(Nindx,fcindx,trialindx,1) = exitflag;
                maxfvalmat(Nindx,fcindx,trialindx,1) = max(abs(fval));
                
                % Calculate stability of steady state:
                [JacobianMat,meandlogfdcj,~] = CreateJacMat_method2(c_ss,...
                    IntParamsMat,phivec,mu);
                lambdavec = eig(JacobianMat);
                JtildeMat = JacobianMat./(-mu*c_ss(end))-eye(N);
                lambdavec_Jtilde = lambdavec./(-mu.*c_ss(end))-1;
                
                meanlambdaMat_Jtilde(Nindx,fcindx,trialindx,1) = mean(lambdavec_Jtilde);
                minlambdaMat_Jtilde(Nindx,fcindx,trialindx,1) = min(real(lambdavec_Jtilde));
                Jtilde_varMat(Nindx,fcindx,trialindx,1) = var(JtildeMat(:));
                var_meandlogfdcj(Nindx,fcindx,trialindx,1) = var(meandlogfdcj);
                cmeandlogfdcj = repmat(c_ss,1,N).*repmat(meandlogfdcj,N,1);
                var_cmeandlogfdcj(Nindx,fcindx,trialindx,1) = var(cmeandlogfdcj(:));

                if exitflag>0
                    c_guess = c_ss;
                end
            end
            
            if ifodesolver == true
                cguess2 = phivec; % rand(N,1);
                cguess2 = cguess2./sum(cguess2);
                opts = odeset('Events',@(t,y) eventfun(t,y,IntPairs,IntParamsMat,phivec,fijintfunc,eps));
                [ttraj_kk,ctraj_kk,te,c_ss,ie] = ode45(@(t,y) SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,mu), ...
                    tspan,cguess2,opts);
                if ~isempty(ie)
                    exitflag = 1;            
                else 
                    exitflag = 0;
                    c_ss = ctraj_kk(end,:);
                    warning('dynamics has not converged!');  
                    te = ttraj_kk(end);
                end
                c_ss = c_ss';
                fval = SetofEqns_v2(c_ss,IntParamsMat,phivec,fijintfunc);
                temat(Nindx,fcindx,trialindx) = te;  
                var_css(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = var(c_ss);
                exitflagmat(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = exitflag;
                maxfvalmat(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = max(abs(fval));
                
                % Calculate stability of steady state:
                [JacobianMat,meandlogfdcj,~] = CreateJacMat_method2(c_ss,...
                    IntParamsMat,phivec,mu);
                lambdavec = eig(JacobianMat);
                JtildeMat = JacobianMat./(-mu*c_ss(end))-eye(N);
                lambdavec_Jtilde = lambdavec./(-mu.*c_ss(end))-1;
                
                meanlambdaMat_Jtilde(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = mean(lambdavec_Jtilde);
                minlambdaMat_Jtilde(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = min(real(lambdavec_Jtilde));
                Jtilde_varMat(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = var(JtildeMat(:));
                var_meandlogfdcj(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = var(meandlogfdcj);
                cmeandlogfdcj = repmat(c_ss,1,N).*repmat(meandlogfdcj,N,1);
                var_cmeandlogfdcj(Nindx,fcindx,trialindx,iffsolve+ifodesolver) = var(cmeandlogfdcj(:));
            end            
        end                
    end
    save(namestr);
end
toc


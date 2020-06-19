% This is a script that looks that how max eigenvalue of Jacobian matrix 
% scales with fraction of negative interactions for bipartitelike networks using the 
% function 'CreateNetwork_bipartitelike.m'.

close all; clear variables;
rng('shuffle');

%% Define parameters
% Network parameters
Nscan = 200; % number of genes (including ribosomes)
k = 0.1;
density = 0.01;
iffixNumSelfInts = true;
if iffixNumSelfInts == true
    numSelfInts = 0;
    NetworkProp.numSelfInts = numSelfInts;
end
numInts_TFTF = 0;
NetworkProp.numInts_TFTF = numInts_TFTF;
Pnegscan = 0:0.1:1;
Pneg_selfint = 0.5;
NetworkProp.Pneg_selfint = Pneg_selfint;
% Pneg_nonselfint = 0.5;
% NetworkProp.Pneg_nonselfint = Pneg_nonselfint;


% interaction/biological parameters
k1 = 1;
phidist = 1;
n = 1;
maxfc = 1.5;

% fold change distribution
fcdist = 'uniform'; % either 'uniform','1overx' or '1overx2'

                             
% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

% Technical parameters
numtrials = 100; %5;
numsets = 5;

% method of solving for steady-state
method = 'ode45';
if strcmp(method, 'ode45')
    tspan = [0 1e12];
    eps = 1e-5;
    tmax = 600; % time limit in seconds for each attempt to find css using odesolver
end

%% Specify filename (where data will be stored)
runindx = 1;
if iffixNumSelfInts == false
    namestr = sprintf('BipartitePnegscaling_Nmin%d_Nmax%d_k%.2f_rho%.2f_n%d_minPneg%.1f_maxPneg%.1f_maxfc%.1f_numtrials%d_numsets%d_run%d',...
        [min(Nscan),max(Nscan),k,density,n,min(Pnegscan),max(Pnegscan),maxfc,numtrials,numsets,runindx]);    
else
    namestr = sprintf('BipartitePnegData_Nmin%d_Nmax%d_k%.2f_rho%.2f_numSelfInts%d_n%d_minPneg%.1f_maxPneg%.1f_maxfc%.1f_numtrials%d_numsets%d_run%d',...
        [min(Nscan),max(Nscan),k,density,numSelfInts,n,min(Pnegscan),max(Pnegscan),maxfc,numtrials,numsets,runindx]);    
end
namestr = strrep(namestr,'.','pt');

%% Define storage arrays
% interaction parameters
gamma_cell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);

% css
if strcmp(method, 'ode45')
    ttraj_cell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);
    ctraj_cell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);
end
cssCell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);
phiCell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);
fiCell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);

% Jacobian matrix elements
Tmat_cell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);

% eigenvalues
maxlambdaMat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets);
lambdaJtilde_cell = cell(length(Nscan),length(Pnegscan),numtrials,numsets);

% For checking convergence:
maxfvalmat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets);
exitflagmat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets);
if strcmp(method, 'ode45')
    temat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets);
end


%% Main simulation
tic
for setindx = 1:numsets
    fprintf('setindx: %d \n',setindx);
    
    for kk = 1:length(Nscan)
        if iffixNumSelfInts == false
            numSelfInts = round(numTFTFints/(q-1));
            NetworkProp.numSelfInts = numSelfInts;
        end
        
        N = Nscan(kk);
        q = k*N;
        numInts = density*N^2;
        NetworkProp.N = N;
        NetworkProp.numTFs = q;
        NetworkProp.numInts = numInts;
        
        fprintf('N: %d \n',N);
        
        IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                             % 2nd col: Kd, 3rd col: Hill's coefficient n
        IntParamsMat(:,3) = n;    
        
        for Pnegindx = 1:length(Pnegscan)
            
            Pneg = Pnegscan(Pnegindx);
            
            fprintf('pneg: %d \n',Pneg);
            NetworkProp.Pneg_nonselfint = Pneg;
    
            for trialindx = 1:numtrials
                fprintf('trialindx: %d \n',trialindx);

                % Draw new network topology 
                IntTypeVec = CreateNetwork_bipartitelike(NetworkProp);

                % Draw interaction parameters
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

                % store interaction parameters
                [intInds,~,gammaVals] = find(gammaVec);
                gamma_cell{kk,Pnegindx,trialindx,setindx} = [intInds,gammaVals];

                % original effective gene copy number
                if phidist == 1
                    %phivec = rand(N-1,1);
                    phivec = rand(N,1);
                elseif phidist == 2
                    logphi = rand(N,1).*10;
                    %logphi = rand(N-1,1).*10;
                    phivec = exp(logphi);
                else % specify custome values
                    %phivec = [1;2;3];
                    %phivec = ones(N-1,1)./(N-1);
                    phivec = ones(N,1)./(N);
                end
                %phivec = [sort(phivec,'descend')./sum(phivec).*(1-phir);phir];
                phivec = phivec./sum(phivec);
                phiCell{kk,trialindx,setindx} = phivec;
                c_guess = phivec; 

                % Interaction parameters 
                IntParamsMat(:,1:2) = [gammaVec,repmat(phivec,N,1)];
                %IntParamsMat(:,1) = gammaVec;
                %IntParamsMat(:,2) = repmat(phivec,N,1);
                %IntParamsMat(:,2) = ones(N*N,1)./N;

                % Find css
                if strcmp(method, 'fsolve')
                    [c_ss,fval,exitflag,~] = ...
                       fsolve(@(x) SetofEqns_v2(x,IntParamsMat,phivec,fijintfunc),c_guess);
                    [~,fivec] = SetofEqns_v2(c_ss,IntParamsMat,phivec,fijintfunc);
                elseif strcmp(method, 'ode45')
                    tstart = tic;
                    %opts = odeset('Events',@(t,y) eventfun(t,y,IntPairs,IntParamsMat,phivec,fijintfunc,eps));
                    opts = odeset('Events',@(t,y) eventfun(t,y,IntParamsMat,phivec,fijintfunc,eps,tstart,tmax));
                    [ttraj_kk,ctraj_kk,te,c_ss,ie] = ode45(@(t,y) SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,k1), ...
                        tspan,c_guess,opts);
                    if isempty(ie)
                        exitflag = 0;
                        c_ss = ctraj_kk(end,:);
                        warning('dynamics has not converged!');  
                        te = ttraj_kk(end);                    
                    end
                    c_ss = c_ss';
                    [fval,fivec] = SetofEqns_v2(c_ss,IntParamsMat,phivec,fijintfunc);
                    if ~isempty(ie)
                        if max(abs(fval)) <= eps
                            exitflag = 1;
                        else
                            exitflag = 0;
                            warning('max time exceeded but dynamics has not converged!');          
                        end
                    end
                    ttraj_cell{kk,Pnegindx,trialindx,setindx} = ttraj_kk;
                    ctraj_cell{kk,Pnegindx,trialindx,setindx} = ctraj_kk;
                    temat(kk,Pnegindx,trialindx,setindx) = te;            
                end
                cssCell{kk,Pnegindx,trialindx,setindx} = c_ss;
                fiCell{kk,Pnegindx,trialindx,setindx} = fivec;
                exitflagmat(kk,Pnegindx,trialindx,setindx) = exitflag;
                maxfvalmat(kk,Pnegindx,trialindx,setindx) = max(abs(fval));

                % Calculate stability of steady state (real network):
                [JacobianMat,~,~] = CreateJacMat_method2(c_ss,IntParamsMat,phivec,k1);
                Tmat = JacobianMat(1:q,1:q);
                Tmat_cell{kk,Pnegindx,trialindx,setindx} = Tmat;
                lambdavec = eig(JacobianMat);
                lambdavec_Jtilde = lambdavec./(k1.*c_ss(end))+1;
                maxlambdaMat(kk,Pnegindx,trialindx,setindx) = max(real(lambdavec_Jtilde));
                lambdaJtilde_cell{kk,Pnegindx,trialindx,setindx} = lambdavec_Jtilde;
                toc
            end
        end
    end
    % Save workspace
    save(namestr);     
    
    if strcmp(method, 'ode45')
        % save ctraj_cell separately since it is sometimes larger than 2gb
        save(strcat(namestr,'_ctraj'),'ctraj_cell','-v7.3');
    end
end                   
toc





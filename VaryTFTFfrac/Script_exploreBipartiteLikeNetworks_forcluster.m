% This is a script to explore how the density of the TF-otherTF
% interactions affect the stability of the system

close all; clear variables;
rng('shuffle');

%% Define parameters
% Network parameters
N = 2275; % total number of different genes/proteins
q = 211;  % number of transcription factors
iffixNumSelfInts = true;
if iffixNumSelfInts == true
    numSelfInts = 0;
    NetworkProp.numSelfInts = numSelfInts;
end
Pneg_selfint = 0.5;
NetworkProp.Pneg_selfint = Pneg_selfint;
numTFTFints_scan = 0:20:400;
conservedqty = 'numInts_TFnonTF'; % can be either 'numInts' or 'numInts_TFnonTF'

NetworkProp.N = N;
NetworkProp.numTFs = q;
if strcmp(conservedqty,'numInts')
    numInts = 5655;
    NetworkProp.numInts = numInts;
elseif strcmp(conservedqty,'numInts_TFnonTF')
    numInts_TFnonTF = 5000;
end

% interaction/biological parameters
k1 = 1;
phidist = 1;
n = 2;
maxfc = 1000;
                             
% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

% Technical parameters
numtrials = 10; %5;
numsets = 1;

% method of solving for steady-state
tspan = [0 1e12];
eps = 1e-4;
tmax = 600; % time limit in seconds for each attempt to find css using odesolver

%% Specify filename (where data will be stored)
runindx = 1;
if iffixNumSelfInts == false
    if strcmp(conservedqty,'numInts')
        namestr = sprintf('BipartiteLikeData_N%d_q%d_numInts%d_minnumTFTF%d_maxnumTFTF%d_n%d_maxfc%d_numtrials%d_run%d',...
            [N,q,numInts,min(numTFTFints_scan),max(numTFTFints_scan),n,maxfc,numtrials,runindx]);
    elseif strcmp(conservedqty,'numInts_TFnonTF')
        namestr = sprintf('BipartiteLikeData_N%d_q%d_numIntsTFnonTF%d_minnumTFTF%d_maxnumTFTF%d_n%d_maxfc%d_numtrials%d_run%d',...
            [N,q,numInts_TFnonTF,min(numTFTFints_scan),max(numTFTFints_scan),n,maxfc,numtrials,runindx]);
    end
else
    if strcmp(conservedqty,'numInts')
        namestr = sprintf('BipartiteLikeData_N%d_q%d_numInts%d_numSelfInts%d_minnumTFTF%d_maxnumTFTF%d_n%d_maxfc%d_numtrials%d_run%d',...
            [N,q,numInts,numSelfInts,min(numTFTFints_scan),max(numTFTFints_scan),n,maxfc,numtrials,runindx]);
    elseif strcmp(conservedqty,'numInts_TFnonTF')
        namestr = sprintf('BipartiteLikeData_N%d_q%d_numSelfInts%d_numIntsTFnonTF%d_minnumTFTF%d_maxnumTFTF%d_n%d_maxfc%d_numtrials%d_run%d',...
            [N,q,numSelfInts,numInts_TFnonTF,min(numTFTFints_scan),max(numTFTFints_scan),n,maxfc,numtrials,runindx]);
    end
end
%% Define storage arrays
% interaction parameters
gamma_cell = cell(length(numTFTFints_scan),numtrials,numsets);

% css
ttraj_cell = cell(length(numTFTFints_scan),numtrials,numsets);
ctraj_cell = cell(length(numTFTFints_scan),numtrials,numsets);
cssMat = zeros(N,length(numTFTFints_scan),numtrials,numsets);
phiMat = zeros(N,length(numTFTFints_scan),numtrials,numsets);
fiMat = zeros(N,length(numTFTFints_scan),numtrials,numsets);

% Jacobian matrix elements
Tmat_cell = cell(length(numTFTFints_scan),numtrials,numsets);

% eigenvalues
maxlambdaMat = zeros(length(numTFTFints_scan),numtrials,numsets);
lambdaJtilde_cell = cell(length(numTFTFints_scan),numtrials,numsets);

% For checking convergence:
maxfvalmat = zeros(length(numTFTFints_scan),numtrials,numsets);
exitflagmat = zeros(length(numTFTFints_scan),numtrials,numsets);
temat = zeros(length(numTFTFints_scan),numtrials,numsets);


%% Main simulation
tic
for setindx = 1:numsets
    fprintf('setindx: %d \n',setindx);
    
    IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                             % 2nd col: Kd, 3rd col: Hill's coefficient n
    IntParamsMat(:,3) = n;
    
    for kk = 1:length(numTFTFints_scan)
        numTFTFints = numTFTFints_scan(kk);
        fprintf('number of TF-TF interactions: %d \n',numTFTFints);    
        NetworkProp.numInts_TFTF = numTFTFints;
        
        if iffixNumSelfInts == false
            numSelfInts = round(numTFTFints/(q-1));
            NetworkProp.numSelfInts = numSelfInts;
        end
        
        if strcmp(conservedqty,'numInts_TFnonTF')
            numInts = numSelfInts + numTFTFints + numInts_TFnonTF;
            NetworkProp.numInts = numInts;
        end
        
        for trialindx = 1:numtrials
            fprintf('trialindx: %d \n',trialindx);
            
            % Draw new network topology 
            IntTypeVec = CreateNetwork_bipartitelike(NetworkProp);
            
            % Draw interaction parameters
            foldchangeVec = (rand(N*N,1).*(maxfc-1)).*(IntTypeVec~=0)+1;
            gammaVec = zeros(N*N,1);
            gammaVec(IntTypeVec==1) = foldchangeVec(IntTypeVec==1)-1;
            gammaVec(IntTypeVec==-1) = 1./foldchangeVec(IntTypeVec==-1)-1;
            
            % store interaction parameters
            [intInds,~,gammaVals] = find(gammaVec);
            gamma_cell{kk,trialindx,setindx} = [intInds,gammaVals];
            
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
            phiMat(:,kk,trialindx,setindx) = phivec;
            c_guess = phivec; 

            % Interaction parameters 
            IntParamsMat(:,1:2) = [gammaVec,repmat(phivec,N,1)];
            %IntParamsMat(:,1) = gammaVec;
            %IntParamsMat(:,2) = repmat(phivec,N,1);
            %IntParamsMat(:,2) = ones(N*N,1)./N;
            
            % Find css
%             [c_ss2,fval,exitflag,~] = ...
%                fsolve(@(x) SetofEqns_v2(x,IntParamsMat,phivec,fijintfunc),c_guess);
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
            cssMat(:,kk,trialindx,setindx) = c_ss;
            fiMat(:,kk,trialindx,setindx) = fivec;
            ttraj_cell{kk,trialindx,setindx} = ttraj_kk;
            ctraj_cell{kk,trialindx,setindx} = ctraj_kk;
            exitflagmat(kk,trialindx,setindx) = exitflag;
            maxfvalmat(kk,trialindx,setindx) = max(abs(fval));
            temat(kk,trialindx,setindx) = te;
            
            % Calculate stability of steady state (real network):
            [JacobianMat,~,~] = CreateJacMat_method2(c_ss,IntParamsMat,phivec,k1);
            Tmat = JacobianMat(1:q,1:q);
            Tmat_cell{kk,trialindx,setindx} = Tmat;
            lambdavec = eig(JacobianMat);
            lambdavec_Jtilde = lambdavec./(k1.*c_ss(end))+1;
            maxlambdaMat(kk,trialindx,setindx) = max(real(lambdavec_Jtilde));
            lambdaJtilde_cell{kk,trialindx,setindx} = lambdavec_Jtilde;
            toc
        end
    end
    % Save workspace
    save(namestr);     
    
    % save ctraj_cell separately since it is sometimes larger than 2gb
    save(strcat(namestr,'_ctraj'),'ctraj_cell','-v7.3');
end                   
toc




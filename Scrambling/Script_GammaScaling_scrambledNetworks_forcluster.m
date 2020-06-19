% This script investigates the scaling of stability with interaction
% strength for scrambled versions of the real network.

close all; clear variables; clc
tic

%% Extract features of real networks

% Import data file interactively
load('PNAS_EcoliTRN_Data');

TFsCell = dataCell(:,2);
genes_regCell = dataCell(:,4);
IntTypeCell = dataCell(:,5);
numints = length(IntTypeCell);

IntTypeCell = strrep(IntTypeCell,'+','1');
IntTypeCell = strrep(IntTypeCell,'-','-1');

TFs = cell2mat(cellfun(@(x) str2double(x(2:end)),TFsCell,'UniformOutput',false));
genes_reg = cell2mat(cellfun(@(x) str2double(x(2:end)),genes_regCell,'UniformOutput',false));
IntType = cell2mat(cellfun(@str2double,IntTypeCell,'UniformOutput',false));

% Find rows with multiple TFs
rows_mulTFs = find(isnan(TFs));
rows_mulreg = find(isnan(genes_reg));
commonrows = intersect(rows_mulTFs,rows_mulreg);
rows_mulTFs = setdiff(rows_mulTFs,commonrows);
rows_mulreg = setdiff(rows_mulreg,commonrows);

% replace nan elements
TFs(rows_mulTFs) = cell2mat(cellfun(@(x) str2double(x(2:5)),TFsCell(rows_mulTFs),...
    'UniformOutput',false));
repeatTFs = cell2mat(cellfun(@(x) str2double(x(8:11)),TFsCell(rows_mulTFs),...
    'UniformOutput',false));
genes_reg(rows_mulreg) = cell2mat(cellfun(@(x) str2double(x(4:7)),genes_regCell(rows_mulreg),...
    'UniformOutput',false));

TFs = [TFs;repeatTFs];
genes_reg = [genes_reg;genes_reg(rows_mulTFs)];
IntType = [IntType;IntType(rows_mulTFs)];
TFs(commonrows) = [];
genes_reg(commonrows) = [];
IntType(commonrows) = [];

% unique protein indices:
allindices = sort(unique([TFs;genes_reg]));
N = length(allindices)+1;
% replacing elements with protein index
[~, TFinds_new] = ismember(TFs, allindices);
[~, genereginds_new] = ismember(genes_reg, allindices);

% create interaction matrix
IntPairs = [TFinds_new,genereginds_new];
[IntPairs,ia] = unique(IntPairs,'rows');
IntType = IntType(ia,:);
IntTypeMat = sparse(IntPairs(:,1),IntPairs(:,2),IntType,N,N);
IntTypeVec = IntTypeMat(:);

disp('interactions extracted.')
toc

%% Properties of the network
numInteractions = size(IntPairs,1);
density = length(IntType)/N^2;
TFinds = unique(TFinds_new);
numTFs = length(TFinds);
numInhibInts = sum(IntTypeVec == -1);
numActInts = sum(IntTypeVec == 1);

fprintf('Number of proteins (including ribosomes) N =  %d \n', N);
fprintf('Number of TFs =  %d \n', numTFs);
fprintf('Number of interactions =  %d, %d of them activating \n', ...
    numInteractions,sum(IntType==1));
fprintf('Interaction density =  %d \n', density);

% Extract self-regulating interactions
selfIntsInds = find(IntPairs(:,1)==IntPairs(:,2));
selfIntType = IntType(selfIntsInds);
numSelfActInts = sum(selfIntType==1);
numSelfInhibInts = sum(selfIntType==-1);
numSelfInts = numSelfActInts + numSelfInhibInts;
fprintf('Number of self interactions =  %d, %d of them activating \n', ...
    numSelfInts,numSelfActInts);

% Extract non-selfregulating interactions:
nonselfIntPairs = IntPairs(IntPairs(:,1)~=IntPairs(:,2),:);
assert(size(nonselfIntPairs,1)==numInteractions-numSelfInts,...
    'Inconsistent number of non-self interacting pairs!');
numNonSelfInts = size(nonselfIntPairs,1);

% For each gene, extract the number of other genes they regulate
% (out-degree), the number of other genes they are regulated by
% (in-degree).
outdeg = accumarray(nonselfIntPairs(:,1),1,[N-1,1]);
indeg = accumarray(nonselfIntPairs(:,2),1,[N-1,1]);

% Out-degree distribution of TFs (since non-TFs must have 0 out-degrees by
% definition)
outdeg_TFs = outdeg(TFinds);

% in-degree distributions of TFs and non-TFs:
indeg_TFs = indeg(TFinds);
indeg_nonTFs = indeg(setdiff(1:N-1,TFinds));

% number of TF-TF interactions that are activating/inhibiting
ifTFTFint = ismember(nonselfIntPairs(:,2),TFinds);
numNonSelf_TFTF_ActInts = sum(IntType(ifTFTFint==1)==1);
numNonSelf_TFTF_InhibInts = sum(IntType(ifTFTFint==1)==-1);
numNonSelf_TFTFints = numNonSelf_TFTF_ActInts+numNonSelf_TFTF_InhibInts;
fprintf('Number of non-self TF-TF interactions =  %d, %d of them activating \n', ...
    numNonSelf_TFTFints,numNonSelf_TFTF_ActInts);

disp('Properties of real network extracted.')
toc

%% Combine properties into structure
NetworkProp.N = N;
NetworkProp.numTFs = numTFs;
NetworkProp.numActInts = numActInts;
NetworkProp.numInhibInts = numInhibInts;
NetworkProp.numSelfActInts = numSelfActInts;
NetworkProp.numSelfInhibInts = numSelfInhibInts;
NetworkProp.outdeg_TFs = outdeg_TFs;
NetworkProp.indeg_TFs = indeg_TFs;
NetworkProp.indeg_nonTFs = indeg_nonTFs;

NetworkProp.numNonSelf_TFTF_ActInts = numNonSelf_TFTF_ActInts;
NetworkProp.numNonSelf_TFTF_InhibInts = numNonSelf_TFTF_InhibInts;
        
%% Define properties of the scrambling and other parameters
ScramblingMethod = 2; % if 1: we keep the same number of TFs and number of regulations
                         % 2: we only keep the number of regulations
                         % 3: keep in and out degrees the same
                         % 4: keep # of TFs and # of autoregulatory vs
                         % non-self interactions.
                         % 5: keep #TFs,#self vs Non-self, #TF-TF vs
                         % TF-nonTF interactions
k1 = 1;
phidist = 1;
maxfoldchangescan = 10.^[0.2,1.2:0.2:2.4];
% maxfoldchangescan = 10.^(2.6:0.2:3.6);
% gammamin = -1;
%gammamaxscan = 10;
n = 2;
numtrials = 10; %5;
numsets = 3;
                             
% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

% method of solving for steady-state
tspan = [0 1e10];
eps = 1e-5;
tmax = 400; % time limit in seconds for each attempt to find css using odesolver

%% Specify filename (where data will be stored)
runindx = 1;
minfc = round(min(maxfoldchangescan));
namestr = sprintf('ScrambledNetwork%d_n%d_minfc%d_numtrials%d_numsets%d_run%d',...
    [ScramblingMethod,n,minfc,numtrials,numsets,runindx]);

%% storage arrays
% css
ttraj_cell = cell(length(maxfoldchangescan),numtrials,numsets);
ctraj_cell = cell(length(maxfoldchangescan),numtrials,numsets);

% eigenvalues
maxlambdaMat = zeros(length(maxfoldchangescan),numtrials,numsets);
lambdaJtilde_cell = cell(length(maxfoldchangescan),numtrials,numsets);

% For checking convergence:
maxfvalmat = zeros(length(maxfoldchangescan),numtrials,numsets);
exitflagmat = zeros(length(maxfoldchangescan),numtrials,numsets);
temat = zeros(length(maxfoldchangescan),numtrials,numsets);

%% Main simulation
tic
for setindx = 1:numsets
    fprintf('setindx: %d \n',setindx);
    
    IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                             % 2nd col: Kd, 3rd col: Hill's coefficient n
    IntParamsMat(:,3) = n;
    
    for fcindx = 1:length(maxfoldchangescan)
        maxfc = maxfoldchangescan(fcindx);
        fprintf('max fold change: %d \n',maxfc);    
        
        for trialindx = 1:numtrials
            fprintf('trialindx: %d \n',trialindx);
            
            % Draw new network topology (scrambling)
            if ScramblingMethod > 0 
                IntTypeVec = ScrambleNetwork(NetworkProp, ScramblingMethod);
            end
            
            % Draw interaction parameters
            foldchangeVec = (rand(N*N,1).*(maxfc-1)).*(IntTypeVec~=0)+1;
            gammaVec = zeros(N*N,1);
            gammaVec(IntTypeVec==1) = foldchangeVec(IntTypeVec==1)-1;
            gammaVec(IntTypeVec==-1) = 1./foldchangeVec(IntTypeVec==-1)-1;
    
        
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
            c_guess = phivec; 

            % Interaction parameters 
            IntParamsMat(:,1:2) = [gammaVec,repmat(phivec,N,1)];
            
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
            fval = SetofEqns_v2(c_ss,IntParamsMat,phivec,fijintfunc);
            if ~isempty(ie)
                if max(abs(fval)) <= eps
                    exitflag = 1;
                else
                    exitflag = 0;
                    warning('max time exceeded but dynamics has not converged!');          
                end
            end
            
            ttraj_cell{fcindx,trialindx,setindx} = ttraj_kk;
            ctraj_cell{fcindx,trialindx,setindx} = ctraj_kk;
            exitflagmat(fcindx,trialindx,setindx) = exitflag;
            maxfvalmat(fcindx,trialindx,setindx) = max(abs(fval));
            temat(fcindx,trialindx,setindx) = te;
            
            % Calculate stability of steady state (real network):
            [JacobianMat,~,~] = CreateJacMat_method2(c_ss,IntParamsMat,phivec,k1);
            lambdavec = eig(JacobianMat);
            lambdavec_Jtilde = lambdavec./(k1.*c_ss(end))+1;
            maxlambdaMat(fcindx,trialindx,setindx) = max(real(lambdavec_Jtilde));
            lambdaJtilde_cell{fcindx,trialindx,setindx} = lambdavec_Jtilde;
            toc
        end
    end
    % Save workspace
    save(namestr);        
end                   
toc


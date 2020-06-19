% Here we investigate what happens to the dyanmics of the system when it is 
% at the edge of stability (e.g. by gradually tuning up the interaction 
% strength.)

close all; clear variables;
addpath('../');
rng('shuffle');

%% Define parameters
maxfcscan_up = 2;
maxfcscan_down = maxfcscan_up;
% maxfcscan_down = ones(1,3).*1e6;
assert(length(maxfcscan_up)==length(maxfcscan_down),...
    'maxfcscan_up and maxfcscan_down must have the same length');
N = 200; % number of genes (including ribosomes)
sparsity = 0.2;
mu = 1;
phidist = 1; % 1= uniform, 2=uniform in log-space
%phir = 0.1;
eps = 1e-5;
n = 1;
pdown = 0;

% fold change distribution
fcdist_up = 'uniform'; % either 'uniform','1overx' or '1overx2'
fcdist_down = 'uniform';

% IntTypeMat = full(sprandn(N,N,sparsity)); 
%IntTypeMat = randn(N,N); 
IntTypeMat = full(sprand(N,N,sparsity)); 
IntTypeMat(1:N+1:end) = 0; % no self-regulation
IntTypeMat(N,:) = 0; % ribosomes don't regulate other genes
%IntTypeMat = IntTypeMat~=0;
IntTypeMat(IntTypeMat~=0) = IntTypeMat(IntTypeMat~=0) - pdown;
IntTypeMat = sign(IntTypeMat);
IntTypeVec = IntTypeMat(:);
% IntType = 1 for upregulation; -1  for downregulation
[IntPairs_i,IntPairs_j] = find(IntTypeMat);
IntPairs = [IntPairs_i,IntPairs_j];

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
IntParamsMat(:,3) = n;
% kdMat = rand(N,N);
% kdMat = kdMat./repmat(sum(kdMat,1),N,1);
% IntParamsMat(:,2) = kdMat(:);


% functional form for interaction
fijintfunc = @(gamma,Kd,n,c) 1 + gamma.*c.^n./(Kd.^n+c.^n); 
dfijdcfunc = @(gamma,Kd,n,c) gamma.*n.*c.^(n-1).*Kd.^n./(Kd.^n+c.^n).^2; 

%% Storage matrices
% interaction parameters
gammaMat = zeros(N^2,length(maxfcscan_up));

cssmat = zeros(N,length(maxfcscan_up));
fiMat = zeros(N,length(maxfcscan_up));

% eigenvalues
eigmat = zeros(N,length(maxfcscan_up));
eigM1mat = zeros(N,length(maxfcscan_up));
eigM2mat = zeros(N,length(maxfcscan_up));


exitflagvec = zeros(1,length(maxfcscan_up));
maxfvalvec = zeros(1,length(maxfcscan_up));
tevec = zeros(1,length(maxfcscan_up));

ttraj_cell = cell(1,length(maxfcscan_up));
ctraj_cell = cell(1,length(maxfcscan_up));

%% Main simulation
maxtint = 1e10;
tspan = [0 maxtint];
tmax = 600; % time limit in seconds for each attempt to find css using odesolver
c_guess = phivec;
currtime = 0;
figure;
colormat = rand(N,3);
for fcindx = 1:length(maxfcscan_up)
    
    maxfc_up = maxfcscan_up(fcindx);
    maxfc_down = maxfcscan_down(fcindx);
            
    if max(maxfc_up,maxfc_down) < 100
        tspan = [0 1e8];
    else
        tspan = [0 maxtint];
    end
    
    fprintf('maxfc_up = %d, maxfc_down = %d \n',[maxfc_up,maxfc_down]);
    
    % Draw interaction parameters
    foldchangeVec = ones(N*N,1);
    if strcmp(fcdist_up,'uniform')
        foldchangeVec_up = (rand(N*N,1).*(maxfc_up-1)).*(IntTypeVec>0)+1;                
    elseif strcmp(fcdist_up,'1overx')
        logfcVec = (rand(N*N,1).*log(maxfc_up)).*(IntTypeVec>0);
        foldchangeVec_up = exp(logfcVec);
    elseif strcmp(fcdist_up,'1overx2')
        inversefcVec = 1-(rand(N*N,1).*(1-1/maxfc_up)).*(IntTypeVec>0);
        foldchangeVec_up = 1./inversefcVec;
    end
    foldchangeVec(IntTypeVec>0) = foldchangeVec_up(IntTypeVec>0);
    if strcmp(fcdist_down,'uniform')
        foldchangeVec_down = (rand(N*N,1).*(maxfc_down-1)).*(IntTypeVec<0)+1;
    elseif strcmp(fcdist_down,'1overx')
        logfcVec = (rand(N*N,1).*log(maxfc_down)).*(IntTypeVec<0);
        foldchangeVec_down = exp(logfcVec);
    elseif strcmp(fcdist_down,'1overx2')
        inversefcVec = 1-(rand(N*N,1).*(1-1/maxfc_down)).*(IntTypeVec<0);
        foldchangeVec_down = 1./inversefcVec;
    end
    foldchangeVec(IntTypeVec<0) = foldchangeVec_down(IntTypeVec<0);
    gammaVec = zeros(N*N,1);
    gammaVec(IntTypeVec==1) = foldchangeVec(IntTypeVec==1)-1;
    gammaVec(IntTypeVec==-1) = 1./foldchangeVec(IntTypeVec==-1)-1;
    
    IntParamsMat(:,1) = gammaVec;
    
    % store interactions:
    gammaMat(:,fcindx) = gammaVec;
    tstart = tic;
    opts = odeset('Events',@(t,y) eventfun(t,y,IntPairs,IntParamsMat,phivec,fijintfunc,eps,tstart,tmax));
    [tvec,ctraj,te,c_ss,ie] = ode45(@(t,y) SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,mu), ...
        tspan, c_guess, opts);
    
    if ~isempty(ie)
        exitflag = 1;
        
    else 
        exitflag = 0;
        te = tvec(end);
        c_ss = ctraj(end,:);
        warning('dynamics has not converged!');    
    end
    c_ss = c_ss.';
    [fval,fivec] = SetofEqns_v2(c_ss,IntParamsMat,phivec,fijintfunc);
    
%     if tvec(end) < maxtint
%         tvec = [tvec;maxtint];
%         ctraj = [ctraj;ctraj(end,:)];
%     end

    % Calculate eigenvalues
    [JacobianMat,meandlogfdcj,dlogfijdcj_mat] = CreateJacMat_method2(c_ss,...
        IntParamsMat,phivec,mu);
    M1 = dlogfijdcj_mat;
    M2 = repmat(meandlogfdcj,N,1);
    cM1 = repmat(c_ss,1,N).*M1;
    cM2 = repmat(c_ss,1,N).*M2;

    lambdavec = eig(JacobianMat);
    JtildeMat = JacobianMat./(-mu*c_ss(end))-eye(N);
    lambdavec_Jtilde = lambdavec./(-mu.*c_ss(end))-1;
    eigmat(:,fcindx) = lambdavec_Jtilde;
    eigM1mat(:,fcindx) = eig(cM1);
    eigM2mat(:,fcindx) = eig(cM2);
    
    ctraj_cell{fcindx} = ctraj;
    ttraj_cell{fcindx} = tvec+currtime;
    
    % Plot
    for kk =1:N
%         plot(log10(tvec+currtime),ctraj(:,kk),'color',colormat(kk,:),'DisplayName',...
%             strcat('p_{',num2str(kk),'}'));
        plot(tvec+currtime,ctraj(:,kk),'color',colormat(kk,:),'DisplayName',...
            strcat('p_{',num2str(kk),'}'));
        hold on
    end
%     plot(log10(tvec+currtime),ctraj(:,N),'k-.');
    plot(tvec+currtime,ctraj(:,N),'k-.');
    hold on
%     plot(log10(ones(1,2).*(tvec(end)+currtime)),[0,max(ctraj(end,:))],'r-.');
    plot(ones(1,2).*(tvec(end)+currtime),[0,max(ctraj(end,:))],'r-.');
    hold on
        
    % update and store variables
    exitflagvec(fcindx) = exitflag;
    maxfvalvec(fcindx) = max(abs(fval));
    tevec(fcindx) = te;
    cssmat(:,fcindx) = c_ss;
    fiMat(:,fcindx) = fivec;
    currtime = currtime+tvec(end);
    c_guess = c_ss;
    
end
    
xlabel('t');
ylabel('c');

%%

% tspan = [0 1e9];
% %[tvec_full,ctraj_full,te,c_ss,ie] = ode45(@(t,y) SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,mu), ...
% %        tspan, c_guess, opts);
% [tvec_full,ctraj_full] = ode45(@(t,y) SetofEqns_forodesolver_v2(t,y,IntParamsMat,phivec,fijintfunc,mu), ...
%     tspan, c_guess);
% figure;
% for kk =1:N
%     plot(tvec_full,ctraj_full(:,kk));
%     hold on
% end
% plot(tvec_full,ctraj_full(:,N),'k-.');
% hold on
% xlabel('t');
% ylabel('c');
% 
% 

%% plot eigenvalue spectrum
for fcindx = 1:min(length(maxfcscan_up),4)
    figure;
    subplot(2,2,1)
    plot(real(eigmat(:,fcindx)),imag(eigmat(:,fcindx)),'x');
    xlabel('real'); ylabel('imag');
    title('spectrum of M matrix');
    
    subplot(2,2,2)
    plot(real(eigM1mat(:,fcindx)),imag(eigM1mat(:,fcindx)),'x');
    xlabel('real'); ylabel('imag');
    title('spectrum of M1 matrix');
    
    subplot(2,2,4)
    plot(real(-eigM2mat(:,fcindx)),imag(-eigM2mat(:,fcindx)),'x');
    xlabel('real'); ylabel('imag');
    title('spectrum of M2 matrix');

end

%% replot concentrations
figure;
for fcindx = 1:min(length(maxfcscan_up),4)
    ttraj = ttraj_cell{fcindx};
    ctraj = ctraj_cell{fcindx};
    subplot(2,2,fcindx)
    for kk =1:N
        plot(ttraj,ctraj(:,kk),'color',colormat(kk,:),'DisplayName',...
            strcat('p_{',num2str(kk),'}'));
        hold on
    end
end
xlim([0 3e4])

%% Plot distributions of concentrations
figure;
for fcindx = 1:min(length(maxfcscan_up),4)
    c_ss = cssmat(:,fcindx);
    subplot(2,2,fcindx);
    histogram(log10(c_ss),'Normalization','pdf');
    xlabel('log10(c)');
    ylabel('P');
end

%% Save
% EgIndx = 1;
% save(strcat('FigChaoticEx',num2str(EgIndx)));
% save(strcat('FigOscEx',num2str(EgIndx)));

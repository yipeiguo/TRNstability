% This is a script for analyzing data from the exploring bipartite-like
% network runs.
% In this version 2, we play with more information like the Jacobian
% elements, looking that the Tmatrix structure, reduced Tmatrix etc.

close all; 
clear variables;
addpath('../../../'); % for using 'CreateJacMat_method2.m'

%% Define parameters
fvalthres_new = 5e-3;
numexts = 2; %2; %3;
numrunsvec = [9,3]; % number of runs for each extension
fncell = cell(1,numexts);

pathname = './';
fncell{1} = 'BipartiteLikeData_N2275_q211_numSelfInts0_numIntsTFnonTF5000_minnumTFTF0_maxnumTFTF400_n2_maxfc1000_numtrials10_run';

fncell{2} = 'BipartiteLikeData_N2275_q211_numSelfInts0_numIntsTFnonTF5000_minnumTFTF500_maxnumTFTF5000_n2_maxfc1000_numtrials10_run';

trialdim = 2; setdim = 3;

%% Load and combine data sets
gammacell_allexts = cell(1,numexts);
cssMatCell_allexts = cell(1,numexts);
phiMatCell_allexts = cell(1,numexts);
fiMatCell_allexts = cell(1,numexts);
maxlambdaMatCell_allexts = cell(1,numexts);
maxfvalMatCell_allexts = cell(1,numexts);
TmatCell_allexts = cell(1,numexts);

% ctraj_cell_allsets = cell(1,numexts);
% ttraj_cell_allsets = cell(1,numexts);

M1_maxeigMat_allexts = cell(1,numexts);
M2_maxeigMat_allexts = cell(1,numexts);
M2_mineigMat_allexts = cell(1,numexts);
cM1_maxeigMat_allexts = cell(1,numexts);
cM2_maxeigMat_allexts = cell(1,numexts);
cM2_mineigMat_allexts = cell(1,numexts);

numsetsvec = zeros(1,numexts);
numtrialsvec = zeros(1,numexts);
toaccept_allexts = cell(1,numexts);
tic
for extindx = 1:numexts
    fn = fncell{extindx};
%     load(strcat(pathname,fn,'1'),'gamma_cell','cssMat','phiMat','fiMat','numTFTFints_scan','maxlambdaMat',...
%         'maxfvalmat','ttraj_cell','ctraj_cell','Tmat_cell','NetworkProp','k1','n','N','q','exitflagmat');
    load(strcat(pathname,fn,'1'),'gamma_cell','cssMat','phiMat','fiMat','numTFTFints_scan','maxlambdaMat',...
        'maxfvalmat','Tmat_cell','NetworkProp','k1','n','N','q','exitflagmat');
    if ~exist('numTFTFints_scan_all','var')
        numTFTFints_scan_all = numTFTFints_scan;
    else
        numTFTFints_scan_all = cat(2,numTFTFints_scan_all,numTFTFints_scan);
    end
    % check number of valid sets 
    ifvalidsets = squeeze(sum(exitflagmat~=0,setdiff(1:ndims(exitflagmat),setdim)))~=0;
    numvalidsets = sum(ifvalidsets);
    if numvalidsets < size(exitflagmat,setdim) % trim variables
        if setdim == 3
            gamma_cell(:,:,ifvalidsets==0) = [];
            cssMat(:,:,:,ifvalidsets==0) = [];
            phiMat(:,:,:,ifvalidsets==0) = [];
            fiMat(:,:,:,ifvalidsets==0) = [];
            maxlambdaMat(:,:,ifvalidsets==0) = [];
            maxfvalmat(:,:,ifvalidsets==0) = [];
%             if exist('ctraj_cell','var')
%                 ctraj_cell(:,:,ifvalidsets==0) = [];
%                 ttraj_cell(:,:,ifvalidsets==0) = [];
%             end
            Tmat_cell(:,:,ifvalidsets==0) = [];
        end
    end
                
    gammacell_all = gamma_cell;
    cssMat_all = cssMat;
    phiMat_all = phiMat;
    fiMat_all = fiMat;
    maxlambdaMat_all = maxlambdaMat;
    maxfvalmat_all = maxfvalmat;
%     if exist('ctraj_cell','var')
%         ctraj_cell_all = ctraj_cell;
%         ttraj_cell_all = ttraj_cell;
%         clear ctraj_cell ttraj_cell
%     end
    Tmatcell_all = Tmat_cell;
    
    for runindx = 2:numrunsvec(extindx)
%         load(strcat(pathname,fn,num2str(runindx)),'gamma_cell','cssMat','phiMat','fiMat',...
%             'maxlambdaMat','maxfvalmat','ttraj_cell','ctraj_cell','Tmat_cell','exitflagmat');
        load(strcat(pathname,fn,num2str(runindx)),'gamma_cell','cssMat','phiMat','fiMat',...
            'maxlambdaMat','maxfvalmat','Tmat_cell','exitflagmat');
        
        % check number of valid sets 
        ifvalidsets = squeeze(sum(exitflagmat~=0,setdiff(1:ndims(exitflagmat),setdim)))~=0;
        numvalidsets = sum(ifvalidsets);
        if numvalidsets < size(exitflagmat,setdim) % trim variables
            if setdim == 3
                gamma_cell(:,:,ifvalidsets==0) = [];
                cssMat(:,:,:,ifvalidsets==0) = [];
                phiMat(:,:,:,ifvalidsets==0) = [];
                fiMat(:,:,:,ifvalidsets==0) = [];
                maxlambdaMat(:,:,ifvalidsets==0) = [];
                maxfvalmat(:,:,ifvalidsets==0) = [];
%                 if exist('ctraj_cell','var')
%                     ctraj_cell(:,:,ifvalidsets==0) = [];
%                     ttraj_cell(:,:,ifvalidsets==0) = [];
%                 end
                Tmat_cell(:,:,ifvalidsets==0) = [];
            end
        end
        
        gammacell_all = cat(setdim,gammacell_all,gamma_cell);
        cssMat_all = cat(setdim+1,cssMat_all,cssMat);
        phiMat_all = cat(setdim+1,phiMat_all,phiMat);
        fiMat_all = cat(setdim+1,fiMat_all,fiMat);
        maxlambdaMat_all = cat(setdim,maxlambdaMat_all,maxlambdaMat);
        maxfvalmat_all = cat(setdim,maxfvalmat_all,maxfvalmat);
        Tmatcell_all = cat(setdim,Tmatcell_all,Tmat_cell);
%         if exist('ctraj_cell','var') && exist('ctraj_cell_all','var') 
%             ctraj_cell_all = cat(setdim,ctraj_cell_all,ctraj_cell);
%             ttraj_cell_all = cat(setdim,ttraj_cell_all,ttraj_cell);
%             clear ctraj_cell ttraj_cell
%         end
    end
    
    toaccept_all = maxfvalmat_all<fvalthres_new;
    numsets = size(toaccept_all,setdim);
    numtrials = size(toaccept_all,trialdim);
    toaccept_allexts{extindx} = toaccept_all;

    numsetsvec(extindx) = numsets;
    numtrialsvec(extindx) = numtrials;
    gammacell_allexts{extindx} = gammacell_all;
    cssMatCell_allexts{extindx} = cssMat_all;
    phiMatCell_allexts{extindx} = phiMat_all;
    fiMatCell_allexts{extindx} = fiMat_all;
    maxlambdaMatCell_allexts{extindx} = maxlambdaMat_all;
    maxfvalMatCell_allexts{extindx} = maxfvalmat_all;
    TmatCell_allexts{extindx} = Tmatcell_all;
    
%     if exist('ctraj_cell_all','var')
%         ctraj_cell_allsets{extindx} = ctraj_cell_all;
%         ttraj_cell_allsets{extindx} = ttraj_cell_all;
%         
%         clear ctraj_cell_all ttraj_cell_all
%     end    
    toc
    
    % This section is for calculating the eigenvalues of M1,M2,c*M1,c*M2 separately
    M1maxeigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    M2maxeigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    M2mineigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    cM1maxeigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    cM2maxeigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    cM2mineigMat = zeros(length(numTFTFints_scan),numtrials,numsets);
    for setindx = 1:numsets
        fprintf('setindx = %d \n',setindx);
        for kk = 1:length(numTFTFints_scan)
%             fprintf('kk = %d \n',kk);
            for trialindx = 1:numtrials
%                 fprintf('trialindx = %d \n',trialindx);
                
                c_ss = cssMat_all(:,kk,trialindx,setindx);
                phivec = phiMat_all(:,kk,trialindx,setindx);
                gammaVals = gammacell_all{kk,trialindx,setindx};
                
                gammaVec = zeros(N*N,1);
                gammaVec(gammaVals(:,1)) = gammaVals(:,2);
                IntParamsMat = zeros(N*N,3); % 1st col: strength of interaction, 
                                         % 2nd col: Kd, 3rd col: Hill's coefficient n
                IntParamsMat(:,3) = n;
                IntParamsMat(:,1:2) = [gammaVec,repmat(phivec,N,1)];

                [~,meandlogfdcj,dlogfijdcj_mat] = CreateJacMat_method2(c_ss,IntParamsMat,phivec,k1);
                M1 = dlogfijdcj_mat;
                M2 = repmat(meandlogfdcj,N,1);
                cM1 = repmat(c_ss,1,N).*M1;
                cM2 = repmat(c_ss,1,N).*M2;
                M1maxeigMat(kk,trialindx,setindx) = max(real(eig(M1(1:q,1:q))));
                M2maxeigMat(kk,trialindx,setindx) = max(real(eig(M2(1:q,1:q))));
                M2mineigMat(kk,trialindx,setindx) = min(real(eig(M2(1:q,1:q))));
                cM1maxeigMat(kk,trialindx,setindx) = max(real(eig(cM1(1:q,1:q))));
                cM2maxeigMat(kk,trialindx,setindx) = max(real(eig(cM2(1:q,1:q))));
                cM2mineigMat(kk,trialindx,setindx) = min(real(eig(cM2(1:q,1:q))));
            end
        end
    end
    M1_maxeigMat_allexts{extindx} = M1maxeigMat;
    M2_maxeigMat_allexts{extindx} = M2maxeigMat;
    M2_mineigMat_allexts{extindx} = M2mineigMat;
    cM1_maxeigMat_allexts{extindx} = cM1maxeigMat;
    cM2_maxeigMat_allexts{extindx} = cM2maxeigMat;
    cM2_mineigMat_allexts{extindx} = cM2mineigMat;
    toc
    
end
toc

disp('All data imported');

%% Plot distribution of matrix elements
figure;
subplot(2,2,1);
histogram(log10(M1(M1>0)),'Normalization','pdf');
hold on
histogram(log10(abs(M1(M1<0))),'Normalization','pdf');
xlabel('log_{10}(M1 elements)');
ylabel('P');
legend('positive','negative');

subplot(2,2,2);
histogram(log10(M2(M2>0)),'Normalization','pdf');
hold on
histogram(log10(abs(M2(M2<0))),'Normalization','pdf');
xlabel('log_{10}(M2 elements)');
ylabel('P');

subplot(2,2,3);
histogram(log10(cM1(cM1>0)),'Normalization','pdf');
hold on
histogram(log10(abs(cM1(cM1<0))),'Normalization','pdf');
xlabel('log_{10}(cM1 elements)');
ylabel('P');

subplot(2,2,4);
histogram(log10(cM2(cM2>0)),'Normalization','pdf');
hold on
histogram(log10(abs(cM2(cM2<0))),'Normalization','pdf');
xlabel('log_{10}(cM2 elements)');
ylabel('P');

%% mean max eigenvalue
for extindx = 1:numexts
    toaccept_all = toaccept_allexts{extindx};
    maxlambdaMat_all = maxlambdaMatCell_allexts{extindx};
    
    maxlambdaMat_mean = sum(maxlambdaMat_all.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    maxlambdaMat_grandmean = mean(maxlambdaMat_mean,setdim);
    maxlambdaMat_std = std(maxlambdaMat_mean,[],setdim);
    
    if ~exist('maxlambdaMat_grandmean_all','var')
        maxlambdaMat_grandmean_all = maxlambdaMat_grandmean;
        maxlambdaMat_std_all = maxlambdaMat_std;
    else
        maxlambdaMat_grandmean_all = cat(1,maxlambdaMat_grandmean_all,maxlambdaMat_grandmean);
        maxlambdaMat_std_all = cat(1,maxlambdaMat_std_all,maxlambdaMat_std);
    end
end

%% mean max eigenvalues of M1,M2,cM1,cM2
for extindx = 1:numexts
    toaccept_all = toaccept_allexts{extindx};
    M1maxeigMat = M1_maxeigMat_allexts{extindx};
    M2maxeigMat = M2_maxeigMat_allexts{extindx};
    M2mineigMat = M2_mineigMat_allexts{extindx};
    cM1maxeigMat = cM1_maxeigMat_allexts{extindx};
    cM2maxeigMat = cM2_maxeigMat_allexts{extindx};
    cM2mineigMat = cM2_mineigMat_allexts{extindx};
    
    M1maxeig_mean = sum(M1maxeigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    M1maxeig_grandmean = mean(M1maxeig_mean,setdim);
    M1maxeig_std = std(M1maxeig_mean,[],setdim);
    
    M2maxeig_mean = sum(M2maxeigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    M2maxeig_grandmean = mean(M2maxeig_mean,setdim);
    M2maxeig_std = std(M2maxeig_mean,[],setdim);
    
    M2mineig_mean = sum(M2mineigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    M2mineig_grandmean = mean(M2mineig_mean,setdim);
    M2mineig_std = std(M2mineig_mean,[],setdim);
    
    cM1maxeig_mean = sum(cM1maxeigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    cM1maxeig_grandmean = mean(cM1maxeig_mean,setdim);
    cM1maxeig_std = std(cM1maxeig_mean,[],setdim);
    
    cM2maxeig_mean = sum(cM2maxeigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    cM2maxeig_grandmean = mean(cM2maxeig_mean,setdim);
    cM2maxeig_std = std(cM2maxeig_mean,[],setdim);
    
    cM2mineig_mean = sum(cM2mineigMat.*toaccept_all,trialdim)./...
        sum(toaccept_all,trialdim);
    cM2mineig_grandmean = mean(cM2mineig_mean,setdim);
    cM2mineig_std = std(cM2mineig_mean,[],setdim);
    
    if ~exist('M1maxeig_grandmean_all','var')
        M1maxeig_grandmean_all = M1maxeig_grandmean;
        M1maxeig_std_all = M1maxeig_std;
        
        M2maxeig_grandmean_all = M2maxeig_grandmean;
        M2maxeig_std_all = M2maxeig_std;
        
        M2mineig_grandmean_all = M2mineig_grandmean;
        M2mineig_std_all = M2mineig_std;
        
        cM1maxeig_grandmean_all = cM1maxeig_grandmean;
        cM1maxeig_std_all = cM1maxeig_std;
        
        cM2maxeig_grandmean_all = cM2maxeig_grandmean;
        cM2maxeig_std_all = cM2maxeig_std;
        
        cM2mineig_grandmean_all = cM2mineig_grandmean;
        cM2mineig_std_all = cM2mineig_std;
    else
        M1maxeig_grandmean_all = cat(1,M1maxeig_grandmean_all,M1maxeig_grandmean);
        M1maxeig_std_all = cat(1,M1maxeig_std_all,M1maxeig_std);
        
        M2maxeig_grandmean_all = cat(1,M2maxeig_grandmean_all,M2maxeig_grandmean);
        M2maxeig_std_all = cat(1,M2maxeig_std_all,M2maxeig_std);
        
        M2mineig_grandmean_all = cat(1,M2mineig_grandmean_all,M2mineig_grandmean);
        M2mineig_std_all = cat(1,M2mineig_std_all,M2mineig_std);
        
        cM1maxeig_grandmean_all = cat(1,cM1maxeig_grandmean_all,cM1maxeig_grandmean);
        cM1maxeig_std_all = cat(1,cM1maxeig_std_all,cM1maxeig_std);
        
        cM2maxeig_grandmean_all = cat(1,cM2maxeig_grandmean_all,cM2maxeig_grandmean);
        cM2maxeig_std_all = cat(1,cM2maxeig_std_all,cM2maxeig_std);
        
        cM2mineig_grandmean_all = cat(1,cM2mineig_grandmean_all,cM2mineig_grandmean);
        cM2mineig_std_all = cat(1,cM2mineig_std_all,cM2mineig_std);
    end
end


%% Probability of being stable
pvec = [0.25,0.75];

for extindx = 1:numexts
    toaccept_all = toaccept_allexts{extindx};

    Pstable = sum(toaccept_all,trialdim)./size(toaccept_all,trialdim);
    meanPstable = mean(Pstable,setdim);
    sdPstable = std(Pstable,[],setdim);
    ulbPstable = quantile(Pstable,pvec,setdim);
    if ~exist('meanPstable_all','var')
        meanPstable_all = meanPstable;
        sdPstable_all = sdPstable;
        ulbPstable_all = ulbPstable;
    else
        meanPstable_all = cat(1,meanPstable_all,meanPstable);
        sdPstable_all = cat(1,sdPstable_all,sdPstable);
        ulbPstable_all = cat(1,ulbPstable_all,ulbPstable);
    end
end

%% Plot
q = NetworkProp.numTFs;
rhoscan = numTFTFints_scan_all./(q*(q-1));
figure;
subplot(2,2,1);
errorbar(rhoscan,maxlambdaMat_grandmean_all,maxlambdaMat_std_all,'x');
hold on
errorbar(rhoscan,cM1maxeig_grandmean_all,cM1maxeig_std_all,'x');
xlabel('\rho_q');
% xlabel('# TF-TF interactions');
ylabel('max lambda');
legend('M','cM1','location','best')

subplot(2,2,3);
errorbar(log10(rhoscan),log10(maxlambdaMat_grandmean_all),...
    log10(maxlambdaMat_grandmean_all) - log10(maxlambdaMat_grandmean_all-maxlambdaMat_std_all),...
    log10(maxlambdaMat_grandmean_all+maxlambdaMat_std_all) - log10(maxlambdaMat_grandmean_all),'x');
hold on
% errorbar(log10(rhoscan),log10(M1maxeig_grandmean_all),...
%     log10(M1maxeig_grandmean_all) - log10(M1maxeig_grandmean_all-M1maxeig_std_all),...
%     log10(M1maxeig_grandmean_all+M1maxeig_std_all) - log10(M1maxeig_grandmean_all),'x');
% hold on
% errorbar(log10(rhoscan),log10(M2maxeig_grandmean_all),...
%     log10(M2maxeig_grandmean_all) - log10(M2maxeig_grandmean_all-M2maxeig_std_all),...
%     log10(M2maxeig_grandmean_all+M2maxeig_std_all) - log10(M2maxeig_grandmean_all),'x');
% hold on
errorbar(log10(rhoscan),log10(cM1maxeig_grandmean_all),...
    log10(cM1maxeig_grandmean_all) - log10(cM1maxeig_grandmean_all-cM1maxeig_std_all),...
    log10(cM1maxeig_grandmean_all+cM1maxeig_std_all) - log10(cM1maxeig_grandmean_all),'x');
hold on
% errorbar(log10(rhoscan),log10(cM2maxeig_grandmean_all),...
%     log10(cM2maxeig_grandmean_all) - log10(cM2maxeig_grandmean_all-cM2maxeig_std_all),...
%     log10(cM2maxeig_grandmean_all+cM2maxeig_std_all) - log10(cM2maxeig_grandmean_all),'x');

xlabel('log10(rho)');
ylabel('log10(max lambda)');
legend('M','cM1','location','best')

subplot(2,2,2);
errorbar(numTFTFints_scan_all, meanPstable_all, ...
    meanPstable_all-ulbPstable_all(:,:,1),...
    ulbPstable_all(:,:,2)-meanPstable_all,'x');
xlabel('# TF-TF interactions');
ylabel('P(stable)');

subplot(2,2,4);
% errorbar(log10(rhoscan),log10(meanPstable_all),...
%     log10(meanPstable_all) - log10(meanPstable_all-ulbPstable_all(:,:,1)),...
%     log10(meanPstable_all+ulbPstable_all(:,:,2)) - log10(meanPstable_all),'x');
errorbar(log10(rhoscan), meanPstable_all, ...
    meanPstable_all-ulbPstable_all(:,:,1),...
    ulbPstable_all(:,:,2)-meanPstable_all,'x');
% xlabel('\rho');
xlabel('log10(\rho)');
ylabel('P(stable)');

%% Exploring eigenvalues of cM1,cM2,M1,M2
figure;
subplot(2,2,1);
errorbar(rhoscan,cM1maxeig_grandmean_all,cM1maxeig_std_all,'x');
hold on
errorbar(rhoscan,cM2maxeig_grandmean_all,cM2maxeig_std_all,'x');
hold on
errorbar(rhoscan,-cM2mineig_grandmean_all,cM2mineig_std_all,'o');
xlabel('\rho_q');
ylabel('\lambda_{max}');
title('cM matrix');
legend('cM1','cM2(max)','cM2(min)');
% legend('cM1','cM2');

subplot(2,2,2);
errorbar(rhoscan,M1maxeig_grandmean_all,M1maxeig_std_all,'x');
hold on
errorbar(rhoscan,M2maxeig_grandmean_all,M2maxeig_std_all,'x');
hold on
errorbar(rhoscan,-M2mineig_grandmean_all,M2mineig_std_all,'o');
xlabel('\rho_q');
ylabel('\lambda_{max}');
title('M matrix');
legend('M1','M2(max)','M2(min)');

subplot(2,2,3);
errorbar(log10(rhoscan),log10(cM1maxeig_grandmean_all),...
    log10(cM1maxeig_grandmean_all) - log10(cM1maxeig_grandmean_all-cM1maxeig_std_all),...
    log10(cM1maxeig_grandmean_all+cM1maxeig_std_all) - log10(cM1maxeig_grandmean_all),'x');
hold on
errorbar(log10(rhoscan),log10(cM2maxeig_grandmean_all),...
    log10(cM2maxeig_grandmean_all) - log10(cM2maxeig_grandmean_all-cM2maxeig_std_all),...
    log10(cM2maxeig_grandmean_all+cM2maxeig_std_all) - log10(cM2maxeig_grandmean_all),'x');
hold on
errorbar(log10(rhoscan),log10(abs(cM2mineig_grandmean_all)),...
    log10(-cM2mineig_grandmean_all) - log10(-cM2mineig_grandmean_all-cM2mineig_std_all),...
    log10(-cM2mineig_grandmean_all+cM2mineig_std_all) - log10(-cM2mineig_grandmean_all),'o');
xlabel('log10(\rho_q)');
ylabel('log10(\lambda_{max})');
legend('cM1','cM2(max)','cM2(min)');

subplot(2,2,4);
errorbar(log10(rhoscan),log10(M1maxeig_grandmean_all),...
    log10(M1maxeig_grandmean_all) - log10(M1maxeig_grandmean_all-M1maxeig_std_all),...
    log10(M1maxeig_grandmean_all+M1maxeig_std_all) - log10(M1maxeig_grandmean_all),'x');
hold on
errorbar(log10(rhoscan),log10(M2maxeig_grandmean_all),...
    log10(M2maxeig_grandmean_all) - log10(M2maxeig_grandmean_all-M2maxeig_std_all),...
    log10(M2maxeig_grandmean_all+M2maxeig_std_all) - log10(M2maxeig_grandmean_all),'x');
hold on
errorbar(log10(rhoscan),log10(-M2mineig_grandmean_all),...
    log10(-M2mineig_grandmean_all) - log10(-M2mineig_grandmean_all-M2mineig_std_all),...
    log10(-M2mineig_grandmean_all+M2mineig_std_all) - log10(-M2mineig_grandmean_all),'x');
xlabel('log10(\rho_q)');
ylabel('log10(\lambda_{max})');
legend('M1','M2 (max)','M2 (min)');

%% Visualize trajectories of interest
% extindxOI = 1;
% toaccept_allOI = toaccept_allexts{extindxOI};
% ttraj_cell_allOI = ttraj_cell_allsets{extindxOI};
% ctraj_cell_allOI = ctraj_cell_allsets{extindxOI};
% if ~isempty(ctraj_cell_allOI)
%     rhoIndxOI = 21; trialIndxOI = 4; setIndxOI = 3;
%     ttrajOI = ttraj_cell_allOI{rhoIndxOI,trialIndxOI,setIndxOI};
%     ctrajOI = ctraj_cell_allOI{rhoIndxOI,trialIndxOI,setIndxOI};
% 
%     figure;
%     if ~exist('N','var')
%         N = size(ctrajOI,2);
%     end
%     proteinind = 1:N;
%     for kk =1:length(proteinind)
%         plot(ttrajOI,ctrajOI(:,proteinind(kk)),'DisplayName',...
%             strcat('p_{',num2str(proteinind(kk)),'}'));
%         hold on
%     end
%     plot(ttrajOI,ctrajOI(:,proteinind(N)),'k-','DisplayName','r');
%     xlabel('t');
%     ylabel('c');
% end

%% save data
% newfn = sprintf('AnalyzedData_N%d_q%d_numIntsTFnonTF%d_minnumTFTF%d_maxnumTFTF%d_n%d_maxfc%d',...
%             [N,q,4500,min(numTFTFints_scan),max(numTFTFints_scan),n,1000]);
% save(strcat(pathname,newfn));


% save('Fig6Data_numSelfInts0_numIntsTFnonTF5000_v3_moreruns');
% save('WorkspaceWithM1M2_numSelfInts0_numIntsTFnonTF5000_repeat2');
% save('Fig6Data_numSelfInts134_numInts5655_v2');









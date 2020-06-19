% This is a script for analyzing data from the scrambling
% network runs.

close all; 
clear variables;
% addpath('../../../'); 

%% Define parameters
fvalthres_new = 5e-3;
numexts = 2; % number of 'foldchange' regions
numrunsvec = [7,7]; % number of runs for each extension
fncell = cell(1,numexts);

pathname = './';
% fncell{1} = 'ScrambledNetwork2_n2_minfc2_numtrials10_numsets2_run';
fncell{1} = 'ScrambledNetwork1_n2_minfc2_numtrials10_numsets2_run';

% fncell{2} = 'ScrambledNetwork2_n2_minfc398_numtrials10_numsets2_run';
fncell{2} = 'ScrambledNetwork1_n2_minfc398_numtrials10_numsets2_run';

trialdim = 2; setdim = 3;

%% Load and combine data sets
% maxlambdaMatCell_allexts = cell(1,numexts);
maxfvalMatCell_allexts = cell(1,numexts);
TmatCell_allexts = cell(1,numexts);

% ctraj_cell_allsets = cell(1,numexts);
% ttraj_cell_allsets = cell(1,numexts);

numsetsvec = zeros(1,numexts);
numtrialsvec = zeros(1,numexts);
toaccept_allexts = cell(1,numexts);
tic
for extindx = 1:numexts
    fn = fncell{extindx};
%     load(strcat(pathname,fn,'1'),'maxfoldchangescan','maxfvalmat',...
%         'ttraj_cell','ctraj_cell','exitflagmat');
    load(strcat(pathname,fn,'1'),'maxfoldchangescan','maxfvalmat','exitflagmat');
    if ~exist('log10fc_scan_all','var')
        log10fc_scan_all = log10(maxfoldchangescan);
    else
        log10fc_scan_all = cat(2,log10fc_scan_all,log10(maxfoldchangescan));
    end
    % check number of valid sets 
    ifvalidsets = squeeze(sum(exitflagmat~=0,setdiff(1:ndims(exitflagmat),setdim)))~=0;
    numvalidsets = sum(ifvalidsets);
    if numvalidsets < size(exitflagmat,setdim) % trim variables
        if setdim == 3
            maxfvalmat(:,:,ifvalidsets==0) = [];
%             exitflagmat(:,:,ifvalidsets==0) = [];
%             if exist('ctraj_cell','var')
%                 ctraj_cell(:,:,ifvalidsets==0) = [];
%                 ttraj_cell(:,:,ifvalidsets==0) = [];
%             end
        end
    end
                
    maxfvalmat_all = maxfvalmat;
%     if exist('ctraj_cell','var')
%         ctraj_cell_all = ctraj_cell;
%         ttraj_cell_all = ttraj_cell;
%         clear ctraj_cell ttraj_cell
%     end
     clear maxfvalmat
    
    for runindx = 2:numrunsvec(extindx)
%         load(strcat(pathname,fn,num2str(runindx)),'maxfoldchangescan',...
%             'maxfvalmat','ttraj_cell','ctraj_cell','exitflagmat');
        load(strcat(pathname,fn,num2str(runindx)),'maxfoldchangescan',...
            'maxfvalmat','exitflagmat');
        
        % check number of valid sets 
        ifvalidsets = squeeze(sum(exitflagmat~=0,setdiff(1:ndims(exitflagmat),setdim)))~=0;
        numvalidsets = sum(ifvalidsets);
        if numvalidsets < size(exitflagmat,setdim) % trim variables
            if setdim == 3
                maxfvalmat(:,:,ifvalidsets==0) = [];
%                 if exist('ctraj_cell','var')   
%                     ctraj_cell(:,:,ifvalidsets==0) = [];
%                     ttraj_cell(:,:,ifvalidsets==0) = [];
%                 end
            end
        end
        
        maxfvalmat_all = cat(setdim,maxfvalmat_all,maxfvalmat);
%         if exist('ctraj_cell','var') && exist('ctraj_cell_all','var')    
%             ctraj_cell_all = cat(setdim,ctraj_cell_all,ctraj_cell);
%             ttraj_cell_all = cat(setdim,ttraj_cell_all,ttraj_cell);
%             clear ctraj_cell ttraj_cell
%         end
        clear maxfvalmat
    end
    
    toaccept_all = maxfvalmat_all<fvalthres_new;
    numsets = size(toaccept_all,setdim);
    numtrials = size(toaccept_all,trialdim);
    toaccept_allexts{extindx} = toaccept_all;

    numsetsvec(extindx) = numsets;
    numtrialsvec(extindx) = numtrials;
    maxfvalMatCell_allexts{extindx} = maxfvalmat_all;
    
%     if exist('ctraj_cell_all','var')    
%         ctraj_cell_allsets{extindx} = ctraj_cell_all;
%         ttraj_cell_allsets{extindx} = ttraj_cell_all;
%         clear ctraj_cell_all ttraj_cell_all
%     end
%     
    
    toc
        
end
toc

disp('All data imported');


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
figure;
errorbar(log10fc_scan_all, meanPstable_all, ...
    meanPstable_all-ulbPstable_all(:,:,1),...
    ulbPstable_all(:,:,2)-meanPstable_all,'x');
% errorbar(log10fc_scan_all, meanPstable_all, sdPstable_all,'x');
% ylim([0.2,1])
xlabel('log10(\Omega_{max})');
ylabel('P(stable)');


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
% save('FigData_Pstable_scrambled1_v4_15sets');










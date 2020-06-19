% This is a script to analyze cluster data for eigenvalue scaling

close all; 
clear variables;

%% load data
pathname = './ClusterData/';
fn = 'EigScaling_randomDAG_phidist1_rho0pt01_n1_Nmin100_Nmax1000_fcmin2_fcmax2_numtrials10_run1';

load(strcat(pathname,fn));

%% Plot variation with N
if length(Nscan)> 1
    figure; 
    subplot(2,2,1);
    names = cell(1,length(maxfcscan));
    storeindx = 1;
    for fcindx = [1]
        avminlambdavec = zeros(1,length(Nscan));
        sdminlambdavec = zeros(1,length(Nscan));    
        for Nindx = 1:length(Nscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avminlambdavec(Nindx) = mean(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
            sdminlambdavec(Nindx) = std(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
        end
        errorbar(Nscan,avminlambdavec,...
            sdminlambdavec,'x','LineStyle','none')
        hold on
        names{storeindx} = strcat('\gamma: ',num2str(maxfcscan(fcindx)));
        storeindx = storeindx+1;
    end
    xlabel('N');
    ylabel('lambda min');
    % ylim([-0.2,0]);
    %legend(names,'location','southwest');
    legend(names,'location','best');

    subplot(2,2,2);
    for fcindx = [1]
        avminlambdavec = zeros(1,length(Nscan));
        sdminlambdavec = zeros(1,length(Nscan));    
        for Nindx = 1:length(Nscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avminlambdavec(Nindx) = mean(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
            sdminlambdavec(Nindx) = std(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
        end
        plot(log10(Nscan),log10(abs(avminlambdavec)),'x',...
            'DisplayName',strcat('max \gamma = ',num2str(maxfcscan(fcindx))));
        errorbar(log10(Nscan),log10(abs(avminlambdavec)),...
            log10(abs(avminlambdavec)) - log10(abs(avminlambdavec)-sdminlambdavec),...
            log10(abs(avminlambdavec)+sdminlambdavec) - log10(abs(avminlambdavec)),'x');

        hold on
    end
    xlabel('log10(N)');
    ylabel('log10(|lambda min|)');
    %ylim([-1.4,-0.6]);
    %legend(names,'location','southwest');


    subplot(2,2,3);
    names = cell(1,length(maxfcscan));
    storeindx = 1;
    for fcindx = 1:length(maxfcscan)
        avJtildeSDvec = zeros(1,length(Nscan));
        sdJtildeSDvec = zeros(1,length(Nscan));
        for Nindx = 1:length(Nscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avJtildeSDvec(Nindx) = mean(sqrt(Jtilde_varMat(Nindx,fcindx,reltrials)));
            sdJtildeSDvec(Nindx) = std(sqrt(Jtilde_varMat(Nindx,fcindx,reltrials)));
        end
        errorbar(Nscan,avJtildeSDvec,...
            sdJtildeSDvec,'x','LineStyle','none')
        hold on
        names{storeindx} = strcat('fc: ',num2str(maxfcscan(fcindx)));
        storeindx = storeindx+1;
    end
    xlabel('N');
    ylabel('\sigma Jtilde');
    legend(names,'location','best');

    subplot(2,2,4);
    hJtilde_pos = histogram(log10(JtildeMat(JtildeMat>0)),100,'Normalization','pdf');
    hold on
    hJtilde_neg = histogram(log10(abs(JtildeMat(JtildeMat<0))),100,'Normalization','pdf');
    xlabel('log10(Jtilde elements)');
    ylabel('P');
    legend('positive','negative');

end

%% Plot variation with maxfc
if length(maxfcscan)> 1
    figure; 
    subplot(2,2,1);
    names = cell(1,length(Nscan));
    storeindx = 1;
    for Nindx = [1]
        avminlambdavec = zeros(1,length(maxfcscan));
        sdminlambdavec = zeros(1,length(maxfcscan));    
        for fcindx = 1:length(maxfcscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avminlambdavec(fcindx) = mean(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
            sdminlambdavec(fcindx) = std(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
        end
        %plot(Nscan,avminlambdavec,'x',...
        %    'DisplayName',strcat('max \gamma = ',num2str(maxfcscan(gammaindx))));
        errorbar(maxfcscan,abs(avminlambdavec),...
            sdminlambdavec,'x','LineStyle','none')
        hold on
        names{storeindx} = strcat('N: ',num2str(Nscan(Nindx)));
        storeindx = storeindx+1;
    end
    xlabel('\Omega_{max}');
    ylabel('lambda {max}');
    % ylim([-0.2,0]);
    %legend(names,'location','southwest');
    legend(names,'location','best');

    subplot(2,2,2);
    for Nindx = [1]
        avminlambdavec = zeros(1,length(maxfcscan));
        sdminlambdavec = zeros(1,length(maxfcscan));    
        for fcindx = 1:length(maxfcscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avminlambdavec(fcindx) = mean(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
            sdminlambdavec(fcindx) = std(minlambdaMat_Jtilde(Nindx,fcindx,reltrials));
        end
        plot(log10(maxfcscan),log10(abs(avminlambdavec)),'x',...
            'DisplayName',strcat('N = ',num2str(Nscan(Nindx))));
        errorbar(log10(maxfcscan),log10(abs(avminlambdavec)),...
            log10(abs(avminlambdavec)) - log10(abs(avminlambdavec)-sdminlambdavec),...
            log10(abs(avminlambdavec)+sdminlambdavec) - log10(abs(avminlambdavec)),'x');

        hold on
    end
    xlabel('log10(\Omega_{max})');
    ylabel('log10(\lambda_{max})');
    %ylim([-1.4,-0.6]);
    %legend(names,'location','southwest');


    subplot(2,2,3);
    names = cell(1,length(Nscan));
    storeindx = 1;
    for Nindx = 1:length(Nscan)
        avJtildeSDvec = zeros(1,length(maxfcscan));
        sdJtildeSDvec = zeros(1,length(maxfcscan));
        for fcindx = 1:length(maxfcscan)      
            reltrials = find(exitflagmat(Nindx,fcindx,:)>0);
            avJtildeSDvec(fcindx) = mean(sqrt(Jtilde_varMat(Nindx,fcindx,reltrials)));
            sdJtildeSDvec(fcindx) = std(sqrt(Jtilde_varMat(Nindx,fcindx,reltrials)));
        end
        %plot(Nscan,avJtildeSDvec,'x',...
        %    'DisplayName',strcat('max \gamma = ',num2str(maxfcscan(gammaindx))));
        errorbar(maxfcscan,avJtildeSDvec,...
            sdJtildeSDvec,'x','LineStyle','none')
        hold on
        names{storeindx} = strcat('N: ',num2str(Nscan(Nindx)));
        storeindx = storeindx+1;
    end
    xlabel('\Omega_{max}');
    ylabel('\sigma Jtilde');
    legend(names,'location','best');

    subplot(2,2,4);
    hJtilde_pos = histogram(log10(JtildeMat(JtildeMat>0)),100,'Normalization','pdf');
    hold on
    hJtilde_neg = histogram(log10(abs(JtildeMat(JtildeMat<0))),100,'Normalization','pdf');
    xlabel('log10(Jtilde elements)');
    ylabel('P');
    legend('positive','negative');

end

%% Save data
% save('SIfig2part_randDAG_maxfc1pt5_pdown0pt5_fcdist1overx');
% save('FigNew2part_randDAG_maxfc1pt5_pdown0pt5');
% save('FigNew2part_randDAG_maxfc2_pdown0pt5');
% save('FigNew2part_randDAG_maxfc2');
% save('FigNew2part_randDAG_maxfc1pt5');

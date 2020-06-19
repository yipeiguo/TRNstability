% This is a script to analyze cluster data for eigenvalue scaling

close all; 
clear variables;

%% load data
pathname = './';
fn = 'EigScaling_fullrandom_phidist1_rho0pt01_n1_Nmin100_Nmax1000_fcmin2_fcmax2_numtrials10_fsolve1_ode1_run1';


load(strcat(pathname,fn));

%% Plot variation with N
trialdim = 3; methodOI = 1;
pvec = [0.25,0.75];
avmeanlambda_Jtilde = mean(meanlambdaMat_Jtilde(:,:,:,methodOI),trialdim);
avminlambda_Jtilde = mean(minlambdaMat_Jtilde(:,:,:,methodOI),trialdim);
sdminlambda_Jtilde = std(minlambdaMat_Jtilde(:,:,:,methodOI),0,trialdim);
ulbminlambda_Jtilde = zeros([size(meanlambdaMat_Jtilde(:,:,:,methodOI),1:2),2]);
for Nindx = 1:length(Nscan)
    for fcindx = 1:length(maxfcscan)
        ulbminlambda_Jtilde(Nindx,fcindx,:) = quantile(abs(minlambdaMat_Jtilde(Nindx,fcindx,:,methodOI)),pvec);
    end
end
avJtilde_sd = mean(sqrt(Jtilde_varMat),trialdim);
sdJtilde_sd = mean(sqrt(Jtilde_varMat),trialdim);

figure; 
subplot(2,2,1);
for fcindx = 1:length(maxfcscan)
%     plot(log10(Nscan),log10(abs(avminlambda_Jtilde(:,1,1,methodOI))),'x');
    % errorbar(Nscan,avminlambda_Jtilde(:,1,1,methodOI),...
    %     sdminlambda_Jtilde(:,1,1,methodOI),'x','LineStyle','none');
    errorbar(log10(Nscan),log10(abs(avminlambda_Jtilde(:,fcindx))),...
        log10(abs(avminlambda_Jtilde(:,fcindx)))-log10(ulbminlambda_Jtilde(:,fcindx,1)),...
        log10(ulbminlambda_Jtilde(:,fcindx,2))-log10(abs(avminlambda_Jtilde(:,fcindx))),...
        'x','LineStyle','none');
    hold on
end
xlabel('N');
ylabel('\lambda_{min}');
% xlim([0 1200]);
% ylim([-1 0]);

subplot(2,2,2);
plot(log10(Nscan),log10(avJtilde_sd(:,1,1,methodOI)),'x');
% errorbar(Nscan,avJtilde_sd(:,1,1,methodOI),...
%     sdJtilde_sd(:,1,1,methodOI),'x','LineStyle','none')    
xlabel('N');
ylabel('\sigma_{Jtilde}');
% xlim([0 1200]);
%ylim([0 15]);

subplot(2,2,4);
hJtilde_pos = histogram(log10(JtildeMat(JtildeMat>0)),100,'Normalization','pdf');
hold on
hJtilde_neg = histogram(log10(abs(JtildeMat(JtildeMat<0))),100,'Normalization','pdf');
xlabel('log10(Jtilde elements)');
ylabel('P');
legend('positive','negative');

%% Save data
% save('SIfig2part_fullyrandom_maxfc1pt5_pdown0pt5_fcdist1overx');
% save('FigNew2part_fullyrandom_maxfc2_pdown1');
% save('FigNew2part_fullyrandom_maxfc2');
% save('FigNew2part_fullyrandom_maxfc1pt2');

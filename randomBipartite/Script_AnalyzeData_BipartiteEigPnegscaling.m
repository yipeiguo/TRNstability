% This is a script for analyzing data from
% 'Script_EigPnegscaling_forcluster.m'

close all; 
clear variables;

%% load data
pathname = './';
fn = 'BipartitePnegData_Nmin200_Nmax200_k0pt10_rho0pt01_numSelfInts0_n1_minPneg0pt0_maxPneg1pt0_maxfc2pt0_numtrials100_numsets5_run1';

load(strcat(pathname,fn));
% load(strcat(pathname,fn),'exitflagmat','maxlambdaMat','Nscan');

%% Calculate eigenvalue using analytical expression (bipartite)
lambdaQ_Mat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets); % analytical expression 
for Nindx = 1:length(Nscan)
    for Pnegindx = 1:length(Pnegscan)
        for trialindx = 1:numtrials
            for setindx = 1:numsets
                Tmat = Tmat_cell{Nindx,Pnegindx,trialindx,setindx};
                c_ss = cssCell{Nindx,Pnegindx,trialindx,setindx};
                TtildeMat = Tmat./(k1.*c_ss(end))+eye(q);
                lambdaQ_Mat(Nindx,Pnegindx,trialindx,setindx) = sum(diag(TtildeMat));
            end
        end
    end
end               

lambdaM_Mat = max(lambdaQ_Mat,0);


%% Maximum eigenvalue
trialdim = 3; setdim = 4;
toaccept_all = exitflagmat>0;
maxlambdaMat_mean = median(lambdaM_Mat,trialdim);
maxlambdaMat_grandmean = nanmean(maxlambdaMat_mean,setdim);
maxlambdaMat_std = nanstd(maxlambdaMat_mean,[],setdim);

pvec = [0.25,0.75];
ulbmaxlambda = zeros(length(Pnegscan),2,length(Nscan));
for Nindx = 1:length(Nscan)
    for pnegindx = 1:length(Pnegscan)
        ulbmaxlambda(pnegindx,:,Nindx) = quantile(lambdaM_Mat(Nindx,pnegindx,:,:),pvec,'all');
    end
end

%% sd of jacobian matrix elements
Jtilde_varMat = zeros(length(Nscan),length(Pnegscan),numtrials,numsets);
for Nindx = 1:length(Nscan)
    for pnegindx = 1:length(Pnegscan)
        for trialindx = 1:numtrials
            for setindx = 1:numsets
                Tmat = Tmat_cell{Nindx,pnegindx,trialindx,setindx};
                css = cssCell{Nindx,pnegindx,trialindx,setindx};
                JtildeVec = Tmat(:)./(k1.*c_ss(end))+1;
                Jtilde_varMat(Nindx,pnegindx,trialindx,setindx) = var(JtildeVec);
            end
        end
    end
end
Jtilde_sdMat = sqrt(Jtilde_varMat);
meanJtildeSD = mean(Jtilde_sdMat,trialdim);
avJtildeSDvec = mean(meanJtildeSD,setdim);
sdJtildeSDvec = std(meanJtildeSD,[],setdim);


%% Plot

figure;
subplot(2,2,1);
errorbar(Pnegscan,maxlambdaMat_grandmean(1,:),maxlambdaMat_std(1,:),'x');
xlabel('Pneg');
ylabel('\lambda_{max}');

subplot(2,2,3);
errorbar(log10(Pnegscan),log10(maxlambdaMat_grandmean(1,:)),...
    log10(maxlambdaMat_grandmean(1,:)) - log10(maxlambdaMat_grandmean(1,:)-maxlambdaMat_std(1,:)),...
    log10(maxlambdaMat_grandmean(1,:)+maxlambdaMat_std(1,:)) - log10(maxlambdaMat_grandmean(1,:)),'x');
xlabel('log10(Pneg)');
ylabel('log10(\lambda_{max})');
% ylim([-16 -14])

subplot(2,2,2);
JtildeMat = JacobianMat./(k1*c_ss(end))+eye(N);       
% JtildeMat = JacobianMat./(-k1*c_ss(end))-eye(N);       
hJtilde_pos = histogram(log10(JtildeMat(JtildeMat>0)),100,'Normalization','pdf');
hold on
hJtilde_neg = histogram(log10(abs(JtildeMat(JtildeMat<0))),100,'Normalization','pdf');
xlabel('log10(Jtilde elements)');
ylabel('P');
legend('positive','negative');

subplot(2,2,4);
Ttilde = JtildeMat(1:q,1:q);
hJtilde_pos2 = histogram(log10(JtildeMat(JtildeMat>0)),100,'Normalization','pdf');
hold on
hJtilde_neg2 = histogram(log10(abs(JtildeMat(JtildeMat<0))),100,'Normalization','pdf');
hold on
hT_pos = histogram(log10(Ttilde(Ttilde>0)),100,'Normalization','pdf');
hold on
if min(Ttilde(:)<0)
    hT_neg = histogram(log10(abs(Ttilde(Ttilde<0))),100,'Normalization','pdf');
end
xlabel('log10(Jtilde elements)');
ylabel('P');
legend('positive','negative','Tpos','Tneg');

%% Visualize jacobian matrix
% JnegInds = find(JtildeMat<=0);
% JposInds = find(JtildeMat>=0);
% Jtoplot = log10(abs(JtildeMat));
% 
% figure;
% imagesc(Jtoplot)

%% Variation of variance of jacobian matrix elements with N
figure;
%plot(Nscan,avJtildeSDvec,'x',...
%    'DisplayName',strcat('max \gamma = ',num2str(Pnegscan(gammaindx))));
errorbar(Pnegscan,avJtildeSDvec(1,:),...
    sdJtildeSDvec(1,:),'x','LineStyle','none')
xlabel('Pneg');
ylabel('\sigma Jtilde');


%% save data
% save('FigBipartiteEigPnegscaling_maxfc1pt5_N500');



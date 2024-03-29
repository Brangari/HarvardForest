% % % % Scripts for 'Shifts in microbial temperature relationships under warming can mitigate heat-induced losses in carbon stocks (Paper: Brangarí et al, 2024, submitted to Nature Climate Change)
% % % % By Albert C Brangarí. University of Amsterdam, Institute for Biodiversity and Ecosystem Dynamics (IBED), 2023.
% % % % To run: download this file and HarvardForest_ExperimentalData.xls, and save them in the same folder. Then, run this Matlab script from there

close all; 
clearvars;
parentdir = pwd;

%% Read experimental data from file and generate variables for simulations

% Import sensor data (Field temperature)
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',4,'Range','E270:AS9030','Format','auto');
Temp = valsExt.Variables;
TimeH = linspace(1,length(Temp),length(Temp));
TimeD = TimeH/24;
TimeY = TimeD/365;

% Import fungal growth rates
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',1,'Range','M4:P45','Format','auto','VariableNamingRule','preserve');
depTg_vec = str2double(erase(valsExt.Properties.VariableNames,"°C"));
depF_vec = valsExt.Variables;
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',1,'Range','R5:U45');
vals = valsExt.Variables;
parF_a = vals(:,1);
parF_Tm = -vals(:,2)./vals(:,1);
fconF = vals(1,4);
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',1,'Range','C5:C45');
SOC = valsExt.Variables/100;

% Import bacterial growth rates
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',2,'Range','M4:P45','Format','auto','VariableNamingRule','preserve');
depB_vec = valsExt.Variables;
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',2,'Range','R5:U45');
vals = valsExt.Variables;
parB_a = vals(:,1);
parB_Tm = -vals(:,2)./vals(:,1);
fconB = vals(1,4);

% Import respiration rates
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',3,'Range','M4:P45','Format','auto','VariableNamingRule','preserve');
depTr_vec = str2double(erase(valsExt.Properties.VariableNames,"°C"));
depR_vec = valsExt.Variables;
valsExt = readtable('HarvardForest_ExperimentalData.xlsx','Sheet',3,'Range','R5:U45');
vals = valsExt.Variables;
parR_a = vals(:,1);
parR_Tm = -vals(:,2)./vals(:,1);
fconR = vals(1,4);

% Define T vector and functions
vecT = linspace(-10,40,100);
df_vec = vecT(2)-vecT(1);
funGsq = @(x_,P) max(P(1).*(x_-P(2)).*(1-exp(P(3).*(x_-P(4)))),0);
funRsq = @(x_,P) max(P(1)*(x_-P(2)),0);
normV25 = 25;

posMinC = [1:5 31:35];
posOrgC = [6:11 36:41];
posMinH = [12:15 21:24];
posOrgH = [16:20 25:30];

%% ANALYSE LAB CURVES - TEMP DEPENDENCES

for ii = 1:41
        funFsq_cal = @(x_,P) max(parF_a(ii).*(x_-parF_Tm(ii)).*(1-exp(P(1).*(x_-P(2)))),0);
    [regYF,parF(ii,:),~] = Reg_Ratkowsky(depTg_vec,depF_vec(ii,:),funFsq_cal,[0.2 50],[0 40],[0.2 53],vecT);
        funBsq_cal = @(x_,P) max(parB_a(ii).*(x_-parB_Tm(ii)).*(1-exp(P(1).*(x_-P(2)))),0);
    [regYB,parB(ii,:),~] = Reg_Ratkowsky(depTg_vec,depB_vec(ii,:),funBsq_cal,[0.2 50],[0 40],[0.2 53],vecT);
    FgrowthT(ii,:) = SOC(ii)*fconF*(regYF).^2;
    BgrowthT(ii,:) = SOC(ii)*fconB*(regYB).^2;
    regYR = max(funRsq(vecT,[parR_a(ii) parR_Tm(ii)]),0);
    TrespirT(ii,:) = SOC(ii)*fconR*(regYR).^2;
    CUET(ii,:) = (FgrowthT(ii,:)+BgrowthT(ii,:))./((FgrowthT(ii,:)+BgrowthT(ii,:))+TrespirT(ii,:));

    if any(posMinC == ii)
        ptreat = 1;
    elseif any(posOrgC == ii)
        ptreat = 2;
    elseif any(posMinH == ii)
        ptreat = 3;
    elseif any(posOrgH == ii)
        ptreat = 4;
    else
        ptreat = 0;
    end
    if any(1:4 == ptreat)
        norm1 = funGsq(normV25,[parF_a(ii) parF_Tm(ii) parF(ii,:)]); norm_depF_vec(ii,:) = depF_vec(ii,:)/norm1;
        Frate_gSOC(ii) = fconF*(norm1).^2;
        Frate_gds(ii) = Frate_gSOC(ii)*SOC(ii);
        norm2 = funGsq(normV25,[parB_a(ii) parB_Tm(ii) parB(ii,:)]); norm_depB_vec(ii,:) = depB_vec(ii,:)/norm2;
        Brate_gSOC(ii) = fconB*(norm2).^2;
        Brate_gds(ii) = Brate_gSOC(ii)*SOC(ii);
        norm3 = funRsq(normV25,[parR_a(ii) parR_Tm(ii)]); norm_depR_vec(ii,:) = depR_vec(ii,:)/norm3;
        Rrate_gSOC(ii) = fconR*(norm3).^2;
        Rrate_gds(ii) = Rrate_gSOC(ii)*SOC(ii);
    end
end

minTempC = min(Temp(:,[posOrgC posMinC]),[],2)';
maxTempC = max(Temp(:,[posOrgC posMinC]),[],2)';
meanTempC = mean(Temp(:,[posOrgC posMinC]),2);
minTempH = min(Temp(:,[posOrgH posMinH]),[],2)';
maxTempH = max(Temp(:,[posOrgH posMinH]),[],2)';
meanTempH = mean(Temp(:,[posOrgH posMinH]),2);

norm_depF_vec_treat(1,:) = mean(norm_depF_vec(posMinC,:));
norm_depF_vec_treat(2,:) = mean(norm_depF_vec(posOrgC,:));
norm_depF_vec_treat(3,:) = mean(norm_depF_vec(posMinH,:));
norm_depF_vec_treat(4,:) = mean(norm_depF_vec(posOrgH,:));
se_norm_depF_vec_treat(1,:) = std(norm_depF_vec(posMinC,:))/sqrt(length(posMinC));
se_norm_depF_vec_treat(2,:) = std(norm_depF_vec(posOrgC,:))/sqrt(length(posOrgC));
se_norm_depF_vec_treat(3,:) = std(norm_depF_vec(posMinH,:))/sqrt(length(posMinH));
se_norm_depF_vec_treat(4,:) = std(norm_depF_vec(posOrgH,:))/sqrt(length(posOrgH));
norm_depB_vec_treat(1,:) = mean(norm_depB_vec(posMinC,:));
norm_depB_vec_treat(2,:) = mean(norm_depB_vec(posOrgC,:));
norm_depB_vec_treat(3,:) = mean(norm_depB_vec(posMinH,:));
norm_depB_vec_treat(4,:) = mean(norm_depB_vec(posOrgH,:));
se_norm_depB_vec_treat(1,:) = std(norm_depB_vec(posMinC,:))/sqrt(length(posMinC));
se_norm_depB_vec_treat(2,:) = std(norm_depB_vec(posOrgC,:))/sqrt(length(posOrgC));
se_norm_depB_vec_treat(3,:) = std(norm_depB_vec(posMinH,:))/sqrt(length(posMinH));
se_norm_depB_vec_treat(4,:) = std(norm_depB_vec(posOrgH,:))/sqrt(length(posOrgH));
norm_depR_vec_treat(1,:) = mean(norm_depR_vec(posMinC,:));
norm_depR_vec_treat(2,:) = mean(norm_depR_vec(posOrgC,:));
norm_depR_vec_treat(3,:) = mean(norm_depR_vec(posMinH,:));
norm_depR_vec_treat(4,:) = mean(norm_depR_vec(posOrgH,:));
se_norm_depR_vec_treat(1,:) = std(norm_depR_vec(posMinC,:))/sqrt(length(posMinC));
se_norm_depR_vec_treat(2,:) = std(norm_depR_vec(posOrgC,:))/sqrt(length(posOrgC));
se_norm_depR_vec_treat(3,:) = std(norm_depR_vec(posMinH,:))/sqrt(length(posMinH));
se_norm_depR_vec_treat(4,:) = std(norm_depR_vec(posOrgH,:))/sqrt(length(posOrgH));
CUET_treat(1,:) = mean(CUET(posMinC,:));
CUET_treat(2,:) = mean(CUET(posOrgC,:));
CUET_treat(3,:) = mean(CUET(posMinH,:));
CUET_treat(4,:) = mean(CUET(posOrgH,:));

Frate_gSOC_treat(1,:) = mean(Frate_gSOC(posMinC));
Frate_gSOC_treat(2,:) = mean(Frate_gSOC(posOrgC));
Frate_gSOC_treat(3,:) = mean(Frate_gSOC(posMinH));
Frate_gSOC_treat(4,:) = mean(Frate_gSOC(posOrgH));
se_Frate_gSOC_treat(1,:) = std(Frate_gSOC(posMinC))/sqrt(length(posMinC));
se_Frate_gSOC_treat(2,:) = std(Frate_gSOC(posOrgC))/sqrt(length(posOrgC));
se_Frate_gSOC_treat(3,:) = std(Frate_gSOC(posMinH))/sqrt(length(posMinH));
se_Frate_gSOC_treat(4,:) = std(Frate_gSOC(posOrgH))/sqrt(length(posOrgH));
Brate_gSOC_treat(1,:) = mean(Brate_gSOC(posMinC));
Brate_gSOC_treat(2,:) = mean(Brate_gSOC(posOrgC));
Brate_gSOC_treat(3,:) = mean(Brate_gSOC(posMinH));
Brate_gSOC_treat(4,:) = mean(Brate_gSOC(posOrgH));
se_Brate_gSOC_treat(1,:) = std(Brate_gSOC(posMinC))/sqrt(length(posMinC));
se_Brate_gSOC_treat(2,:) = std(Brate_gSOC(posOrgC))/sqrt(length(posOrgC));
se_Brate_gSOC_treat(3,:) = std(Brate_gSOC(posMinH))/sqrt(length(posMinH));
se_Brate_gSOC_treat(4,:) = std(Brate_gSOC(posOrgH))/sqrt(length(posOrgH));
Rrate_gSOC_treat(1,:) = mean(Rrate_gSOC(posMinC));
Rrate_gSOC_treat(2,:) = mean(Rrate_gSOC(posOrgC));
Rrate_gSOC_treat(3,:) = mean(Rrate_gSOC(posMinH));
Rrate_gSOC_treat(4,:) = mean(Rrate_gSOC(posOrgH));
se_Rrate_gSOC_treat(1,:) = std(Rrate_gSOC(posMinC))/sqrt(length(posMinC));
se_Rrate_gSOC_treat(2,:) = std(Rrate_gSOC(posOrgC))/sqrt(length(posOrgC));
se_Rrate_gSOC_treat(3,:) = std(Rrate_gSOC(posMinH))/sqrt(length(posMinH));
se_Rrate_gSOC_treat(4,:) = std(Rrate_gSOC(posOrgH))/sqrt(length(posOrgH));

Frate_gds_treat(1,:) = mean(Frate_gds(posMinC));
Frate_gds_treat(2,:) = mean(Frate_gds(posOrgC));
Frate_gds_treat(3,:) = mean(Frate_gds(posMinH));
Frate_gds_treat(4,:) = mean(Frate_gds(posOrgH));
se_Frate_gds_treat(1,:) = std(Frate_gds(posMinC))/sqrt(length(posMinC));
se_Frate_gds_treat(2,:) = std(Frate_gds(posOrgC))/sqrt(length(posOrgC));
se_Frate_gds_treat(3,:) = std(Frate_gds(posMinH))/sqrt(length(posMinH));
se_Frate_gds_treat(4,:) = std(Frate_gds(posOrgH))/sqrt(length(posOrgH));
Brate_gds_treat(1,:) = mean(Brate_gds(posMinC));
Brate_gds_treat(2,:) = mean(Brate_gds(posOrgC));
Brate_gds_treat(3,:) = mean(Brate_gds(posMinH));
Brate_gds_treat(4,:) = mean(Brate_gds(posOrgH));
se_Brate_gds_treat(1,:) = std(Brate_gds(posMinC))/sqrt(length(posMinC));
se_Brate_gds_treat(2,:) = std(Brate_gds(posOrgC))/sqrt(length(posOrgC));
se_Brate_gds_treat(3,:) = std(Brate_gds(posMinH))/sqrt(length(posMinH));
se_Brate_gds_treat(4,:) = std(Brate_gds(posOrgH))/sqrt(length(posOrgH));
Rrate_gds_treat(1,:) = mean(Rrate_gds(posMinC));
Rrate_gds_treat(2,:) = mean(Rrate_gds(posOrgC));
Rrate_gds_treat(3,:) = mean(Rrate_gds(posMinH));
Rrate_gds_treat(4,:) = mean(Rrate_gds(posOrgH));
se_Rrate_gds_treat(1,:) = std(Rrate_gds(posMinC))/sqrt(length(posMinC));
se_Rrate_gds_treat(2,:) = std(Rrate_gds(posOrgC))/sqrt(length(posOrgC));
se_Rrate_gds_treat(3,:) = std(Rrate_gds(posMinH))/sqrt(length(posMinH));
se_Rrate_gds_treat(4,:) = std(Rrate_gds(posOrgH))/sqrt(length(posOrgH));

for ii = 1:4
    funLin_cal = @(x_,P) P(1).*(x_-P(2));

    [~,norm_parF_All(ii,1:2),~] = Reg_Ratkowsky(depTg_vec(1:3),norm_depF_vec_treat(ii,1:3),funLin_cal,[1/30 0],[0 -20],[1 10],vecT);
    funFsq_cal = @(x_,P) norm_parF_All(ii,1).*(x_-norm_parF_All(ii,2)).*(1-exp(P(1).*(x_-P(2))));
    [norm_regYF(ii,:),norm_parF_All(ii,3:4),SqEf] = Reg_Ratkowsky(depTg_vec,norm_depF_vec_treat(ii,:),funFsq_cal,[0.1 45],[0 40],[2 50],vecT);

    [~,norm_parB_All(ii,1:2),~] = Reg_Ratkowsky(depTg_vec(1:3),norm_depB_vec_treat(ii,1:3),funLin_cal,[1/30 0],[0 -20],[1 10],vecT);
    funBsq_cal = @(x_,P) norm_parB_All(ii,1).*(x_-norm_parB_All(ii,2)).*(1-exp(P(1).*(x_-P(2))));
    [norm_regYB(ii,:),norm_parB_All(ii,3:4),SqEb] = Reg_Ratkowsky(depTg_vec,norm_depB_vec_treat(ii,:),funBsq_cal,[0.1 45],[0 40],[0.6 50],vecT);

    [norm_regYR(ii,:),norm_parR_All(ii,:),SqEr] = Reg_Ratkowsky(depTr_vec(1:3),norm_depR_vec_treat(ii,1:3),funLin_cal,[1/30 0],[0 -20],[1 10],vecT); 
end

% Define plotting properties
sizeF = 14;
sizeFm = 14;
ngroups = 2;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    xposC(:,i) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
end
colorMinC = [0.3010 0.7450 0.9330];
colorOrgC = [0 0.4470 0.7410];
colorMinH = [0.8500 0.3250 0.0980];
colorOrgH = [0.6350 0.0780 0.1840];
colorOrgT_p5 = [0.6 0 0.6];
colorMinT_p5 = [1 0.6 0.8];
colorOrgT_m5 = [0.5 0 1];
colorMinT_m5 = [0.6 0.6 1];

% plot Figure S1
figure(10); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    patch('XData',[TimeD fliplr(TimeD)],'YData',[minTempC fliplr(maxTempC)],'FaceColor',colorOrgC,'FaceAlpha',.3,'EdgeColor','none'); hold on;
    patch('XData',[TimeD fliplr(TimeD)],'YData',[minTempH fliplr(maxTempH)],'FaceColor',colorOrgH,'FaceAlpha',.3,'EdgeColor','none');
    h(1)=plot(TimeD,meanTempC,'Color',colorOrgC,'LineWidth',2);
    h(2)=plot(TimeD,meanTempH,'Color',colorOrgH,'LineWidth',2);
    xticks(linspace(0,TimeD(end),13)); set(gca,'Box','on'); xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}); set(gca,'fontsize',sizeF); ylabel('Temperature [°C]'); axis([0 TimeD(end) -2 27]); legend(h,{'Control','Warming'});

% plot Figure 1
figure(1); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    subplot(3,3,[1 4]); 
        plot(depTg_vec,norm_depB_vec_treat(2,:),'o','MarkerEdgeColor',colorOrgC); hold on; errorbar(depTg_vec,norm_depB_vec_treat(2,:),se_norm_depB_vec_treat(2,:),'Color',colorOrgC,'linestyle','none'); plot(vecT,norm_regYB(2,:),'Color',colorOrgC,'LineWidth',2);
        plot(depTg_vec,norm_depB_vec_treat(4,:),'o','MarkerEdgeColor',colorOrgH); hold on; errorbar(depTg_vec,norm_depB_vec_treat(4,:),se_norm_depB_vec_treat(4,:),'Color',colorOrgH,'linestyle','none'); plot(vecT,norm_regYB(4,:),'Color',colorOrgH,'LineWidth',2); 
        plot(depTg_vec,norm_depB_vec_treat(1,:),'o','MarkerEdgeColor',colorMinC); hold on; errorbar(depTg_vec,norm_depB_vec_treat(1,:),se_norm_depB_vec_treat(1,:),'Color',colorMinC,'linestyle','none'); plot(vecT,norm_regYB(1,:),'Color',colorMinC,'LineWidth',2);
        plot(depTg_vec,norm_depB_vec_treat(3,:),'o','MarkerEdgeColor',colorMinH); hold on; errorbar(depTg_vec,norm_depB_vec_treat(3,:),se_norm_depB_vec_treat(3,:),'Color',colorMinH,'linestyle','none'); plot(vecT,norm_regYB(3,:),'Color',colorMinH,'LineWidth',2);
        b(1) = bar(1,NaN,'FaceColor',colorOrgC);
        b(2) = bar(2,NaN,'FaceColor',colorOrgH);
        b(3) = bar(3,NaN,'FaceColor',colorMinC);
        b(4) = bar(4,NaN,'FaceColor',colorMinH);
        set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylabel('Norm \surd Bacterial growth'); axis([-10 40 0 1.5]); legend(b,{'Control ORG','Warming ORG','Control MIN','Warming MIN',},'Orientation','horizontal');
    subplot(3,3,[2 5]); 
        plot(depTg_vec,norm_depF_vec_treat(2,:),'o','MarkerEdgeColor',colorOrgC); hold on; errorbar(depTg_vec,norm_depF_vec_treat(2,:),se_norm_depF_vec_treat(2,:),'Color',colorOrgC,'linestyle','none'); plot(vecT,norm_regYF(2,:),'Color',colorOrgC,'LineWidth',2); 
        plot(depTg_vec,norm_depF_vec_treat(4,:),'o','MarkerEdgeColor',colorOrgH); hold on; errorbar(depTg_vec,norm_depF_vec_treat(4,:),se_norm_depF_vec_treat(4,:),'Color',colorOrgH,'linestyle','none'); plot(vecT,norm_regYF(4,:),'Color',colorOrgH,'LineWidth',2); 
        plot(depTg_vec,norm_depF_vec_treat(1,:),'o','MarkerEdgeColor',colorMinC); hold on; errorbar(depTg_vec,norm_depF_vec_treat(1,:),se_norm_depF_vec_treat(1,:),'Color',colorMinC,'linestyle','none'); plot(vecT,norm_regYF(1,:),'Color',colorMinC,'LineWidth',2);
        plot(depTg_vec,norm_depF_vec_treat(3,:),'o','MarkerEdgeColor',colorMinH); hold on; errorbar(depTg_vec,norm_depF_vec_treat(3,:),se_norm_depF_vec_treat(3,:),'Color',colorMinH,'linestyle','none'); plot(vecT,norm_regYF(3,:),'Color',colorMinH,'LineWidth',2);
        set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylabel('Norm \surd Fungal growth'); axis([-10 40 0 1.5]);
   subplot(3,3,[3 6]); 
        plot(depTr_vec,norm_depR_vec_treat(2,:),'o','MarkerEdgeColor',colorOrgC); hold on; errorbar(depTr_vec,norm_depR_vec_treat(2,:),se_norm_depR_vec_treat(2,:),'Color',colorOrgC,'linestyle','none'); plot(vecT,norm_regYR(2,:),'Color',colorOrgC,'LineWidth',2);
        plot(depTr_vec,norm_depR_vec_treat(4,:),'o','MarkerEdgeColor',colorOrgH); hold on; errorbar(depTr_vec,norm_depR_vec_treat(4,:),se_norm_depR_vec_treat(4,:),'Color',colorOrgH,'linestyle','none'); plot(vecT,norm_regYR(4,:),'Color',colorOrgH,'LineWidth',2);
        plot(depTr_vec,norm_depR_vec_treat(1,:),'o','MarkerEdgeColor',colorMinC); hold on; errorbar(depTr_vec,norm_depR_vec_treat(1,:),se_norm_depR_vec_treat(1,:),'Color',colorMinC,'linestyle','none'); plot(vecT,norm_regYR(1,:),'Color',colorMinC,'LineWidth',2);
        plot(depTr_vec,norm_depR_vec_treat(3,:),'o','MarkerEdgeColor',colorMinH); hold on; errorbar(depTr_vec,norm_depR_vec_treat(3,:),se_norm_depR_vec_treat(3,:),'Color',colorMinH,'linestyle','none'); plot(vecT,norm_regYR(3,:),'Color',colorMinH,'LineWidth',2);
        set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylabel('Norm \surd Respiration'); axis([-10 40 0 2]); 
    subplot(3,3,7)
        plot(xposC(1,1),Brate_gSOC_treat(2),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgC);  hold on;
        plot(xposC(1,2),Brate_gSOC_treat(4),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgH);
        plot(xposC(2,1),Brate_gSOC_treat(1),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinC);
        plot(xposC(2,2),Brate_gSOC_treat(3),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinH);
        errorbar([xposC(1,:) xposC(2,:)],Brate_gSOC_treat([2 4 1 3]),se_Brate_gSOC_treat([2 4 1 3]),'k','linestyle','none');
        set(gca,'fontsize',sizeF); ylabel({'Bacterial growth';'[\mugC/gSOC/h]'}); set(gca,'xtick',[]); set(gca,'xticklabels',[]); ylim([0 6]); xlim([0.5 2.5]); yticks([0 2 4 6]);
    subplot(3,3,8)
        plot(xposC(1,1),Frate_gSOC_treat(2),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgC);  hold on;
        plot(xposC(1,2),Frate_gSOC_treat(4),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgH);
        plot(xposC(2,1),Frate_gSOC_treat(1),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinC);
        plot(xposC(2,2),Frate_gSOC_treat(3),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinH);
        errorbar([xposC(1,:) xposC(2,:)],Frate_gSOC_treat([2 4 1 3]),se_Frate_gSOC_treat([2 4 1 3]),'k','linestyle','none');
        set(gca,'fontsize',sizeF); ylabel({'Fungal growth';'[\mugC/gSOC/h]'}); set(gca,'xtick',[]); set(gca,'xticklabels',[]); ylim([0 1]); xlim([0.5 2.5]);
    subplot(3,3,9)
        plot(xposC(1,1),Rrate_gSOC_treat(2),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgC);  hold on;
        plot(xposC(1,2),Rrate_gSOC_treat(4),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorOrgH);
        plot(xposC(2,1),Rrate_gSOC_treat(1),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinC);
        plot(xposC(2,2),Rrate_gSOC_treat(3),'LineStyle','none','Marker','o','MarkerSize',20,'MarkerFaceColor',colorMinH);
        errorbar([xposC(1,:) xposC(2,:)],Rrate_gSOC_treat([2 4 1 3]),se_Rrate_gSOC_treat([2 4 1 3]),'k','linestyle','none');
        set(gca,'fontsize',sizeF); ylabel({'Respiration';'[\mugC/gSOC/h]'}); set(gca,'xtick',[]); set(gca,'xticklabels',[]); ylim([0 45]); xlim([0.5 2.5]);

%% UPSCALING FIELD SIMULATIONS

parF_All = [parF_a parF_Tm parF];
parB_All = [parB_a parB_Tm parB];
parR_All = [parR_a parR_Tm];

for ii = 1:41
    Fgrowth_gSOC(ii,:) = fconF*(funGsq(Temp(:,ii),parF_All(ii,:))).^2;
    Bgrowth_gSOC(ii,:) = fconB*(funGsq(Temp(:,ii),parB_All(ii,:))).^2;
    Trespir_gSOC(ii,:) = fconR*(funRsq(Temp(:,ii),parR_All(ii,:))).^2;
    Fgrowth_gds(ii,:) = Fgrowth_gSOC(ii,:)*SOC(ii);
    Bgrowth_gds(ii,:) = Bgrowth_gSOC(ii,:)*SOC(ii);
    Mgrowth_gds(ii,:) = (Fgrowth_gSOC(ii,:)+Bgrowth_gSOC(ii,:))*SOC(ii);
    Trespir_gds(ii,:) = Trespir_gSOC(ii,:)*SOC(ii);

    yrFgrowth_gSOC(ii) = trapz(Fgrowth_gSOC(ii,:))/TimeY(end);
    yrBgrowth_gSOC(ii) = trapz(Bgrowth_gSOC(ii,:))/TimeY(end);
    yrTrespir_gSOC(ii) = trapz(Trespir_gSOC(ii,:))/TimeY(end);
    yrFgrowth_gds(ii) = yrFgrowth_gSOC(ii)*SOC(ii);
    yrBgrowth_gds(ii) = yrBgrowth_gSOC(ii)*SOC(ii);
    yrMgrowth_gds(ii) = yrFgrowth_gds(ii) + yrBgrowth_gds(ii);
    yrTrespir_gds(ii) = yrTrespir_gSOC(ii)*SOC(ii);
    yrMgrowth_gSOC(ii) = yrMgrowth_gds(ii)/SOC(ii);

    CUE(ii) = (yrFgrowth_gSOC(ii)+yrBgrowth_gSOC(ii))/(yrFgrowth_gSOC(ii)+yrBgrowth_gSOC(ii)+yrTrespir_gSOC(ii));
end

SOC_treat(1) = mean(SOC(posMinC));
SOC_treat(2) = mean(SOC(posOrgC));
SOC_treat(3) = mean(SOC(posMinH));
SOC_treat(4) = mean(SOC(posOrgH));
                      
Fgrowth_gds_treat(1,:) = mean(Fgrowth_gds(posMinC,:));
Bgrowth_gds_treat(1,:) = mean(Bgrowth_gds(posMinC,:));
Mgrowth_gds_treat(1,:) = mean(Mgrowth_gds(posMinC,:));
Trespir_gds_treat(1,:) = mean(Trespir_gds(posMinC,:));
yrFgrowth_gds_treat(1) = mean(yrFgrowth_gds(posMinC)); se_yrFgrowth_gds_treat(1) = std(yrFgrowth_gds(posMinC))/sqrt(length(posMinC));
yrBgrowth_gds_treat(1) = mean(yrBgrowth_gds(posMinC)); se_yrBgrowth_gds_treat(1) = std(yrBgrowth_gds(posMinC))/sqrt(length(posMinC));
yrMgrowth_gds_treat(1) = mean(yrFgrowth_gds(posMinC)+yrBgrowth_gds(posMinC)); se_yrMgrowth_gds_treat(1) = std(yrFgrowth_gds(posMinC)+yrBgrowth_gds(posMinC))/sqrt(length(posMinC));
yrTrespir_gds_treat(1) = mean(yrTrespir_gds(posMinC)); se_yrTrespir_gds_treat(1) = std(yrTrespir_gds(posMinC))/sqrt(length(posMinC));
CUE_treat(1) = mean(CUE(posMinC)); seCUE_treat(1) = std(CUE(posMinC))/sqrt(length(posMinC));
Fgrowth_gds_treat(2,:) = mean(Fgrowth_gds(posOrgC,:));
Bgrowth_gds_treat(2,:) = mean(Bgrowth_gds(posOrgC,:));
Mgrowth_gds_treat(2,:) = mean(Mgrowth_gds(posOrgC,:));
Trespir_gds_treat(2,:) = mean(Trespir_gds(posOrgC,:));
yrFgrowth_gds_treat(2) = mean(yrFgrowth_gds(posOrgC)); se_yrFgrowth_gds_treat(2) = std(yrFgrowth_gds(posOrgC))/sqrt(length(posOrgC));
yrBgrowth_gds_treat(2) = mean(yrBgrowth_gds(posOrgC)); se_yrBgrowth_gds_treat(2) = std(yrBgrowth_gds(posOrgC))/sqrt(length(posOrgC));
yrMgrowth_gds_treat(2) = mean(yrFgrowth_gds(posOrgC)+yrBgrowth_gds(posOrgC)); se_yrMgrowth_gds_treat(2) = std(yrFgrowth_gds(posOrgC)+yrBgrowth_gds(posOrgC))/sqrt(length(posOrgC));
yrTrespir_gds_treat(2) = mean(yrTrespir_gds(posOrgC)); se_yrTrespir_gds_treat(2) = std(yrTrespir_gds(posOrgC))/sqrt(length(posOrgC));
CUE_treat(2) = mean(CUE(posOrgC)); seCUE_treat(2) = std(CUE(posOrgC))/sqrt(length(posOrgC));
Fgrowth_gds_treat(3,:) = mean(Fgrowth_gds(posMinH,:));
Bgrowth_gds_treat(3,:) = mean(Bgrowth_gds(posMinH,:));
Mgrowth_gds_treat(3,:) = mean(Mgrowth_gds(posMinH,:));
Trespir_gds_treat(3,:) = mean(Trespir_gds(posMinH,:));
yrFgrowth_gds_treat(3) = mean(yrFgrowth_gds(posMinH)); se_yrFgrowth_gds_treat(3) = std(yrFgrowth_gds(posMinH))/sqrt(length(posMinH));
yrBgrowth_gds_treat(3) = mean(yrBgrowth_gds(posMinH)); se_yrBgrowth_gds_treat(3) = std(yrBgrowth_gds(posMinH))/sqrt(length(posMinH));
yrMgrowth_gds_treat(3) = mean(yrFgrowth_gds(posMinH)+yrBgrowth_gds(posMinH)); se_yrMgrowth_gds_treat(3) = std(yrFgrowth_gds(posMinH)+yrBgrowth_gds(posMinH))/sqrt(length(posMinH));
yrTrespir_gds_treat(3) = mean(yrTrespir_gds(posMinH)); se_yrTrespir_gds_treat(3) = std(yrTrespir_gds(posMinH))/sqrt(length(posMinH));
CUE_treat(3) = mean(CUE(posMinH)); seCUE_treat(3) = std(CUE(posMinH))/sqrt(length(posMinH));
Fgrowth_gds_treat(4,:) = mean(Fgrowth_gds(posOrgH,:));
Bgrowth_gds_treat(4,:) = mean(Bgrowth_gds(posOrgH,:));
Mgrowth_gds_treat(4,:) = mean(Mgrowth_gds(posOrgH,:));
Trespir_gds_treat(4,:) = mean(Trespir_gds(posOrgH,:));
yrFgrowth_gds_treat(4) = mean(yrFgrowth_gds(posOrgH)); se_yrFgrowth_gds_treat(4) = std(yrFgrowth_gds(posOrgH))/sqrt(length(posOrgH));
yrBgrowth_gds_treat(4) = mean(yrBgrowth_gds(posOrgH)); se_yrBgrowth_gds_treat(4) = std(yrBgrowth_gds(posOrgH))/sqrt(length(posOrgH));
yrMgrowth_gds_treat(4) = mean(yrFgrowth_gds(posOrgH)+yrBgrowth_gds(posOrgH)); se_yrMgrowth_gds_treat(4) = std(yrFgrowth_gds(posOrgH)+yrBgrowth_gds(posOrgH))/sqrt(length(posOrgH));
yrTrespir_gds_treat(4) = mean(yrTrespir_gds(posOrgH)); se_yrTrespir_gds_treat(4) = std(yrTrespir_gds(posOrgH))/sqrt(length(posOrgH));
CUE_treat(4) = mean(CUE(posOrgH)); seCUE_treat(4) = std(CUE(posOrgH))/sqrt(length(posOrgH));

densMinC = 0.54; %soil sample characteristics
densOrgC = 0.22;
densMinH = 0.54;
densOrgH = 0.27;
thickMinC = 6.1;
thickOrgC = 3.9;
thickMinH = 6.6;
thickOrgH = 3.4;
xMinC = thickMinC/(thickOrgC+thickMinC);       %fraction core
xOrgC = thickOrgC/(thickOrgC+thickMinC); 
xMinH = thickMinH/(thickOrgH+thickMinH); 
xOrgH = thickOrgH/(thickOrgH+thickMinH);

convF = 10000/1000000;      % to units of gC/m2/y
yrFgrowth_m2_treat = yrFgrowth_gds_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF; % total C in soil (depth factored in)
se_yrFgrowth_m2_treat(1) = std(yrFgrowth_gds(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrFgrowth_m2_treat(2) = std(yrFgrowth_gds(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrFgrowth_m2_treat(3) = std(yrFgrowth_gds(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrFgrowth_m2_treat(4) = std(yrFgrowth_gds(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrBgrowth_m2_treat = yrBgrowth_gds_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrBgrowth_m2_treat(1) = std(yrBgrowth_gds(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrBgrowth_m2_treat(2) = std(yrBgrowth_gds(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrBgrowth_m2_treat(3) = std(yrBgrowth_gds(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrBgrowth_m2_treat(4) = std(yrBgrowth_gds(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrMgrowth_m2_treat = (yrFgrowth_gds_treat+yrBgrowth_gds_treat).*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrMgrowth_m2_treat(1) = std((yrFgrowth_gds(posMinC)+yrBgrowth_gds(posMinC))*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrMgrowth_m2_treat(2) = std((yrFgrowth_gds(posOrgC)+yrBgrowth_gds(posOrgC))*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrMgrowth_m2_treat(3) = std((yrFgrowth_gds(posMinH)+yrBgrowth_gds(posMinH))*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrMgrowth_m2_treat(4) = std((yrFgrowth_gds(posOrgH)+yrBgrowth_gds(posOrgH))*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrTrespir_m2_treat = yrTrespir_gds_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrTrespir_m2_treat(1) = std(yrTrespir_gds(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrTrespir_m2_treat(2) = std(yrTrespir_gds(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrTrespir_m2_treat(3) = std(yrTrespir_gds(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrTrespir_m2_treat(4) = std(yrTrespir_gds(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));

Fgrowth_gSOC_treat(1,:) = Fgrowth_gds_treat(1,:)/SOC_treat(1);
Bgrowth_gSOC_treat(1,:) = Bgrowth_gds_treat(1,:)/SOC_treat(1);
Mgrowth_gSOC_treat(1,:) = Mgrowth_gds_treat(1,:)/SOC_treat(1);
Trespir_gSOC_treat(1,:) = Trespir_gds_treat(1,:)/SOC_treat(1);
yrFgrowth_gSOC_treat(1) = yrFgrowth_gds_treat(1)/SOC_treat(1);
    se_yrFgrowth_gSOC_treat(1) = se_yrFgrowth_gds_treat(1)/SOC_treat(1);
yrBgrowth_gSOC_treat(1) = yrBgrowth_gds_treat(1)/SOC_treat(1);
    se_yrBgrowth_gSOC_treat(1) = se_yrBgrowth_gds_treat(1)/SOC_treat(1);
yrMgrowth_gSOC_treat(1) = yrMgrowth_gds_treat(1)/SOC_treat(1);
    se_yrMgrowth_gSOC_treat(1) = se_yrMgrowth_gds_treat(1)/SOC_treat(1);
yrTrespir_gSOC_treat(1) = yrTrespir_gds_treat(1)/SOC_treat(1);
    se_yrTrespir_gSOC_treat(1) = se_yrTrespir_gds_treat(1)/SOC_treat(1);
Fgrowth_gSOC_treat(2,:) = Fgrowth_gds_treat(2,:)/SOC_treat(2);
Bgrowth_gSOC_treat(2,:) = Bgrowth_gds_treat(2,:)/SOC_treat(2);
Mgrowth_gSOC_treat(2,:) = Mgrowth_gds_treat(2,:)/SOC_treat(2);
Trespir_gSOC_treat(2,:) = Trespir_gds_treat(2,:)/SOC_treat(2);
yrFgrowth_gSOC_treat(2) = yrFgrowth_gds_treat(2)/SOC_treat(2);
    se_yrFgrowth_gSOC_treat(2) = se_yrFgrowth_gds_treat(2)/SOC_treat(2);
yrBgrowth_gSOC_treat(2) = yrBgrowth_gds_treat(2)/SOC_treat(2);
    se_yrBgrowth_gSOC_treat(2) = se_yrBgrowth_gds_treat(2)/SOC_treat(2);
yrMgrowth_gSOC_treat(2) = yrMgrowth_gds_treat(2)/SOC_treat(2);
    se_yrMgrowth_gSOC_treat(2) = se_yrMgrowth_gds_treat(2)/SOC_treat(2);
yrTrespir_gSOC_treat(2) = yrTrespir_gds_treat(2)/SOC_treat(2);
    se_yrTrespir_gSOC_treat(2) = se_yrTrespir_gds_treat(2)/SOC_treat(2);
Fgrowth_gSOC_treat(3,:) = Fgrowth_gds_treat(3,:)/SOC_treat(3);
Bgrowth_gSOC_treat(3,:) = Bgrowth_gds_treat(3,:)/SOC_treat(3);
Mgrowth_gSOC_treat(3,:) = Mgrowth_gds_treat(3,:)/SOC_treat(3);
Trespir_gSOC_treat(3,:) = Trespir_gds_treat(3,:)/SOC_treat(3);
yrFgrowth_gSOC_treat(3) = yrFgrowth_gds_treat(3)/SOC_treat(3);
    se_yrFgrowth_gSOC_treat(3) = se_yrFgrowth_gds_treat(3)/SOC_treat(3);
yrBgrowth_gSOC_treat(3) = yrBgrowth_gds_treat(3)/SOC_treat(3);
    se_yrBgrowth_gSOC_treat(3) = se_yrBgrowth_gds_treat(3)/SOC_treat(3);
yrMgrowth_gSOC_treat(3) = yrMgrowth_gds_treat(3)/SOC_treat(3);
    se_yrMgrowth_gSOC_treat(3) = se_yrMgrowth_gds_treat(3)/SOC_treat(3);
yrTrespir_gSOC_treat(3) = yrTrespir_gds_treat(3)/SOC_treat(3);
    se_yrTrespir_gSOC_treat(3) = se_yrTrespir_gds_treat(3)/SOC_treat(3);
Fgrowth_gSOC_treat(4,:) = Fgrowth_gds_treat(4,:)/SOC_treat(4);
Bgrowth_gSOC_treat(4,:) = Bgrowth_gds_treat(4,:)/SOC_treat(4);
Mgrowth_gSOC_treat(4,:) = Mgrowth_gds_treat(4,:)/SOC_treat(4);
Trespir_gSOC_treat(4,:) = Trespir_gds_treat(4,:)/SOC_treat(4);
yrFgrowth_gSOC_treat(4) = yrFgrowth_gds_treat(4)/SOC_treat(4);
    se_yrFgrowth_gSOC_treat(4) = se_yrFgrowth_gds_treat(4)/SOC_treat(4);
yrBgrowth_gSOC_treat(4) = yrBgrowth_gds_treat(4)/SOC_treat(4);
    se_yrBgrowth_gSOC_treat(4) = se_yrBgrowth_gds_treat(4)/SOC_treat(4);
yrMgrowth_gSOC_treat(4) = yrMgrowth_gds_treat(4)/SOC_treat(4);
    se_yrMgrowth_gSOC_treat(4) = se_yrMgrowth_gds_treat(4)/SOC_treat(4);
yrTrespir_gSOC_treat(4) = yrTrespir_gds_treat(4)/SOC_treat(4);
    se_yrTrespir_gSOC_treat(4) = se_yrTrespir_gds_treat(4)/SOC_treat(4);

% plot Figure 3
figure(3); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    subplot(2,2,1); plot(TimeD,Fgrowth_gSOC_treat(2,:)+Bgrowth_gSOC_treat(2,:),'Color',colorOrgC,'LineWidth',2); hold on;
                      plot(TimeD,Fgrowth_gSOC_treat(4,:)+Bgrowth_gSOC_treat(4,:),'Color',colorOrgH,'LineWidth',2);
                      plot(TimeD,Fgrowth_gSOC_treat(1,:)+Bgrowth_gSOC_treat(1,:),'Color',colorMinC,'LineWidth',2);
                      plot(TimeD,Fgrowth_gSOC_treat(3,:)+Bgrowth_gSOC_treat(3,:),'Color',colorMinH,'LineWidth',2); 
                      set(gca,'fontsize',sizeF);  ylabel({'Microbial growth';'[\mugC/gSOC/h]'}); xticks(linspace(0,TimeD(end),13)); set(gca,'Box','on'); xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}); xlim([0 max(TimeD)]);
                      b(1) = bar(1,NaN,'FaceColor',colorOrgC);
                      b(2) = bar(2,NaN,'FaceColor',colorOrgH);
                      b(3) = bar(3,NaN,'FaceColor',colorMinC);
                      b(4) = bar(4,NaN,'FaceColor',colorMinH);
                      legend(b,{'Control ORG','Warming ORG','Control MIN','Warming MIN',},'Orientation','horizontal');
    subplot(2,2,2); plot(TimeD,Trespir_gSOC_treat(2,:),'Color',colorOrgC,'LineWidth',2); hold on;
                      plot(TimeD,Trespir_gSOC_treat(4,:),'Color',colorOrgH,'LineWidth',2);
                      plot(TimeD,Trespir_gSOC_treat(1,:),'Color',colorMinC,'LineWidth',2);
                      plot(TimeD,Trespir_gSOC_treat(3,:),'Color',colorMinH,'LineWidth',2);
                      set(gca,'fontsize',sizeF); xticks(linspace(0,TimeD(end),13)); set(gca,'Box','on'); xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}); ylabel({'Microbial respiration';'[\mugC/gSOC/h]'}); xlim([0 max(TimeD)]);         

yrMgrowth_m2_treat_x(1) = yrMgrowth_m2_treat(1)+yrMgrowth_m2_treat(2);
yrMgrowth_m2_treat_x(2) = yrMgrowth_m2_treat(3)+yrMgrowth_m2_treat(4);
se_yrMgrowth_m2_treat_x(1) = (se_yrMgrowth_m2_treat(1)*sqrt(length(posMinC)) + se_yrMgrowth_m2_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
se_yrMgrowth_m2_treat_x(2) = (se_yrMgrowth_m2_treat(3)*sqrt(length(posMinH)) + se_yrMgrowth_m2_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));
yrTrespir_m2_treat_x(1) = yrTrespir_m2_treat(1) + yrTrespir_m2_treat(2);
yrTrespir_m2_treat_x(2) = yrTrespir_m2_treat(3) + yrTrespir_m2_treat(4);
se_yrTrespir_m2_treat_x(1) = (se_yrTrespir_m2_treat(1)*sqrt(length(posMinC)) + se_yrTrespir_m2_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
se_yrTrespir_m2_treat_x(2) = (se_yrTrespir_m2_treat(3)*sqrt(length(posMinH)) + se_yrTrespir_m2_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));
CUE_treat_x(1) = xMinC*CUE_treat(1) + xOrgC*CUE_treat(2);
CUE_treat_x(2) = xMinH*CUE_treat(3) + xOrgH*CUE_treat(4);
seCUE_treat_x(1) = (xMinC*seCUE_treat(1)*sqrt(length(posMinC)) + xOrgC*seCUE_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
seCUE_treat_x(2) = (xMinH*seCUE_treat(3)*sqrt(length(posMinH)) + xOrgH*seCUE_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));

CUET_treat_x(1,:) = xMinC*CUET_treat(1,:) + xOrgC*CUET_treat(2,:);
CUET_treat_x(2,:) = xMinH*CUET_treat(3,:) + xOrgH*CUET_treat(4,:);

for ii = 1:41
    meanG(ii) = mean(CUET(ii,29:89));
    stdG(ii) = std(CUET(ii,29:89));
    auxpos = find(CUET(ii,60:100)<meanG(ii)-stdG(ii),1);
    if isempty(auxpos)
        posG(ii) = 40;
    else
        posG(ii) = auxpos;
    end
    posCUE(ii) = posG(ii)+59;
    tempCUE(ii) = vecT(posCUE(ii));
end
mean_tempCUE(1) = mean(tempCUE([posMinC posOrgC]));
mean_tempCUE(2) = mean(tempCUE([posMinH posOrgH]));
se_tempCUE(1) = std(tempCUE([posMinC posOrgC]))/sqrt(length([posMinC posOrgC]));
se_tempCUE(2) = std(tempCUE([posMinH posOrgH]))/sqrt(length([posMinH posOrgH]));
mean_CUEpl(1) = interp1(vecT,CUET_treat_x(1,:)*100,mean_tempCUE(1));
mean_CUEpl(2) = interp1(vecT,CUET_treat_x(2,:)*100,mean_tempCUE(2));

% plot Figure 2
figure(2);
    set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    plot(vecT,CUET_treat_x(1,:)*100,'Color',colorOrgC,'LineStyle',':','LineWidth',2); hold on;
    h(1)=plot(vecT(21:97),CUET_treat_x(1,21:97)*100,'Color',colorOrgC,'LineWidth',2);
    plot(mean_tempCUE(1),mean_CUEpl(1),'Marker','o','MarkerEdgeColor',colorOrgC,'MarkerSize',5,'LineWidth',2);
    errorbar(mean_tempCUE(1),mean_CUEpl(1),se_tempCUE(1),'horizontal','k','linestyle','none');
    
    plot(vecT,CUET_treat_x(2,:)*100,'Color',colorOrgH,'LineStyle',':','LineWidth',2); 
    h(2)=plot(vecT(21:97),CUET_treat_x(2,21:97)*100,'Color',colorOrgH,'LineWidth',2);
    plot(mean_tempCUE(2),mean_CUEpl(2),'Marker','o','MarkerEdgeColor',colorOrgH,'MarkerSize',5,'LineWidth',2);
    errorbar(mean_tempCUE(2),mean_CUEpl(2),se_tempCUE(2),'horizontal','k','linestyle','none');
    set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylabel('CUE [%]'); axis([-10 40 0 30]); legend(h,{'Control','Warming'});

block = [1 2 3 4 5 6 7 8 9 12 14 16 17 18 19 20 21 23];
minpos = [5 13 24 22 12 35 23 15 32 2 4 31 21 14 3 1 33 34];
orgpos = [7 16 26 25 17 36 27 20 39 8 6 38 29 19 9 10 37 41];
thickorg_i = [3.38 3.25 3.80 3.73 2.58 3.75 2.58 3.13 3.63 4.90 3.50 4.88 3.33 2.98 2.94 4.06 4.00 4.00];
thickmin_i = 10-thickorg_i;
densorg_i = [0.22 0.26 0.15 0.28 0.22 0.22 0.30 0.32 0.22 0.16 0.21 0.09 0.24 0.26 0.21 0.28 0.19 0.21];
densmin_i = [0.62 0.50 0.43 0.73 0.53 0.40 0.52 0.53 0.45 0.45 0.54 0.84 0.63 0.43 0.66 0.50 0.39 0.54];

for ii = 1:length(block)
    BgrowthT_block(ii,:) = BgrowthT(minpos(ii),:)*thickmin_i(ii)*densmin_i(ii) + BgrowthT(orgpos(ii),:)*thickorg_i(ii)*densorg_i(ii);
    FgrowthT_block(ii,:) = FgrowthT(minpos(ii),:)*thickmin_i(ii)*densmin_i(ii) + FgrowthT(orgpos(ii),:)*thickorg_i(ii)*densorg_i(ii);
    TrespirT_block(ii,:) = TrespirT(minpos(ii),:)*thickmin_i(ii)*densmin_i(ii) + TrespirT(orgpos(ii),:)*thickorg_i(ii)*densorg_i(ii);
    CUET_block(ii,:) = (BgrowthT_block(ii,:)+FgrowthT_block(ii,:))./(BgrowthT_block(ii,:)+FgrowthT_block(ii,:)+TrespirT_block(ii,:));
       
    meanG_block(ii) = mean(CUET_block(ii,29:89));
    stdG_block(ii) = std(CUET_block(ii,29:89));
    auxpos = find(CUET_block(ii,60:100)<meanG_block(ii)-stdG_block(ii),1);
    if isempty(auxpos)
        posG_block(ii) = 40;
    else
        posG_block(ii) = auxpos;
    end
    posCUE_block(ii) = posG_block(ii)+59;
    tempCUE_block(ii) = vecT(posCUE_block(ii));
end

%% TRANSPLANTATION

Temp_transp(:,[posOrgH posMinH]) = Temp(:,[posOrgH posMinH])-5;
Temp_transp(:,[posOrgC posMinC]) = Temp(:,[posOrgC posMinC])+5;

for ii = 1:41
    Fgrowth_gSOC_transp(ii,:) = fconF*(funGsq(Temp_transp(:,ii),parF_All(ii,:))).^2;
    Bgrowth_gSOC_transp(ii,:) = fconB*(funGsq(Temp_transp(:,ii),parB_All(ii,:))).^2;
    Trespir_gSOC_transp(ii,:) = fconR*(funRsq(Temp_transp(:,ii),parR_All(ii,:))).^2;
    Fgrowth_gds_transp(ii,:) = Fgrowth_gSOC_transp(ii,:)*SOC(ii);
    Bgrowth_gds_transp(ii,:) = Bgrowth_gSOC_transp(ii,:)*SOC(ii);
    Trespir_gds_transp(ii,:) = Trespir_gSOC_transp(ii,:)*SOC(ii);

    yrFgrowth_gSOC_transp(ii) = trapz(Fgrowth_gSOC_transp(ii,:))/TimeY(end);
    yrBgrowth_gSOC_transp(ii) = trapz(Bgrowth_gSOC_transp(ii,:))/TimeY(end);
    yrTrespir_gSOC_transp(ii) = trapz(Trespir_gSOC_transp(ii,:))/TimeY(end);
    yrFgrowth_gds_transp(ii) = yrFgrowth_gSOC_transp(ii)*SOC(ii);
    yrBgrowth_gds_transp(ii) = yrBgrowth_gSOC_transp(ii)*SOC(ii);
    yrTrespir_gds_transp(ii) = yrTrespir_gSOC_transp(ii)*SOC(ii);

    CUE_transp(ii) = (yrFgrowth_gSOC_transp(ii)+yrBgrowth_gSOC_transp(ii))/(yrFgrowth_gSOC_transp(ii)+yrBgrowth_gSOC_transp(ii)+yrTrespir_gSOC_transp(ii));
end
               
Fgrowth_gds_transp_treat(1,:) = mean(Fgrowth_gds_transp(posMinC,:));
Bgrowth_gds_transp_treat(1,:) = mean(Bgrowth_gds_transp(posMinC,:));
Trespir_gds_transp_treat(1,:) = mean(Trespir_gds_transp(posMinC,:));
yrFgrowth_gds_transp_treat(1) = mean(yrFgrowth_gds_transp(posMinC)); se_yrFgrowth_gds_transp_treat(1) = std(yrFgrowth_gds_transp(posMinC))/sqrt(length(posMinC));
yrBgrowth_gds_transp_treat(1) = mean(yrBgrowth_gds_transp(posMinC)); se_yrBgrowth_gds_transp_treat(1) = std(yrBgrowth_gds_transp(posMinC))/sqrt(length(posMinC));
yrMgrowth_gds_transp_treat(1) = mean(yrFgrowth_gds_transp(posMinC)+yrBgrowth_gds_transp(posMinC)); se_yrMgrowth_gds_transp_treat(1) = std(yrFgrowth_gds_transp(posMinC)+yrBgrowth_gds_transp(posMinC))/sqrt(length(posMinC));
yrTrespir_gds_transp_treat(1) = mean(yrTrespir_gds_transp(posMinC)); se_yrTrespir_gds_transp_treat(1) = std(yrTrespir_gds_transp(posMinC))/sqrt(length(posMinC));
CUE_transp_treat(1) = mean(CUE_transp(posMinC)); seCUE_transp_treat(1) = std(CUE_transp(posMinC))/sqrt(length(posMinC));
Fgrowth_gds_transp_treat(2,:) = mean(Fgrowth_gds_transp(posOrgC,:));
Bgrowth_gds_transp_treat(2,:) = mean(Bgrowth_gds_transp(posOrgC,:));
Trespir_gds_transp_treat(2,:) = mean(Trespir_gds_transp(posOrgC,:));
yrFgrowth_gds_transp_treat(2) = mean(yrFgrowth_gds_transp(posOrgC)); se_yrFgrowth_gds_transp_treat(2) = std(yrFgrowth_gds_transp(posOrgC))/sqrt(length(posOrgC));
yrBgrowth_gds_transp_treat(2) = mean(yrBgrowth_gds_transp(posOrgC)); se_yrBgrowth_gds_transp_treat(2) = std(yrBgrowth_gds_transp(posOrgC))/sqrt(length(posOrgC));
yrMgrowth_gds_transp_treat(2) = mean(yrFgrowth_gds_transp(posOrgC)+yrBgrowth_gds_transp(posOrgC)); se_yrMgrowth_gds_transp_treat(2) = std(yrFgrowth_gds_transp(posOrgC)+yrBgrowth_gds_transp(posOrgC))/sqrt(length(posOrgC));
yrTrespir_gds_transp_treat(2) = mean(yrTrespir_gds_transp(posOrgC)); se_yrTrespir_gds_transp_treat(2) = std(yrTrespir_gds_transp(posOrgC))/sqrt(length(posOrgC));
CUE_transp_treat(2) = mean(CUE_transp(posOrgC)); seCUE_transp_treat(2) = std(CUE_transp(posOrgC))/sqrt(length(posOrgC));
Fgrowth_gds_transp_treat(3,:) = mean(Fgrowth_gds_transp(posMinH,:));
Bgrowth_gds_transp_treat(3,:) = mean(Bgrowth_gds_transp(posMinH,:));
Trespir_gds_transp_treat(3,:) = mean(Trespir_gds_transp(posMinH,:));
yrFgrowth_gds_transp_treat(3) = mean(yrFgrowth_gds_transp(posMinH)); se_yrFgrowth_gds_transp_treat(3) = std(yrFgrowth_gds_transp(posMinH))/sqrt(length(posMinH));
yrBgrowth_gds_transp_treat(3) = mean(yrBgrowth_gds_transp(posMinH)); se_yrBgrowth_gds_transp_treat(3) = std(yrBgrowth_gds_transp(posMinH))/sqrt(length(posMinH));
yrMgrowth_gds_transp_treat(3) = mean(yrFgrowth_gds_transp(posMinH)+yrBgrowth_gds_transp(posMinH)); se_yrMgrowth_gds_transp_treat(3) = std(yrFgrowth_gds_transp(posMinH)+yrBgrowth_gds_transp(posMinH))/sqrt(length(posMinH));
yrTrespir_gds_transp_treat(3) = mean(yrTrespir_gds_transp(posMinH)); se_yrTrespir_gds_transp_treat(3) = std(yrTrespir_gds_transp(posMinH))/sqrt(length(posMinH));
CUE_transp_treat(3) = mean(CUE_transp(posMinH)); seCUE_transp_treat(3) = std(CUE_transp(posMinH))/sqrt(length(posMinH));
Fgrowth_gds_transp_treat(4,:) = mean(Fgrowth_gds_transp(posOrgH,:));
Bgrowth_gds_transp_treat(4,:) = mean(Bgrowth_gds_transp(posOrgH,:));
Trespir_gds_transp_treat(4,:) = mean(Trespir_gds_transp(posOrgH,:));
yrFgrowth_gds_transp_treat(4) = mean(yrFgrowth_gds_transp(posOrgH)); se_yrFgrowth_gds_transp_treat(4) = std(yrFgrowth_gds_transp(posOrgH))/sqrt(length(posOrgH));
yrBgrowth_gds_transp_treat(4) = mean(yrBgrowth_gds_transp(posOrgH)); se_yrBgrowth_gds_transp_treat(4) = std(yrBgrowth_gds_transp(posOrgH))/sqrt(length(posOrgH));
yrMgrowth_gds_transp_treat(4) = mean(yrFgrowth_gds_transp(posOrgH)+yrBgrowth_gds_transp(posOrgH)); se_yrMgrowth_gds_transp_treat(4) = std(yrFgrowth_gds_transp(posOrgH)+yrBgrowth_gds_transp(posOrgH))/sqrt(length(posOrgH));
yrTrespir_gds_transp_treat(4) = mean(yrTrespir_gds_transp(posOrgH)); se_yrTrespir_gds_transp_treat(4) = std(yrTrespir_gds_transp(posOrgH))/sqrt(length(posOrgH));
CUE_transp_treat(4) = mean(CUE_transp(posOrgH)); seCUE_transp_treat(4) = std(CUE_transp(posOrgH))/sqrt(length(posOrgH));

yrFgrowth_m2_transp_treat = yrFgrowth_gds_transp_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrFgrowth_m2_transp_treat(1) = std(yrFgrowth_gds_transp(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrFgrowth_m2_transp_treat(2) = std(yrFgrowth_gds_transp(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrFgrowth_m2_transp_treat(3) = std(yrFgrowth_gds_transp(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrFgrowth_m2_transp_treat(4) = std(yrFgrowth_gds_transp(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrBgrowth_m2_transp_treat = yrBgrowth_gds_transp_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrBgrowth_m2_transp_treat(1) = std(yrBgrowth_gds_transp(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrBgrowth_m2_transp_treat(2) = std(yrBgrowth_gds_transp(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrBgrowth_m2_transp_treat(3) = std(yrBgrowth_gds_transp(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrBgrowth_m2_transp_treat(4) = std(yrBgrowth_gds_transp(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrMgrowth_m2_transp_treat = (yrFgrowth_gds_transp_treat+yrBgrowth_gds_transp_treat).*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrMgrowth_m2_transp_treat(1) = std((yrFgrowth_gds_transp(posMinC)+yrBgrowth_gds_transp(posMinC))*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrMgrowth_m2_transp_treat(2) = std((yrFgrowth_gds_transp(posOrgC)+yrBgrowth_gds_transp(posOrgC))*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrMgrowth_m2_transp_treat(3) = std((yrFgrowth_gds_transp(posMinH)+yrBgrowth_gds_transp(posMinH))*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrMgrowth_m2_transp_treat(4) = std((yrFgrowth_gds_transp(posOrgH)+yrBgrowth_gds_transp(posOrgH))*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));
yrTrespir_m2_transp_treat = yrTrespir_gds_transp_treat.*[densMinC densOrgC densMinH densOrgH].*[thickMinC thickOrgC thickMinH thickOrgH]*convF;
se_yrTrespir_m2_transp_treat(1) = std(yrTrespir_gds_transp(posMinC)*densMinC*thickMinC*convF)/sqrt(length(posMinC));
se_yrTrespir_m2_transp_treat(2) = std(yrTrespir_gds_transp(posOrgC)*densOrgC*thickOrgC*convF)/sqrt(length(posOrgC));
se_yrTrespir_m2_transp_treat(3) = std(yrTrespir_gds_transp(posMinH)*densMinH*thickMinH*convF)/sqrt(length(posMinH));
se_yrTrespir_m2_transp_treat(4) = std(yrTrespir_gds_transp(posOrgH)*densOrgH*thickOrgH*convF)/sqrt(length(posOrgH));

yrMgrowth_m2_transp_treat_x(1) = yrMgrowth_m2_transp_treat(1)+yrMgrowth_m2_transp_treat(2);
yrMgrowth_m2_transp_treat_x(2) = yrMgrowth_m2_transp_treat(3)+yrMgrowth_m2_transp_treat(4);
se_yrMgrowth_m2_transp_treat_x(1) = (se_yrMgrowth_m2_transp_treat(1)*sqrt(length(posMinC)) + se_yrMgrowth_m2_transp_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
se_yrMgrowth_m2_transp_treat_x(2) = (se_yrMgrowth_m2_transp_treat(3)*sqrt(length(posMinH)) + se_yrMgrowth_m2_transp_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));
yrTrespir_m2_transp_treat_x(1) = yrTrespir_m2_transp_treat(1) + yrTrespir_m2_transp_treat(2);
yrTrespir_m2_transp_treat_x(2) = yrTrespir_m2_transp_treat(3) + yrTrespir_m2_transp_treat(4);
se_yrTrespir_m2_transp_treat_x(1) = (se_yrTrespir_m2_transp_treat(1)*sqrt(length(posMinC)) + se_yrTrespir_m2_transp_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
se_yrTrespir_m2_transp_treat_x(2) = (se_yrTrespir_m2_transp_treat(3)*sqrt(length(posMinH)) + se_yrTrespir_m2_transp_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));
CUE_transp_treat_x(1) = xMinC*CUE_transp_treat(1) + xOrgC*CUE_transp_treat(2);
CUE_transp_treat_x(2) = xMinH*CUE_transp_treat(3) + xOrgH*CUE_transp_treat(4);
seCUE_transp_treat_x(1) = (xMinC*seCUE_transp_treat(1)*sqrt(length(posMinC)) + xOrgC*seCUE_transp_treat(2)*sqrt(length(posOrgC)))/sqrt(length([posMinC posOrgC]));
seCUE_transp_treat_x(2) = (xMinH*seCUE_transp_treat(3)*sqrt(length(posMinH)) + xOrgH*seCUE_transp_treat(4)*sqrt(length(posOrgH)))/sqrt(length([posMinH posOrgH]));

% plot Figure 4
figure(4); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    subplot(2,3,1); b = bar(yrMgrowth_m2_treat_x,'FaceColor','flat'); hold on; 
                    b.CData(1,:) = colorOrgC; b.CData(2,:) = colorOrgH;
                    errorbar([1 2],yrMgrowth_m2_treat_x,se_yrMgrowth_m2_treat_x,'k','linestyle','none');
                    bar(1,yrMgrowth_m2_transp_treat_x(1),'FaceAlpha',0,'EdgeColor',colorOrgT_p5,'LineWidth',2,'LineStyle','--'); 
                    errorbar(1,yrMgrowth_m2_transp_treat_x(1),se_yrMgrowth_m2_transp_treat_x(1),'k','linestyle','none');
                    c = bar((yrFgrowth_gds_treat([1 3])+yrBgrowth_gds_treat([1 3])).*[densMinC densMinH].*[thickMinC thickMinH]*convF,'FaceColor','flat'); 
                    c.CData(1,:) = colorMinC; c.CData(2,:) = colorMinH;
                    ylabel({'Total microbial growth';'[gC/m^2/y]'}); set(gca,'fontsize',sizeF); set(gca,'xtick',[]); set(gca,'xticklabels',[]); 
                    d(1) = bar(1,NaN,'FaceColor',colorOrgC);
                    d(2) = bar(2,NaN,'FaceColor',colorOrgH);
                    d(3) = bar(1,NaN,'FaceColor',colorMinC);
                    d(4) = bar(2,NaN,'FaceColor',colorMinH);
                    legend(d,{'Control ORG','Warming ORG','Control MIN','Warming MIN',},'Orientation','horizontal');
    subplot(2,3,2); b = bar(yrTrespir_m2_treat_x,'FaceColor','flat'); hold on;
                    b.CData(1,:) = colorOrgC; b.CData(2,:) = colorOrgH;
                    c = bar(yrTrespir_gds_treat([1 3]).*[densMinC densMinH].*[thickMinC thickMinH]*convF,'FaceColor','flat'); 
                    c.CData(1,:) = colorMinC; c.CData(2,:) = colorMinH;
                    errorbar([1 2],yrTrespir_m2_treat_x,se_yrTrespir_m2_treat_x,'k','linestyle','none');
                    d = bar(1,yrTrespir_m2_transp_treat_x(1),'FaceAlpha',0,'EdgeColor',colorOrgT_p5,'LineWidth',2,'LineStyle','--'); 
                    errorbar(1,yrTrespir_m2_transp_treat_x(1),se_yrTrespir_m2_transp_treat_x(1),'k','linestyle','none');
                    ylabel({'Total respiration';'[gC/m^2/y]'}); set(gca,'fontsize',sizeF); set(gca,'xtick',[]); set(gca,'xticklabels',[]); 
                    legend(d,{'Control +5°C',},'Orientation','horizontal');
    subplot(2,3,3); b = bar(CUE_treat_x*100,'FaceColor','flat'); hold on;  
                    b.CData(1,:) = colorOrgC; b.CData(2,:) = colorOrgH;
                    errorbar([0.95 2],CUE_treat_x*100,seCUE_treat_x*100,'k','linestyle','none');
                    bar(1,CUE_transp_treat_x(1)*100,'FaceAlpha',0,'EdgeColor',colorOrgT_p5,'LineWidth',2,'LineStyle','--'); 
                    errorbar(1.1,CUE_transp_treat_x(1)*100,seCUE_transp_treat_x(1)*100,'k','linestyle','none');
                    ylabel('CUE [%]'); set(gca,'fontsize',sizeF); set(gca,'xtick',[]); set(gca,'xticklabels',[]);

%% SOC balance
xx = 0:1:9;
xxL = 0:0.1:9;
yy = [100 100 - cumsum([linspace(629,411,80) 411*ones(1,10)]-412)/10/5879*100];
yyM = [100 100 - cumsum(629*ones(1,9)-412)/5879*100];
yym1 = [100 100 - cumsum(411*ones(1,9)-412)/5879*100];
yym2 = [100 100 - cumsum(412*ones(1,9)-412)/5879*100];

% plot Figure 5
figure(5);
    plot(xxL,yy,'Color','g','LineWidth',4); hold on;
    plot(xx,yyM,'--m','LineWidth',2);
    plot(xx,yym1,'--r','LineWidth',2);
    plot(xx,yym2,'-.b','LineWidth',2);
    xlabel('Time [y]');
    ylabel('SOC [%]');
    set(gca,'fontsize',14);  set(gcf,'color','w'); ylim([0 105]); xlim([0 9]); 
    set(gca,'XTick',[1 2 3 4 5 6 7 8 9]); set(gca,'YTick',[0 25 50 75 100]);

beep

%% FUNCTIONS CALLED

function [regY,P,Rsq] = Reg_Ratkowsky(t,y,eq,P0,minP,maxP,regT)

    opts = optimset('MaxIter',100000,'MaxFunEval',100000,'TolFun',1e-20,'TolX',1e-20,'TolConSQP',1e-20,'Display','off'); % solver options definition
    errfh = @(P,x,z) sum((z(:)-eq(x(:),P)).^2);

%     [P,err,a,b] = fminsearch(errfh,P0,opts,t,y);
    [P,err,a,b] = fmincon(errfh,P0,[],[],[],[],minP,maxP,[],opts,t,y);
    
    Rsq = 1-errfh(P,t,y)/sum((y-mean(y)).^2);

    regY = eq(regT,P);

end

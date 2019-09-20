clc; clear all; dbstop if error; close all;

n_exp = 4;
n_rep = 10;
% experiments = [0, 1];

% GHEXIN
psnr_ghexin = zeros(n_rep, n_exp);
nrs_ghexin = zeros(n_rep, n_exp);
silhouette_ghexin = zeros(n_rep, n_exp);
dbindex_ghexin = zeros(n_rep, n_exp);
ET_ghexin = zeros(n_rep, n_exp);

% DGSOT
psnr_dgsot = zeros(n_rep, n_exp);
nrs_dgsot = zeros(n_rep, n_exp);
silhouette_dgsot = zeros(n_rep, n_exp);
dbindex_dgsot = zeros(n_rep, n_exp);
ET_dgsot = zeros(n_rep, n_exp);

%GHNG
psnr_ghng = zeros(n_rep, n_exp);
nrs_ghng = zeros(n_rep, n_exp);
silhouette_ghng = zeros(n_rep, n_exp);
dbindex_ghng = zeros(n_rep, n_exp);
ET_ghng = zeros(n_rep, n_exp);

j = 1;
for i = 1:n_rep
    rng(i)
    [psnr_ghexin(i, j), nrs_ghexin(i, j), silhouette_ghexin(i, j), dbindex_ghexin(i, j), ET_ghexin(i, j), ...
        psnr_dgsot(i, j), nrs_dgsot(i, j), silhouette_dgsot(i, j), dbindex_dgsot(i, j), ET_dgsot(i, j), ...
        psnr_ghng(i, j), nrs_ghng(i, j), silhouette_ghng(i, j), dbindex_ghng(i, j), ET_ghng(i, j)] = test_models_X(i);
    close all
end

j = 2;
for i = 1:n_rep
    rng(i)
    [psnr_ghexin(i, j), nrs_ghexin(i, j), silhouette_ghexin(i, j), dbindex_ghexin(i, j), ET_ghexin(i, j), ...
        psnr_dgsot(i, j), nrs_dgsot(i, j), silhouette_dgsot(i, j), dbindex_dgsot(i, j), ET_dgsot(i, j), ...
        psnr_ghng(i, j), nrs_ghng(i, j), silhouette_ghng(i, j), dbindex_ghng(i, j), ET_ghng(i, j)] = test_models_density(i);
    close all
end

j = 3;
for i = 1:n_rep
    rng(i)
    [psnr_ghexin(i, j), nrs_ghexin(i, j), silhouette_ghexin(i, j), dbindex_ghexin(i, j), ET_ghexin(i, j), ...
        psnr_dgsot(i, j), nrs_dgsot(i, j), silhouette_dgsot(i, j), dbindex_dgsot(i, j), ET_dgsot(i, j), ...
        psnr_ghng(i, j), nrs_ghng(i, j), silhouette_ghng(i, j), dbindex_ghng(i, j), ET_ghng(i, j)] = test_models_gaussian(i);
    close all
end

j = 4;
for i = 1:n_rep
    rng(i)
    [psnr_ghexin(i, j), nrs_ghexin(i, j), silhouette_ghexin(i, j), dbindex_ghexin(i, j), ET_ghexin(i, j), ...
        psnr_dgsot(i, j), nrs_dgsot(i, j), silhouette_dgsot(i, j), dbindex_dgsot(i, j), ET_dgsot(i, j), ...
        psnr_ghng(i, j), nrs_ghng(i, j), silhouette_ghng(i, j), dbindex_ghng(i, j), ET_ghng(i, j)] = test_models_videos(i,3);
    close all
end

% PSNR
psnr_ghexin_avg = mean(psnr_ghexin);
psnr_ghexin_sem = std(psnr_ghexin)/sqrt(n_rep);
psnr_dgsot_avg = mean(psnr_dgsot);
psnr_dgsot_sem = std(psnr_dgsot)/sqrt(n_rep);
psnr_ghng_avg = mean(psnr_ghng);
psnr_ghng_sem = std(psnr_ghng)/sqrt(n_rep);

% NRS
nrs_ghexin_avg = mean(nrs_ghexin);
nrs_ghexin_sem = std(nrs_ghexin)/sqrt(n_rep);
nrs_dgsot_avg = mean(nrs_dgsot);
nrs_dgsot_sem = std(nrs_dgsot)/sqrt(n_rep);
nrs_ghng_avg = mean(nrs_ghng);
nrs_ghng_sem = std(nrs_ghng)/sqrt(n_rep);

% SILHOUETTE
silhouette_ghexin_avg = mean(silhouette_ghexin);
silhouette_ghexin_sem = std(silhouette_ghexin)/sqrt(n_rep);
silhouette_dgsot_avg = mean(silhouette_dgsot);
silhouette_dgsot_sem = std(silhouette_dgsot)/sqrt(n_rep);
silhouette_ghng_avg = mean(silhouette_ghng);
silhouette_ghng_sem = std(silhouette_ghng)/sqrt(n_rep);

% DBINDEX
dbindex_ghexin_avg = mean(dbindex_ghexin);
dbindex_ghexin_sem = std(dbindex_ghexin)/sqrt(n_rep);
dbindex_dgsot_avg = mean(dbindex_dgsot);
dbindex_dgsot_sem = std(dbindex_dgsot)/sqrt(n_rep);
dbindex_ghng_avg = mean(dbindex_ghng);
dbindex_ghng_sem = std(dbindex_ghng)/sqrt(n_rep);

% ET
ET_ghexin_avg = mean(ET_ghexin);
ET_ghexin_sem = std(ET_ghexin)/sqrt(n_rep);
ET_dgsot_avg = mean(ET_dgsot);
ET_dgsot_sem = std(ET_dgsot)/sqrt(n_rep);
ET_ghng_avg = mean(ET_ghng);
ET_ghng_sem = std(ET_ghng)/sqrt(n_rep);

ngroups = 4;
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = zeros(ngroups, nbars);
for i = 1:nbars
    x(:, i) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
end
x=x';

% PSNR
f = figure();
hold on
bar(1:n_exp, [psnr_ghexin_avg; psnr_dgsot_avg; psnr_ghng_avg]', 'EdgeColor', 'none')
errorbar(x(1,:), psnr_ghexin_avg, psnr_ghexin_sem,'k.')
errorbar(x(2,:), psnr_dgsot_avg, psnr_dgsot_sem,'k.')
errorbar(x(3,:), psnr_ghng_avg, psnr_ghng_sem,'k.')
hold off
xticks(x(2,:))
xticklabels({'X letter','square','Gaussians','videos'})
legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3)
legend('boxoff')
% xlabel('experiments')
ylabel('PSNR index')
box off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_psnr','-dpdf') % then print it
print(f,'fig_psnr','-dpng') % then print it

% NRS
f = figure();
hold on
bar(1:n_exp, [nrs_ghexin_avg; nrs_dgsot_avg; nrs_ghng_avg]', 'EdgeColor', 'none')
errorbar(x(1,:), nrs_ghexin_avg, nrs_ghexin_sem,'k.')
errorbar(x(2,:), nrs_dgsot_avg, nrs_dgsot_sem,'k.')
errorbar(x(3,:), nrs_ghng_avg, nrs_ghng_sem,'k.')
hold off
xticks(x(2,:))
xticklabels({'X letter','square','Gaussians','videos'})
legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3)
legend('boxoff')
% xlabel('experiments')
ylabel('number of neurons')
box off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_neurons','-dpdf') % then print it
print(f,'fig_neurons','-dpng') % then print it

% SILHOUETTE
f = figure();
hold on
bar(1:n_exp, [silhouette_ghexin_avg; silhouette_dgsot_avg; silhouette_ghng_avg]', 'EdgeColor', 'none')
errorbar(x(1,:), silhouette_ghexin_avg, silhouette_ghexin_sem,'k.')
errorbar(x(2,:), silhouette_dgsot_avg, silhouette_dgsot_sem,'k.')
errorbar(x(3,:), silhouette_ghng_avg, silhouette_ghng_sem,'k.')
hold off
xticks(x(2,:))
xticklabels({'X letter','square','Gaussians','videos'})
legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3)
legend('boxoff')
% xlabel('experiments')
ylabel('silhouette index')
box off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_silhouette','-dpdf') % then print it
print(f,'fig_silhouette','-dpng') % then print it

% DBINDEX
f = figure();
hold on
bar(1:n_exp, [dbindex_ghexin_avg; dbindex_dgsot_avg; dbindex_ghng_avg]', 'EdgeColor', 'none')
errorbar(x(1,:), dbindex_ghexin_avg, dbindex_ghexin_sem,'k.')
errorbar(x(2,:), dbindex_dgsot_avg, dbindex_dgsot_sem,'k.')
errorbar(x(3,1:3), dbindex_ghng_avg(1:3), dbindex_ghng_sem(1:3),'k.')
% set(gca,'yscale','log')
hold off
xticks(x(2,:))
xticklabels({'X letter','square','Gaussians','videos'})
l = legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3);
legend('boxoff')
% xlabel('experiments')
ylabel('Davies-Bouldin index')
box off
ylim([0 99])
xlim([0.5, 4.5])
breakyaxis([3.75, 97]);
l.TextColor = 'black';
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_dbindex','-dpdf') % then print it
print(f,'fig_dbindex','-dpng') % then print it

% ET
f = figure();
hold on
bar(1:n_exp, [ET_ghexin_avg; ET_dgsot_avg; ET_ghng_avg]', 'EdgeColor', 'none')
errorbar(x(1,:), ET_ghexin_avg, ET_ghexin_sem,'k.')
errorbar(x(2,:), ET_dgsot_avg, ET_dgsot_sem,'k.')
errorbar(x(3,:), ET_ghng_avg, ET_ghng_sem,'k.')
hold off
xticks(x(2,:))
xticklabels({'X letter','square','Gaussians','videos'})
legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3)
legend('boxoff')
% xlabel('experiments')
ylabel('training time')
box off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_training_time','-dpdf') % then print it
print(f,'fig_training_time','-dpng') % then print it








% for i = 1:n_rep
%     [psnr_ghexin(i, 2), nrs_ghexin(i, 2), ...
%         psnr_dgsot(i, 2), nrs_dgsot(i, 2), ...
%         psnr_ghng(i, 2), nrs_ghng(i, 2)] = test_models_spiral();
% end
%
% for i = 1:1
%     [psnr_ghexin(i, 4), nrs_ghexin(i, 4), ...
%         psnr_dgsot(i, 4), nrs_dgsot(i, 4), ...
%         psnr_ghng(i, 4), nrs_ghng(i, 4)] = test_models_colors();
% end
% 
% for i = 1:n_rep
%     [psnr_ghexin(i, 5), nrs_ghexin(i, 5), ...
%         psnr_dgsot(i, 5), nrs_dgsot(i, 5), ...
%         psnr_ghng(i, 5), nrs_ghng(i, 5)] = test_models_twinp;
% end

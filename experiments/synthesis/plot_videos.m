clc; clear all; dbstop if error; close all;

levels = 2:4;
n_exp = 1;
n_rep = 2;
% experiments = [0, 1];

dbindex_ghexin_avg = zeros(length(levels), 1);
dbindex_ghexin_sem = zeros(length(levels), 1);
dbindex_dgsot_avg = zeros(length(levels), 1);
dbindex_dgsot_sem = zeros(length(levels), 1);
dbindex_ghng_avg = zeros(length(levels), 1);
dbindex_ghng_sem = zeros(length(levels), 1);

lev = 1;
for k = levels

    j = 1;    
    dbindex_ghexin = zeros(n_rep, n_exp);
    dbindex_dgsot = zeros(n_rep, n_exp);
    dbindex_ghng = zeros(n_rep, n_exp);
    
    for i = 1:n_rep
        rng(i)
        [psnr_ghexin(i, j), nrs_ghexin(i, j), silhouette_ghexin(i, j), dbindex_ghexin(i, j), ET_ghexin(i, j), ...
            psnr_dgsot(i, j), nrs_dgsot(i, j), silhouette_dgsot(i, j), dbindex_dgsot(i, j), ET_dgsot(i, j), ...
            psnr_ghng(i, j), nrs_ghng(i, j), silhouette_ghng(i, j), dbindex_ghng(i, j), ET_ghng(i, j)] = test_models_videos(i, k);
        close all
    end

    % DBINDEX
    dbindex_ghexin_avg(lev) = mean(dbindex_ghexin);
    dbindex_ghexin_sem(lev) = std(dbindex_ghexin)/sqrt(n_rep);
    dbindex_dgsot_avg(lev) = mean(dbindex_dgsot);
    dbindex_dgsot_sem(lev) = std(dbindex_dgsot)/sqrt(n_rep);
    dbindex_ghng_avg(lev) = mean(dbindex_ghng);
    dbindex_ghng_sem(lev) = std(dbindex_ghng)/sqrt(n_rep);
    
    lev = lev + 1;
end

ngroups = 3;
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = zeros(ngroups, nbars);
for i = 1:nbars
    x(i, :) = (1:nbars) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
end

% DBINDEX
f = figure();
hold on
% set(gca,'yscale','log')
bar(1:3, [dbindex_ghexin_avg, dbindex_dgsot_avg, dbindex_ghng_avg], 'EdgeColor', 'none')
errorbar(x(1,:), dbindex_ghexin_avg, dbindex_ghexin_sem,'k.')
errorbar(x(2,:), dbindex_dgsot_avg, dbindex_dgsot_sem,'k.')
errorbar(x(3,:), dbindex_ghng_avg, dbindex_ghng_sem,'k.')
breakyaxis([8, 48]);
% breakyaxis([55, 60]);
hold off
xticks(x(2,:))
xticklabels({'level 2','level 3','level 4'})
legend({'gh-exin','dgsot','ghng'},'Location', 'southoutside', 'NumColumns', 3)
legend('boxoff')
% xlabel('experiments')
ylabel('Davies-Bouldin index')
box off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(f,'PaperSize',[5.5 4]); %set the paper size to what you want
print(f,'fig_dbindex_videos','-dpdf') % then print it
print(f,'fig_dbindex_videos','-dpng') % then print it

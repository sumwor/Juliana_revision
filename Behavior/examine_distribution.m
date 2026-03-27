function [cutoff95_overall, cutoff95] = examine_distribution(data_session, genotypes, label, color_list, savefigpath)

% use fitted parameters to estimate actual behavior


all_data = vertcat(data_session{:});
nbins = 40;
edges = logspace(log10(min(all_data)), log10(max(all_data)), nbins+1);
bin_centers = edges(1:end-1) + diff(edges)/2;

nsess = numel(data_session);
pdf_mat = nan(nsess, nbins);
cutoff95 = nan(nsess,1);
for i = 1:nsess
    iti = data_session{i};
    
    counts = histcounts(iti, edges, 'Normalization', 'pdf');
    pdf_mat(i,:) = counts;
    cutoff95(i) = prctile(iti,95);
end
cutoff95_overall = prctile(all_data,95);
isWT  = strcmp(genotypes, 'WT');
isHet = strcmp(genotypes, 'HET');

pdf_WT  = pdf_mat(isWT,:);
pdf_Het = pdf_mat(isHet,:);

mean_WT  = mean(pdf_WT,1,'omitnan');
sem_WT   = std(pdf_WT,0,1,'omitnan') ./ sqrt(sum(isWT));

mean_Het = mean(pdf_Het,1,'omitnan');
sem_Het  = std(pdf_Het,0,1,'omitnan') ./ sqrt(sum(isHet));

figure; hold on
sgtitle(label)
% WT
subplot(2,2,1)
title('WT')
hold on;
for i = find(isWT)
    plot(bin_centers, pdf_mat(i,:), 'Color', color_list{2},'LineWidth', 0.5,'LineStyle',':')
    
end
errorbar(bin_centers, mean_WT, sem_WT,  'Color', color_list{2}, 'LineWidth', 2)
xline(cutoff95_overall, 'r--', 'LineWidth', 2);
set(gca,'XScale','log')
xlabel('Time (s)')
ylabel('Probability')


% Het
subplot(2,2,2)
title('Mutant')
hold on;
for i = find(isHet)
     plot(bin_centers, pdf_mat(i,:), 'Color', color_list{1},'LineWidth', 0.5,'LineStyle',':')
end
errorbar(bin_centers, mean_Het, sem_Het,'Color', color_list{1}, 'LineWidth', 2)
xline(cutoff95_overall, 'r--', 'LineWidth', 2);
set(gca,'XScale','log')
xlabel('Time (s)')

% test bimidol 
deltaBIC = zeros(1, nsess);
for i = 1:nsess

    iti = data_session{i};
    logITI = log(iti(:));
    
    % Fit 1-component
    gm1 = fitgmdist(logITI,1);
    
    % Fit 2-component
    gm2 = fitgmdist(logITI,2);
    
    % Compare BIC
    bic1 = gm1.BIC;
    bic2 = gm2.BIC;
    deltaBIC(i) = bic2-bic1;
end


deltaBIC_WT  = deltaBIC(isWT);
deltaBIC_Het = deltaBIC(isHet);

subplot(2,2,3)
histogram(deltaBIC_WT, 'FaceColor', color_list{2}, 'FaceAlpha', 0.5)
xlabel('BIC2 - BIC1')
ylabel('Number of sessions')
set(gca,'box','off')

subplot(2,2,4)
histogram(deltaBIC_Het, 'FaceColor', color_list{1}, 'FaceAlpha', 0.5)
xlabel('BIC2 - BIC1')
set(gca,'box','off')

%% save them in a csv file
nSub = size(pdf_mat,1);
nBins = size(pdf_mat,2);

% Repeat subject/genotype for each bin
Subject  = repelem((1:nSub)', nBins);
Genotype = repelem(genotypes(:), nBins);

% Repeat bin centers for each subject
BinCenter = repmat(bin_centers(:), nSub, 1);

% Flatten pdf values
PDF_temp = pdf_mat';
PDF = PDF_temp(:);

% Add cutoff column (same value for all rows)
Cutoff95 = repmat(cutoff95_overall, nSub*nBins, 1);

% Create table
T = table(Subject, Genotype, BinCenter, PDF, Cutoff95);

% Write to CSV
savecsvfile = fullfile(savefigpath(1:end-8),'Results', [label,'ITL_distribution.csv']);
writetable(T,savecsvfile);

    print(gcf,'-dpng',fullfile(savefigpath, label));    %png format
    saveas(gcf, fullfile(savefigpath, label), 'fig');
    saveas(gcf, fullfile(savefigpath, label),'svg');
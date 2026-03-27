function plot_performance_quantile(data, strain, genotype, savefigpath, tlabel)

if contains(tlabel,'AB') || contains(tlabel, 'CD')
    nPlot = 3;
elseif contains(tlabel, 'DC')
    nPlot = 6;
end

figure;
sgtitle(['Performance in the ', tlabel, ' sessions in quantile'])

for pp = 1:nPlot
    subplot(2,3,pp)

    genotype_list = unique(genotype);
    
   if strcmp(strain, 'Cntnap2')  % only look at WT and KO
        genotype_list = {'KO', 'WT'};
   elseif strcmp(strain, 'Shank3B')
       genotype_list={'HET','WT'};
    end
    %color_list = {'red','blue', 'magenta'};

    % temporary color for JC revision
    colors = struct;
    colors.TSC2_det_male = [72, 161, 217] / 255;
    colors.TSC2_det_female = [22, 111, 56] / 255;
    colors.Shank3_det_male = [50, 133, 141] / 255;
    colors.Shank3_det_female = [248, 151, 29] / 255;
    colors.TSC2_prob_male =  [125, 83, 162] / 255;
    colors.Shnk3_prob_male = [237, 28, 36] / 255;
    color_list = {colors.Shank3_det_male, 'black'}
    for gg = 1:length(genotype_list)
        mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),:,pp),1);
        if sum(strcmp(genotype, genotype_list{gg}))>1
            ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),:,pp),1)/sum(strcmp(genotype, genotype_list{gg}));
        else
            ste_perf = [0, 0, 0, 0];
        end
        %plot(mean_perf)
        hold on;
        errorbar([1:4], mean_perf, ste_perf, 'LineWidth', 2, 'Color',color_list{gg});

        % plot each session in thinner lines
        tempData= data(strcmp(genotype, genotype_list{gg}),:,pp);
        for tt =1:size(tempData,1)
            plot([1:4], tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
        end
    end
    ylim([0, 1]);
    xticks([1,2,3,4]);
    %legend(genotype_list)
    title(['Session ', num2str(pp)])
    ylabel('Performance');

    % two-way anova
    if ismember('HEM', genotype_list)
        Hetdata = data(strcmp(genotype, 'HEM'),:,pp);
        nHet = sum(strcmp(genotype, 'HEM'));
    elseif ismember('KO', genotype_list)
        Hetdata = data(strcmp(genotype, 'KO'),:,pp);
        nHet = sum(strcmp(genotype, 'KO'));
    else
        Hetdata = data(strcmp(genotype, 'HET'),:,pp);
        nHet = sum(strcmp(genotype, 'HET'));
    end
    WTdata = data(strcmp(genotype, 'WT'),:,pp);
    nWT = sum(strcmp(genotype, 'WT'));
    
    dat = [WTdata(:); Hetdata(:)];
    group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
    quantile = [repmat((1:4)', size(WTdata, 1), 1); repmat((1:4)', size(Hetdata, 1), 1)];   % Quantile labels
    nQuantile = 4;
    subject = [repelem((1:nWT)', nQuantile, 1);
                repelem((1:nHet)' + nWT, nQuantile, 1)];
 
    tblData = table(dat, categorical(group), quantile, categorical(subject), ...
    'VariableNames', {'Y', 'Group', 'Block', 'Subject'});

    % Ensure WT is the reference genotype
    tblData.Group = reordercats(tblData.Group, {'WT', 'HET'});
    tblData.Block2 = tblData.Block.^2;
    tblData.lBlock = log(tblData.Block);
    % --- Fit linear mixed-effects model ---
    % Block as random intercept (repeated measure within Subject)
    %lme = fitlme(tblData, 'Y ~ Group*Block + Group * Block2 + (Block|Subject)');
    
    % generalized linear model
    tblData.NumCorrect = round(tblData.Y * 100);   % Convert proportion → count
    TotalTrials = 100;  
    glme = fitglme(tblData, ...
        'NumCorrect ~ Group*Block +  (1+Block|Subject)', ...
        'Distribution','Binomial','Link','logit','BinomialSize',100);
            % 


           % --- Display results ---
        disp(anova(glme)); % Mixed ANOVA table
        pTable = anova(glme);
        pGroup = pTable.pValue(strcmp(pTable.Term, 'Group'));
        pLearn = pTable.pValue(strcmp(pTable.Term, 'Block'));
        pInt   = pTable.pValue(strcmp(pTable.Term, 'Group:Block'));

        % --- Show p-values on figure ---
        text(1, 0.25, ['P(geno): ', num2str(pGroup, '%.3g')], 'FontSize', 15)
        text(1, 0.18, ['P(learn): ', num2str(pLearn, '%.3g')], 'FontSize', 15)
        text(1, 0.11, ['P(geno/timeslearn): ', num2str(pInt, '%.3g')], 'FontSize', 15)
    
    
    if pp==nPlot
        legend(genotype_list, 'Location','southeast')
        legend('Box','off')
        text(3, 0.1, ['WT ', num2str(nWT)], 'FontSize', 10 );
        text(3, 0.2, ['HET ', num2str(nHet)], 'FontSize', 10 );
    end
end
%% savefigs


print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in quantile']));    %png format
savefig(fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile.fig']));
% saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile']),'svg');


% close all including anova windows
delete(findall(0, 'Type', 'figure'))

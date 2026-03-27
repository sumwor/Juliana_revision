% generate csv files quickly
root_dir = 'Y:\HongliWang\Juvi_ASD Deterministic\HitStrikeCheck';

matfiles = dir(fullfile(root_dir, '*.mat'));

% disable warning
warning('off', 'all'); % Disable all warnings
for ii = 1:length(matfiles)
    matFile = fullfile(matfiles(ii).folder, matfiles(ii).name);
    outfname = fullfile(matfiles(ii).folder, [matfiles(ii).name(1:end-4),'.csv']);
    %if ~exist(outfname)
        resultdf = extract_behavior_df(matFile);
        writetable(resultdf, outfname);
    
        % check if there is rewardsize 1 exist in the log
        if sum(abs(resultdf.trial_types-0.01) < 0.0005) > 0
            display(['Reward size 1 detected', matfiles(ii).name, num2str(sum((resultdf.trial_types-0.01) < 0.0005))])
        end
    %end
end

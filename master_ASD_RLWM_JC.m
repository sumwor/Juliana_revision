%script master_ASD_RLWM
% for Juliana revision 
clear all;
%% to do:
% check the timestamp of trial start to make sure there is no carrying over
% from last session under same protocol


%% masterfile to process the behavior log files of ASD_RLWM behavior

root_dir = 'Z:\HongliWang\Juliana_revision_data';

strain_list = {'TSC2_male_det'; 'Shank3_male_det'; 'Shank3_female_det'; 'TSC2_female_det'; 'Shank3_male_prob'; 'TSC2_male_prob'};

%% add for loop later
colors = struct;
colors.TSC2_det_male = [72, 161, 217] / 255;
colors.TSC2_det_female = [22, 111, 56] / 255;
colors.Shank3_det_male = [50, 133, 141] / 255;
colors.Shank3_det_female = [248, 151, 29] / 255;
colors.TSC2_prob_male =  [125, 83, 162] / 255;
colors.Shnk3_prob_male = [237, 28, 36] / 255;

strainNum =1;

% save the dataIndex data
dataIndexPath = fullfile(root_dir, strain_list{strainNum},'dataIndex.csv');

%if ~exist(dataIndexPath)
animalList = readtable(fullfile(root_dir, strain_list{strainNum},'Data','AnimalList.csv'));

%dataIndex = makeDataIndex_ASD(logfilepath);

logfilepath = fullfile(root_dir, strain_list{strainNum},'Data');
analysispath = fullfile(root_dir, strain_list{strainNum},'Analysis');
dataIndex = makeDataIndex_ASD(logfilepath, analysispath);

nFiles = size(dataIndex,1);
% parse every .mat file, generate .csv files

ErrorList = table([],[],'VariableNames',{'Session','ErrorMessage'});

odors = cell(nFiles,1);
for ii = 1:nFiles

    if ii==1
        dataIndex.BehCSV = cell(nFiles,1);
        dataIndex.OdorPresented = cell(nFiles,1);
        dataIndex.RRABL = nan(nFiles,1);
        dataIndex.RRABR = nan(nFiles,1);
        dataIndex.RRCDL = nan(nFiles,1);
        dataIndex.RRCDR = nan(nFiles,1);
        dataIndex.RRDCL = nan(nFiles,1);
        dataIndex.RRDCR = nan(nFiles,1);
        dataIndex.EntropyAB = nan(nFiles,1);
        dataIndex.EntropyCD = nan(nFiles,1);
        dataIndex.EntropyDC = nan(nFiles,1);
        dataIndex.nAB = nan(nFiles,1);
        dataIndex.nCD = nan(nFiles, 1);
        dataIndex.nDC = nan(nFiles,1);
    end
    if ~exist(dataIndex.BehPath{ii})
        mkdir(dataIndex.BehPath{ii});
    end
    outfname = fullfile(dataIndex.BehPath{ii}, sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii}));
    dataIndex.BehCSV{ii} = sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii});

    try
        if ~exist(outfname)
            resultdf = extract_behavior_df(fullfile(dataIndex.LogFilePath{ii}, dataIndex.LogFileName{ii}));
            writetable(resultdf, outfname);
        else
            resultdf = readtable(outfname);
        end
        %   end

        % parse through the result to check what odor presented in the behavior
        % file
        odor_presented = unique(resultdf.schedule(~isnan(resultdf.schedule)));
        dataIndex.OdorPresented{ii} = odor_presented;

        % check the if the protocol is correct
        odor_plan.AB = [1,2];
        odor_plan.ABCD = [1,2,3,4];
        odor_plan.ABCDDC = [1,2,3,4,5,6];
        odor_plan.ABDC = [1,2,5,6];

        if strcmp(dataIndex.Protocol{ii},'AB')
            odor = odor_plan.AB;
        elseif strcmp(dataIndex.Protocol{ii},'AB-CD')
            odor = odor_plan.ABCD;
        elseif strcmp(dataIndex.Protocol{ii},'AB-CD-DC')
            odor = odor_plan.ABCDDC;
        elseif strcmp(dataIndex.Protocol{ii}, 'AB-DC')
            odor = odor_plan.ABDC;
        end
        if ~isequal(odor_presented,odor')
            display(['Unmatching odor and file name:',dataIndex.LogFileName{ii}])
        end

        %% calculate average reward rate - to compare with the notebook
        % double check to make sure correct .mat file is logged here
        if strcmp(dataIndex.Protocol{ii},'AB')
            ABEnd = size(resultdf,1);
            CDEnd = NaN;
            DCEnd = NaN;
            dataIndex.RRABL(ii) = sum((resultdf.schedule(1:ABEnd) == 1) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 1);
            dataIndex.RRABR(ii) = sum((resultdf.schedule(1:ABEnd) == 2) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 2);
            %dataIndex.EntropyAB(ii) = get_stimulus_entropy(resultdf.schedule(1:ABEnd));
            dataIndex.nAB(ii) = size(resultdf,1);

        elseif strcmp(dataIndex.Protocol{ii},'AB-CD')
            ABEnd = find(resultdf.schedule == 3 | resultdf.schedule == 4, 1, 'first')-1;
            CDEnd = size(resultdf,1);
            DCEnd = NaN;
            dataIndex.RRABL(ii) = sum((resultdf.schedule(1:ABEnd) == 1) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 1);
            dataIndex.RRABR(ii) = sum((resultdf.schedule(1:ABEnd) == 2) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 2);
            dataIndex.RRCDL(ii) = sum((resultdf.schedule(ABEnd+1:CDEnd) == 3) & ~isnan(resultdf.reward(ABEnd+1:CDEnd)))/sum(resultdf.schedule(ABEnd+1:CDEnd) == 3);
            dataIndex.RRCDR(ii) = sum((resultdf.schedule(ABEnd+1:CDEnd) == 4) & ~isnan(resultdf.reward(ABEnd+1:CDEnd)))/sum(resultdf.schedule(ABEnd+1:CDEnd) == 4);
            dataIndex.EntropyAB(ii) = get_stimulus_entropy(resultdf.schedule(1:ABEnd));
            dataIndex.EntropyCD(ii) = get_stimulus_entropy(resultdf.schedule(ABEnd+1:CDEnd));
            dataIndex.nAB(ii) = ABEnd;
            dataIndex.nCD(ii) = CDEnd-ABEnd;
        elseif strcmp(dataIndex.Protocol{ii}, 'AB-CD-DC')
            ABEnd = find(resultdf.schedule == 3 | resultdf.schedule == 4, 1, 'first')-1;
            CDEnd = find(resultdf.schedule == 5 | resultdf.schedule == 6, 1, 'first')-1;
            DCEnd = size(resultdf,1);
            dataIndex.RRABL(ii) = sum((resultdf.schedule(1:ABEnd) == 1) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 1);
            dataIndex.RRABR(ii) = sum((resultdf.schedule(1:ABEnd) == 2) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 2);
            dataIndex.RRCDL(ii) = sum((resultdf.schedule(ABEnd+1:CDEnd) == 3) & ~isnan(resultdf.reward(ABEnd+1:CDEnd)))/sum(resultdf.schedule(ABEnd+1:CDEnd) == 3);
            dataIndex.RRCDR(ii) = sum((resultdf.schedule(ABEnd+1:CDEnd) == 4) & ~isnan(resultdf.reward(ABEnd+1:CDEnd)))/sum(resultdf.schedule(ABEnd+1:CDEnd) == 4);
            dataIndex.RRDCL(ii) = sum((resultdf.schedule(CDEnd+1:DCEnd) == 6) & ~isnan(resultdf.reward(CDEnd+1:DCEnd)))/sum(resultdf.schedule(CDEnd+1:DCEnd) == 6);
            dataIndex.RRDCR(ii)= sum((resultdf.schedule(CDEnd+1:DCEnd) == 5) & ~isnan(resultdf.reward(CDEnd+1:DCEnd)))/sum(resultdf.schedule(CDEnd+1:DCEnd) == 5);
            %dataIndex.EntropyAB(ii) = get_stimulus_entropy(resultdf.schedule(1:ABEnd));
            %dataIndex.EntropyCD(ii) = get_stimulus_entropy(resultdf.schedule(ABEnd+1:CDEnd));
            %dataIndex.EntropyDC(ii) = get_stimulus_entropy(resultdf.schedule(CDEnd+1:DCEnd));
            dataIndex.nAB(ii) = ABEnd;
            dataIndex.nCD(ii) = CDEnd-ABEnd;
            dataIndex.nDC(ii) = DCEnd-CDEnd;
        elseif strcmp(dataIndex.Protocol{ii}, 'AB-DC')
            ABEnd = find(resultdf.schedule == 5 | resultdf.schedule == 6, 1, 'first')-1;
            CDEnd = NaN;
            DCEnd = size(resultdf,1);
            dataIndex.RRABL(ii) = sum((resultdf.schedule(1:ABEnd) == 1) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 1);
            dataIndex.RRABR(ii) = sum((resultdf.schedule(1:ABEnd) == 2) & ~isnan(resultdf.reward(1:ABEnd)))/sum(resultdf.schedule(1:ABEnd) == 2);
            dataIndex.RRDCL(ii) = sum((resultdf.schedule(ABEnd+1:DCEnd) == 6) & ~isnan(resultdf.reward(ABEnd+1:DCEnd)))/sum(resultdf.schedule(ABEnd+1:DCEnd) == 6);
            dataIndex.RRDCR(ii) = sum((resultdf.schedule(ABEnd+1:DCEnd) == 5) & ~isnan(resultdf.reward(ABEnd+1:DCEnd)))/sum(resultdf.schedule(ABEnd+1:DCEnd) == 5);
            dataIndex.EntropyAB(ii) = get_stimulus_entropy(resultdf.schedule(1:ABEnd));
            dataIndex.EntropyDC(ii) = get_stimulus_entropy(resultdf.schedule(ABEnd+1:DCEnd));
            dataIndex.nAB(ii) = ABEnd;
            dataIndex.nDC(ii) = DCEnd-ABEnd;
        end




    catch ME
        newEntry = {[dataIndex.Animal{ii},'_', dataIndex.Session{ii}], ME.message};
        ErrorList = [ErrorList; newEntry];
    end

end

% save the error table
writetable(ErrorList, fullfile(root_dir,strain_list{strainNum}, 'BuggedSessions.csv'));


writetable(dataIndex, dataIndexPath);
%else
%dataIndex=readtable(dataIndexPath);
%nFiles = size(dataIndex,1);
%end

savesummaryfolder = fullfile(root_dir, strain_list{strainNum},'Summary');
if ~exist(savesummaryfolder)
    mkdir(savesummaryfolder)
end
savefigpath = fullfile(root_dir,strain_list{strainNum},'Summary','BehPlot');
if ~exist(savefigpath)
    mkdir(savefigpath)
end
savedatapath = fullfile(root_dir,strain_list{strainNum},'Summary','Results');
if ~exist(savedatapath)
    mkdir(savedatapath)
end


%% RL hybrid model
% Jing-jing model
% prepare data for model fit first
files = ASD_hybrid_dataPrep(dataIndex,savedatapath);

% fit the model
ASD_hybrid_modelFit(dataIndex, files, savedatapath, savefigpath);

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB



ASD_odor_summary(dataIndex, strain_list{strainNum}, savefigpath, savedatapath);


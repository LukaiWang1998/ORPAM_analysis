clear;
close all;
% TablePath = '\\10.229.121.108\Workspace\Lukai\EndometrialCancer\imagingRecord.xlsx';
TablePath = '/Volumes/Workspace/Lukai/EndometrialCancer/imagingRecord.xlsx';
imagingRecordTable = readtable(TablePath, 'Sheet', 'Apr18_processing', 'VariableNamingRule', 'preserve');
% filepath = '\\10.229.121.108\DataArchive\ORPAM_Data\EndometrialCancer';
filepath = '/Volumes/DataArchive/ORPAM_Data/EndometrialCancer';

Cases = getFolderNames(filepath);

i = 1;

aqdate = Cases{i};
parts = split(Cases{i}, '_');
ID = parts{4};

Dir = fullfile(filepath, aqdate);
S_temp = dir(fullfile(Dir,'*.bin'));
N_bins = numel(S_temp);
fprintf('>>>>>>>> FOUND %d SCANS\n', N_bins)
%% SELECT SCAN AND LOAD DATA
sample_idx = 1;
F = fullfile(Dir,S_temp(sample_idx).name);
tokens = regexp(S_temp(sample_idx).name, '^(.*)\.bin$', 'tokens');
filename = tokens{1}{1};





% Find the matching row
row_idx = find(imagingRecordTable.("Patient ID") == str2double(ID) & ...
               strcmp(imagingRecordTable.Position, filename));

% Display result
if isempty(row_idx)
    disp("No match found.");
elseif length(row_idx) > 1
    warning("Multiple matches found:");
    disp(imagingRecordTable(row_idx, :));
else
    disp("Matching row:");
    disp(imagingRecordTable(row_idx, :));
end






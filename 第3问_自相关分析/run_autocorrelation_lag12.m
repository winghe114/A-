clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataFile = fullfile(projectDir, 'data.xlsx');
outCsv = fullfile(scriptDir, 'autocorrelation_lag12_results.csv');
outMat = fullfile(scriptDir, 'autocorrelation_lag12_results.mat');

opts = detectImportOptions(dataFile, 'Sheet', 2, 'VariableNamingRule', 'preserve');
T = readtable(dataFile, opts);

t = T{:, 1};
x = T{:, 2};
valid = isfinite(t) & isfinite(x);
t = t(valid);
x = x(valid);

xCentered = x - mean(x);
denominator = sum(xCentered .^ 2);
lags = [1; 2];
autocorrValues = zeros(numel(lags), 1);

for i = 1:numel(lags)
    k = lags(i);
    autocorrValues(i) = sum(xCentered(1:end-k) .* xCentered(1+k:end)) / denominator;
end

resultTable = table(lags, autocorrValues, 'VariableNames', {'lag', 'autocorrelation'});
writetable(resultTable, outCsv);
save(outMat, 't', 'x', 'lags', 'autocorrValues', 'resultTable');

disp(resultTable);

%% Question 3 add-on: sub-bin two-frequency super-resolution experiment
% This script tests whether joint nonlinear frequency refinement can resolve
% two close sinusoidal fault components beyond the full-length FFT bin spacing.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q3Dir = fileparts(scriptDir);
rootDir = fileparts(q3Dir);
dataFile = fullfile(rootDir, 'data.xlsx');
resultDir = scriptDir;
figureDir = scriptDir;

fprintf('Question 3 super-resolution experiment started...\n');

[t, ~] = readMultiSourceData(dataFile);
t = t(:);
dt = median(diff(t));
fs = 1 / dt;
duration = t(end) - t(1);
fftBinHz = 1 / duration;

baseFreq = 10.0;
A1 = 0.030;
A2 = 0.026;
phi1 = 0.35;
phi2 = -1.20;

deltaList = [0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0050]';
snrList = [-20, -16, -12, -8]';
trials = 3;
rng(20260621);

rows = numel(deltaList) * numel(snrList);
delta_f_Hz = zeros(rows, 1);
SNR_dB = zeros(rows, 1);
fft_success_rate = zeros(rows, 1);
varpro_success_rate = zeros(rows, 1);
varpro_mean_abs_error_Hz = NaN(rows, 1);

fftMatrix = zeros(numel(snrList), numel(deltaList));
varproMatrix = zeros(numel(snrList), numel(deltaList));

row = 0;
for iS = 1:numel(snrList)
    for iD = 1:numel(deltaList)
        delta = deltaList(iD);
        trueFreqs = [baseFreq; baseFreq + delta];
        clean = A1 * sin(2*pi*trueFreqs(1)*t + phi1) + ...
            A2 * sin(2*pi*trueFreqs(2)*t + phi2);
        cleanPower = mean(clean .^ 2);
        noiseSigma = sqrt(cleanPower / (10^(snrList(iS) / 10)));

        fftOk = 0;
        varproOk = 0;
        errSum = 0;
        errCount = 0;
        tolerance = max(0.00012, min(0.0012, 0.35 * delta));

        for tr = 1:trials
            x = clean + noiseSigma * randn(size(t));
            x = x - mean(x, 'omitnan');

            fftFreqs = estimateTwoPeaksByFullLengthFft(x, fs, baseFreq - 0.02, baseFreq + 0.04);
            [okFft, ~] = checkTwoFrequencySuccess(fftFreqs, trueFreqs, tolerance);
            fftOk = fftOk + okFft;

            center0 = estimateSingleFrequencyByProjection(t, x, baseFreq - 0.015, baseFreq + 0.025);
            varproFreqs = refineTwoCloseFrequenciesVarpro(t, x, center0, ...
                baseFreq - 0.02, baseFreq + 0.04);
            [okVarpro, errVarpro] = checkTwoFrequencySuccess(varproFreqs, trueFreqs, tolerance);
            varproOk = varproOk + okVarpro;
            if okVarpro
                errSum = errSum + errVarpro;
                errCount = errCount + 1;
            end
        end

        row = row + 1;
        delta_f_Hz(row) = delta;
        SNR_dB(row) = snrList(iS);
        fft_success_rate(row) = fftOk / trials;
        varpro_success_rate(row) = varproOk / trials;
        if errCount > 0
            varpro_mean_abs_error_Hz(row) = errSum / errCount;
        end
        fftMatrix(iS, iD) = fft_success_rate(row);
        varproMatrix(iS, iD) = varpro_success_rate(row);
    end
end

resultTable = table(delta_f_Hz, SNR_dB, fft_success_rate, ...
    varpro_success_rate, varpro_mean_abs_error_Hz);
writetable(resultTable, fullfile(resultDir, 'q3_superresolution_comparison.csv'));

fftTable = matrixToTable(fftMatrix, deltaList, snrList);
varproTable = matrixToTable(varproMatrix, deltaList, snrList);
writetable(fftTable, fullfile(resultDir, 'q3_superresolution_fft_success_matrix.csv'));
writetable(varproTable, fullfile(resultDir, 'q3_superresolution_varpro_success_matrix.csv'));

makeComparisonHeatmap(figureDir, fftMatrix, varproMatrix, deltaList, snrList, fftBinHz);
writeSuperresolutionSummary(resultDir, fftBinHz, resultTable);

fprintf('Question 3 super-resolution experiment finished. FFT bin spacing = %.6g Hz\n', fftBinHz);

%% Local functions
function [t, x] = readMultiSourceData(dataFile)
    try
        tbl = readtable(dataFile, 'Sheet', 2, 'VariableNamingRule', 'preserve');
        t = tbl{:, 1};
        x = tbl{:, 2};
    catch
        num = readmatrix(dataFile, 'Sheet', 2);
        t = num(:, 1);
        x = num(:, 2);
    end
    good = isfinite(t) & isfinite(x);
    t = t(good);
    x = x(good);
end

function fHat = estimateSingleFrequencyByProjection(t, x, lower, upper)
    objective = @(f) -projectionEnergyAtFrequency(t, x, f);
    opts = optimset('TolX', 1e-13, 'Display', 'off');
    fHat = fminbnd(objective, lower, upper, opts);
end

function freqs = estimateTwoPeaksByFullLengthFft(x, fs, lower, upper)
    x = x(:) - mean(x, 'omitnan');
    n = numel(x);
    freq = (0:floor(n / 2))' * fs / n;
    X = fft(x);
    power = abs(X(1:numel(freq))).^2;
    valid = freq >= lower & freq <= upper;
    idx = find(valid);
    loc = zeros(0, 1);
    for ii = 2:numel(idx)-1
        k = idx(ii);
        if power(k) > power(k-1) && power(k) > power(k+1)
            loc(end+1, 1) = k; %#ok<AGROW>
        end
    end
    if numel(loc) < 2
        [~, order] = sort(power(idx), 'descend');
        loc = idx(order(1:min(2, numel(order))));
    else
        [~, order] = sort(power(loc), 'descend');
        loc = loc(order(1:min(2, numel(order))));
    end
    freqs = sort(freq(loc(:)));
end

function freqs = refineTwoCloseFrequenciesVarpro(t, x, center0, lower, upper)
    centerGrid = linspace(center0 - 0.0025, center0 + 0.0025, 41);
    deltaGrid = [0.00025, 0.00040, 0.00060, 0.00080, 0.00100, 0.00125, ...
        0.00150, 0.00175, 0.00200, 0.00250, 0.00300, 0.00400, 0.00500, 0.00750];
    [bestCenter, bestDelta] = gridSearchTwoFrequency(t, x, centerGrid, deltaGrid, lower, upper);

    centerFine = linspace(bestCenter - 0.00045, bestCenter + 0.00045, 25);
    deltaFine = linspace(max(0.00012, bestDelta - 0.00045), bestDelta + 0.00045, 31);
    [bestCenter, bestDelta] = gridSearchTwoFrequency(t, x, centerFine, deltaFine, lower, upper);
    freqs = sort([bestCenter - bestDelta / 2; bestCenter + bestDelta / 2]);
end

function [bestCenter, bestDelta] = gridSearchTwoFrequency(t, x, centerGrid, deltaGrid, lower, upper)
    bestSse = inf;
    bestCenter = centerGrid(1);
    bestDelta = deltaGrid(1);
    for c = centerGrid(:)'
        for d = deltaGrid(:)'
            freqs = sort([c - d/2; c + d/2]);
            if freqs(1) < lower || freqs(2) > upper
                continue;
            end
            [~, ~, residual] = fitMultiSinusoid(t, x, freqs);
            sse = sum(residual .^ 2);
            if sse < bestSse
                bestSse = sse;
                bestCenter = c;
                bestDelta = d;
            end
        end
    end
end

function energy = projectionEnergyAtFrequency(t, x, f)
    [~, fit, ~] = fitSinusoidAtFrequency(t, x, f);
    energy = sum(fit .^ 2);
end

function [theta, fit, residual] = fitSinusoidAtFrequency(t, x, f)
    t = t(:);
    x = x(:) - mean(x, 'omitnan');
    H = [cos(2*pi*f*t), sin(2*pi*f*t)];
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function [theta, fit, residual] = fitMultiSinusoid(t, x, freqs)
    t = t(:);
    x = x(:) - mean(x, 'omitnan');
    freqs = sort(freqs(:));
    H = zeros(numel(t), 2 * numel(freqs));
    for k = 1:numel(freqs)
        H(:, 2*k - 1) = cos(2*pi*freqs(k)*t);
        H(:, 2*k) = sin(2*pi*freqs(k)*t);
    end
    theta = H \ x;
    fit = H * theta;
    residual = x - fit;
end

function [ok, meanErr] = checkTwoFrequencySuccess(estFreqs, truthFreqs, tolerance)
    estFreqs = sort(estFreqs(:));
    truthFreqs = sort(truthFreqs(:));
    ok = false;
    meanErr = NaN;
    if numel(estFreqs) < 2
        return;
    end
    d1 = abs(estFreqs - truthFreqs(1));
    [err1, idx1] = min(d1);
    d2 = abs(estFreqs - truthFreqs(2));
    [err2, idx2] = min(d2);
    if idx1 ~= idx2 && err1 <= tolerance && err2 <= tolerance
        ok = true;
        meanErr = mean([err1, err2]);
    end
end

function out = matrixToTable(matrix, deltaList, snrList)
    out = array2table(matrix, ...
        'VariableNames', matlab.lang.makeValidName("df_" + string(deltaList(:)')));
    out.SNR_dB = snrList(:);
    out = movevars(out, 'SNR_dB', 'Before', 1);
end

function makeComparisonHeatmap(outDir, fftMatrix, varproMatrix, deltaList, snrList, fftBinHz)
    deltaMilliHz = 1000 * deltaList;
    fftBinMilliHz = 1000 * fftBinHz;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1260, 560]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    imagesc(deltaMilliHz, snrList, fftMatrix);
    set(gca, 'YDir', 'normal');
    caxis([0, 1]);
    colorbar;
    xline(fftBinMilliHz, 'w--', 'LineWidth', 1.5);
    title('Full-length FFT two-peak success rate');
    xlabel('\Delta f / mHz');
    ylabel('SNR / dB');
    xticks(deltaMilliHz);
    nexttile;
    imagesc(deltaMilliHz, snrList, varproMatrix);
    set(gca, 'YDir', 'normal');
    caxis([0, 1]);
    colorbar;
    xline(fftBinMilliHz, 'w--', 'LineWidth', 1.5);
    title('Variable projection two-frequency success rate');
    xlabel('\Delta f / mHz');
    ylabel('SNR / dB');
    xticks(deltaMilliHz);
    colormap(parula);
    sgtitle(sprintf('Sub-bin frequency resolution comparison, FFT bin spacing = %.1f mHz (%.4f Hz)', ...
        fftBinMilliHz, fftBinHz));
    exportgraphics(fig, fullfile(outDir, 'fig_q3_superresolution_fft_vs_varpro.png'), 'Resolution', 220);
    close(fig);
end

function writeSuperresolutionSummary(outDir, fftBinHz, resultTable)
    fid = fopen(fullfile(outDir, 'q3_superresolution_summary.txt'), 'w', 'n', 'UTF-8');
    fprintf(fid, 'Question 3 super-resolution comparison\n');
    fprintf(fid, 'FFT bin spacing 1/T: %.12g Hz\n\n', fftBinHz);
    fprintf(fid, 'Rows: delta_f_Hz, SNR_dB, FFT success, variable-projection success\n');
    for i = 1:height(resultTable)
        fprintf(fid, 'df=%.6g Hz, SNR=%g dB, FFT=%.3f, VarPro=%.3f, VarPro mean error=%.6g Hz\n', ...
            resultTable.delta_f_Hz(i), resultTable.SNR_dB(i), ...
            resultTable.fft_success_rate(i), resultTable.varpro_success_rate(i), ...
            resultTable.varpro_mean_abs_error_Hz(i));
    end
    fprintf(fid, ['\nInterpretation: variable projection assumes a local two-frequency model and ', ...
        'jointly optimizes both frequencies. It can resolve sub-bin spacing when SNR is sufficient, ', ...
        'whereas direct FFT two-peak picking remains limited by grid spacing and spectral leakage.\n']);
    fclose(fid);
end

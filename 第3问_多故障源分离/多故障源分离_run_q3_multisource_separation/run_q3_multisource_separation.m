%% Question 3: multi-source weak periodic signal separation
% Full-length sinusoidal projection OMP with robust residual thresholding.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q3Dir = fileparts(scriptDir);
rootDir = fileparts(q3Dir);
dataFile = fullfile(rootDir, 'data.xlsx');
resultDir = scriptDir;
figureDir = scriptDir;

fprintf('Question 3 multi-source separation started...\n');
fprintf('Data file: %s\n', dataFile);

[t, xRaw] = readMultiSourceData(dataFile);
t = t(:);
xRaw = xRaw(:);
y = xRaw - mean(xRaw, 'omitnan');
dt = median(diff(t));
fs = 1 / dt;
duration = t(end) - t(1);
dfTheory = 1 / duration;
totalEnergy = sum(y .^ 2);

settings = defaultOmpSettings(fs, duration);
model = runProjectionOmp(t, y, settings, true);

[theta, fitAll, residual, componentMatrix] = fitMultiSinusoid(t, y, model.frequencies_Hz);
componentTable = buildComponentTable(model.frequencies_Hz, theta, componentMatrix, totalEnergy);

sse = sum(residual .^ 2);
rSquared = 1 - sse / totalEnergy;
reconstructedEnergy = sum(fitAll .^ 2);
reconstructedEnergyFraction = reconstructedEnergy / totalEnergy;
residualRms = rmsLocal(residual);
observedRms = rmsLocal(y);
rmsReductionDb = 20 * log10(observedRms / max(residualRms, eps));
snrReconstructedVsResidualDb = 10 * log10(reconstructedEnergy / max(sse, eps));

writetable(componentTable, fullfile(resultDir, 'q3_multisource_results.csv'));
writetable(model.history, fullfile(resultDir, 'q3_omp_history.csv'));
writetable(model.snapshotTable, fullfile(resultDir, 'q3_residual_projection_spectra.csv'));

summaryTable = table(numel(model.frequencies_Hz), fs, duration, dfTheory, totalEnergy, ...
    reconstructedEnergy, reconstructedEnergyFraction, sse, rSquared, ...
    observedRms, residualRms, rmsReductionDb, snrReconstructedVsResidualDb, ...
    settings.thresholdMultiplier, settings.minPeakToMedianRatio, settings.minSeparationHz, ...
    'VariableNames', {'num_sources', 'fs_Hz', 'duration_s', 'theoretical_resolution_Hz', ...
    'total_energy', 'reconstructed_energy', 'reconstructed_energy_fraction', ...
    'residual_energy', 'R_squared', 'observed_RMS', 'residual_RMS', ...
    'RMS_reduction_dB', 'SNR_reconstructed_vs_residual_dB', ...
    'threshold_multiplier', 'min_peak_to_median_ratio', 'min_separation_Hz'});
writetable(summaryTable, fullfile(resultDir, 'q3_summary_metrics.csv'));

writeSummaryText(fullfile(resultDir, 'q3_summary.txt'), componentTable, summaryTable, model);

makeOriginalSpectrumFigure(figureDir, model);
makeResidualProjectionFigure(figureDir, model, model.frequencies_Hz);
makeTimeRecoveryFigure(figureDir, t, y, fitAll, residual);
makeComponentsFigure(figureDir, t, componentMatrix, model.frequencies_Hz);
makeFrequencyComparisonFigure(figureDir, y, fitAll, residual, fs, model.frequencies_Hz);

fprintf('Running threshold sensitivity check...\n');
thresholdSensitivity = runThresholdSensitivity(t, y, settings);
writetable(thresholdSensitivity, fullfile(resultDir, 'q3_threshold_sensitivity.csv'));

fprintf('Running near-frequency resolution simulation...\n');
[simLongTable, simMatrixTable, successMatrix, deltaList, snrList] = runResolutionSimulation(t, settings);
writetable(simLongTable, fullfile(resultDir, 'q3_simulation_resolution_results.csv'));
writetable(simMatrixTable, fullfile(resultDir, 'q3_resolution_success_matrix.csv'));
makeSimulationHeatmap(figureDir, successMatrix, deltaList, snrList);

fprintf('Question 3 finished. Outputs saved in %s\n', q3Dir);

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

function settings = defaultOmpSettings(fs, duration)
    settings.fMinHz = 0.05;
    settings.fMaxHz = min(49.5, fs / 2 - 0.1);
    settings.maxSources = 10;
    settings.thresholdMultiplier = 18;
    settings.minPeakToMedianRatio = 20;
    settings.minSeparationHz = max(0.01, 4 / duration);
    settings.refineHalfWidthHz = max(0.015, 8 / duration);
    settings.coordinateRefineHalfWidthHz = max(0.006, 3 / duration);
    settings.coordinateRefineSweeps = 2;
    settings.snapshotMaxHz = 20;
end

function model = runProjectionOmp(t, y, settings, keepSnapshots)
    t = t(:);
    y = y(:) - mean(y, 'omitnan');
    fs = 1 / median(diff(t));
    totalEnergy = sum(y .^ 2);
    selectedFreqs = zeros(0, 1);
    residual = y;
    history = emptyHistoryTable();
    snapshotTable = emptySnapshotTable();
    stopReason = "Reached maximum source count";

    for iter = 1:settings.maxSources
        [gridFreq, gridEnergy] = gridProjectionSpectrum(residual, fs);
        valid = gridFreq >= settings.fMinHz & gridFreq <= settings.fMaxHz;
        for j = 1:numel(selectedFreqs)
            valid = valid & abs(gridFreq - selectedFreqs(j)) > settings.minSeparationHz;
        end

        validEnergy = gridEnergy(valid);
        validFreq = gridFreq(valid);
        if isempty(validEnergy)
            stopReason = "No valid frequency bins remain";
            break;
        end

        [noiseMedian, noiseMadSigma, threshold] = robustNoiseThreshold(validEnergy, settings.thresholdMultiplier);
        [peakEnergy, peakIdx] = max(validEnergy);
        coarseFreq = validFreq(peakIdx);
        peakToMedian = peakEnergy / max(noiseMedian, eps);

        if keepSnapshots
            snapshotTable = appendSnapshot(snapshotTable, iter, gridFreq, gridEnergy, ...
                threshold, noiseMedian, settings.snapshotMaxHz);
        end

        accept = peakEnergy >= threshold && peakToMedian >= settings.minPeakToMedianRatio;
        if ~accept
            stopReason = sprintf('Peak %.6g below robust threshold %.6g or ratio %.3f below %.3f', ...
                peakEnergy, threshold, peakToMedian, settings.minPeakToMedianRatio);
            break;
        end

        refinedFreq = refineFrequency(t, residual, coarseFreq, settings.refineHalfWidthHz, settings);
        selectedFreqs = sort([selectedFreqs; refinedFreq]);
        selectedFreqs = coordinateRefineFrequencies(t, y, selectedFreqs, settings);
        [~, fit, residual, ~] = fitMultiSinusoid(t, y, selectedFreqs);

        history = appendHistory(history, iter, coarseFreq, refinedFreq, peakEnergy, ...
            noiseMedian, noiseMadSigma, threshold, peakToMedian, numel(selectedFreqs), ...
            sum(fit .^ 2) / totalEnergy, rmsLocal(residual));
    end

    model.frequencies_Hz = sort(selectedFreqs(:));
    model.history = history;
    model.snapshotTable = snapshotTable;
    model.stopReason = stopReason;
end

function history = emptyHistoryTable()
    history = table('Size', [0, 11], ...
        'VariableTypes', {'double','double','double','double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'iteration','coarse_frequency_Hz','refined_frequency_Hz', ...
        'peak_grid_projection_energy','noise_median','noise_MAD_sigma', ...
        'robust_threshold','peak_to_median_ratio','num_sources_after', ...
        'cumulative_energy_fraction','residual_RMS_after'});
end

function snapshotTable = emptySnapshotTable()
    snapshotTable = table('Size', [0, 5], ...
        'VariableTypes', {'double','double','double','double','double'}, ...
        'VariableNames', {'iteration','frequency_Hz','projection_energy', ...
        'robust_threshold','noise_median'});
end

function history = appendHistory(history, iteration, coarseFreq, refinedFreq, peakEnergy, ...
    noiseMedian, noiseMadSigma, threshold, peakToMedian, numSourcesAfter, energyFraction, residualRms)
    row = table(iteration, coarseFreq, refinedFreq, peakEnergy, noiseMedian, noiseMadSigma, ...
        threshold, peakToMedian, numSourcesAfter, energyFraction, residualRms, ...
        'VariableNames', history.Properties.VariableNames);
    history = [history; row]; %#ok<AGROW>
end

function snapshotTable = appendSnapshot(snapshotTable, iteration, freq, energy, threshold, noiseMedian, maxHz)
    keep = freq <= maxHz;
    freq = freq(keep);
    energy = energy(keep);
    n = numel(freq);
    row = table(repmat(iteration, n, 1), freq(:), energy(:), ...
        repmat(threshold, n, 1), repmat(noiseMedian, n, 1), ...
        'VariableNames', snapshotTable.Properties.VariableNames);
    snapshotTable = [snapshotTable; row]; %#ok<AGROW>
end

function [freq, energy] = gridProjectionSpectrum(x, fs)
    x = x(:) - mean(x, 'omitnan');
    n = numel(x);
    X = fft(x);
    freq = (0:floor(n / 2))' * fs / n;
    energy = 2 / n * abs(X(1:numel(freq))).^2;
    energy(1) = 0;
end

function [noiseMedian, noiseMadSigma, threshold] = robustNoiseThreshold(values, multiplier)
    values = values(:);
    noiseMedian = median(values, 'omitnan');
    madRaw = median(abs(values - noiseMedian), 'omitnan');
    noiseMadSigma = 1.4826 * madRaw;
    threshold = noiseMedian + multiplier * max(noiseMadSigma, eps);
end

function fHat = refineFrequency(t, x, coarseFreq, halfWidth, settings)
    lower = max(settings.fMinHz, coarseFreq - halfWidth);
    upper = min(settings.fMaxHz, coarseFreq + halfWidth);
    if upper <= lower
        fHat = coarseFreq;
        return;
    end
    fHat = maximizeProjectionOnInterval(t, x, lower, upper, 161);
end

function freqs = coordinateRefineFrequencies(t, y, freqs, settings)
    freqs = sort(freqs(:));
    for sweep = 1:settings.coordinateRefineSweeps %#ok<NASGU>
        [~, ~, ~, componentMatrix] = fitMultiSinusoid(t, y, freqs);
        for k = 1:numel(freqs)
            if numel(freqs) == 1
                partial = y;
            else
                others = setdiff(1:numel(freqs), k);
                partial = y - sum(componentMatrix(:, others), 2);
            end
            lower = max(settings.fMinHz, freqs(k) - settings.coordinateRefineHalfWidthHz);
            upper = min(settings.fMaxHz, freqs(k) + settings.coordinateRefineHalfWidthHz);
            if k > 1
                lower = max(lower, freqs(k - 1) + settings.minSeparationHz);
            end
            if k < numel(freqs)
                upper = min(upper, freqs(k + 1) - settings.minSeparationHz);
            end
            if upper > lower
                freqs(k) = maximizeProjectionOnInterval(t, partial, lower, upper, 101);
            end
        end
        freqs = sort(freqs(:));
    end
end

function fHat = maximizeProjectionOnInterval(t, x, lower, upper, gridCount)
    if upper <= lower
        fHat = 0.5 * (lower + upper);
        return;
    end
    fGrid = linspace(lower, upper, gridCount)';
    energy = zeros(size(fGrid));
    for i = 1:numel(fGrid)
        energy(i) = projectionEnergyAtFrequency(t, x, fGrid(i));
    end
    [~, idx] = max(energy);
    coarse = fGrid(idx);
    if numel(fGrid) > 1
        step = fGrid(2) - fGrid(1);
    else
        step = upper - lower;
    end
    fineLower = max(lower, coarse - step);
    fineUpper = min(upper, coarse + step);
    if fineUpper <= fineLower
        fHat = coarse;
        return;
    end
    objective = @(f) -projectionEnergyAtFrequency(t, x, f);
    opts = optimset('TolX', 1e-12, 'Display', 'off');
    fHat = fminbnd(objective, fineLower, fineUpper, opts);
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

function [theta, fit, residual, componentMatrix] = fitMultiSinusoid(t, y, freqs)
    t = t(:);
    y = y(:) - mean(y, 'omitnan');
    freqs = sort(freqs(:));
    n = numel(y);
    kNum = numel(freqs);
    if kNum == 0
        theta = zeros(0, 1);
        fit = zeros(n, 1);
        residual = y;
        componentMatrix = zeros(n, 0);
        return;
    end
    H = zeros(n, 2 * kNum);
    for k = 1:kNum
        w = 2 * pi * freqs(k);
        H(:, 2*k - 1) = cos(w * t);
        H(:, 2*k) = sin(w * t);
    end
    theta = H \ y;
    fit = H * theta;
    residual = y - fit;
    componentMatrix = zeros(n, kNum);
    for k = 1:kNum
        cols = 2*k - 1 : 2*k;
        componentMatrix(:, k) = H(:, cols) * theta(cols);
    end
end

function componentTable = buildComponentTable(freqs, theta, componentMatrix, totalEnergy)
    kNum = numel(freqs);
    source_id = (1:kNum)';
    frequency_Hz = freqs(:);
    cos_coefficient_a = zeros(kNum, 1);
    sin_coefficient_b = zeros(kNum, 1);
    A_hat = zeros(kNum, 1);
    phi_hat_rad = zeros(kNum, 1);
    component_energy = zeros(kNum, 1);
    component_energy_fraction = zeros(kNum, 1);
    for k = 1:kNum
        a = theta(2*k - 1);
        b = theta(2*k);
        cos_coefficient_a(k) = a;
        sin_coefficient_b(k) = b;
        A_hat(k) = sqrt(a^2 + b^2);
        phi_hat_rad(k) = atan2(a, b);
        component_energy(k) = sum(componentMatrix(:, k) .^ 2);
        component_energy_fraction(k) = component_energy(k) / totalEnergy;
    end
    componentTable = table(source_id, frequency_Hz, A_hat, phi_hat_rad, ...
        cos_coefficient_a, sin_coefficient_b, component_energy, component_energy_fraction);
end

function writeSummaryText(summaryFile, componentTable, summaryTable, model)
    fid = fopen(summaryFile, 'w', 'n', 'UTF-8');
    fprintf(fid, 'Question 3: multi-source weak periodic signal separation\n');
    fprintf(fid, 'Data: data.xlsx, second worksheet\n');
    fprintf(fid, 'Method: full-length sinusoidal projection OMP with robust residual thresholding\n\n');
    fprintf(fid, 'Automatically selected source count K: %d\n', summaryTable.num_sources(1));
    fprintf(fid, 'Stop reason: %s\n\n', model.stopReason);
    fprintf(fid, 'Detected components sorted by frequency:\n');
    for k = 1:height(componentTable)
        fprintf(fid, ['  source %d: f = %.12f Hz, A = %.12g, phi = %.12g rad, ', ...
            'energy fraction = %.12g\n'], componentTable.source_id(k), ...
            componentTable.frequency_Hz(k), componentTable.A_hat(k), ...
            componentTable.phi_hat_rad(k), componentTable.component_energy_fraction(k));
    end
    fprintf(fid, '\nOverall reconstruction metrics on observed real data:\n');
    fprintf(fid, '  reconstructed energy fraction: %.12g\n', summaryTable.reconstructed_energy_fraction(1));
    fprintf(fid, '  R_squared: %.12g\n', summaryTable.R_squared(1));
    fprintf(fid, '  residual RMS: %.12g\n', summaryTable.residual_RMS(1));
    fprintf(fid, '  RMS reduction: %.12g dB\n', summaryTable.RMS_reduction_dB(1));
    fprintf(fid, '  SNR reconstructed vs residual: %.12g dB\n\n', summaryTable.SNR_reconstructed_vs_residual_dB(1));
    fprintf(fid, 'Threshold settings:\n');
    fprintf(fid, '  robust threshold = median + %.6g * MAD_sigma\n', summaryTable.threshold_multiplier(1));
    fprintf(fid, '  minimum peak-to-median ratio = %.6g\n', summaryTable.min_peak_to_median_ratio(1));
    fprintf(fid, '  minimum frequency separation = %.6g Hz\n\n', summaryTable.min_separation_Hz(1));
    fprintf(fid, ['Note: only a single channel is provided in Question 3. Therefore, localization is ', ...
        'interpreted as frequency localization of fault sources, not physical spatial localization.\n']);
    fclose(fid);
end

function makeOriginalSpectrumFigure(outDir, model)
    rows = model.snapshotTable.iteration == 1;
    freq = model.snapshotTable.frequency_Hz(rows);
    energy = model.snapshotTable.projection_energy(rows);
    freqs = model.frequencies_Hz(:);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 620]);
    ax = axes(fig);
    semilogy(ax, freq, energy + eps, 'Color', [0.18, 0.20, 0.22], 'LineWidth', 0.95);
    hold on;
    selectedEnergy = interpolateSpectrumValues(freq, energy, freqs);
    colors = componentColors(numel(freqs));
    for k = 1:numel(freqs)
        plot(ax, freqs(k), selectedEnergy(k) + eps, 'o', ...
            'MarkerSize', 5.5, 'MarkerFaceColor', colors(k, :), ...
            'MarkerEdgeColor', 'w', 'LineWidth', 0.7, 'HandleVisibility', 'off');
    end
    xlim(ax, [0, 20]);
    setProjectionAxisLimits(ax, energy, []);
    addTopFrequencyMarkers(ax, freqs, ...
        arrayfun(@(f) sprintf('%.4f Hz', f), freqs(:), 'UniformOutput', false), colors);
    grid(ax, 'on');
    xlabel(ax, 'Frequency / Hz');
    ylabel(ax, 'Projection energy J_1(f)');
    title(ax, 'Full-length projection energy spectrum and identified frequencies');
    legend(ax, {'Projection energy'}, 'Location', 'northeast');
    exportgraphics(fig, fullfile(outDir, 'fig_q3_original_spectrum_identified.png'), 'Resolution', 220);
    close(fig);
end

function makeResidualProjectionFigure(outDir, model, freqs)
    if isempty(model.snapshotTable)
        return;
    end
    iterations = unique(model.snapshotTable.iteration);
    numPanels = numel(iterations);
    figHeight = max(560, 205 * numPanels);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1120, figHeight]);
    tiledlayout(numPanels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    colors = componentColors(numel(freqs));
    for p = 1:numPanels
        iter = iterations(p);
        rows = model.snapshotTable.iteration == iter;
        f = model.snapshotTable.frequency_Hz(rows);
        e = model.snapshotTable.projection_energy(rows);
        threshold = model.snapshotTable.robust_threshold(find(rows, 1));
        ax = nexttile;
        semilogy(ax, f, e + eps, 'Color', [0.18, 0.20, 0.22], 'LineWidth', 0.82);
        hold on;
        yline(ax, threshold + eps, '--', 'Color', [0.92, 0.55, 0.55], ...
            'LineWidth', 0.8);
        xlim(ax, [0, 20]);
        setProjectionAxisLimits(ax, e, threshold);
        grid(ax, 'on');
        ylabel(ax, 'J_m(f)');
        if p <= height(model.history)
            historyRow = model.history(model.history.iteration == iter, :);
            if ~isempty(historyRow)
                selectedFreq = historyRow.refined_frequency_Hz(1);
                selectedEnergy = interpolateSpectrumValues(f, e, selectedFreq);
                colorIdx = nearestFrequencyIndex(freqs, selectedFreq);
                plot(ax, selectedFreq, selectedEnergy + eps, 'o', ...
                    'MarkerSize', 5.8, 'MarkerFaceColor', colors(colorIdx, :), ...
                    'MarkerEdgeColor', 'w', 'LineWidth', 0.7, 'HandleVisibility', 'off');
                addTopFrequencyMarkers(ax, selectedFreq, ...
                    {sprintf('%.4f Hz', selectedFreq)}, colors(colorIdx, :), 'horizontal');
                title(ax, sprintf('Residual projection spectrum before selecting source %d', iter));
            else
                title(ax, sprintf('Residual projection spectrum at iteration %d', iter));
            end
        else
            title(ax, 'Residual projection spectrum at stopping step');
        end
        if p == 1
            legend(ax, {'Residual projection energy', 'Robust threshold'}, 'Location', 'northeast');
        end
        if p == numPanels
            xlabel(ax, 'Frequency / Hz');
        end
    end
    exportgraphics(fig, fullfile(outDir, 'fig_q3_omp_residual_projection_spectra.png'), 'Resolution', 220);
    close(fig);
end

function addTopFrequencyMarkers(ax, freqs, labels, colors, labelMode)
    if isempty(freqs)
        return;
    end
    if nargin < 5
        labelMode = 'vertical';
    end
    freqs = freqs(:);
    if size(colors, 1) == 1 && numel(freqs) > 1
        colors = repmat(colors, numel(freqs), 1);
    end
    yl = ylim(ax);
    if any(~isfinite(yl)) || yl(1) <= 0 || yl(2) <= yl(1)
        return;
    end
    logRange = log10(yl(2)) - log10(yl(1));
    markerTop = 10 ^ (log10(yl(2)) - 0.035 * logRange);
    markerBottom = 10 ^ (log10(yl(2)) - 0.075 * logRange);
    textY = 10 ^ (log10(yl(2)) - 0.105 * logRange);
    for k = 1:numel(freqs)
        line(ax, [freqs(k), freqs(k)], [markerBottom, markerTop], ...
            'Color', colors(k, :), 'LineWidth', 1.8, 'HandleVisibility', 'off');
        plot(ax, freqs(k), markerTop, 'v', 'MarkerSize', 5.2, ...
            'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', colors(k, :), ...
            'HandleVisibility', 'off');
        if strcmpi(labelMode, 'horizontal')
            text(ax, freqs(k), textY, labels{k}, 'Color', colors(k, :), ...
                'FontSize', 8.2, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', 'Margin', 1.0, 'Clipping', 'on');
        else
            text(ax, freqs(k), textY, labels{k}, 'Color', colors(k, :), ...
                'FontSize', 8.2, 'FontWeight', 'bold', 'Rotation', 90, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'Clipping', 'on');
        end
    end
end

function setProjectionAxisLimits(ax, energy, threshold)
    values = energy(:);
    values = values(isfinite(values) & values > 0);
    if nargin >= 3 && ~isempty(threshold) && isfinite(threshold) && threshold > 0
        values = [values; threshold];
    end
    if isempty(values)
        return;
    end
    lower = max(percentileLocal(values, 0.5) * 0.55, min(values) * 0.85);
    upper = max(values) * 5.0;
    if ~isfinite(lower) || lower <= 0 || lower >= upper
        lower = max(upper * 1e-6, eps);
    end
    ylim(ax, [lower, upper]);
end

function valuesAtTargets = interpolateSpectrumValues(freq, energy, targets)
    valuesAtTargets = zeros(size(targets));
    for i = 1:numel(targets)
        [~, idx] = min(abs(freq(:) - targets(i)));
        valuesAtTargets(i) = energy(idx);
    end
end

function idx = nearestFrequencyIndex(freqs, target)
    [~, idx] = min(abs(freqs(:) - target));
end

function q = percentileLocal(x, p)
    x = sort(x(:));
    if isempty(x)
        q = NaN;
        return;
    end
    pos = 1 + (numel(x) - 1) * p / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (x(hi) - x(lo)) * (pos - lo);
    end
end

function makeTimeRecoveryFigure(outDir, t, y, fitAll, residual)
    tSpan = 4.0;
    idx = t <= t(1) + tSpan;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 700]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(t(idx), y(idx), 'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.8);
    hold on;
    plot(t(idx), fitAll(idx), 'r-', 'LineWidth', 1.4);
    grid on;
    xlabel('Time / s');
    ylabel('Acceleration');
    title('Observed multi-source signal and reconstructed periodic components');
    legend({'Observed signal', 'Multi-frequency reconstruction'}, 'Location', 'best');
    nexttile;
    plot(t(idx), residual(idx), 'Color', [0.1, 0.35, 0.68], 'LineWidth', 0.8);
    grid on;
    xlabel('Time / s');
    ylabel('Residual');
    title('Residual after removing reconstructed components');
    exportgraphics(fig, fullfile(outDir, 'fig_q3_time_reconstruction_short.png'), 'Resolution', 220);
    close(fig);
end

function makeComponentsFigure(outDir, t, componentMatrix, freqs)
    if isempty(componentMatrix)
        return;
    end
    tSpan = 4.0;
    idx = t <= t(1) + tSpan;
    kNum = numel(freqs);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 80, 1120, max(560, 165*kNum)]);
    tiledlayout(kNum, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    colors = componentColors(kNum);
    for k = 1:kNum
        nexttile;
        plot(t(idx), componentMatrix(idx, k), 'Color', colors(k, :), 'LineWidth', 1.1);
        grid on;
        ylabel(sprintf('s_%d(t)', k));
        title(sprintf('Separated component %d: f = %.6f Hz', k, freqs(k)));
        if k == kNum
            xlabel('Time / s');
        end
    end
    exportgraphics(fig, fullfile(outDir, 'fig_q3_separated_components.png'), 'Resolution', 220);
    close(fig);
end

function makeFrequencyComparisonFigure(outDir, y, fitAll, residual, fs, freqs)
    [freq, ampY] = singleSidedAmplitude(y, fs);
    [~, ampFit] = singleSidedAmplitude(fitAll, fs);
    [~, ampResidual] = singleSidedAmplitude(residual, fs);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 660]);
    semilogy(freq, ampY + eps, 'Color', [0.50, 0.50, 0.50], 'LineWidth', 0.9);
    hold on;
    semilogy(freq, ampFit + eps, 'r-', 'LineWidth', 1.2);
    semilogy(freq, ampResidual + eps, 'Color', [0.1, 0.35, 0.72], 'LineWidth', 0.9);
    colors = componentColors(numel(freqs));
    for k = 1:numel(freqs)
        xline(freqs(k), '--', 'Color', colors(k, :), 'LineWidth', 1.1);
    end
    xlim([0, 20]);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Single-sided amplitude');
    title('Frequency-domain comparison before and after multi-source separation');
    legend({'Observed signal', 'Reconstructed components', 'Residual'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outDir, 'fig_q3_frequency_before_after.png'), 'Resolution', 220);
    close(fig);
end

function sensitivity = runThresholdSensitivity(t, y, baseSettings)
    multipliers = [12; 14; 16; 18; 20; 24];
    num_sources = zeros(numel(multipliers), 1);
    frequency_list_Hz = strings(numel(multipliers), 1);
    reconstructed_energy_fraction = zeros(numel(multipliers), 1);
    residual_RMS = zeros(numel(multipliers), 1);
    totalEnergy = sum((y - mean(y, 'omitnan')) .^ 2);
    for i = 1:numel(multipliers)
        settings = baseSettings;
        settings.thresholdMultiplier = multipliers(i);
        model = runProjectionOmp(t, y, settings, false);
        [~, fit, residual, ~] = fitMultiSinusoid(t, y, model.frequencies_Hz);
        num_sources(i) = numel(model.frequencies_Hz);
        frequency_list_Hz(i) = join(string(round(model.frequencies_Hz(:), 9)), ", ");
        reconstructed_energy_fraction(i) = sum(fit .^ 2) / totalEnergy;
        residual_RMS(i) = rmsLocal(residual);
    end
    sensitivity = table(multipliers, num_sources, frequency_list_Hz, ...
        reconstructed_energy_fraction, residual_RMS, ...
        'VariableNames', {'threshold_multiplier','num_sources','frequency_list_Hz', ...
        'reconstructed_energy_fraction','residual_RMS'});
end

function [longTable, matrixTable, successMatrix, deltaList, snrList] = runResolutionSimulation(t, baseSettings)
    rng(20260620);
    t = t(:);
    deltaList = [0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0050, 0.0100, 0.0200]';
    snrList = [-40, -36, -32, -28, -24]';
    trials = 4;
    baseFreq = 10.0;
    A1 = 0.030;
    A2 = 0.025;
    phi1 = 0.40;
    phi2 = -1.10;
    simSettings = baseSettings;
    simSettings.maxSources = 2;
    simSettings.minSeparationHz = 0.0002;
    simSettings.thresholdMultiplier = baseSettings.thresholdMultiplier;
    simSettings.minPeakToMedianRatio = baseSettings.minPeakToMedianRatio;
    simSettings.refineHalfWidthHz = max(0.006, baseSettings.refineHalfWidthHz);
    simSettings.coordinateRefineHalfWidthHz = max(0.003, baseSettings.coordinateRefineHalfWidthHz);

    totalRows = numel(deltaList) * numel(snrList);
    delta_f_Hz = zeros(totalRows, 1);
    SNR_dB = zeros(totalRows, 1);
    success_rate = zeros(totalRows, 1);
    mean_abs_frequency_error_Hz = zeros(totalRows, 1);
    row = 0;
    successMatrix = zeros(numel(snrList), numel(deltaList));

    for iS = 1:numel(snrList)
        for iD = 1:numel(deltaList)
            delta = deltaList(iD);
            truth = [baseFreq; baseFreq + delta];
            clean = A1 * sin(2*pi*truth(1)*t + phi1) + A2 * sin(2*pi*truth(2)*t + phi2);
            cleanPower = mean(clean .^ 2);
            noiseSigma = sqrt(cleanPower / (10^(snrList(iS) / 10)));
            okCount = 0;
            errSum = 0;
            errCount = 0;
            tolerance = max(0.00015, min(0.003, 0.30 * delta));
            for tr = 1:trials
                x = clean + noiseSigma * randn(size(t));
                x = x - mean(x, 'omitnan');
                simModel = runProjectionOmp(t, x, simSettings, false);
                [ok, err] = checkTwoFrequencySuccess(simModel.frequencies_Hz, truth, tolerance);
                okCount = okCount + ok;
                if ok
                    errSum = errSum + err;
                    errCount = errCount + 1;
                end
            end
            row = row + 1;
            delta_f_Hz(row) = delta;
            SNR_dB(row) = snrList(iS);
            success_rate(row) = okCount / trials;
            if errCount > 0
                mean_abs_frequency_error_Hz(row) = errSum / errCount;
            else
                mean_abs_frequency_error_Hz(row) = NaN;
            end
            successMatrix(iS, iD) = success_rate(row);
        end
    end

    longTable = table(delta_f_Hz, SNR_dB, success_rate, mean_abs_frequency_error_Hz);
    matrixTable = array2table(successMatrix, ...
        'VariableNames', matlab.lang.makeValidName("df_" + string(deltaList(:)')));
    matrixTable.SNR_dB = snrList;
    matrixTable = movevars(matrixTable, 'SNR_dB', 'Before', 1);
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

function makeSimulationHeatmap(outDir, successMatrix, deltaList, snrList)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 980, 620]);
    imagesc(deltaList, snrList, successMatrix);
    set(gca, 'YDir', 'normal');
    colormap(parula);
    c = colorbar;
    c.Label.String = 'Identification success rate';
    caxis([0, 1]);
    xlabel('Frequency spacing \Delta f / Hz');
    ylabel('SNR / dB');
    title('Near-frequency resolution experiment under full-duration sampling');
    xticks(deltaList);
    grid on;
    exportgraphics(fig, fullfile(outDir, 'fig_q3_near_frequency_resolution_heatmap.png'), 'Resolution', 220);
    close(fig);
end

function [freq, amp] = singleSidedAmplitude(x, fs)
    x = x(:) - mean(x, 'omitnan');
    n = numel(x);
    X = fft(x);
    freq = (0:floor(n / 2))' * fs / n;
    amp = abs(X(1:numel(freq))) / n * 2;
    amp(1) = amp(1) / 2;
end

function y = rmsLocal(x)
    y = sqrt(mean(x(:) .^ 2));
end

function colors = componentColors(kNum)
    base = [0.82, 0.18, 0.18;
            0.10, 0.42, 0.72;
            0.16, 0.58, 0.32;
            0.58, 0.28, 0.70;
            0.85, 0.48, 0.08;
            0.20, 0.55, 0.58;
            0.45, 0.45, 0.45;
            0.70, 0.20, 0.45];
    if kNum <= size(base, 1)
        colors = base(1:kNum, :);
    else
        colors = lines(kNum);
    end
end

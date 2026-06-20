%% Question 3: redraw key spectrum figures only
% This script refreshes two presentation figures from existing CSV results.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
q3Dir = fileparts(scriptDir);

mainDir = fullfile(q3Dir, '多故障源分离_run_q3_multisource_separation');
resultDir = mainDir;
figureDir = mainDir;

snapshotTable = readtable(fullfile(resultDir, 'q3_residual_projection_spectra.csv'), ...
    'VariableNamingRule', 'preserve');
history = readtable(fullfile(resultDir, 'q3_omp_history.csv'), ...
    'VariableNamingRule', 'preserve');
componentTable = readtable(fullfile(resultDir, 'q3_multisource_results.csv'), ...
    'VariableNamingRule', 'preserve');

model.snapshotTable = snapshotTable;
model.history = history;
model.frequencies_Hz = componentTable.frequency_Hz(:);

makeOriginalProjectionFigure(figureDir, model);
makeResidualProjectionFigure(figureDir, model, model.frequencies_Hz);

fprintf('Refreshed key Question 3 figures in %s\n', figureDir);

function makeOriginalProjectionFigure(outDir, model)
    rows = model.snapshotTable.iteration == 1;
    freq = model.snapshotTable.frequency_Hz(rows);
    energy = model.snapshotTable.projection_energy(rows);
    freqs = model.frequencies_Hz(:);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 620]);
    ax = axes(fig);
    semilogy(ax, freq, energy + eps, 'Color', [0.18, 0.20, 0.22], 'LineWidth', 0.95);
    hold on;

    colors = componentColors(numel(freqs));
    selectedEnergy = interpolateSpectrumValues(freq, energy, freqs);
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
        yline(ax, threshold + eps, '--', 'Color', [0.92, 0.55, 0.55], 'LineWidth', 0.8);
        xlim(ax, [0, 20]);
        setProjectionAxisLimits(ax, e, threshold);
        grid(ax, 'on');
        ylabel(ax, 'J_m(f)');

        if p <= height(model.history)
            historyRow = model.history(model.history.iteration == iter, :);
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

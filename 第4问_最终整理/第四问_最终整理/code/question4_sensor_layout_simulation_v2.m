%% Question 4 full simulation framework
% Purpose:
%   This script follows the six-step model framework more strictly than the
%   first simple validation script.
%
% What it validates:
%   1. Different candidate sensor structures: M=12, M=16, M=20.
%   2. Different numbers of potential fault sources: K=1, K=3, K=5.
%   3. Repeated random fault-source locations.
%   4. Sensitivity to alpha, beta, and omega.
%   5. Reliability under one-sensor failure.
%   6. Performance gain under different SNR levels.
%
% Notes:
%   - This is still a normalized gearbox simulation, not a real gearbox test.
%   - kappa and nu are AHP results or extensions by the same AHP mapping.
%   - SUBJECTIVE / BASELINE marks parameters that are modeling assumptions.

clear; clc; close all;
rng(2026);

%% =========================
% 1. Baseline settings
% =========================

% SUBJECTIVE / BASELINE: normalized attenuation parameters.
baseAlpha = 0.33;
baseBeta = 1.00;

% SUBJECTIVE / BASELINE: omega is not a physical constant. We test multiple
% values and use baseOmega only as the main reporting value.
baseOmega = 0.70;
omegaGrid = [0.3 0.5 0.7 0.9];

% SUBJECTIVE / BASELINE: false alarm constraint.
pfaMax = 0.01;

% SUBJECTIVE / BASELINE: signal and sampling settings.
fs = 100;
windowSeconds = 5;
L = fs * windowSeconds;
t = (0:L-1)' / fs;
f0 = 2.0;
sigma0 = 1.0;
baseAmp = 0.25 * sigma0;

% SUBJECTIVE / BASELINE: Monte Carlo counts. Increase for final smoother
% values. These defaults are chosen to keep runtime acceptable.
nMCDetect = 500;
nMCFalse = 3000;
nRandomSourceCases = 25;
nRandomLayouts = 40;

structureModes = {'M12_baseline', 'M16_extended', 'M20_dense'};
faultCounts = [1 3 5];
alphaGrid = [0.20 0.33 0.50];
betaGrid = [0.80 1.00 1.20];
snrScales = [0.25 0.50 1.00 1.50 2.00];

% SUBJECTIVE / BASELINE: perturbation levels for AHP-derived kappa and nu.
% These are used to test whether subjective AHP coefficients strongly affect
% the final layout. 0.10 means +/-10% multiplicative perturbation.
ahpPerturbLevel = 0.10;
nAHPParamCases = 40;

% SUBJECTIVE / BASELINE: source parameter perturbation for extra robustness.
% Frequency varies within +/-10%, amplitude varies within +/-20%.
nSourceParamCases = 25;
sourceFreqRelRange = 0.10;
sourceAmpRelRange = 0.20;

resultDir = fullfile(fileparts(mfilename('fullpath')), 'question4_simulation_v2_outputs');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

%% =========================
% 2. Main baseline structure and fixed sources
% =========================

sensors = make_candidate_sensors('M12_baseline');
sourcesFixed = make_fixed_sources(3, f0);
layouts = enumerate_layouts(numel(sensors), 3);

bestLayout = optimize_layout_by_lambda(sensors, sourcesFixed, layouts, baseAlpha, baseBeta, baseOmega, baseAmp, sigma0, L);
eta3 = chi2_threshold_even(2*numel(bestLayout), pfaMax);
eta1 = chi2_threshold_even(2, pfaMax);

fprintf('\n=== Baseline optimal layout ===\n');
disp(layout_to_text(bestLayout, sensors));

%% =========================
% 3. Baseline Monte Carlo detection and false alarm
% =========================

conLayout = concentrated_layout(sensors, layouts);
singleBest = best_single_sensor_by_lambda(sensors, sourcesFixed, baseAlpha, baseBeta, baseOmega, baseAmp, sigma0, L);
randLayouts = sample_random_layouts(layouts, nRandomLayouts);

bestSim = simulate_layout_pd(sensors, sourcesFixed, bestLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, nMCDetect);
conSim = simulate_layout_pd(sensors, sourcesFixed, conLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, nMCDetect);
singleSim = simulate_layout_pd(sensors, sourcesFixed, singleBest, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta1, nMCDetect);

randMean = zeros(numel(randLayouts), 1);
randMin = zeros(numel(randLayouts), 1);
for k = 1:numel(randLayouts)
    simR = simulate_layout_pd(sensors, sourcesFixed, randLayouts{k}, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, nMCDetect);
    randMean(k) = simR.meanPd;
    randMin(k) = simR.minPd;
end

pfaEmp = simulate_false_alarm(sensors, bestLayout, sigma0, f0, t, eta3, nMCFalse);

baselineTable = table( ...
    {'Optimal'; 'Concentrated'; 'Random average'; 'Single best'}, ...
    [bestSim.meanPd; conSim.meanPd; mean(randMean); singleSim.meanPd], ...
    [bestSim.minPd; conSim.minPd; mean(randMin); singleSim.minPd], ...
    [bestSim.meanPd-mean(randMean); conSim.meanPd-mean(randMean); 0; singleSim.meanPd-mean(randMean)], ...
    'VariableNames', {'Layout','MeanPd','MinPd','AbsGainVsRandom'});

disp('=== Baseline Monte Carlo result ===');
disp(baselineTable);
fprintf('Empirical Pfa of optimal layout = %.4f\n', pfaEmp);

%% =========================
% 4. Omega sensitivity
% =========================

omegaRows = table();
for w = omegaGrid
    layoutW = optimize_layout_by_lambda(sensors, sourcesFixed, layouts, baseAlpha, baseBeta, w, baseAmp, sigma0, L);
    simW = simulate_layout_pd(sensors, sourcesFixed, layoutW, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, nMCDetect);
    omegaRows = [omegaRows; table(w, {layout_ids(layoutW, sensors)}, simW.meanPd, simW.minPd, ...
        'VariableNames', {'Omega','BestLayout','MeanPd','MinPd'})]; %#ok<AGROW>
end

disp('=== Omega sensitivity ===');
disp(omegaRows);

%% =========================
% 5. Alpha/Beta sensitivity
% =========================

paramRows = table();
for a = alphaGrid
    for b = betaGrid
        layoutAB = optimize_layout_by_lambda(sensors, sourcesFixed, layouts, a, b, baseOmega, baseAmp, sigma0, L);
        simAB = simulate_layout_pd(sensors, sourcesFixed, layoutAB, a, b, baseAmp, sigma0, fs, t, eta3, nMCDetect);
        paramRows = [paramRows; table(a, b, {layout_ids(layoutAB, sensors)}, simAB.meanPd, simAB.minPd, ...
            'VariableNames', {'Alpha','Beta','BestLayout','MeanPd','MinPd'})]; %#ok<AGROW>
    end
end

disp('=== Alpha/Beta sensitivity ===');
disp(paramRows);

%% =========================
% 6. AHP parameter stability: kappa and nu perturbation
% =========================

ahpRows = table();
layoutTexts = cell(nAHPParamCases, 1);
meanVals = zeros(nAHPParamCases, 1);
minVals = zeros(nAHPParamCases, 1);
for r = 1:nAHPParamCases
    sensorsPert = perturb_ahp_parameters(sensors, ahpPerturbLevel);
    layoutPert = optimize_layout_by_lambda(sensorsPert, sourcesFixed, layouts, baseAlpha, baseBeta, baseOmega, baseAmp, sigma0, L);
    simPert = simulate_layout_pd(sensorsPert, sourcesFixed, layoutPert, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, max(200, round(nMCDetect/2)));
    layoutTexts{r} = layout_ids(layoutPert, sensorsPert);
    meanVals(r) = simPert.meanPd;
    minVals(r) = simPert.minPd;
end
ahpRows = [ahpRows; table(ahpPerturbLevel, nAHPParamCases, {most_frequent_layout(layoutTexts)}, ...
    mean(meanVals), min(meanVals), std(meanVals), mean(minVals), ...
    'VariableNames', {'PerturbLevel','CaseCount','MostFrequentLayout','MeanPdAvg','MeanPdMin','MeanPdStd','MinPdAvg'})];

disp('=== AHP parameter stability: kappa/nu perturbation ===');
disp(ahpRows);

%% =========================
% 7. Structure size and fault count sweep
% =========================

scenarioRows = table();
for m = 1:numel(structureModes)
    mode = structureModes{m};
    sensorsM = make_candidate_sensors(mode);
    layoutsM = enumerate_layouts(numel(sensorsM), 3);
    for fc = faultCounts
        meanPdCases = zeros(nRandomSourceCases, 1);
        minPdCases = zeros(nRandomSourceCases, 1);
        layoutText = cell(nRandomSourceCases, 1);
        for r = 1:nRandomSourceCases
            srcRand = make_random_sources(fc, f0);
            layoutR = optimize_layout_by_lambda(sensorsM, srcRand, layoutsM, baseAlpha, baseBeta, baseOmega, baseAmp, sigma0, L);
            etaR = chi2_threshold_even(2*numel(layoutR), pfaMax);
            simR = simulate_layout_pd(sensorsM, srcRand, layoutR, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, etaR, max(200, round(nMCDetect/2)));
            meanPdCases(r) = simR.meanPd;
            minPdCases(r) = simR.minPd;
            layoutText{r} = layout_ids(layoutR, sensorsM);
        end
        mostLayout = most_frequent_layout(layoutText);
        scenarioRows = [scenarioRows; table({mode}, numel(sensorsM), fc, {mostLayout}, ...
            mean(meanPdCases), min(meanPdCases), std(meanPdCases), mean(minPdCases), ...
            'VariableNames', {'StructureMode','CandidateCount','FaultCount','MostFrequentLayout','MeanPdAvg','MeanPdMin','MeanPdStd','MinPdAvg'})]; %#ok<AGROW>
    end
end

disp('=== Structure/fault-count sweep ===');
disp(scenarioRows);

%% =========================
% 8. Source-position robustness with fixed optimal layout
% =========================

robustRows = table();
for fc = faultCounts
    meanVals = zeros(nRandomSourceCases, 1);
    minVals = zeros(nRandomSourceCases, 1);
    for r = 1:nRandomSourceCases
        srcRand = make_random_sources(fc, f0);
        simFix = simulate_layout_pd(sensors, srcRand, bestLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, max(200, round(nMCDetect/2)));
        meanVals(r) = simFix.meanPd;
        minVals(r) = simFix.minPd;
    end
    robustRows = [robustRows; table(fc, mean(meanVals), min(meanVals), std(meanVals), mean(minVals), ...
        'VariableNames', {'FaultCount','MeanPdAvg','MeanPdMin','MeanPdStd','MinPdAvg'})]; %#ok<AGROW>
end

disp('=== Random source-position robustness of baseline optimal layout ===');
disp(robustRows);

%% =========================
% 9. Source frequency/amplitude robustness
% =========================

sourceParamRows = table();
bestParamMean = zeros(nSourceParamCases, 1);
randParamMean = zeros(nSourceParamCases, 1);
for r = 1:nSourceParamCases
    srcParam = perturb_source_parameters(sourcesFixed, sourceFreqRelRange, sourceAmpRelRange);
    % Keep the baseline optimal layout fixed. This tests whether the chosen
    % layout still works when fault frequencies and amplitudes are not identical.
    simBestParam = simulate_layout_pd(sensors, srcParam, bestLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, max(200, round(nMCDetect/2)));
    bestParamMean(r) = simBestParam.meanPd;

    tempRand = zeros(numel(randLayouts), 1);
    for k = 1:numel(randLayouts)
        simRandParam = simulate_layout_pd(sensors, srcParam, randLayouts{k}, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta3, max(150, round(nMCDetect/3)));
        tempRand(k) = simRandParam.meanPd;
    end
    randParamMean(r) = mean(tempRand);
end
sourceParamRows = [sourceParamRows; table(nSourceParamCases, sourceFreqRelRange, sourceAmpRelRange, ...
    mean(bestParamMean), min(bestParamMean), std(bestParamMean), mean(randParamMean), ...
    mean(bestParamMean-randParamMean), ...
    'VariableNames', {'CaseCount','FreqRelRange','AmpRelRange','OptimalMeanPdAvg','OptimalMeanPdMin','OptimalMeanPdStd','RandomMeanPdAvg','AbsGainVsRandomAvg'})];

disp('=== Source frequency/amplitude robustness ===');
disp(sourceParamRows);

%% =========================
% 10. Single-sensor failure reliability
% =========================

failureRows = table();
failureRows = [failureRows; make_failure_row('No failure', bestLayout, bestSim)];
for k = 1:numel(bestLayout)
    remain = bestLayout;
    failedId = sensors(bestLayout(k)).id;
    remain(k) = [];
    etaRemain = chi2_threshold_even(2*numel(remain), pfaMax);
    simRemain = simulate_layout_pd(sensors, sourcesFixed, remain, baseAlpha, baseBeta, baseAmp, sigma0, fs, t, etaRemain, nMCDetect);
    failureRows = [failureRows; make_failure_row(['Sensor ', failedId, ' failed'], remain, simRemain)]; %#ok<AGROW>
end

disp('=== Single-sensor failure reliability ===');
disp(failureRows);

%% =========================
% 11. SNR performance gain
% =========================

snrRows = table();
for s = snrScales
    ampS = baseAmp * sqrt(s);
    simBestS = simulate_layout_pd(sensors, sourcesFixed, bestLayout, baseAlpha, baseBeta, ampS, sigma0, fs, t, eta3, nMCDetect);
    simConS = simulate_layout_pd(sensors, sourcesFixed, conLayout, baseAlpha, baseBeta, ampS, sigma0, fs, t, eta3, nMCDetect);
    simSingleS = simulate_layout_pd(sensors, sourcesFixed, singleBest, baseAlpha, baseBeta, ampS, sigma0, fs, t, eta1, nMCDetect);
    randS = zeros(numel(randLayouts), 1);
    for k = 1:numel(randLayouts)
        simRS = simulate_layout_pd(sensors, sourcesFixed, randLayouts{k}, baseAlpha, baseBeta, ampS, sigma0, fs, t, eta3, max(200, round(nMCDetect/2)));
        randS(k) = simRS.meanPd;
    end
    randAvg = mean(randS);
    snrRows = [snrRows; table(s, simBestS.meanPd, randAvg, simConS.meanPd, simSingleS.meanPd, ...
        simBestS.meanPd-randAvg, (simBestS.meanPd-randAvg)/max(randAvg, eps), ...
        'VariableNames', {'SNRScale','OptimalMeanPd','RandomMeanPd','ConcentratedMeanPd','SingleMeanPd','AbsGainVsRandom','RelGainVsRandom'})]; %#ok<AGROW>
end

disp('=== SNR performance gain ===');
disp(snrRows);

%% =========================
% 12. Fusion spectrum visualization
% =========================

[fAxis, fusionBest] = one_fusion_spectrum(sensors, sourcesFixed(2), bestLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t);
[~, fusionCon] = one_fusion_spectrum(sensors, sourcesFixed(2), conLayout, baseAlpha, baseBeta, baseAmp, sigma0, fs, t);

fig1 = figure('Color','w');
plot(fAxis, fusionBest, 'LineWidth', 1.6); hold on;
plot(fAxis, fusionCon, '--', 'LineWidth', 1.4);
xlim([0 10]); grid on;
xlabel('Frequency / Hz');
ylabel('Weighted fusion amplitude');
legend('Optimal layout', 'Concentrated layout', 'Location', 'northeast');
title('Fusion spectrum comparison');
saveas(fig1, fullfile(resultDir, 'fusion_spectrum_comparison.png'));

fig2 = figure('Color','w');
plot(snrRows.SNRScale, snrRows.OptimalMeanPd, '-o', 'LineWidth', 1.6); hold on;
plot(snrRows.SNRScale, snrRows.RandomMeanPd, '-s', 'LineWidth', 1.3);
plot(snrRows.SNRScale, snrRows.ConcentratedMeanPd, '-^', 'LineWidth', 1.3);
plot(snrRows.SNRScale, snrRows.SingleMeanPd, '-d', 'LineWidth', 1.3);
grid on;
xlabel('SNR scale');
ylabel('Mean detection probability');
legend('Optimal', 'Random average', 'Concentrated', 'Single best', 'Location', 'southeast');
title('Performance under different SNR levels');
saveas(fig2, fullfile(resultDir, 'snr_performance_gain.png'));

fig3 = figure('Color','w');
scatter(robustRows.FaultCount, robustRows.MeanPdAvg, 80, 'filled'); hold on;
errorbar(robustRows.FaultCount, robustRows.MeanPdAvg, robustRows.MeanPdStd, 'LineStyle', 'none', 'LineWidth', 1.2);
grid on;
xlabel('Number of fault sources');
ylabel('Mean detection probability');
title('Robustness under random fault-source positions');
saveas(fig3, fullfile(resultDir, 'random_source_robustness.png'));

fig4 = figure('Color','w');
plot(omegaRows.Omega, omegaRows.MeanPd, '-o', 'LineWidth', 1.6); hold on;
plot(omegaRows.Omega, omegaRows.MinPd, '-s', 'LineWidth', 1.4);
grid on;
xlabel('\omega');
ylabel('Detection probability');
legend('Mean Pd', 'Min Pd', 'Location', 'best');
title('Omega sensitivity');
saveas(fig4, fullfile(resultDir, 'omega_sensitivity.png'));

fig5 = figure('Color','w');
histogram(meanVals, 10);
grid on;
xlabel('Mean detection probability');
ylabel('Count');
title('\kappa and \nu perturbation stability');
saveas(fig5, fullfile(resultDir, 'ahp_parameter_stability.png'));

fig6 = figure('Color','w');
bar(categorical({'Optimal fixed layout','Random average'}), [mean(bestParamMean), mean(randParamMean)]);
grid on;
ylabel('Mean detection probability');
title('Robustness to source frequency/amplitude perturbation');
saveas(fig6, fullfile(resultDir, 'source_parameter_robustness.png'));

%% =========================
% 13. Save outputs
% =========================

writetable(baselineTable, fullfile(resultDir, 'baseline_detection_table.csv'));
writetable(omegaRows, fullfile(resultDir, 'omega_sensitivity_table.csv'));
writetable(paramRows, fullfile(resultDir, 'alpha_beta_sensitivity_table.csv'));
writetable(ahpRows, fullfile(resultDir, 'ahp_parameter_stability_table.csv'));
writetable(scenarioRows, fullfile(resultDir, 'structure_faultcount_sweep_table.csv'));
writetable(robustRows, fullfile(resultDir, 'random_source_robustness_table.csv'));
writetable(sourceParamRows, fullfile(resultDir, 'source_parameter_robustness_table.csv'));
writetable(failureRows, fullfile(resultDir, 'failure_reliability_table.csv'));
writetable(snrRows, fullfile(resultDir, 'snr_performance_table.csv'));

save(fullfile(resultDir, 'question4_simulation_v2_results.mat'), ...
    'sensors', 'sourcesFixed', 'bestLayout', 'conLayout', 'singleBest', ...
    'baselineTable', 'omegaRows', 'paramRows', 'ahpRows', 'scenarioRows', 'robustRows', 'sourceParamRows', ...
    'failureRows', 'snrRows', 'pfaEmp', 'baseAlpha', 'baseBeta', 'baseOmega', ...
    'pfaMax', 'fs', 'windowSeconds', 'baseAmp', 'sigma0', 'ahpPerturbLevel', 'nAHPParamCases', ...
    'nSourceParamCases', 'sourceFreqRelRange', 'sourceAmpRelRange');

fprintf('\nAll v2 simulation outputs saved to:\n%s\n', resultDir);

%% ========================================================================
% Local functions
% ========================================================================

function sensors = make_candidate_sensors(mode)
    base = repmat(make_sensor('', '', [0 0 0], 0, 0), 1, 20);
    base(1)  = make_sensor('P1',  'Input shaft left bearing seat',   [0.20 0.25 0.35], 1.0000, 1.3376);
    base(2)  = make_sensor('P2',  'Input shaft right bearing seat',  [0.20 0.75 0.35], 1.0000, 1.3376);
    base(3)  = make_sensor('P3',  'Middle shaft left bearing seat',  [0.50 0.25 0.35], 1.0000, 1.3376);
    base(4)  = make_sensor('P4',  'Middle shaft right bearing seat', [0.50 0.75 0.35], 1.0000, 1.3376);
    base(5)  = make_sensor('P5',  'Output shaft left bearing seat',  [0.80 0.25 0.35], 1.0000, 1.3376);
    base(6)  = make_sensor('P6',  'Output shaft right bearing seat', [0.80 0.75 0.35], 1.0000, 1.3376);
    base(7)  = make_sensor('P7',  'Front gearbox shell',             [0.50 0.00 0.55], 0.6500, 1.0211);
    base(8)  = make_sensor('P8',  'Back gearbox shell',              [0.50 1.00 0.55], 0.6500, 1.0211);
    base(9)  = make_sensor('P9',  'Top gearbox shell',               [0.50 0.50 1.00], 0.6500, 0.8000);
    base(10) = make_sensor('P10', 'Bottom gearbox shell',            [0.50 0.50 0.00], 0.6500, 1.3493);
    base(11) = make_sensor('P11', 'Motor input end',                 [0.00 0.50 0.45], 0.7420, 1.5000);
    base(12) = make_sensor('P12', 'Load output end',                 [1.00 0.50 0.45], 0.7420, 1.1835);

    % SUBJECTIVE / BASELINE extensions: extra shell/stiffener points used
    % only to test whether the method works when candidate count changes.
    base(13) = make_sensor('P13', 'Input upper shell auxiliary point',  [0.20 0.50 0.80], 0.6500, 0.9000);
    base(14) = make_sensor('P14', 'Output upper shell auxiliary point', [0.80 0.50 0.80], 0.6500, 0.9000);
    base(15) = make_sensor('P15', 'Input lower shell auxiliary point',  [0.20 0.50 0.08], 0.6500, 1.3000);
    base(16) = make_sensor('P16', 'Output lower shell auxiliary point', [0.80 0.50 0.08], 0.6500, 1.3000);
    base(17) = make_sensor('P17', 'Front-left shell stiffener',         [0.35 0.00 0.35], 0.7420, 1.0800);
    base(18) = make_sensor('P18', 'Back-right shell stiffener',         [0.65 1.00 0.35], 0.7420, 1.0800);
    base(19) = make_sensor('P19', 'Input side cover',                  [0.10 0.50 0.30], 0.7420, 1.3000);
    base(20) = make_sensor('P20', 'Output side cover',                 [0.90 0.50 0.30], 0.7420, 1.1500);

    if strcmp(mode, 'M12_baseline')
        sensors = base(1:12);
    elseif strcmp(mode, 'M16_extended')
        sensors = base(1:16);
    elseif strcmp(mode, 'M20_dense')
        sensors = base(1:20);
    else
        error('Unknown structure mode.');
    end
end

function sensorsPert = perturb_ahp_parameters(sensors, level)
    sensorsPert = sensors;
    for i = 1:numel(sensors)
        % Multiplicative perturbation around AHP baseline.
        kFactor = 1 + level * (2*rand() - 1);
        vFactor = 1 + level * (2*rand() - 1);
        sensorsPert(i).kappa = min(max(sensors(i).kappa * kFactor, 0.05), 1.20);
        sensorsPert(i).nu = max(sensors(i).nu * vFactor, 0.10);
    end
end

function s = make_sensor(id, name, point, kappa, nu)
    s.id = id;
    s.name = name;
    s.point = point;
    s.kappa = kappa;
    s.nu = nu;
end

function sources = make_fixed_sources(k, f0)
    if k == 1
        points = [0.50 0.50 0.50];
    elseif k == 3
        points = [0.25 0.50 0.50; 0.50 0.50 0.50; 0.75 0.50 0.50];
    elseif k == 5
        points = [0.20 0.50 0.50; 0.35 0.50 0.50; 0.50 0.50 0.50; 0.65 0.50 0.50; 0.80 0.50 0.50];
    else
        error('Only K=1,3,5 fixed sources are supported.');
    end
    sources = repmat(make_source('', '', [0 0 0], f0, 1.0), 1, size(points,1));
    for i = 1:size(points,1)
        sources(i) = make_source(['R', num2str(i)], ['Fault source ', num2str(i)], points(i,:), f0, 1.0);
    end
end

function sources = make_random_sources(k, f0)
    % SUBJECTIVE / BASELINE: random fault sources are sampled from the
    % internal gear transmission region, not the whole gearbox shell.
    sources = repmat(make_source('', '', [0 0 0], f0, 1.0), 1, k);
    for i = 1:k
        point = [0.15 + 0.70*rand(), 0.35 + 0.30*rand(), 0.35 + 0.30*rand()];
        sources(i) = make_source(['R', num2str(i)], ['Random fault ', num2str(i)], point, f0, 1.0);
    end
end

function sourcesNew = perturb_source_parameters(sources, freqRelRange, ampRelRange)
    sourcesNew = sources;
    for i = 1:numel(sources)
        fFactor = 1 + freqRelRange * (2*rand() - 1);
        aFactor = 1 + ampRelRange * (2*rand() - 1);
        sourcesNew(i).f = sources(i).f * fFactor;
        sourcesNew(i).ampScale = sources(i).ampScale * aFactor;
    end
end

function r = make_source(id, name, point, f, ampScale)
    r.id = id;
    r.name = name;
    r.point = point;
    r.f = f;
    r.ampScale = ampScale;
end

function layouts = enumerate_layouts(n, k)
    C = nchoosek(1:n, k);
    layouts = mat2cell(C, ones(size(C,1),1), k);
end

function bestLayout = optimize_layout_by_lambda(sensors, sources, layouts, alpha, beta, omega, amp, sigma0, L)
    bestObj = -inf;
    bestMin = -inf;
    bestSpacing = -inf;
    bestLayout = layouts{1};
    for m = 1:numel(layouts)
        layout = layouts{m};
        lambdas = zeros(numel(sources), 1);
        for j = 1:numel(sources)
            lambdas(j) = layout_lambda(sensors, sources(j), layout, alpha, beta, amp, sigma0, L);
        end
        obj = omega * min(lambdas) + (1-omega) * mean(lambdas);
        spacing = mean_pair_spacing(sensors, layout);
        if obj > bestObj + 1e-12 || (abs(obj-bestObj) <= 1e-12 && min(lambdas) > bestMin) || ...
                (abs(obj-bestObj) <= 1e-12 && abs(min(lambdas)-bestMin) <= 1e-12 && spacing > bestSpacing)
            bestObj = obj;
            bestMin = min(lambdas);
            bestSpacing = spacing;
            bestLayout = layout;
        end
    end
end

function lam = layout_lambda(sensors, source, layout, alpha, beta, amp, sigma0, L)
    lam = 0;
    for idx = layout
        d = norm(sensors(idx).point - source.point);
        g = attenuation(d, alpha, beta);
        sigmai = sigma0 * sqrt(sensors(idx).nu);
        ai = sensors(idx).kappa * g * amp * source.ampScale;
        lam = lam + L * ai^2 / (2 * sigmai^2);
    end
end

function g = attenuation(d, alpha, beta)
    g = exp(-alpha*d) / (1 + beta*d);
end

function eta = chi2_threshold_even(df, pfa)
    lo = 0; hi = 1;
    while chi2sf_even(hi, df) > pfa
        hi = hi * 2;
    end
    for iter = 1:80
        mid = (lo + hi) / 2;
        if chi2sf_even(mid, df) > pfa
            lo = mid;
        else
            hi = mid;
        end
    end
    eta = (lo + hi) / 2;
end

function sf = chi2sf_even(x, df)
    m = df / 2;
    z = x / 2;
    s = 0;
    for k = 0:(m-1)
        s = s + z^k / factorial(k);
    end
    sf = exp(-z) * s;
end

function sim = simulate_layout_pd(sensors, sources, layout, alpha, beta, amp, sigma0, fs, t, eta, nMC)
    sourcePd = zeros(numel(sources), 1);
    for j = 1:numel(sources)
        hit = 0;
        for r = 1:nMC
            T = detection_statistic_trial(sensors, sources(j), layout, alpha, beta, amp, sigma0, t, true);
            hit = hit + (T > eta);
        end
        sourcePd(j) = hit / nMC;
    end
    sim.sourcePd = sourcePd;
    sim.meanPd = mean(sourcePd);
    sim.minPd = min(sourcePd);
end

function pfa = simulate_false_alarm(sensors, layout, sigma0, f0, t, eta, nMC)
    falseCount = 0;
    dummySource = make_source('H0', 'No fault', [0 0 0], f0, 0);
    for r = 1:nMC
        T = detection_statistic_trial(sensors, dummySource, layout, 0, 0, 0, sigma0, t, false);
        falseCount = falseCount + (T > eta);
    end
    pfa = falseCount / nMC;
end

function Tsum = detection_statistic_trial(sensors, source, layout, alpha, beta, amp, sigma0, t, hasSignal)
    Tsum = 0;
    L = numel(t);
    phi = 2*pi*rand();
    u = sin(2*pi*source.f*t);
    v = cos(2*pi*source.f*t);
    for idx = layout
        sigmai = sigma0 * sqrt(sensors(idx).nu);
        noise = sigmai * randn(L, 1);
        if hasSignal
            d = norm(sensors(idx).point - source.point);
            g = attenuation(d, alpha, beta);
            ai = sensors(idx).kappa * g * amp * source.ampScale;
            x = ai * sin(2*pi*source.f*t + phi) + noise;
        else
            x = noise;
        end
        I = sum(x .* u);
        Q = sum(x .* v);
        Ti = 2 / (L * sigmai^2) * (I^2 + Q^2);
        Tsum = Tsum + Ti;
    end
end

function layout = concentrated_layout(sensors, layouts)
    best = inf; layout = layouts{1};
    for k = 1:numel(layouts)
        ids = layouts{k};
        s = 0;
        for a = 1:numel(ids)
            for b = (a+1):numel(ids)
                s = s + norm(sensors(ids(a)).point - sensors(ids(b)).point);
            end
        end
        if s < best
            best = s;
            layout = ids;
        end
    end
end

function randLayouts = sample_random_layouts(layouts, n)
    n = min(n, numel(layouts));
    idx = randperm(numel(layouts), n);
    randLayouts = layouts(idx);
end

function best = best_single_sensor_by_lambda(sensors, sources, alpha, beta, omega, amp, sigma0, L)
    bestObj = -inf; best = 1;
    for i = 1:numel(sensors)
        lambdas = zeros(numel(sources),1);
        for j = 1:numel(sources)
            lambdas(j) = layout_lambda(sensors, sources(j), i, alpha, beta, amp, sigma0, L);
        end
        obj = omega*min(lambdas) + (1-omega)*mean(lambdas);
        if obj > bestObj
            bestObj = obj;
            best = i;
        end
    end
end

function spacing = mean_pair_spacing(sensors, layout)
    ds = [];
    for a = 1:numel(layout)
        for b = (a+1):numel(layout)
            ds(end+1) = norm(sensors(layout(a)).point - sensors(layout(b)).point); %#ok<AGROW>
        end
    end
    spacing = mean(ds);
end

function txt = layout_to_text(layout, sensors)
    txt = cell(numel(layout), 1);
    for k = 1:numel(layout)
        txt{k} = [sensors(layout(k)).id, ' - ', sensors(layout(k)).name];
    end
end

function ids = layout_ids(layout, sensors)
    parts = cell(1, numel(layout));
    for k = 1:numel(layout)
        parts{k} = sensors(layout(k)).id;
    end
    ids = strjoin(parts, ',');
end

function most = most_frequent_layout(layoutText)
    uniqueLayouts = unique(layoutText);
    counts = zeros(numel(uniqueLayouts), 1);
    for i = 1:numel(uniqueLayouts)
        counts(i) = sum(strcmp(layoutText, uniqueLayouts{i}));
    end
    [~, idx] = max(counts);
    most = uniqueLayouts{idx};
end

function row = make_failure_row(caseName, layout, sim)
    ids = cell(1, numel(layout));
    for k = 1:numel(layout)
        ids{k} = ['P', num2str(layout(k))];
    end
    row = table({caseName}, {strjoin(ids, ',')}, sim.meanPd, sim.minPd, ...
        'VariableNames', {'Case','RemainingSensors','MeanPd','MinPd'});
end

function [fAxis, fusionAmp] = one_fusion_spectrum(sensors, source, layout, alpha, beta, amp, sigma0, fs, t)
    L = numel(t);
    spectra = zeros(floor(L/2)+1, numel(layout));
    weights = zeros(numel(layout), 1);
    phi = 2*pi*rand();
    for k = 1:numel(layout)
        idx = layout(k);
        d = norm(sensors(idx).point - source.point);
        g = attenuation(d, alpha, beta);
        ai = sensors(idx).kappa * g * amp * source.ampScale;
        sigmai = sigma0 * sqrt(sensors(idx).nu);
        x = ai * sin(2*pi*source.f*t + phi) + sigmai * randn(L,1);
        X = abs(fft(x)) / L * 2;
        spectra(:,k) = X(1:floor(L/2)+1);
        weights(k) = ai^2 / (sigmai^2);
    end
    weights = weights / sum(weights);
    fusionAmp = spectra * weights;
    fAxis = (0:floor(L/2))' * fs / L;
end

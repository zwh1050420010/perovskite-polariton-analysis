%% Fit lower polariton and cavity dispersion
% Fixed exciton energy: 2.35 eV
% Input file format: LP_peak_points.csv with columns:
% angle_deg, LP_eV

clear; clc; close all;

Ex = 2.3964;                         % exciton energy, eV
inputFile = "LP18.xlsx";   % replace with your measured LP peak points
mapInputFile = "angle_energy_intensity_18in.xlsx";
outputFile = "polariton_fit_result_18in";
paramFile = "polariton_fit_parameters_18in.txt";

%% Load LP peak points
if isfile(inputFile)
    data = readmatrix(inputFile);
    thetaData = data(:, 1);
    lpData = data(:, 2);
else
    % Approximate points from the example figure. Replace these with your real
    % peak positions for a publication-quality fit.
    thetaData = (-45:5:45).';
    lpData = [2.286; 2.275; 2.266; 2.259; 2.254; 2.252; ...
              2.251; 2.252; 2.254; 2.259; 2.266; 2.275; 2.286];
end

valid = isfinite(thetaData) & isfinite(lpData);
thetaData = thetaData(valid);
lpData = lpData(valid);

%% Coupled oscillator model
% p = [Ec0, neff, theta0, Omega]
% Ec0   : cavity energy at theta0, eV
% neff  : effective refractive index
% theta0: angular offset, degree
% Omega : Rabi splitting, eV

cavityFun = @(theta, p) p(1) ./ sqrt(1 - (sin((theta - p(3)) * pi / 180) ./ p(2)).^2);
lpFun = @(theta, p) 0.5 .* (cavityFun(theta, p) + Ex ...
    - sqrt((cavityFun(theta, p) - Ex).^2 + p(4).^2));
upFun = @(theta, p) 0.5 .* (cavityFun(theta, p) + Ex ...
    + sqrt((cavityFun(theta, p) - Ex).^2 + p(4).^2));

%% Fit
% The LP-only fit is weakly constrained, so cavity anchor points are included.
% If you have measured cavity points, replace cavityAnchorTheta/EcavityAnchor.
cavityAnchorTheta = [-45; 0; 45];
EcavityAnchor = [2.305; 2.260; 2.305];
cavityWeight = 0.01;

costFun = @(p) sum(((lpFun(thetaData, p) - lpData) ./ 0.002).^2) ...
    + cavityWeight * sum(((cavityFun(cavityAnchorTheta, p) - EcavityAnchor) ./ 0.02).^2) ...
    + boundPenalty(p);

p0 = [2.26, 1.55, 0.0, 0.06];
options = optimset("Display", "final", "MaxFunEvals", 20000, "MaxIter", 10000, "TolX", 1e-10, "TolFun", 1e-10);
pFit = fminsearch(costFun, p0, options);

Ec0 = pFit(1);
neff = pFit(2);
theta0 = pFit(3);
Omega = pFit(4);

%% Generate fitted curves
thetaFit = (-45:0.1:45).';
cavityFit = cavityFun(thetaFit, pFit);
excitonFit = Ex .* ones(size(thetaFit));
lpFit = lpFun(thetaFit, pFit);
upFit = upFun(thetaFit, pFit);

resultTable = table(thetaFit, cavityFit, excitonFit, lpFit, upFit, ...
    'VariableNames', {'angle_deg', 'cavity_eV', 'exciton_eV', 'LP_eV', 'UP_eV'});
writetable(resultTable, outputFile);

fid = fopen(paramFile, "w");
fprintf(fid, "Fixed exciton energy Ex = %.6f eV\n", Ex);
fprintf(fid, "Ec0 = %.6f eV\n", Ec0);
fprintf(fid, "neff = %.6f\n", neff);
fprintf(fid, "theta0 = %.6f deg\n", theta0);
fprintf(fid, "Omega = %.6f eV\n", Omega);
fprintf(fid, "Rabi splitting = %.2f meV\n", Omega * 1000);
fclose(fid);

%% Diagnostic plot
figure("Color", "w");
hold on; box on;
plot(thetaFit, cavityFit, "-", "Color", [1 1 1] * 0.15, "LineWidth", 1.8);
plot(thetaFit, excitonFit, "-", "Color", [0.2 0.45 0.95], "LineWidth", 1.5);
plot(thetaFit, lpFit, "--", "Color", [0.1 0.1 0.1], "LineWidth", 1.8);
plot(thetaData, lpData, "o", "MarkerSize", 5, "MarkerFaceColor", [0.95 0.2 0.1], ...
    "MarkerEdgeColor", "none");
xlabel("Emission angle (degree)");
ylabel("Energy (eV)");
legend("Cavity", "Exciton", "LP fit", "LP peak points", "Location", "northwest");
xlim([-45 45]);
ylim([2.10 2.42]);
set(gca, "FontName", "Arial", "FontSize", 12, "LineWidth", 1);

print(gcf, "polariton_fit_plot.png", "-dpng", "-r300");

%% Final Origin-style polariton map
% If angle_energy_intensity.csv exists, it is used as the real background map.
% Required long-table format:
% angle_deg, energy_eV, intensity
% If the file is absent, a weak LP emission map is simulated only for layout.
fprintf("=== 读取强度背景文件：%s ===\n",mapInputFile);
[thetaMap, energyMap, intensity] = loadIntensityMap(mapInputFile, thetaFit, pFit, lpFun);

fig2 = figure("Color", "w", "Position", [100 100 560 500]);
ax1 = axes(fig2);
plotArea = [0.15 0.14 0.66 0.70];
imagesc(ax1, thetaMap, energyMap, intensity);
set(ax1, "YDir", "normal");
hold(ax1, "on");
colormap(ax1, jet(256));
caxis(ax1, [0 1]);

plot(ax1, thetaFit, cavityFit, "-", "Color", "w", "LineWidth", 1.5);
plot(ax1, thetaFit, excitonFit, "-", "Color", "w", "LineWidth", 1.2);
plot(ax1, thetaFit, lpFit, ":", "Color", "w", "LineWidth", 1.6);

text(ax1, -25, 2.390, "Cavity", ...
    "Color", "w", "FontName", "Arial", "FontSize", 12);
text(ax1, -26, Ex - 0.006, "Exciton", ...
    "Color", "w", "FontName", "Arial", "FontSize", 12);
text(ax1, 22, interp1(thetaFit, lpFit, 22) - 0.018, "LP", ...
    "Color", "w", "FontName", "Arial", "FontSize", 12);
text(ax1, -0.08, 1.04, "D", "Units", "normalized", "Color", "k", ...
    "FontName", "Arial", "FontWeight", "bold", "FontSize", 12);

xlabel(ax1, "Emission angle (degree)");
ylabel(ax1, "Energy (eV)");
xlim(ax1, [-45 45]);
ylim(ax1, [2.10 2.42]);
set(ax1, "FontName", "Arial", "FontSize", 12, "LineWidth", 1.1, ...
    "TickDir", "out", "Layer", "top", "Box", "on");

cb = colorbar(ax1);
cb.Label.String = "Intensity (a.u.)";
cb.Ticks = [0 1];
cb.FontName = "Arial";
cb.FontSize = 11;
ax1.Position = plotArea;
cb.Position = [0.84 0.14 0.035 0.70];

% Top k_parallel axis. For a visual guide, convert k to angle using
% k_parallel = k0 * sin(theta), with k0 from the exciton energy.
hbarc_eV_um = 0.1973269804;
k0 = Ex / hbarc_eV_um;
kTicks = -4:2:4;
thetaTicksTop = asind(kTicks ./ k0);
ax2 = axes(fig2, "Position", ax1.Position, "Color", "none", ...
    "XAxisLocation", "top", "YAxisLocation", "right", ...
    "XLim", ax1.XLim, "YLim", ax1.YLim, "YTick", [], ...
    "XTick", thetaTicksTop, "XTickLabel", string(kTicks), ...
    "FontName", "Arial", "FontSize", 11, "LineWidth", 1.1, ...
    "TickDir", "out", "Box", "off");
xlabel(ax2, "k_{||} (\mum^{-1})");

print(fig2, "polariton_origin_style.png", "-dpng", "-r300");
print(fig2, "polariton_final_overlay.png", "-dpng", "-r300");

fprintf("\nFit finished.\n");
fprintf("Ec0 = %.6f eV\n", Ec0);
fprintf("neff = %.6f\n", neff);
fprintf("theta0 = %.6f deg\n", theta0);
fprintf("Omega = %.6f eV = %.2f meV\n", Omega, Omega * 1000);
fprintf("Wrote: %s\n", outputFile);
fprintf("Wrote: %s\n", paramFile);
fprintf("Wrote: polariton_fit_plot.png\n");
fprintf("Wrote: polariton_origin_style.png\n");
fprintf("Wrote: polariton_final_overlay.png\n");

%% Local function
function penalty = boundPenalty(p)
    Ec0 = p(1);
    neff = p(2);
    theta0 = p(3);
    Omega = p(4);

    penalty = 0;
    penalty = penalty + softBound(Ec0, 2.20, 2.42, 1e6);
    penalty = penalty + softBound(neff, 1.05, 4.00, 1e6);
    penalty = penalty + softBound(theta0, -8.00, 8.00, 1e5);
    penalty = penalty + softBound(Omega, 0.02, 0.50, 1e6);
end

function v = softBound(x, lo, hi, weight)
    v = 0;
    if x < lo
        v = weight * (lo - x)^2;
    elseif x > hi
        v = weight * (x - hi)^2;
    end
end

function [thetaMap, energyMap, intensity] = loadIntensityMap(mapInputFile, thetaFit, pFit, lpFun)
    if isfile(mapInputFile)
        raw = readmatrix(mapInputFile);
        raw = raw(all(isfinite(raw), 2), :);
        if size(raw, 2) < 3
            error("angle_energy_intensity.csv needs three columns: angle_deg, energy_eV, intensity");
        end

        thetaRaw = raw(:, 1);
        energyRaw = raw(:, 2);
        intensityRaw = raw(:, 3);
        thetaMap = sort(unique(thetaRaw));
        energyMap = sort(unique(energyRaw));
        %thetaMap = unique(thetaRaw).';
        %energyMap = unique(energyRaw).';
        [thetaMesh, energyMesh] = meshgrid(thetaMap, energyMap);
        intensity = griddata(thetaRaw, energyRaw, intensityRaw, thetaMesh, energyMesh, "natural");
        intensity(isnan(intensity)) = 0;
    else
        thetaMap = -30:0.2:30;
        energyMap = 2.10:0.001:2.42;
        [thetaMesh, ~] = meshgrid(thetaMap, energyMap);
        lpMap = lpFun(thetaMap, pFit);

        intensity = zeros(size(thetaMesh));
        for ii = 1:numel(thetaMap)
            amp = exp(-(thetaMap(ii) / 4.0)^2);
            intensity(:, ii) = amp .* exp(-((energyMap(:) - lpMap(ii)) / 0.0016).^2);
        end
    end

    intensity = intensity - min(intensity(:));
    if max(intensity(:)) > 0
        intensity = intensity ./ max(intensity(:));
    end
end

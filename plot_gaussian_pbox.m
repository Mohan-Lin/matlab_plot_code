% clear; clc; close all;

% --- 1. Distribution Parameter Interval Settings ---
mu_min = 10;
mu_max = 15;

sigma = 2.5;       % Fixed standard deviation to avoid piecewise breakpoints at the boundaries

% --- 2. Style Control ---
lineWidth = 15;
lineColor = [0, 0.447, 0.698];
fillColor = [0.88, 0.88, 0.88];
fillAlpha = 0.45;

% --- 3. Data Generation ---
x_min = mu_min - 4 * sigma;
x_max = mu_max + 4 * sigma;
x = linspace(x_min, x_max, 1000);

% --- 4. Calculate Smooth P-box Upper and Lower Bounds ---
% For a normal CDF, a smaller mu shifts the CDF left, yielding the upper bound;
% a larger mu shifts the CDF right, yielding the lower bound.
F_upper = 0.5 * (1 + erf((x - mu_min) / (sigma * sqrt(2))));
F_lower = 0.5 * (1 + erf((x - mu_max) / (sigma * sqrt(2))));

% --- 5. Create Transparent Figure and Axes ---
fig = figure('Color', 'none', 'Name', 'Gaussian P-box Pure');
ax = axes('Parent', fig, 'Color', 'none');
hold(ax, 'on');

% --- 6. Fill the P-box Region ---
fill(ax, ...
    [x, fliplr(x)], ...
    [F_upper, fliplr(F_lower)], ...
    fillColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', fillAlpha);

% --- 7. Plot Upper and Lower Bounds ---
plot(ax, x, F_upper, ...
    'Color', lineColor, ...
    'LineWidth', lineWidth);

plot(ax, x, F_lower, ...
    'Color', lineColor, ...
    'LineWidth', lineWidth);

% --- 8. Remove Axes and Grid ---
axis(ax, 'off');
grid(ax, 'off');
box(ax, 'off');

xlim(ax, [x_min, x_max]);
ylim(ax, [0, 1]);

hold(ax, 'off');

% --- 9. Export ---
% exportgraphics(fig, 'Gaussian_Pbox_Pure.svg', 'ContentType', 'vector', 'BackgroundColor', 'none');
% exportgraphics(fig, 'Gaussian_Pbox_Pure.png', 'Resolution', 300, 'BackgroundColor', 'none');
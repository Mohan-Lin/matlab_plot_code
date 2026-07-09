% =========================================================================
% Generate Smooth Weibull P-box:
% Uncertain scale parameter A, fixed shape parameter B
% Transparent background, no axes, ultra-thick blue upper/lower boundaries, light gray shading fill
% =========================================================================
clear; clc; close all;

% --- 1. Weibull Distribution Parameter Interval Settings ---
% Weibull CDF:
% F(x; A, B) = 1 - exp(-(x/A)^B), x >= 0
%
% A: scale parameter
% B: shape parameter
A_min = 6.0;      % Lower bound of the scale parameter
A_max = 15.0;     % Upper bound of the scale parameter
B = 2.0;          % Fixed shape parameter, commonly used to model wind speed Weibull distributions

% --- 2. Style Control ---
lineWidth = 15;                  % Ultra-thick style for blue boundary lines
lineColor = [0, 0.447, 0.698];   % Nature scientific blue
fillColor = [0.88, 0.88, 0.88];  % Light gray shading
fillAlpha = 0.35;                % Transparency of the p-box region

% --- 3. Data Generation ---
% Weibull distribution is defined for x >= 0
% Taking approximately 4 times A_max basically covers the main variation region of the CDF from 0 to 1
x_min = 0;
x_max = 4 * A_max;
x = linspace(x_min, x_max, 1000);

% --- 4. Calculate Weibull P-box Upper and Lower Bounds ---
% A smaller A results in a faster-rising CDF, thus corresponding to the upper bound;
% A larger A results in a slower-rising CDF, thus corresponding to the lower bound.
F_upper = 1 - exp(-(x ./ A_min) .^ B);
F_lower = 1 - exp(-(x ./ A_max) .^ B);

% Numerical protection to prevent extremely small errors from exceeding [0, 1]
F_upper = min(max(F_upper, 0), 1);
F_lower = min(max(F_lower, 0), 1);

% --- 5. Create Transparent Figure and Axes ---
fig = figure('Color', 'none', 'Name', 'Weibull P-box Pure');
ax = axes('Parent', fig, 'Color', 'none');
hold(ax, 'on');

% --- 6. Fill the P-box Region ---
fill(ax, ...
    [x, fliplr(x)], ...
    [F_upper, fliplr(F_lower)], ...
    fillColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', fillAlpha);

% --- 7. Plot P-box Upper Bound ---
plot(ax, x, F_upper, ...
    'Color', lineColor, ...
    'LineWidth', lineWidth);

% --- 8. Plot P-box Lower Bound ---
plot(ax, x, F_lower, ...
    'Color', lineColor, ...
    'LineWidth', lineWidth);

% --- 9. Remove Axes and Grid ---
axis(ax, 'off');
grid(ax, 'off');
box(ax, 'off');

% Fixed display range to prevent thick lines from being clipped
xlim(ax, [x_min, x_max]);
ylim(ax, [0, 1]);

hold(ax, 'off');

% --- 10. Vector / High-Resolution Export Commands ---
% exportgraphics(fig, 'Weibull_Pbox_Pure.svg', 'ContentType', 'vector', 'BackgroundColor', 'none');
% exportgraphics(fig, 'Weibull_Pbox_Pure.png', 'Resolution', 300, 'BackgroundColor', 'none');
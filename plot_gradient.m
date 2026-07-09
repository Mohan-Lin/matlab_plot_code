% =========================================================================
% Nature-style 3D surface plot of gradient descent
% Requirements:
%   1. Surface with clear undulations, but not overly deep
%   2. Gradient descent path closely follows the surface, no visual lift
%   3. No axes, no grid, no borders
%   4. No mesh, no lighting, no specular highlights
%   5. Transparent background
%   6. Nature/NPG style color scheme
% =========================================================================

clear; clc; close all;

% --- 1. Nature/NPG style color scheme ---
npgRed     = [0.902, 0.294, 0.208];   % #E64B35
npgBlue    = [0.235, 0.329, 0.533];   % #3C5488
npgCyan    = [0.302, 0.733, 0.835];   % #4DBBD5
npgGreen   = [0.000, 0.627, 0.529];   % #00A087
npgMint    = [0.569, 0.820, 0.757];   % #91D1C2
npgOrange  = [0.949, 0.639, 0.455];   % #F39B7F
softSand   = [0.930, 0.860, 0.650];
softPurple = [0.520, 0.460, 0.700];
whiteLine  = [1.000, 1.000, 1.000];

% --- 2. Define a non-convex loss function similar to peaks ---
% This function has peaks, valleys, and saddle points, appearing more natural
% than a simple bowl-shaped surface.
J = @(x, y) ...
    3 .* (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) ...
    - 10 .* (x ./ 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) ...
    - 1/3 .* exp(-(x + 1).^2 - y.^2);

gradJ = @(x, y) numericalGradient(J, x, y);

xPlotLim = [-2.25, 2.25];
yPlotLim = [-2.25, 2.25];
% --- 3. Gradient descent initial point ---
% This initial point gradually descends along the surface toward the main valley.
x0 = 0.00;
y0 = 1.20;

maxIter = 90;
alpha0  = 0.16;
tolGrad = 1e-5;

x_path = zeros(maxIter + 1, 1);
y_path = zeros(maxIter + 1, 1);
z_path = zeros(maxIter + 1, 1);

x_path(1) = x0;
y_path(1) = y0;
z_path(1) = J(x0, y0);

% --- 4. Armijo backtracking line search gradient descent ---
for k = 1:maxIter

    xk = x_path(k);
    yk = y_path(k);
    fk = J(xk, yk);

    g = gradJ(xk, yk);
    gnorm2 = g(1)^2 + g(2)^2;

    if sqrt(gnorm2) < tolGrad
        x_path = x_path(1:k);
        y_path = y_path(1:k);
        z_path = z_path(1:k);
        break;
    end

    alpha = alpha0;

    while true
        x_new = xk - alpha * g(1);
        y_new = yk - alpha * g(2);

        x_new = min(max(x_new, xPlotLim(1)), xPlotLim(2));
        y_new = min(max(y_new, yPlotLim(1)), yPlotLim(2));

        f_new = J(x_new, y_new);

        if f_new <= fk - 1e-4 * alpha * gnorm2 || alpha < 1e-7
            break;
        end

        alpha = 0.5 * alpha;
    end

    x_path(k + 1) = x_new;
    y_path(k + 1) = y_new;
    z_path(k + 1) = f_new;
end

% --- 5. Construct surface mesh ---
x = linspace(xPlotLim(1), xPlotLim(2), 800);
y = linspace(yPlotLim(1), yPlotLim(2), 800);
[X, Y] = meshgrid(x, y);
Z = J(X, Y);

% --- 6. Surface height processing: key modification ---
% No longer using zScale = 0.05.
% Use targetZRange to control visual height: maintain undulations without being too deep.
Zmin_raw = min([Z(:); z_path(:)]);
Zmax_raw = max([Z(:); z_path(:)]);

denom = Zmax_raw - Zmin_raw;
if denom < eps
    denom = 1;
end

Z_norm = (Z - Zmin_raw) ./ denom;
z_path_norm = (z_path - Zmin_raw) ./ denom;

Z_norm = min(max(Z_norm, 0), 1);
z_path_norm = min(max(z_path_norm, 0), 1);

targetZRange = 1.35;       % Controls overall surface height; larger values give more pronounced undulations
bottomOffset  = 0.18;      % Raise the lowest point to avoid overly deep valleys

Z_plot = bottomOffset + targetZRange .* Z_norm;
z_path_plot = bottomOffset + targetZRange .* z_path_norm;

% --- 7. Compute global minimum of the surface for verification ---
[Zmin_grid, idxMin] = min(Z(:));
x_min_global = X(idxMin);
y_min_global = Y(idxMin);

% --- 8. Smooth the gradient descent path ---
% Path is smoothed only in the x-y plane; z is recomputed from the surface
% function to ensure it remains on the surface.
ds = sqrt(diff(x_path).^2 + diff(y_path).^2);
s = [0; cumsum(ds)];

valid = [true; diff(s) > 1e-10];

s_valid = s(valid);
x_valid = x_path(valid);
y_valid = y_path(valid);

if numel(s_valid) < 4
    x_smooth = x_valid;
    y_smooth = y_valid;
else
    s_fine = linspace(s_valid(1), s_valid(end), 700);
    x_smooth = interp1(s_valid, x_valid, s_fine, 'pchip');
    y_smooth = interp1(s_valid, y_valid, s_fine, 'pchip');
end

z_smooth_raw = J(x_smooth, y_smooth);
z_smooth_norm = (z_smooth_raw - Zmin_raw) ./ denom;
z_smooth_norm = min(max(z_smooth_norm, 0), 1);
z_smooth = bottomOffset + targetZRange .* z_smooth_norm;

% % --- 9. Create transparent figure ---
% fig = figure( ...
%     'Color', 'none', ...
%     'Name', 'Nature Style Gradient Descent Surface', ...
%     'Position', [80, 120, 1500, 850]);
% 
% ax = axes( ...
%     'Parent', fig, ...
%     'Color', 'none', ...
%     'Position', [0.01, 0.01, 0.98, 0.98]);
% 
% hold(ax, 'on');
% --- 9. Create transparent 4K figure ---
fig = figure( ...
    'Color', 'none', ...
    'Name', 'Nature Style Gradient Descent Surface', ...
    'Units', 'pixels', ...
    'Position', [100, 100, 3840, 2160], ...
    'Renderer', 'opengl');

fig.GraphicsSmoothing = 'on';
set(fig, 'InvertHardcopy', 'off');

ax = axes( ...
    'Parent', fig, ...
    'Color', 'none', ...
    'Position', [0.00, 0.00, 1.00, 1.00]);

hold(ax, 'on');

% --- 10. Nature/NPG style colormap ---
nSeg = 90;

cmap = [
    interpColor(npgBlue,    softPurple, nSeg);
    interpColor(softPurple, npgCyan,    nSeg);
    interpColor(npgCyan,    npgMint,    nSeg);
    interpColor(npgMint,    softSand,   nSeg);
    interpColor(softSand,   npgOrange,  nSeg);
    interpColor(npgOrange,  npgRed,     nSeg)
];

colormap(ax, cmap);

% --- 11. Plot surface: no mesh, no lighting, no specular highlights ---
surf(ax, X, Y, Z_plot, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1.00);

shading(ax, 'interp');

% --- 12. Plot clearer and more concise gradient descent path ---
% Note:
%   trajLift and markerLift are only used for display to prevent
%   the trajectory from being obscured by the surface;
%   they do not alter the actual gradient descent results.
trajLift   = 0.014;
markerLift = 0.022;

z_smooth_show = z_smooth + trajLift;

% Trajectory color: Nature-style vermillion, more prominent than dark blue
routeColor = [0.902, 0.294, 0.208];   % NPG red, #E64B35
haloColor  = [0.985, 0.985, 0.960];   % Soft white outline, not harsh

% Soft white outline: ensures the curve is clearly visible on both dark and light surfaces
plot3(ax, x_smooth, y_smooth, z_smooth_show, ...
    '-', ...
    'Color', haloColor, ...
    'LineWidth', 12.2);

% Main optimization trajectory: vermillion main line
plot3(ax, x_smooth, y_smooth, z_smooth_show + 0.002, ...
    '-', ...
    'Color', routeColor, ...
    'LineWidth', 10.4);

% --- 13. Plot discrete iteration points ---
% Points should not be too dense to avoid breaking the overall curve aesthetics
markerStep = 12;
idx_marker = unique([1:markerStep:length(x_path), length(x_path)]);

z_marker_show = z_path_plot + markerLift;

scatter3(ax, ...
    x_path(idx_marker), ...
    y_path(idx_marker), ...
    z_marker_show(idx_marker), ...
    62, ...
    routeColor, ...
    'filled', ...
    'MarkerEdgeColor', haloColor, ...
    'LineWidth', 1.25);

% --- 14. Start and end points ---
% Start point: Nature green
scatter3(ax, x_path(1), y_path(1), z_path_plot(1) + markerLift + 0.008, ...
    190, npgGreen, ...
    'filled', ...
    'MarkerEdgeColor', haloColor, ...
    'LineWidth', 1.6);

% End point: deep purple-blue, indicating convergence to the minimum
endColor = [0.365, 0.267, 0.557];     % Stable purple-blue

scatter3(ax, x_path(end), y_path(end), z_path_plot(end) + markerLift + 0.008, ...
    235, endColor, ...
    'filled', ...
    'MarkerEdgeColor', haloColor, ...
    'LineWidth', 1.6);

% --- 15. Turn off axes, grid, and borders ---
axis(ax, 'off');
grid(ax, 'off');
box(ax, 'off');

% --- 16. Control view range and aspect ratio ---
xlim(ax, xPlotLim);
ylim(ax, yPlotLim);

zlim(ax, [0, bottomOffset + targetZRange + 0.12]);

% Key: do not compress too much here.
% A value around 0.55 preserves undulations while avoiding the overly deep
% appearance of the original peaks function.
daspect(ax, [1.0, 1.0, 0.55]);

view(ax, [-42, 25]);

hold(ax, 'off');

% --- 17. High-resolution transparent background export ---
% exportgraphics(fig, 'Nature_Gradient_Descent_Surface_Balanced.png', ...
%     'Resolution', 600, ...
%     'BackgroundColor', 'none');

% exportgraphics(fig, 'Nature_Gradient_Descent_Surface_Balanced.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'none');

% --- 18. Terminal output verification information ---
fprintf('Gradient descent endpoint: x = %.4f, y = %.4f, J = %.6f\n', ...
    x_path(end), y_path(end), z_path(end));

fprintf('Grid global minimum:       x = %.4f, y = %.4f, J = %.6f\n', ...
    x_min_global, y_min_global, Zmin_grid);

fprintf('Distance to grid minimum:  %.6f\n', ...
    sqrt((x_path(end) - x_min_global)^2 + (y_path(end) - y_min_global)^2));

% =========================================================================
% Local function 1: Central difference numerical gradient
% =========================================================================
function g = numericalGradient(fun, x, y)

    h = 1e-4;

    dfdx = (fun(x + h, y) - fun(x - h, y)) / (2 * h);
    dfdy = (fun(x, y + h) - fun(x, y - h)) / (2 * h);

    g = [dfdx; dfdy];

end

% =========================================================================
% Local function 2: Linear interpolation between two colors
% =========================================================================
function cmap = interpColor(c1, c2, n)

    cmap = [ ...
        linspace(c1(1), c2(1), n)', ...
        linspace(c1(2), c2(2), n)', ...
        linspace(c1(3), c2(3), n)' ...
    ];

end
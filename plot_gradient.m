% =========================================================================
% Nature 风格梯度下降三维曲面图
% 要求：
%   1. 曲面有明显起伏，但不要过深
%   2. 梯度下降路径贴合曲面，不做视觉抬升
%   3. 无坐标轴、无网格、无边框
%   4. 无网面、无光照、无光泽
%   5. 透明背景
%   6. Nature/NPG 风格配色
% =========================================================================

clear; clc; close all;

% --- 1. Nature/NPG 风格配色 ---
npgRed     = [0.902, 0.294, 0.208];   % #E64B35
npgBlue    = [0.235, 0.329, 0.533];   % #3C5488
npgCyan    = [0.302, 0.733, 0.835];   % #4DBBD5
npgGreen   = [0.000, 0.627, 0.529];   % #00A087
npgMint    = [0.569, 0.820, 0.757];   % #91D1C2
npgOrange  = [0.949, 0.639, 0.455];   % #F39B7F
softSand   = [0.930, 0.860, 0.650];
softPurple = [0.520, 0.460, 0.700];
whiteLine  = [1.000, 1.000, 1.000];

% --- 2. 定义类似 peaks 的非凸损失函数 ---
% 该函数具有峰、谷、鞍点，视觉上比单一碗形曲面更自然。
J = @(x, y) ...
    3 .* (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) ...
    - 10 .* (x ./ 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) ...
    - 1/3 .* exp(-(x + 1).^2 - y.^2);

gradJ = @(x, y) numericalGradient(J, x, y);

xPlotLim = [-2.25, 2.25];
yPlotLim = [-2.25, 2.25];
% --- 3. 梯度下降初始点 ---
% 该初始点能沿曲面逐步下降到主要低谷附近。
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

% --- 4. Armijo 回溯线搜索梯度下降 ---
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

% --- 5. 构造曲面网格 ---
x = linspace(xPlotLim(1), xPlotLim(2), 800);
y = linspace(yPlotLim(1), yPlotLim(2), 800);
[X, Y] = meshgrid(x, y);
Z = J(X, Y);

% --- 6. 曲面高度处理：关键修改 ---
% 不再使用 zScale = 0.05。
% 用 targetZRange 控制视觉高度：既有起伏，又不会太深。
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

targetZRange = 1.35;       % 控制曲面整体高度；越大起伏越明显
bottomOffset  = 0.18;      % 抬高最低点，避免谷底过深

Z_plot = bottomOffset + targetZRange .* Z_norm;
z_path_plot = bottomOffset + targetZRange .* z_path_norm;

% --- 7. 计算曲面全局最低点，用于验证 ---
[Zmin_grid, idxMin] = min(Z(:));
x_min_global = X(idxMin);
y_min_global = Y(idxMin);

% --- 8. 平滑梯度下降路径 ---
% 路径只在 x-y 平面做平滑，z 重新由曲面函数计算，保证贴合曲面。
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

% % --- 9. 创建透明图窗 ---
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
% --- 9. 创建透明 4K 图窗 ---
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

% --- 10. Nature/NPG 风格 colormap ---
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

% --- 11. 绘制普通曲面：无网面、无光照、无光泽 ---
surf(ax, X, Y, Z_plot, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1.00);

shading(ax, 'interp');

% --- 12. 绘制更清晰、更简洁的梯度下降路径 ---
% 说明：
%   trajLift 和 markerLift 只用于显示，避免轨迹被曲面遮挡；
%   不改变真实梯度下降计算结果。
trajLift   = 0.014;
markerLift = 0.022;

z_smooth_show = z_smooth + trajLift;

% 轨迹颜色：Nature 风格朱红色，比深蓝更醒目
routeColor = [0.902, 0.294, 0.208];   % NPG red, #E64B35
haloColor  = [0.985, 0.985, 0.960];   % 柔和白色描边，不刺眼

% 白色柔和描边：保证曲线在深色/浅色曲面上都清楚
plot3(ax, x_smooth, y_smooth, z_smooth_show, ...
    '-', ...
    'Color', haloColor, ...
    'LineWidth', 12.2);

% 主寻优轨迹：朱红色主线
plot3(ax, x_smooth, y_smooth, z_smooth_show + 0.002, ...
    '-', ...
    'Color', routeColor, ...
    'LineWidth', 10.4);

% --- 13. 绘制离散迭代点 ---
% 点不要太密，避免破坏曲线整体感
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

% --- 14. 起点和终点 ---
% 起点：Nature 绿色
scatter3(ax, x_path(1), y_path(1), z_path_plot(1) + markerLift + 0.008, ...
    190, npgGreen, ...
    'filled', ...
    'MarkerEdgeColor', haloColor, ...
    'LineWidth', 1.6);

% 终点：深紫蓝，表示收敛到最低点
endColor = [0.365, 0.267, 0.557];     % 稳重紫蓝色

scatter3(ax, x_path(end), y_path(end), z_path_plot(end) + markerLift + 0.008, ...
    235, endColor, ...
    'filled', ...
    'MarkerEdgeColor', haloColor, ...
    'LineWidth', 1.6);

% --- 15. 关闭坐标轴、网格、边框 ---
axis(ax, 'off');
grid(ax, 'off');
box(ax, 'off');

% --- 16. 控制画面范围和比例 ---
xlim(ax, xPlotLim);
ylim(ax, yPlotLim);

zlim(ax, [0, bottomOffset + targetZRange + 0.12]);

% 关键：这里不要再压得太狠。
% 0.55 左右能保留起伏，同时不会像原始 peaks 那样很深。
daspect(ax, [1.0, 1.0, 0.55]);

view(ax, [-42, 25]);

hold(ax, 'off');

% --- 17. 高清透明背景导出 ---
% exportgraphics(fig, 'Nature_Gradient_Descent_Surface_Balanced.png', ...
%     'Resolution', 600, ...
%     'BackgroundColor', 'none');

% exportgraphics(fig, 'Nature_Gradient_Descent_Surface_Balanced.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'none');

% --- 18. 终端输出验证信息 ---
fprintf('Gradient descent endpoint: x = %.4f, y = %.4f, J = %.6f\n', ...
    x_path(end), y_path(end), z_path(end));

fprintf('Grid global minimum:       x = %.4f, y = %.4f, J = %.6f\n', ...
    x_min_global, y_min_global, Zmin_grid);

fprintf('Distance to grid minimum:  %.6f\n', ...
    sqrt((x_path(end) - x_min_global)^2 + (y_path(end) - y_min_global)^2));

% =========================================================================
% 局部函数 1：中心差分数值梯度
% =========================================================================
function g = numericalGradient(fun, x, y)

    h = 1e-4;

    dfdx = (fun(x + h, y) - fun(x - h, y)) / (2 * h);
    dfdy = (fun(x, y + h) - fun(x, y - h)) / (2 * h);

    g = [dfdx; dfdy];

end

% =========================================================================
% 局部函数 2：两种颜色之间线性插值
% =========================================================================
function cmap = interpColor(c1, c2, n)

    cmap = [ ...
        linspace(c1(1), c2(1), n)', ...
        linspace(c1(2), c2(2), n)', ...
        linspace(c1(3), c2(3), n)' ...
    ];

end
function analyze_dt_statistics(source_id, fre, dt_dir, out_dir)
% ANALYZE_DT_STATISTICS  Per-source dt scatter + linear fit + radar chart
% =========================================================================
%  For one source point and one period band, this function:
%    1) Loads dt picks for all 4 velocity models (yao/zhang/bao/feng)
%    2) Splits by 12 azimuthal sectors
%    3) Scatter plot of dt vs distance per sector, with linear fit (k, σ)
%    4) Radar (polar) chart of per-sector mean dt with ±σ band
%
%  Inputs:
%    source_id - source point ID, e.g. 104
%    fre       - frequency band, one of {'5-10','8-18','15-35','20-45'}
%    dt_dir    - directory containing dt text files:
%                <dt_dir>/<source_id>_<model>_dt_<fre>_v9.txt
%    out_dir   - output directory for PNG figures
%
%  dt file format (8 columns):
%    [sta_meta(1-5), dt(6), distance_m(7), azimuth_deg(8)]
%
%  Example:
%    analyze_dt_statistics(179, '20-45', './output/179dt', './figures/179');
% =========================================================================

if nargin < 4
    error('Usage: analyze_dt_statistics(source_id, fre, dt_dir, out_dir)');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% --- Build filenames for the 4 models ---
models = {'yao', 'zhang', 'bao', 'feng'};
data   = cell(4, 1);
for im = 1:4
    fn = fullfile(dt_dir, ...
        sprintf('%d_%s_dt_%s_v9.txt', source_id, models{im}, fre));
    if ~exist(fn, 'file')
        warning('Missing: %s — will be skipped in plots', fn);
        data{im} = [];
    else
        data{im} = load(fn);
    end
end
yao = data{1}; zhang = data{2}; bao = data{3}; feng = data{4};

% --- 12-sector configuration (matches rtspg_pick_dt.m) ---
dire_vals   = [90, 270, 180, 0, 60, 240, 30, 210, 300, 120, 330, 150];
subplot_pos = [4, 10, 7, 1, 3, 9, 2, 8, 11, 5, 12, 6];
dire_titles = {'90°','270°','180°','0°','60°','240°', ...
               '30°','210°','300°','120°','330°','150°'};

% --- Colors (consistent across all figures in paper) ---
c_yao   = [0.8500 0.3250 0.0980];   % SWChinaCVM-2.0
c_zhang = [0.0000 0.4470 0.7410];   % USTClitho2.0
c_bao   = [0.4660 0.6740 0.1880];   % Bao20
c_feng  = [0.9290 0.6940 0.1250];   % Feng20

% --- Preallocate per-sector mean & std arrays ---
yao_mean   = nan(1,12); yao_std   = nan(1,12);
zhang_mean = nan(1,12); zhang_std = nan(1,12);
bao_mean   = nan(1,12); bao_std   = nan(1,12);
feng_mean  = nan(1,12); feng_std  = nan(1,12);

% =========================================================================
%  Figure 1: dt vs distance scatter + linear fit per sector
% =========================================================================
figure(1);
set(gcf, 'position', [800 800 1350 900]);

for ia = 1:12
    az = dire_vals(ia);

    % --- Filter by azimuth (column 8) ---
    p_yao   = idx_by_az(yao,   az);
    p_zhang = idx_by_az(zhang, az);
    p_bao   = idx_by_az(bao,   az);
    p_feng  = idx_by_az(feng,  az);

    if isempty(p_yao) && isempty(p_zhang) && ...
       isempty(p_bao) && isempty(p_feng)
        continue;
    end

    % --- Per-model mean & std ---
    [yao_mean(ia),   yao_std(ia)]   = mean_std(yao,   p_yao);
    [zhang_mean(ia), zhang_std(ia)] = mean_std(zhang, p_zhang);
    [bao_mean(ia),   bao_std(ia)]   = mean_std(bao,   p_bao);
    [feng_mean(ia),  feng_std(ia)]  = mean_std(feng,  p_feng);

    % --- Subplot ---
    all_dist = [get_col(yao,p_yao,7); get_col(zhang,p_zhang,7); ...
                get_col(bao,p_bao,7); get_col(feng,p_feng,7)];
    xmax = max(all_dist) * 1e-3;

    subplot(4, 3, subplot_pos(ia));
    line([0 xmax], [0 0], 'linestyle','--','Color','black','LineWidth',1.5);
    hold on;

    % --- Per-model scatter + linear fit ---
    models_loop = { p_yao,   yao,   c_yao,   80; ...
                    p_zhang, zhang, c_zhang, 60; ...
                    p_bao,   bao,   c_bao,   50; ...
                    p_feng,  feng,  c_feng,  40 };
    k_arr     = nan(4, 1);
    sigma_arr = nan(4, 1);

    for im = 1:4
        p_m   = models_loop{im, 1};
        dat_m = models_loop{im, 2};
        clr   = models_loop{im, 3};
        msz   = models_loop{im, 4};
        if isempty(p_m), continue; end

        r_m  = dat_m(p_m, 7) * 1e-3;   % km
        dt_m = dat_m(p_m, 6);          % s
        scatter(r_m, dt_m, msz, clr, 'filled');

        valid = isfinite(dt_m);
        r_v = r_m(valid); dt_v = dt_m(valid);
        if numel(r_v) >= 2
            p_fit = polyfit(r_v, dt_v, 1);
            k_arr(im) = p_fit(1);
            resid = dt_v - polyval(p_fit, r_v);
            sigma_arr(im) = sqrt(mean(resid.^2));
        end
    end

    set(gca, 'YDir', 'reverse');
    ylim([-10 10]); xlim([0 xmax]);
    set(gca, 'FontSize', 13);
    yticks(linspace(-10, 10, 5));
    grid on;
    title(dire_titles{ia});

    % --- Annotation: k and sigma per model ---
    clrs_ann = {c_yao, c_zhang, c_bao, c_feng};
    y_ann_start = -9.0; dy_ann = 2.0;
    for im = 1:4
        if isnan(k_arr(im)), continue; end
        ann = sprintf('k=%.3f \\sigma=%.2f', k_arr(im), sigma_arr(im));
        text(xmax * 0.98, y_ann_start + (im-1)*dy_ann, ann, ...
            'Color', clrs_ann{im}, 'FontSize', 7, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end

    if ismember(subplot_pos(ia), [10 11 12])
        xlabel('Distance (km)', 'FontSize', 13);
    end
    if ismember(subplot_pos(ia), [1 4 7 10])
        ylabel('\Deltat (s)', 'FontSize', 13);
    end
end

% --- Common legend ---
subplot(4, 3, 6);
h(1) = scatter(nan, nan, 80, c_yao,   'filled');
h(2) = scatter(nan, nan, 60, c_zhang, 'filled');
h(3) = scatter(nan, nan, 50, c_bao,   'filled');
h(4) = scatter(nan, nan, 40, c_feng,  'filled');
legend(h, 'SWChinaCVM-2.0','USTClitho2.0','Bao20','Feng20', ...
    'Orientation','horizontal','Location',[0.5 0.01 0.03 0.03]);

f1_name = fullfile(out_dir, sprintf('%d_dt%s.png', source_id, fre));
print(gcf, f1_name, '-dpng', '-r300');
fprintf('[Save] dt scatter: %s\n', f1_name);

% =========================================================================
%  Figure 2: Radar chart of per-sector mean dt with ±σ band
% =========================================================================
az_order = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330];

% --- Reorder means/stds by az_order ---
[yao_m, yao_s]     = reorder(yao_mean,   yao_std,   dire_vals, az_order);
[zhang_m, zhang_s] = reorder(zhang_mean, zhang_std, dire_vals, az_order);
[bao_m, bao_s]     = reorder(bao_mean,   bao_std,   dire_vals, az_order);
[feng_m, feng_s]   = reorder(feng_mean,  feng_std,  dire_vals, az_order);

mean_dt = [yao_m; zhang_m; bao_m; feng_m];

figure;
set(gcf, 'position', [800 800 1000 700]);

RC = radarChart(mean_dt);
RC.RLim      = [-5, 5];
RC.RTick     = [-5, 0, 5];
RC.PropName  = {'0°','30°','60°','90°','120°','150°', ...
                '180°','210°','240°','270°','300°','330°'};
RC.ClassName = {'SWChinaCVM-2.0','USTClitho2.0','Bao20','Feng20'};
RC = RC.draw();
RC.legend();
RC.setRTick('LineWidth', 1.5, 'Color', [0,0,0]);
RC.setPropLabel('FontSize', 13);
RC.setRLabel('FontSize', 14);
RC.setPatchN(1, 'Color', c_yao,   'MarkerFaceColor', c_yao);
RC.setPatchN(2, 'Color', c_zhang, 'MarkerFaceColor', c_zhang);
RC.setPatchN(3, 'Color', c_bao,   'MarkerFaceColor', c_bao);
RC.setPatchN(4, 'Color', c_feng,  'MarkerFaceColor', c_feng);
hold on;

% --- Extract plotted coordinates from radarChart for std-band overlay ---
%   radarChart's line objects (XData length=13) are the data polygons,
%   drawn in reverse order: feng, bao, zhang, yao.
ch = get(gca, 'Children');
data_lines = [];
for ii = 1:length(ch)
    if strcmp(get(ch(ii), 'Type'), 'line') ...
       && length(get(ch(ii), 'XData')) == 13
        data_lines = [data_lines, ch(ii)];  %#ok<AGROW>
    end
end
model_xy = cell(4, 1);
model_xy{1} = [get(data_lines(4),'XData'); get(data_lines(4),'YData')]';  % yao
model_xy{2} = [get(data_lines(3),'XData'); get(data_lines(3),'YData')]';  % zhang
model_xy{3} = [get(data_lines(2),'XData'); get(data_lines(2),'YData')]';  % bao
model_xy{4} = [get(data_lines(1),'XData'); get(data_lines(1),'YData')]';  % feng
for im = 1:4
    model_xy{im} = model_xy{im}(1:12, :);   % drop closing duplicate point
end

% --- Compute data-unit -> plot-radius scale by averaging over all points ---
rLim_lo = RC.RLim(1);
scales = [];
all_mean_arr = {yao_m, zhang_m, bao_m, feng_m};
for im = 1:4
    xy = model_xy{im}; mv = all_mean_arr{im};
    for k = 1:12
        if isnan(mv(k)), continue; end
        r_plot = sqrt(xy(k,1)^2 + xy(k,2)^2);
        dv = mv(k) - rLim_lo;
        if abs(dv) > 0.1
            scales = [scales, r_plot / dv];   %#ok<AGROW>
        end
    end
end
r_per_unit = mean(scales);

% --- Draw transparent std bands per sector ---
nDir = 12;
half_sector = pi / nDir;
nArc = 30;
all_std  = {yao_s, zhang_s, bao_s, feng_s};
all_clr  = {c_yao, c_zhang, c_bao, c_feng};
alpha_val = 0.20;
std_scale = 0.3;

for im = 1:4
    xy = model_xy{im};
    s_vals = all_std{im};
    clr = all_clr{im};

    for k = 1:nDir
        if isnan(s_vals(k)), continue; end

        r_mean = sqrt(xy(k,1)^2 + xy(k,2)^2);
        theta_center = atan2(xy(k,2), xy(k,1));
        half_band = s_vals(k) * std_scale * r_per_unit;
        r_upper = r_mean + half_band;
        r_lower = max(0, r_mean - half_band);

        theta_arc = linspace(theta_center - half_sector, ...
                             theta_center + half_sector, nArc);
        x_outer = r_upper * cos(theta_arc);
        y_outer = r_upper * sin(theta_arc);
        x_inner = r_lower * cos(fliplr(theta_arc));
        y_inner = r_lower * sin(fliplr(theta_arc));

        fill([x_outer, x_inner], [y_outer, y_inner], clr, ...
             'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    end
end

% --- σ reference scale bar ---
ref_std_s = 1;
ref_width = ref_std_s * std_scale * r_per_unit * 2;
sb_x = 0.55; sb_y = -0.58; sb_w = 0.03;
fill([sb_x-sb_w, sb_x+sb_w, sb_x+sb_w, sb_x-sb_w], ...
     [sb_y-ref_width/2, sb_y-ref_width/2, sb_y+ref_width/2, sb_y+ref_width/2], ...
     [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'k', 'LineWidth', 1);
text(sb_x + sb_w + 0.02, sb_y, ['\sigma = ' num2str(ref_std_s) ' s'], ...
    'FontSize', 11, 'VerticalAlignment', 'middle');

% --- Per-model summary text (mean dt, mean σ across sectors) ---
text(-0.55, -0.75, sprintf('SWCVM:  %.3f s, %.2f s', ...
    mean(yao_m,'omitnan'), mean(yao_s,'omitnan')), 'FontSize', 10);
text( 0.05, -0.75, sprintf('USTC:   %.3f s, %.2f s', ...
    mean(zhang_m,'omitnan'), mean(zhang_s,'omitnan')), 'FontSize', 10);
text(-0.55, -0.82, sprintf('Bao20:  %.3f s, %.2f s', ...
    mean(bao_m,'omitnan'), mean(bao_s,'omitnan')), 'FontSize', 10);
text( 0.05, -0.82, sprintf('Feng20: %.3f s, %.2f s', ...
    mean(feng_m,'omitnan'), mean(feng_s,'omitnan')), 'FontSize', 10);

f2_name = fullfile(out_dir, sprintf('%d_rida%s.png', source_id, fre));
print(gcf, f2_name, '-dpng', '-r300');
fprintf('[Save] radar plot: %s\n', f2_name);

end


%% ========================================================================
%   Helper functions
% =========================================================================
function p = idx_by_az(data, az)
% Get row indices where column 8 (azimuth) equals az; returns empty if no data.
if isempty(data), p = []; return; end
p = find(data(:, 8) == az);
end

function [m, s] = mean_std(data, p)
% Compute mean and std of column 6 (dt) for selected rows
if isempty(p) || isempty(data)
    m = nan; s = nan;
else
    m = mean(data(p, 6), 'omitnan');
    s = std(data(p, 6),  'omitnan');
end
end

function v = get_col(data, p, col)
% Extract column 'col' for selected rows, empty if no data
if isempty(p) || isempty(data), v = []; else, v = data(p, col); end
end

function [m_out, s_out] = reorder(m_in, s_in, dire_vals, az_order)
% Reorder per-sector statistics to match radar chart azimuth ordering
m_out = nan(1, 12); s_out = nan(1, 12);
for k = 1:12
    idx = find(dire_vals == az_order(k), 1);
    if ~isempty(idx)
        m_out(k) = m_in(idx);
        s_out(k) = s_in(idx);
    end
end
end

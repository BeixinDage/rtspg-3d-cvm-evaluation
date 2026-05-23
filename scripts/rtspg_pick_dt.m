function rtspg_pick_dt(source_id, model_name, paths)
% RTSPG_PICK_DT  Reverse-Time Source-Point Gather dt picking
% =========================================================================
%  Pick surface-wave traveltime residuals (dt) from reverse-time
%  source-point gathers (RTSPG) constructed by convolving 3D synthetic
%  waveforms (CGFD3D outputs) with time-reversed empirical Green's
%  functions (EGFs) from ambient-noise cross-correlations.
%
%  本函数对一个源点 + 一个速度模型，在 4 个周期带内完成：
%    1) 读取已滤波的 RTSPG 道（按 ChinArray 台站）
%    2) 按方位角分 12 个扇区，扇区内按距离排序
%    3) 在零时刻附近的搜索窗内拾取粗峰
%    4) 空间一致性检查（与前面相邻道中值对比）
%    5) 振幅检查（归一化振幅阈值）
%    6) 输出 dt = peak_time - t0 (s)
%
%  Inputs:
%    source_id   - source point index in station list (e.g. 60/104/142/...)
%    model_name  - velocity model tag, one of {'yao','zhang','bao','feng'}
%                  for {SWChinaCVM-2.0, USTClitho2.0, Bao20, Feng20}
%    paths       - struct with the following fields (all absolute paths):
%                  .station_list   path to StaX1.list (台站坐标列表)
%                  .data_root      root dir of filtered RTSPG SAC files
%                                  (data_root/<model>con<band>/*.sac)
%                  .out_dt         output dir for dt text files
%                  .out_daoji      output dir for gather diagnostic PNGs
%                  .out_circle     output dir for circle PNGs
%
%  Outputs (written to disk):
%    <out_dt>/<src>_<model>_dt_<band>_v9.txt
%    <out_daoji>/<src>_<model>_daoji_<band>_v9.png
%    <out_circle>/<src>_<model>_circle_<band>_v9.png
%
%  Example:
%    paths.station_list = './data_example/StaX1.list';
%    paths.data_root    = './data_example/179prefil_con';
%    paths.out_dt       = './output/179dt';
%    paths.out_daoji    = './output/179daoji';
%    paths.out_circle   = './output/179circle';
%    rtspg_pick_dt(179, 'zhang', paths);
%
%  Note on per-source parameters (azimuthal band widths):
%    Sector half-widths in some azimuthal sectors vary between 30-50 km
%    among source points. This is because ChinArray station density is
%    uneven; sectors with sparser coverage were widened to ensure
%    sufficient sampling. See get_source_config() for source-specific
%    settings. This is an empirical choice driven by data availability,
%    not a tunable algorithmic parameter.
%
%  Algorithm reference:
%    [Your manuscript citation, JGR: Solid Earth, 2026]
% =========================================================================

% --- 输入检查 ---
if nargin < 3
    error('Usage: rtspg_pick_dt(source_id, model_name, paths)');
end
valid_models = {'yao','zhang','bao','feng'};
if ~ismember(model_name, valid_models)
    error('model_name must be one of: %s', strjoin(valid_models, ', '));
end

% --- Load source-specific config (azimuthal band widths) ---
cfg = get_source_config(source_id);
fprintf('========================================\n');
fprintf('Source ID: %d   Model: %s\n', source_id, model_name);
fprintf('========================================\n');

% =========================================================================
%  全局参数（所有源点统一）
% =========================================================================
t0 = 299;     % 零时刻采样点（对应 RTSPG 卷积输出的中心）
nt = 597;     % 总采样点数

% --- 频段配置: {name, cHW (coarse-peak half-window), dtmax, scTol} ---
%   cHW   = coarse peak search half-window (samples = seconds @ 1Hz)
%   dtmax = max allowed |dt| (s)
%   scTol = spatial consistency tolerance (s)
%   统一使用 (4, 6, 11, 15) 作为 dtmax —— 见 README §Parameters
bands = {
    %  name      cHW  dtmax  scTol
    {'5-10',    5,    4,    4.0};
    {'8-18',    9,    6,    6.0};
    {'15-35',  18,   11,   11.0};
    {'20-45',  23,   15,   15.0};
};

sc_ncan      = 4;     % number of peak candidates retained per trace
sc_ref_hw    = 6;     % reference window: up to N previous traces
min_peak_amp = 0.2;   % normalized amplitude threshold for peak acceptance

% =========================================================================
%  加载台站坐标并构造距离矩阵
% =========================================================================
station = load(paths.station_list);
NumSta  = size(station, 1);

% --- pairwise distance matrix among all stations ---
dist   = zeros(NumSta, NumSta-1);
dissta = zeros(NumSta, NumSta);
for i = 1:NumSta
    evex_i = station(i,4); evey_i = station(i,5);
    dx = evex_i - station(:,4);
    dy = evey_i - station(:,5);
    for j = 1:NumSta
        dissta(i,j) = sqrt(dx(j)^2 + dy(j)^2);
        if j <= i
            dist(i,j) = dissta(i,j);
        else
            dist(i,j-1) = dissta(i,j);
        end
    end
end
dist(:, NumSta) = [];
fprintf('[Load] %d stations\n', NumSta);

% --- 当前源点位置 ---
evex = station(source_id, 4);
evey = station(source_id, 5);

% 12 个扇区的方位角中心 (degrees, geographic convention: 0°=East)
dire_vals = [90, 270, 180, 0, 60, 240, 30, 210, 300, 120, 330, 150];

% =========================================================================
%  对每个频段执行完整 dt 拾取流程
% =========================================================================
nbands = numel(bands);

for ib = 1:nbands
    band_info = bands{ib};
    fre       = band_info{1};
    coarse_hw = band_info{2};
    dtmax     = band_info{3};
    sc_tol    = band_info{4};

    fprintf('\n========== Band %s s ==========\n', fre);
    fprintf('peakWin=±%ds, dtmax=%d, scTol=%.1f, minAmp=%.1f, refHW=%d\n', ...
            coarse_hw, dtmax, sc_tol, min_peak_amp, sc_ref_hw);

    % --- Read RTSPG SAC files for this model+band ---
    %   directory: <data_root>/<model>con<band>/*.sac
    eventdir = fullfile(paths.data_root, [model_name, 'con', fre]);
    if ~exist(eventdir, 'dir')
        fprintf('[Warn] Directory not found: %s, skipping\n', eventdir);
        continue;
    end
    cd(eventdir);
    sac_files = dir(fullfile(eventdir, '*.sac'));
    Ncon = length(sac_files);
    if Ncon == 0
        fprintf('[Warn] No SAC files in %s, skipping\n', eventdir);
        continue;
    end

    % --- Read and amplitude-normalize each trace ---
    condata_norm = zeros(Ncon, nt);
    for k = 1:Ncon
        sacdata = readsac(sac_files(k).name);
        tr = sacdata.DATA1(:)';
        amp = max(abs(tr));
        if amp > 0
            condata_norm(k,:) = tr / amp;
        end
    end
    fprintf('[Load] %d RTSPG traces from %s\n', Ncon, eventdir);

    % --- 按方位角扇区分组 + 按震中距排序 ---
    %   12 个扇区，每个扇区有自己的几何判据 (见 cfg.sector_widths)
    [stations_sec, seiscon_sec, dist_sec] = sort_into_sectors( ...
        station, condata_norm, dist, source_id, evex, evey, cfg);

    % --- 拾取窗口 ---
    tp1 = t0 - coarse_hw;
    tp2 = t0 + coarse_hw;

    data_dt_all      = cell(12,1);
    seiscon_all      = cell(12,1);
    coarse_peaks_all = cell(12,1);

    for ia = 1:12
        sc   = seiscon_sec{ia};
        st_i = stations_sec{ia};
        do_i = dist_sec{ia};

        if isempty(sc)
            data_dt_all{ia}      = [];
            seiscon_all{ia}      = [];
            coarse_peaks_all{ia} = [];
            continue;
        end
        N = size(sc, 1);

        % ---- 粗峰候选: 每道最多 sc_ncan 个，按幅值降序 ----
        cand_list = cell(N, 1);
        for ii = 1:N
            tr = sc(ii, :);
            if ~any(tr), cand_list{ii} = []; continue; end
            cand_list{ii} = find_positive_peaks_descending( ...
                tr, tp1, tp2, sc_ncan, 0.05);
            if isempty(cand_list{ii})
                fprintf('  [NoPeak] sec%d trc%2d: no positive peak in [%d,%d]\n', ...
                        ia, ii, tp1, tp2);
            end
        end

        % ---- 初始粗峰 = 每道最大幅值峰 ----
        coarse_pk = nan(N, 1);
        for ii = 1:N
            if ~isempty(cand_list{ii})
                coarse_pk(ii) = cand_list{ii}(1);
            end
        end
        n_no_cand = sum(cellfun(@isempty, cand_list));
        if n_no_cand > 0
            fprintf('  [Gate1] sec%d: %d/%d traces have no candidate\n', ...
                    ia, n_no_cand, N);
        end

        % ---- 空间一致性检查 (与前 ref_hw 道中值对比) ----
        n_before = sum(isfinite(coarse_pk));
        coarse_pk = enforce_spatial_consistency_coarse( ...
            coarse_pk, cand_list, sc_tol, t0, sc_ref_hw);
        n_after = sum(isfinite(coarse_pk));
        if n_before - n_after > 0
            fprintf('  [Gate2] sec%d: consistency rejected %d (%d->%d)\n', ...
                    ia, n_before-n_after, n_before, n_after);
        end

        % ---- 振幅检查 ----
        n_before_amp = sum(isfinite(coarse_pk));
        for ii = 1:N
            if ~isfinite(coarse_pk(ii)), continue; end
            pk_samp = round(coarse_pk(ii));
            if abs(sc(ii, pk_samp)) < min_peak_amp
                fprintf('  [AmpRej] sec%d trc%2d: amp=%.3f < %.1f, set NaN\n', ...
                        ia, ii, abs(sc(ii, pk_samp)), min_peak_amp);
                coarse_pk(ii) = nan;
            end
        end
        n_after_amp = sum(isfinite(coarse_pk));
        if n_before_amp - n_after_amp > 0
            fprintf('  [Gate3] sec%d: amplitude rejected %d (%d->%d)\n', ...
                    ia, n_before_amp-n_after_amp, n_before_amp, n_after_amp);
        end

        coarse_peaks_all{ia} = coarse_pk;
        seiscon_all{ia}      = sc;

        % ---- 输出 dt = peak - t0，并应用 dtmax 上限 ----
        dt_arr = nan(N, 1);
        for ii = 1:N
            if ~isfinite(coarse_pk(ii)), continue; end
            dt_val = coarse_pk(ii) - t0;
            if abs(dt_val) <= dtmax
                dt_arr(ii) = dt_val;
            else
                fprintf('  [DtMaxRej] sec%d trc%2d: |dt|=%.1f > %d\n', ...
                        ia, ii, abs(dt_val), dtmax);
            end
        end

        % 8 列输出: [station_info(1-5), dt(6), distance(7), azimuth(8)]
        data_dt_all{ia} = [st_i, dt_arr, do_i(:), ones(N,1)*dire_vals(ia)];

        dtv = data_dt_all{ia}(:,6);
        n_valid = sum(isfinite(dtv));
        fprintf('[sec%2d|%3ddeg] N=%d, valid=%d, dt=[%.2f,%.2f]s\n', ...
            ia, dire_vals(ia), N, n_valid, ...
            min(dtv, 'omitnan'), max(dtv, 'omitnan'));
    end

    data_dt = vertcat(data_dt_all{:});

    % --- 保存 dt 文本文件 ---
    if ~exist(paths.out_dt, 'dir'), mkdir(paths.out_dt); end
    dt_name = fullfile(paths.out_dt, ...
        sprintf('%d_%s_dt_%s_v9.txt', source_id, model_name, fre));
    save(dt_name, 'data_dt', '-ascii');
    fprintf('[Save] dt: %s\n', dt_name);

    % --- 道集诊断图 ---
    if ~exist(paths.out_daoji, 'dir'), mkdir(paths.out_daoji); end
    plot_daoji(100+ib, ...
        sprintf('v9: coarse-pick + consistency + amp check  band %s s', fre), ...
        seiscon_all, data_dt_all, dire_vals, t0, tp1, tp2, ...
        coarse_peaks_all, source_id, model_name, fre, paths.out_daoji);

    % --- circle 图 ---
    if ~exist(paths.out_circle, 'dir'), mkdir(paths.out_circle); end
    draw_circle(200+ib, data_dt, source_id, model_name, fre, ...
        cfg, paths.out_circle);
end

fprintf('\n===== Done: source %d, model %s =====\n', source_id, model_name);
end


%% ========================================================================
%   Per-source configuration (azimuthal sector half-widths)
% =========================================================================
function cfg = get_source_config(source_id)
% Returns the per-source configuration:
%
%   cfg.sector_widths  [N, S, W, E] perpendicular half-widths in meters
%                      for the 4 cardinal azimuthal sectors. These vary
%                      between 30-50 km based on local ChinArray station
%                      density (sparser sectors get a wider band to
%                      ensure sufficient station sampling).
%
%   cfg.diag_width     constant 30 km half-width for the 8 diagonal
%                      sectors (5-12); used as the perpendicular
%                      distance threshold to each sector's bisector.
%
%   cfg.circle_grid    [r_min_m, r_max_m, n_radial] defining the
%                      radial grid for the circle plot. For the 6
%                      source points reported in the manuscript these
%                      are fixed for exact reproducibility. For any
%                      other source_id this field is returned empty,
%                      which triggers automatic grid detection from
%                      the data range inside draw_circle().
%
%   cfg.circle_dtheta  azimuthal step of the circle grid in degrees
%                      (10° for all sources).
%
%  Note: To add a new source point, either add a new case below with
%  custom parameters, or simply let it fall through to the otherwise
%  branch to get automatic defaults.

switch source_id
    % --- Six source points reported in the manuscript ---
    case 60
        cfg.sector_widths = [5e4, 3e4, 4e4, 4e4];
        cfg.circle_grid   = [26000, 562000, 30];
    case 104
        cfg.sector_widths = [3e4, 5e4, 4e4, 3e4];
        cfg.circle_grid   = [28000, 640000, 26];
    case 142
        cfg.sector_widths = [3e4, 5e4, 4e4, 4e4];
        cfg.circle_grid   = [42000, 670000, 25];
    case 179
        cfg.sector_widths = [5e4, 5e4, 4e4, 4e4];
        cfg.circle_grid   = [26000, 562000, 30];
    case 203
        cfg.sector_widths = [5e4, 5e4, 4e4, 4e4];
        cfg.circle_grid   = [28000, 633000, 26];
    case 246
        cfg.sector_widths = [3e4, 5e4, 4e4, 4e4];
        cfg.circle_grid   = [41000, 620000, 23];

    % --- Generic defaults for any other source point ---
    otherwise
        fprintf(['[Info] Source %d not in manuscript set; ', ...
                 'using default sector_widths = [40,40,40,40] km ', ...
                 'and auto-detecting circle_grid from data.\n'], source_id);
        cfg.sector_widths = [4e4, 4e4, 4e4, 4e4];
        cfg.circle_grid   = [];   % empty = trigger auto-detect in draw_circle
end

% Constants shared across all sources
cfg.diag_width    = 3e4;   % diagonal sectors: 30 km half-width
cfg.circle_dtheta = 10;    % azimuthal step of circle plot grid (degrees)
end


%% ========================================================================
%   按方位角扇区分组并按距离排序
% =========================================================================
function [stations_sec, seiscon_sec, dist_sec] = sort_into_sectors( ...
    station, condata_norm, dist, source_id, evex, evey, cfg)

stations_sec = cell(12, 1);
seiscon_sec  = cell(12, 1);
dist_sec     = cell(12, 1);
for ia = 1:12
    stations_sec{ia} = []; seiscon_sec{ia} = []; dist_sec{ia} = [];
end

NumSta = size(station, 1);
w = cfg.sector_widths;   % [N, S, W, E] half-widths
wd = cfg.diag_width;

for j = 1:NumSta
    if j == source_id, continue; end
    x = station(j,4); y = station(j,5);
    sta_j = station(j,:);
    if j < source_id, con_idx = j;     d_j = dist(source_id, j);
    else,             con_idx = j - 1; d_j = dist(source_id, j-1);
    end
    tr_j = condata_norm(con_idx, :);

    % ---- Cardinal sectors: N(1) S(2) W(3) E(4) ----
    if abs(x-evex) <= w(1) && y > evey
        [stations_sec{1}, seiscon_sec{1}, dist_sec{1}] = ...
            append_one(stations_sec{1}, seiscon_sec{1}, dist_sec{1}, sta_j, tr_j, d_j);
    end
    if abs(x-evex) <= w(2) && y < evey
        [stations_sec{2}, seiscon_sec{2}, dist_sec{2}] = ...
            append_one(stations_sec{2}, seiscon_sec{2}, dist_sec{2}, sta_j, tr_j, d_j);
    end
    if x < evex && abs(y-evey) <= w(3)
        [stations_sec{3}, seiscon_sec{3}, dist_sec{3}] = ...
            append_one(stations_sec{3}, seiscon_sec{3}, dist_sec{3}, sta_j, tr_j, d_j);
    end
    if x > evex && abs(y-evey) <= w(4)
        [stations_sec{4}, seiscon_sec{4}, dist_sec{4}] = ...
            append_one(stations_sec{4}, seiscon_sec{4}, dist_sec{4}, sta_j, tr_j, d_j);
    end

    % ---- Diagonal sectors (5-12): perpendicular distance to bisector ----
    % 60°/240° (NE/SW)
    da1 = 0.5 * abs(sqrt(3)*x - y + evey - sqrt(3)*evex);
    if da1 <= wd && x > evex && y > evey
        [stations_sec{5}, seiscon_sec{5}, dist_sec{5}] = ...
            append_one(stations_sec{5}, seiscon_sec{5}, dist_sec{5}, sta_j, tr_j, d_j);
    end
    if da1 <= wd && x < evex && y < evey
        [stations_sec{6}, seiscon_sec{6}, dist_sec{6}] = ...
            append_one(stations_sec{6}, seiscon_sec{6}, dist_sec{6}, sta_j, tr_j, d_j);
    end
    % 30°/210° (NNE/SSW)
    da2 = 0.5 * abs(x - sqrt(3)*y + sqrt(3)*evey - evex);
    if da2 <= wd && x > evex && y > evey
        [stations_sec{7}, seiscon_sec{7}, dist_sec{7}] = ...
            append_one(stations_sec{7}, seiscon_sec{7}, dist_sec{7}, sta_j, tr_j, d_j);
    end
    if da2 <= wd && x < evex && y < evey
        [stations_sec{8}, seiscon_sec{8}, dist_sec{8}] = ...
            append_one(stations_sec{8}, seiscon_sec{8}, dist_sec{8}, sta_j, tr_j, d_j);
    end
    % 300°/120° (SE/NW)
    da3 = 0.5 * abs(-sqrt(3)*x - y + evey + sqrt(3)*evex);
    if da3 <= wd && x > evex && y < evey
        [stations_sec{9}, seiscon_sec{9}, dist_sec{9}] = ...
            append_one(stations_sec{9}, seiscon_sec{9}, dist_sec{9}, sta_j, tr_j, d_j);
    end
    if da3 <= wd && x < evex && y > evey
        [stations_sec{10}, seiscon_sec{10}, dist_sec{10}] = ...
            append_one(stations_sec{10}, seiscon_sec{10}, dist_sec{10}, sta_j, tr_j, d_j);
    end
    % 330°/150° (SSE/NNW)
    da4 = 0.5 * abs(-x - sqrt(3)*y + sqrt(3)*evey + evex);
    if da4 <= wd && x > evex && y < evey
        [stations_sec{11}, seiscon_sec{11}, dist_sec{11}] = ...
            append_one(stations_sec{11}, seiscon_sec{11}, dist_sec{11}, sta_j, tr_j, d_j);
    end
    if da4 <= wd && x < evex && y > evey
        [stations_sec{12}, seiscon_sec{12}, dist_sec{12}] = ...
            append_one(stations_sec{12}, seiscon_sec{12}, dist_sec{12}, sta_j, tr_j, d_j);
    end
end

% --- Sort each sector by epicentral distance (ascending) ---
for ia = 1:12
    [stations_sec{ia}, seiscon_sec{ia}, dist_sec{ia}] = ...
        sort_by_dist(stations_sec{ia}, seiscon_sec{ia}, dist_sec{ia});
end
end


%% ========================================================================
%   Find top-N positive peaks in [t1,t2], sorted by amplitude descending
% =========================================================================
function peak_locs = find_positive_peaks_descending(tr, t1, t2, ncan, min_amp)
peak_locs = [];
nt = numel(tr);
t1c = max(2, t1);
t2c = min(nt-1, t2);
if t2c < t1c, return; end

seg = tr(t1c:t2c);
is_peak = false(1, numel(seg));
for k = 2:numel(seg)-1
    if seg(k) > min_amp && seg(k) >= seg(k-1) && seg(k) >= seg(k+1)
        is_peak(k) = true;
    end
end
peak_idx = find(is_peak);
if isempty(peak_idx), return; end

[~, order] = sort(seg(peak_idx), 'descend');
peak_idx = peak_idx(order);
peak_idx = peak_idx(1:min(ncan, numel(peak_idx)));
peak_locs = (t1c - 1) + peak_idx(:);
end


%% ========================================================================
%   Spatial consistency check on coarse peaks
%
%   For each trace i (sorted by distance), compute ref_dt = median of
%   previous up-to-ref_hw valid dt values. If |dt_i - ref_dt| > sc_tol,
%   try alternative candidates (amplitude-descending order); if none
%   satisfies, set NaN.
% =========================================================================
function coarse_out = enforce_spatial_consistency_coarse( ...
    coarse_in, cand_list, sc_tol, t0, ref_hw)

coarse_out = coarse_in;
N = numel(coarse_in);
if N < 1, return; end
dtv = coarse_out - t0;

for i = 1:N
    if ~isfinite(coarse_out(i)), continue; end

    % Reference dt: 第1道与0比；其余道用前面最多 ref_hw 道的中值
    if i == 1
        ref_dt = 0;
    else
        j_start = max(1, i - ref_hw);
        prev_valid = dtv(j_start : i-1);
        prev_valid = prev_valid(isfinite(prev_valid));
        if isempty(prev_valid), ref_dt = 0;
        else,                   ref_dt = median(prev_valid);
        end
    end

    cur_dt = dtv(i);
    if abs(cur_dt - ref_dt) <= sc_tol, continue; end

    % Try alternative candidates
    cands = cand_list{i};
    if isempty(cands)
        coarse_out(i) = nan; dtv(i) = nan; continue;
    end

    found = false;
    for c = 1:numel(cands)
        cand_dt = cands(c) - t0;
        if abs(cand_dt - ref_dt) <= sc_tol
            coarse_out(i) = cands(c);
            dtv(i) = cand_dt;
            found = true;
            fprintf('  [Fix] trc%2d: dt %.2f->%.2f (ref=%.2f, cand#%d)\n', ...
                    i, cur_dt, cand_dt, ref_dt, c);
            break;
        end
    end

    if ~found
        coarse_out(i) = nan; dtv(i) = nan;
        fprintf('  [SetNaN] trc%2d: dt %.2f no candidate satisfies (ref=%.2f, tol=%.1f)\n', ...
                i, cur_dt, ref_dt, sc_tol);
    end
end
end


%% ========================================================================
%   Plot 12-sector gather diagnostic figure
% =========================================================================
function plot_daoji(figno, titlestr, seiscon_cell, data_dt_cell, ...
                    dire_vals, t0, tp1, tp2, coarse_peaks, ...
                    source_id, model_name, fre, out_dir)

plim = 269:329;
figure(figno);
set(gcf, 'position', [100 100 1400 1000], 'Color', 'w');

for ii = 1:12
    sc = seiscon_cell{ii};
    if isempty(sc), continue; end
    N  = size(sc, 1);
    dd = data_dt_cell{ii};

    subplot(4, 3, ii);
    wiggle(1:N, plim, (sc(:, plim))', 'VA');
    hold on;

    % Zero-time reference line (red dashed)
    line([0.5 N+0.5], [t0 t0], 'linestyle','--', 'Color','r', 'LineWidth',1.2);

    % Coarse peak search window bounds (blue dotted)
    line([0.5 N+0.5], [tp1 tp1], 'linestyle',':', 'Color',[0.2 0.4 0.9], 'LineWidth',0.8);
    line([0.5 N+0.5], [tp2 tp2], 'linestyle',':', 'Color',[0.2 0.4 0.9], 'LineWidth',0.8);

    % Coarse peaks (orange circles)
    if ~isempty(coarse_peaks)
        pk = coarse_peaks{ii};
        for kk = 1:N
            if isfinite(pk(kk)) && pk(kk) >= plim(1) && pk(kk) <= plim(end)
                plot(kk, pk(kk), 'o', 'MarkerSize', 4, ...
                     'MarkerFaceColor', [1 0.6 0], ...
                     'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
            end
        end
    end

    % Final dt markers (green triangles + numeric label)
    if ~isempty(dd) && size(dd,1) == N
        for kk = 1:N
            dt_kk = dd(kk, 6);
            if ~isfinite(dt_kk), continue; end
            peak_samp = t0 + dt_kk;
            if peak_samp < plim(1) || peak_samp > plim(end), continue; end
            plot(kk, peak_samp, 'v', 'MarkerSize', 6, ...
                 'MarkerFaceColor', [0.1 0.75 0.2], ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
            if dt_kk > 0.2,      clr = [0.85 0.1 0.1];
            elseif dt_kk < -0.2, clr = [0.1 0.2 0.85];
            else,                clr = [0.45 0.45 0.45]; end
            text(kk+0.25, peak_samp, sprintf('%.1f', dt_kk), ...
                 'FontSize', 6, 'Color', clr, ...
                 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        end
    end

    set(gca, 'XTick', 1:2:N);
    yticks(linspace(plim(1), plim(end), 5));
    set(gca, 'yticklabel', linspace(plim(1)-t0, plim(end)-t0, 5), ...
             'FontSize', 9, 'LineWidth', 1);
    xlim([0.5 N+0.5]);
    ylim([plim(1) plim(end)]);

    if ii >= 10, xlabel('Trace #', 'FontSize', 12); end
    if any(ii == [1 4 7 10])
        ylabel('Time relative to t_0 (s)', 'FontSize', 12);
    end
    title(sprintf('Sec%d %d° N=%d', ii, dire_vals(ii), N), 'FontSize', 10);
    hold off;
end

sgtitle(titlestr, 'FontSize', 13, 'FontWeight', 'bold');
outfile = fullfile(out_dir, ...
    sprintf('%d_%s_daoji_%s_v9.png', source_id, model_name, fre));
exportgraphics(gcf, outfile, 'Resolution', 200);
fprintf('[Save] gather plot: %s\n', outfile);
end


%% ========================================================================
%   Circle (polar) plot of dt distribution
%
%   Logic preserved verbatim from selecdata_drawdaoji179_newpipeline_v9_test.m,
%   with the only generalization being that the hardcoded Nr=30 for the
%   periodic-azimuth replication is replaced by cfg.circle_grid(3) so the
%   routine works for all source points.
%
%   If cfg.circle_grid is empty, the radial grid is auto-detected from
%   the data range — this is the fallback path for any source not in
%   the manuscript set.
% =========================================================================
function draw_circle(figno, data_dt, source_id, model_name, fre, cfg, out_dir)

datathe = data_dt(:, 8);
datar   = data_dt(:, 7);
datadt  = data_dt(:, 6);

% --- Build polar grid ---
dtheta = cfg.circle_dtheta;
the    = 0 : dtheta : 360;

if isempty(cfg.circle_grid)
    % Auto-detect from data: round to 1 km, use 25 km radial step
    valid_r = datar(isfinite(datadt));
    if isempty(valid_r)
        warning('No valid dt to determine circle grid; skipping plot');
        return;
    end
    r_min = floor(min(valid_r) / 1000) * 1000;
    r_max = ceil( max(valid_r) / 1000) * 1000;
    Nr    = max(10, round((r_max - r_min) / 25000));
    fprintf('[Info] Auto-detected circle grid: r=[%.0f,%.0f] m, Nr=%d\n', ...
            r_min, r_max, Nr);
else
    r_min = cfg.circle_grid(1);
    r_max = cfg.circle_grid(2);
    Nr    = cfg.circle_grid(3);
end

r   = linspace(r_min, r_max, Nr) * 0.001;   % km
len = numel(the) * numel(r);

% --- k-index: 2D nearest-neighbor in (theta, r_km) space ---
[The, R] = meshgrid(the, r);
A  = reshape(The, len, 1);
B  = reshape(R,   len, 1);
P  = [A, B];
PQ = [datathe, datar * 0.001];
k  = dsearchn(P, PQ);

% --- Manual accumulation and averaging (matches verified test version) ---
inter_dt = nan(len, 1);
count    = zeros(len, 1);
for i = 1:numel(k)
    if isfinite(datadt(i))
        if isnan(inter_dt(k(i)))
            inter_dt(k(i)) = datadt(i);
            count(k(i))    = 1;
        else
            inter_dt(k(i)) = inter_dt(k(i)) + datadt(i);
            count(k(i))    = count(k(i)) + 1;
        end
    end
end
multi = count > 1;
inter_dt(multi) = inter_dt(multi) ./ count(multi);

% --- Replicate index pattern (periodic-azimuth filling) ---
%     NOTE: use k (all stations) not k_g (only valid); this is what the
%     verified test version does. Replacing the hardcoded 30 (Nr for src 179)
%     with Nr to generalize across source points.
inter_dt(k + Nr)       = inter_dt(k);
inter_dt(k + 2 * Nr)   = inter_dt(k);
inter_dt = reshape(inter_dt, numel(r), numel(the));

% --- Render polar pseudocolor ---
figure(figno);
polarPcolor(r, the, inter_dt);
col = colormap('jetwr');
caxis([-10, 10]);
set(gcf, 'position', [600 600 600 600], ...
         'Color', 0.75*[1,1,1], 'InvertHardCopy', 'off');
colormap(flipud(col));
set(gca, 'color', 0.75*[1,1,1]);

filename = fullfile(out_dir, ...
    sprintf('%d_%s_circle_%s_v9.png', source_id, model_name, fre));
print(gcf, filename, '-dpng', '-r300');
fprintf('[Save] circle plot: %s\n', filename);
end


%% ========================================================================
%   Helper: append a row to (station, trace, distance) cell triplet
% =========================================================================
function [stK, scK, dK] = append_one(stK, scK, dK, sta, tr, d)
stK = [stK; sta];
scK = [scK; tr];
dK  = [dK;  d];
end

function [stK, scK, dK] = sort_by_dist(stK, scK, dK)
if isempty(dK), return; end
[dK, idx] = sort(dK, 'ascend');
stK = stK(idx, :);
scK = scK(idx, :);
end

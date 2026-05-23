% =========================================================================
%  run_analysis_example.m
%
%  Example driver for analyze_dt_statistics.m
%  Generates dt scatter + linear fit + radar chart for one source/band.
% =========================================================================
clear; clc;

repo_root = pwd;
addpath(fullfile(repo_root, 'scripts'));

% --- For one source + one band ---
source_id = 179;
fre       = '20-45';
dt_dir    = fullfile(repo_root, 'output', sprintf('%ddt', source_id));
out_dir   = fullfile(repo_root, 'figures', sprintf('%d', source_id));

analyze_dt_statistics(source_id, fre, dt_dir, out_dir);


% --- All sources + all bands (uncomment to run full set) ---
% all_sources = [60, 104, 142, 179, 203, 246];
% all_bands   = {'5-10','8-18','15-35','20-45'};
%
% for src = all_sources
%     dt_dir  = fullfile(repo_root, 'output', sprintf('%ddt', src));
%     out_dir = fullfile(repo_root, 'figures', sprintf('%d', src));
%     for ib = 1:numel(all_bands)
%         try
%             analyze_dt_statistics(src, all_bands{ib}, dt_dir, out_dir);
%         catch ME
%             fprintf('[Error] src=%d band=%s: %s\n', src, all_bands{ib}, ME.message);
%         end
%     end
% end

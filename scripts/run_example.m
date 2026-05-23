% =========================================================================
%  run_example.m
%
%  Example driver script for rtspg_pick_dt.m
%
%  This script demonstrates how to invoke the dt picking routine for
%  a single source point and a single velocity model. To process all
%  source-model combinations, wrap the call in nested loops as shown
%  in the second example block below.
% =========================================================================
clear; clc;

% ---- Set up paths (EDIT THESE TO MATCH YOUR ENVIRONMENT) ----
% Root directory of this repository
repo_root = pwd;   % or set explicitly, e.g., '/home/user/rtspg-3d-cvm-evaluation'

% Paths to data and outputs
paths.station_list = fullfile(repo_root, 'data_example', 'StaX1.list');
paths.data_root    = fullfile(repo_root, 'data_example', '179prefil_con');
paths.out_dt       = fullfile(repo_root, 'output', '179dt');
paths.out_daoji    = fullfile(repo_root, 'output', '179daoji');
paths.out_circle   = fullfile(repo_root, 'output', '179circle');

% ---- Add scripts to MATLAB path ----
addpath(fullfile(repo_root, 'scripts'));

% =========================================================================
%  Example 1: Single source + single model
% =========================================================================
source_id  = 179;
model_name = 'zhang';   % options: 'yao', 'zhang', 'bao', 'feng'

rtspg_pick_dt(source_id, model_name, paths);


% =========================================================================
%  Example 2: All sources x all models (uncomment to run full set)
% =========================================================================
% all_sources = [60, 104, 142, 179, 203, 246];
% all_models  = {'yao', 'zhang', 'bao', 'feng'};
%
% for src = all_sources
%     % Update paths for this source
%     paths.data_root  = fullfile(repo_root, 'data_example', sprintf('%dprefil_con', src));
%     paths.out_dt     = fullfile(repo_root, 'output', sprintf('%ddt', src));
%     paths.out_daoji  = fullfile(repo_root, 'output', sprintf('%ddaoji', src));
%     paths.out_circle = fullfile(repo_root, 'output', sprintf('%dcircle', src));
%
%     for im = 1:numel(all_models)
%         try
%             rtspg_pick_dt(src, all_models{im}, paths);
%         catch ME
%             fprintf('[Error] src=%d model=%s: %s\n', src, all_models{im}, ME.message);
%         end
%     end
% end

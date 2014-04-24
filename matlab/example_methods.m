% This script contains examples for reading .gct, .gctx, and .lxb files
% in MATLAB, and performing z-scoring for a data data set.
% For this script to run properly, run the script from its directory.

% add required path
addpath('lib')

% read a .gct file
ds = parse_gct('../data/log_ybio_epsilon.gct');

% read a .gctx file and inspect its contents
ds = parse_gctx('../data/modzs_n272x978.gctx');
disp(ds.mat(1:10,1:10))
disp(ds.rid(1:10))
disp(ds.cid(1:10))
disp(ds.rdesc(1:10,:))
disp(ds.cdesc(1:10,:))

% z-score a .gctx matrix
zs = robust_zscore(ds.mat);

% read an .lxb file
lxb = l1kt_parse_lxb('../data/A10.lxb');
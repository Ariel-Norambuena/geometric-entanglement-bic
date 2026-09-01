%% test_previous_paper_regression.m
% Regression test for the exact finite-mode dynamics of the previous paper.
%
% This test recomputes selected production-grid observables from the old
% Fig. 3 and compares them with data/baseline_previous_paper.json. It is
% intentionally compact: the full figure uses 100000 output times, while the
% regression only checks representative times that capture the ideal branch
% and the lambda-detuned decay.

clear; clc;

testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
addpath(fullfile(projectRoot, 'scripts'));

reproduce_previous_paper_baseline("verify");

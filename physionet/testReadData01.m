% Sample test plots
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% May 2019
%
% Note: Please use the most recent updates of the Open-Source
% Electrophysiological Toolbos (OSET) online available at: https://gitlab.com/rsameni/OSET
clear;
close all;
clc;

D = dir('D:\Users\samenir\RezaSAMENI\SourceCodes\DaniloFetalECGData\*.mat');

k = 12;
fs = 2048;

load(D(k).name, 'PORTI');
abdominal = PORTI(1:24, :);
thoracic = PORTI(25:27, :);

data = [thoracic ; abdominal];

x = data - LPFilter(data, 15.0/fs);
x = LPFilter(x, 80.0/fs);

w = jadeR(x);
s = w*x;

PlotECG(abdominal, 4, 'b', fs, 'raw abdominal signals');
PlotECG(thoracic, 3, 'r', fs, 'raw thoracic signals');
PlotECG(x, 4, 'k', fs, 'preprocessed signals');
PlotECG(s, 4, 'm', fs, 'signals after JADE');
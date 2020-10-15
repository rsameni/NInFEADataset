% Batch file converter from MAT to Binary
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% Jan 2020
%
% Note: Please use the most recent updates of the Open-Source
% Electrophysiological Toolbos (OSET) online available at: https://gitlab.com/rsameni/OSET

clear;
close all;
clc;

pth = 'D:\Dropbox\Ninfea dataset (shared)\physionet mat';
path(path, pth);
D = dir([pth '*.mat']);
fs = 2048.0;
for subject = 1 : length(D)
    ifname = D(subject).name; 
    ofname = [pth 'BinaryFormatData\' ifname(1:end-4) '.bin'];

    load(ifname, 'PORTI');
    % remove some channels if needed
    %     PORTI = PORTI(1: 27, :);
    
    % convert from MAT to Binary
    MatToBinaryFileConverter(PORTI, fs, ofname);
    
    % read-back from Binary to MAT
    [PORTI_readback, fs_readback] = ReadBinaryFile(ofname);
    
    % Check equality of read and write
    if(isequal(fs, fs_readback) && isequal(PORTI, PORTI_readback))
        disp(['write success: ' ofname]);
    else
        disp(['write failure: ' ofname]);
    end
end

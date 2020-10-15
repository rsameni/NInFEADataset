% Batch file converter from Binary to WFDB
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% Oct 2020
%
% Note: Please use the most recent updates of the Open-Source
% Electrophysiological Toolbos (OSET) online available at: https://gitlab.com/rsameni/OSET

clear;
close all;
clc;

in_folder = './BinaryFormatData';
out_folder = 'WFDBFormatData';
D = dir([in_folder '/*.bin']);
for subject = 1 : length(D)
    in_fname = D(subject).name;
    in_fname_full = [in_folder '/' in_fname];
    
    out_fname = D(subject).name;
    out_fname_full = [out_folder '/' out_fname(end-4:end)];

    % read-back from Binary to MAT
    [data, fs] = ReadBinaryFile(in_fname_full);
    [xbit] = mat2wfdb(data', ['/Users/rsameni/Documents/GitRepository/NInFEADataset/physionet/' out_fname_full], fs);

    disp(num2str(subject));
end

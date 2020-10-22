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
    out_fname_full = [out_folder '/' out_fname(1 : end-4)];
    
    % read-back from Binary to MAT
    [data, fs] = ReadBinaryFile(in_fname_full);
    %     data = data / 1000.0; % convert data into mV
    
    bit_res = 32;
    adu = 'uV';
    isquant = 0;
    isdigital = 0;
    info = {' NInFEADataset'};
    info = cat(1, info, {'The following device voltage offsets (in uV) have been removed from the data:'});
    baseline = round(median(data, 2));
    % % %     mx = max(data, [], 2);
    % % %     mn = min(data, [], 2);
    % % %     exceeds_limit = mx > 2^(bit_res-1)-1 | mn < -2^(bit_res-1)
    % % %     baseline = baseline .* exceeds_limit;
    for kk = 1 : length(baseline)
        info = cat(1, info, {['Ch' num2str(kk) '=' num2str(baseline(kk), '%d')]});
    end
    data_baseline_removed = data - baseline*ones(1, size(data, 2));
    
    %     [xbit] = mat2wfdb(data_baseline_removed', ['/Users/rsameni/Documents/GitRepository/NInFEADataset/physionet/' out_fname_full], fs, bit_res, adu, [], [], [], baseline,isquant, isdigital);
    %    [xbit] = mat2wfdb(data', out_fname(1 : end-4), fs, bit_res, adu, info);
    
    [xbit] = mat2wfdb(data_baseline_removed', out_fname(1 : end-4), fs, bit_res, adu, info);
    
    %     x = dlmread(['foo' out_fname(1 : end-4)]);
    %
    %     er = 10*log10(mean(data_baseline_removed.^2, 2)./mean((data_baseline_removed - x(:, 2:end)').^2, 2))
    
    disp(num2str(subject));
end

% % % figure
% % % hold on
% % % plot(data(1, :));
% % % plot(x(:, 2));
% % % grid
% % % % !rdsamp -r 60 > foo
% % % % x = dlmread('bar');
% % % % subplot(211)
% % % % plot(sig)
% % % % subplot(212)
% % % % plot(x(:,1),x(:,2));hold on;plot(x(:,1),x(:,3),'k');plot(x(:,1),x(:,4),'r')

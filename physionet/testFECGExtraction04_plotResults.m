% Test plot the fetal and maternal heart rate results
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% May 2019
%
% Note: Please use the most recent updates of the Open-Source
% Electrophysiological Toolbos (OSET) online available at: https://gitlab.com/rsameni/OSET

clear;
close all;
clc;

pth = 'D:\Users\samenir\RezaSAMENI\SourceCodes\DaniloFetalECGData';
D = dir([pth '\*.mat']);

% algorithm parameters
fs = 2048; % sampling frequency in Hz
f0 = 50; % notch frequency
Q = 15; % notch filter Q-factor
fm = 1.3; % approximate maternal heart rate in Hz
ITR = 4; % number of deflation iterations used for mECG cancellation
MM = 1; % number of top channels denoised by WDEN per deflation iteration
TPTR = 'rigrsure'; % WDEN threshold selection rule
SORH = 's'; % WDEN thresholding approach
SCAL = 'one'; % WDEN  multiplicative threshold rescaling
NDEN = 1; % Number of WDEN levels decrease to follow the average as it is, increase to make smoother
WNAME = 'coif5'; % WDEN mother wavelat
% mpeak_detection_smoothing_len = round(fs*0.030); % maternal HR detector window length (not always used)
fpeak_detection_smoothing_len = round(fs*0.025); % fetal HR detector window length
MultiNotchNum = 1; % number of channels to remove using multichannel notch filter
wlenL = round(0.030*fs); % before R-peak
wlenR = round(0.030*fs); % after R-peak

% subject = 5; % test data index in the data folder (1...68)
for subject = 1 : 68,
    % for subject = 38 : 68,
    ff = 2.4; % approximate fetal heart rate in Hz
    
    load([pth, '\results\results01_' D(subject).name], 'mref', 'fref0', 'fref1', 'mpeaks', 'fpeaks0', 'fpeaks1', 'fHR0', 'fHR1', 'fHR1_med', 'fHR1smoothed');

    I_mpeaks = find(mpeaks);
    mHR = 60*fs./diff(I_mpeaks);
    I_fpeaks0 = find(fpeaks0);
    I_fpeaks1 = find(fpeaks1);
    T = length(mpeaks);
    t = (0:T-1)/fs;
    
    % L1 = size(y, 1);
    % L = 4;
    % for i = 1:L1
    %     if(mod(i,L)==1 || L==1)
    %         figure;
    %     end
    %     subplot(L,1,mod(i-1,L)+1);
    %     plot(t, y(i, :),'b');
    %     hold on
    %     plot(t, yy(i, :),'r');
    %     ylabel(num2str(i));
    %     grid;
    %     if(mod(i,L)==1 || L==1)
    %         title(D(subject).name, 'Interpreter', 'none');
    %     end
    %     if(mod(i,L)==0 || L==1)
    %         xlabel('time(s)');
    %     end
    % end
    
    %     figure
    %     hold on;
    %     plot(t, mref);
    %     plot(t(I_mpeaks), mref(I_mpeaks), 'ro');
    %     grid
    %     title('maternal reference channel and extracted R-peaks');
    
    %     figure
    %     hold on;
    %     plot(t, fref0/std(fref0), 'b');
    %     plot(t, fref1/std(fref1), 'r');
    %     plot(t(I_fpeaks0), fref0(I_fpeaks0)/std(fref0), 'co');
    %     plot(t(I_fpeaks1), fref1(I_fpeaks1)/std(fref1), 'ko');
    %     legend('fref0', 'fref1', 'fpeaks0', 'fpeaks1');
    %     % plot(t, fref/std(fref));
    %     grid
    %     title(['fetal reference channels for ' D(subject).name], 'Interpreter', 'none');
    
    % PlotECG(s0, 4, 'b', fs, 'JADE on Raw Data');
    % PlotECG(abdominal, 4, 'b', fs, 'raw abdominal signals');
    % PlotECG(thoracic, 3, 'r', fs, 'raw thoracic signals');
    % PlotECG(x, 4, 'k', fs, 'preprocessed signals');
    % PlotECG(s, 4, 'm', fs, 'signals after JADE');
    % PlotECG(s_sorted, 4, 'g', fs, 'signals after sorting');
    % PlotECG(s, 4, 'm', fs, 'signals after PiCA');
    
    figure
    hold on
    plot(t(I_mpeaks(1:end-1)), mHR, 'b', 'linewidth', 2);
    plot(t(I_fpeaks0(1:end-1)), fHR0, 'r');
    plot(t(I_fpeaks1(1:end-1)), fHR1, 'g');
    plot(t(I_fpeaks1(1:end-1)), fHR1_med, 'm');
    plot(t(I_fpeaks1(1:end-1)), fHR1smoothed, 'k', 'linewidth', 2);
    grid
    legend('mHR', 'fHR0', 'fHR1', 'fHR1med', 'fHR1smoothed');
    xlabel('time(s)');
    ylabel('HR (BPM)');
    title(['HR for subject:' num2str(subject) ', fname ' D(subject).name], 'Interpreter', 'none');
end

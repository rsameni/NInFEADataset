% Fetal and maternal ECG extraction and heart rate calculation
% adapted for the NInFEA Dataset (https://physionet.org/content/ninfea/1.0.0/)
%
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% May 2019
%
% Note: Use the most recent updates of the Open-Source
% Electrophysiological Toolbos (OSET) online available at: https://gitlab.com/rsameni/OSET

clear;
close all;
clc;

pth = ''; % Replace with the correct path of the .mat files

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
mpeak_detection_smoothing_len = round(fs*0.040); % maternal HR detector window length (not always used)
fpeak_detection_smoothing_len = round(fs*0.025); % fetal HR detector window length
MultiNotchNum = 1; % number of channels to remove using multichannel notch filter
wlenL = round(0.030*fs); % before R-peak
wlenR = round(0.030*fs); % after R-peak

% subject = 5; % test data index in the data folder (1...68)
% for subject = 38 : 68
for subject = 1 : 68
    ff = 2.4; % approximate fetal heart rate in Hz
    
    % load data
    load(D(subject).name, 'PORTI');
    abdominal = PORTI(1:24, :); % 1:10*fs
    thoracic = PORTI(25:27, :);
    
    % remove DC channels (if any). DC channels are later problematic, resulting in zer eigenvalues
    abdominal = abdominal(std(abdominal, 1, 2)~=0, :);
    thoracic = thoracic(std(thoracic, 1, 2)~=0, :);
    
    data = [thoracic ; abdominal];
    T = size(data, 2);
    
    % preprocess the data (wide-band mode)
    x_abd = abdominal - LPFilter(abdominal, 7.0/fs);
    x_abd = LPFilter(x_abd, 150.0/fs);
    
    % preprocess the data (narrow-band mode)
    x = data - LPFilter(data, 15.0/fs);
    x = LPFilter(x, 80.0/fs);
    mref = x(1, :);
    
    % powerline notch filter
    
    % option 1: bypass the notch filter
    %     y_abd = x_abd;
    %     y = x;

    % option 2: A simple Notch filter:
    Wc = f0/(fs/2);
    BW = Wc/Q;
    [b_notch_filter, a_notch_filter] = iirnotch(Wc, BW);
    y = zeros(size(x));
    for i = 1:size(x, 1)
        y(i, :) = filtfilt(b_notch_filter, a_notch_filter, x(i, :));
    end
    y_abd = zeros(size(x_abd));
    for i = 1:size(x, 1)
        y_abd(i, :) = filtfilt(b_notch_filter, a_notch_filter, x_abd(i, :));
    end
    
    % option 3: a multichannel notch filter
    % y_abd = MultichannelNotchFilter(x_abd, f0, Q, fs); y_abd = y_abd(1:end-MultiNotchNum, :);
    % y = MultichannelNotchFilter(x, f0, Q, fs); y = y(1:end-MultiNotchNum, :);
    
    % maternal R-peak detection
    % option 1: simple method
    mpeaks = PeakDetection(mref, fm/fs);
    % option 2: advanced method
    % [mHR mpeaks] = HRCalculation2(mref, ff, fs, 1, mpeak_detection_smoothing_len, 'trmean');
    
    I_mpeaks = find(mpeaks);
    mHR = 60*fs./diff(I_mpeaks);
    
    % mECG cancellation by deflation
    [tt0, tt1] = SynchPhaseTimes2(mpeaks);
    yy = PeriodicDeflDecompositionWDEN(y, ITR, MM, tt0, tt1, TPTR, SORH, SCAL, NDEN, WNAME);
    
    % abdominal_pre_processed = LPFilter(abdominal-LPFilter(abdominal, 5.0/fs), 300.0/fs);
    % W = jadeR(abdominal_pre_processed);
    % s0 = W*abdominal_pre_processed;
    % s0 = LPFilter(s0-LPFilter(s0, 15.0/fs), 40.0/fs);
    
    % BSS
    % option 1: JADE
    W = jadeR(yy);
    s = W*yy;
    % option 2: FastICA
    % [s, A, w] = fastica (xx, 'approach', 'defl', 'displayMode', 'off', 'lastEig', 24, 'numOfIC', 20);
    % [s, A, w] = fastica (xx, 'approach', 'defl', 'displayMode', 'off', 'lastEig', 20, 'numOfIC', 10);
    
    % A two-run fetal ECG detector
    for run = 1:2 % repeat twice to make the fetal heart rate guess more accurate
        % ind = 1./ChannelIndex9(s, ff, fs); % Channel selection based on variance around average beat
        ind = ChannelIndex10(s, ff, fs); % Channel selection based on fixed template matched filter
        [sorted_indexes, II] = sort(ind, 1, 'ascend');
        s_sorted = s(II, :);
        fref0 = s_sorted(1, :);
        % fpeaks = PeakDetection(fref, ff/fs);
        [fhr0, fpeaks0] = HRCalculation2(fref0, ff, fs, 1, fpeak_detection_smoothing_len, 'trmean');
        %     [fhr0 fpeaks0] = HRCalculation4(fref0, ff, fs, 1, fpeak_detection_smoothing_len, 'trmean');
        I_fpeaks0 = find(fpeaks0);
        ff = fs/median(diff(I_fpeaks0));
    end
    disp([D(subject).name ', ff = ' num2str(ff)]);
    
    % fetal ECG template extraction
    segment_width = wlenL + wlenR + 1;
    X = zeros(length(I_fpeaks0), segment_width);
    for k = 1 : length(I_fpeaks0)
        start = max(I_fpeaks0(k)-wlenL, 1);
        stop = min(I_fpeaks0(k)+wlenR, T);
        xx = [zeros(1, start - (I_fpeaks0(k)-wlenL)), fref0(start : stop), zeros(1, I_fpeaks0(k)+wlenR - T)];
        X(k, :) = xx;% - fref0(I_fpeaks0(k));
    end
    avg_fECG_beat = RWAverage(reshape(X, length(I_fpeaks0), segment_width));
    fpeaks0 = PeakDetection4(fref0, fs, avg_fECG_beat, ff);
    I_fpeaks0 = find(fpeaks0);
    
    % PiCA using fetal R-peaks
    % [s, W, A] = PiCA(x_abd, fpeaks, mpeaks);
    % [s, W, A] = PiCA(xx, fpeaks0);
    [s, W, A] = PiCA(y_abd, fpeaks0);
    fref1 = s(1, :);
    c1 = corrcoef(fref0, fref1);
    fref1 = sign(c1(1,2))*fref1;
    
    % % % [ind fpeaks1] = ChannelIndex10(fref1, ff, fs); % Channel selection based on fixed template matched filter
    %     [fhr1 fpeaks1] = HRCalculation2(fref1, ff, fs, 1, fpeak_detection_smoothing_len, 'trmean'); %fpeaks1 = squeeze(fpeaks1)';
    % % % [fhr1 fpeaks1] = HRCalculation3(s(1:2, :), ff, fs, 1, fpeak_detection_smoothing_len, 'trmean');
    % % % [fhr1 fpeaks1] = HRCalculation4(fref1, ff, fs, 1, fpeak_detection_smoothing_len, 'trmean');
    
    % fetal ECG template extraction
    segment_width = wlenL + wlenR + 1;
    X = zeros(length(I_fpeaks0), segment_width);
    t_segment = (0 : segment_width-1)/fs;
    for k = 1 : length(I_fpeaks0)
        start = max(I_fpeaks0(k)-wlenL, 1);
        stop = min(I_fpeaks0(k)+wlenR, T);
        xx = [zeros(1, start - (I_fpeaks0(k)-wlenL)), fref1(start : stop), zeros(1, I_fpeaks0(k)+wlenR - T)];
        X(k, :) = xx;% - fref1(I_fpeaks0(k));
    end
    %     [avg_fECG_beat vr] = RWAverage(reshape(X(ch, :, :), length(peak_indexes), segment_width));
    avg_fECG_beat = RWAverage(reshape(X, length(I_fpeaks0), segment_width));
    fpeaks1 = PeakDetection4(fref1, fs, avg_fECG_beat, ff);
    I_fpeaks1 = find(fpeaks1);
    
    % HR refinement
    wlen = 2; % # of beats before and after the current beat used for smoothing
    th = 10; % HR error detection threshold in BPM
    fHR0 = 60*fs./diff(I_fpeaks0);
    fHR1 = 60*fs./diff(I_fpeaks1);
    fHR1_med = zeros(1, length(fHR1));
    fHR1smoothed = fHR1;
    hrlen = length(fHR1);
    for i = 1 : hrlen
        index = max(i-wlen, 1) : min(i+wlen, hrlen);
        fHR1_med(i) = median(fHR1(index));
        if(abs(fHR1(i) - fHR1_med(i)) > th)
            k = find(I_fpeaks0 >= I_fpeaks1(i), 1);
            if(~isempty(k) && k <= length(fHR0))
                fHR1smoothed(i) = fHR0(k);
            else
                fHR1smoothed(i) = fHR1_med(i);
            end
        end
    end
    
    save([pth, '\results\results01_' D(subject).name], 'mref', 'fref0', 'fref1', 'mpeaks', 'fpeaks0', 'fpeaks1', 'fHR0', 'fHR1', 'fHR1_med', 'fHR1smoothed');
    
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
    title(['HR for ' D(subject).name], 'Interpreter', 'none');
end

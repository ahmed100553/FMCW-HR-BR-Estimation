%% Initial
close all
clear
clc

%% Select Participant Number ( 1 - 50 )
%Participant#5 have a high accuracy
%Participant#9 have high accuracy
ParticipantNum = 2; 
if ( ParticipantNum < 1 || ParticipantNum > 50 )
    disp('Wrong Participant Number')
    return;
end
fprintf('Processing Participant %d\n', ParticipantNum);

%% Set Parameters
datastart = 21;
dataend = 260;
para = f_Parameter(); 
% Define expected heart rate and breathing rate physiological limits (BPM)
para.hr_bpm_limits = [40 160]; 
para.br_bpm_limits = [15 35];   

% Sliding average filter window size
para.sliding_avg_window_size = 4; 
para.sliding_avg_passes = 1; 

% Bandpass filter order
para.bpf_order = 4;

para.expected_dist_min_m = 0.30; 
para.expected_dist_max_m = 1.5; 

% Option to handle movement artifacts
para.handle_movement_artifacts = true; % Set to true to NaN estimates during movement


%% Read CSV Data
disp('Reading CSV data...');
% read radar rawdata
try
    if ~isdeployed 
        script_path = fileparts(mfilename('fullpath')); 
        addpath(fullfile(script_path, '../FMCW Radar/Rawdata/')); % Adjust if your path is different
        addpath(fullfile(script_path, '../Nihon Kohden/Heart Rate & Breathing Rate/'));
    end
    
    fid = fopen(['Rawdata_' num2str(ParticipantNum) '.csv']);
    if fid == -1
        error(['Could not open Rawdata file for Participant ' num2str(ParticipantNum) '. Check path and filename. Current directory: ' pwd]);
    end
    tstream = textscan(fid, '%d', 'Delimiter', ',');
    tstream = tstream{1,1};
    streamsize = size(tstream,1)/para.AntNum;
    Rawdata = zeros(para.AntNum, streamsize);
    for i = 1:para.AntNum
        Rawdata(i,:) = tstream((i-1)*streamsize+1:i*streamsize);
    end
    fclose(fid);
    disp('Rawdata loaded.');
catch ME
    disp('Error reading radar raw data:');
    disp(ME.message);
    return;
end

% read Nihon Kohden results (Reference)
Ref_BR_processed = []; Ref_HR_processed = []; 
Ref_HR_data_raw = []; 
try
    Ref_BR_data = readmatrix(['Ref_Breath_' num2str(ParticipantNum) '.csv']);
    Ref_HR_data_raw = readmatrix(['Ref_Heart_' num2str(ParticipantNum) '.csv']); 
    disp('Reference BR/HR data loaded.');
    disp('Reference data loaded. Processing will be done after STFT time vector is known.');
catch ME
    disp('Warning: Could not read or process reference BR/HR data.');
    disp(ME.message);
end

%% Reshape Radar RawData & Calculate Duration

RawData_reshaped = reshape(Rawdata, para.AntNum, para.adcsample*para.chirploops, para.datalength);
% Calculate total duration
para.total_duration_sec = para.datalength / para.fps;
disp(['Radar data reshaped. Final para.datalength = ' num2str(para.datalength) ...
      ', Total Duration = ' num2str(para.total_duration_sec) ' seconds.']);

%% 1. Range FFT
disp('Performing Range FFT...');
Range_FFT_output = f_RangeFFT(RawData_reshaped, para); 
disp('Range FFT complete.');

%% 2. Eliminate Static Interference (MTI)
disp('Performing MTI...');
MTI_output = f_MTI(Range_FFT_output, para); 
disp('MTI complete.');
figure('Name', ['Participant ' num2str(ParticipantNum) ' - MTI Output (User Func)']);
imagesc(1:para.datalength, (0:para.FFTSize/2-1)*para.RangeBin, squeeze(abs(MTI_output(1,:,:))) ); 
title('Range vs. Slow Time (Magnitude - After MTI)');
xlabel('Frame Index (Slow Time)'); ylabel('Distance (m)'); colorbar;
drawnow;

%% 3. Movement Measure (Needed for Artifact Handling)
disp('Calculating Movement Measure...');
[Total_Movement_Sum, Total_MovementFlag] = f_MovementMeasure(MTI_output); 
plot(Total_Movement_Sum(datastart*para.fps:dataend*para.fps),'color','b')
datasize = size(Total_Movement_Sum(datastart*para.fps:dataend*para.fps), 1);
k = 0 : 50*para.fps : datasize;
set(gca,'xtick', k)
set(gca, 'xticklabel', round(k./ para.fps))
legend('Movement','fontsize',14)
xlabel('Time (sec)', 'fontsize', 14);
ylabel('Movement Quantity', 'fontsize', 14);
xlim([0 dataend-datastart+1] .* para.fps)
disp('Movement Measure calculated.');

%% --- ADDED: Visualize Range Information ---
disp('Plotting Range Information for Bin Selection Aid...');
num_range_bins_mti = size(MTI_output, 2); % Number of bins in MTI output (FFTSize/2)
range_axis_m = (0:num_range_bins_mti-1) * para.RangeBin;

% Calculate Average Range FFT Magnitude (using first FFTSize/2 bins)
avg_range_fft_mag = mean(abs(squeeze(Range_FFT_output(1, 1:num_range_bins_mti, :))), 2);

% Calculate MTI Magnitude Variance
mti_magnitudes_ant1 = abs(squeeze(MTI_output(1,:,:))); 
variance_per_bin_mti = var(mti_magnitudes_ant1, 0, 2); % Variance along time

figure('Name', ['Participant ' num2str(ParticipantNum) ' - Range Bin Visualization']);
% Plot Average FFT Magnitude
subplot(2,1,1);
plot(range_axis_m, avg_range_fft_mag);
title('Average Range FFT Magnitude (Ant 1)');
xlabel('Distance (m)'); ylabel('Average Magnitude'); grid on;
xlim([range_axis_m(1) range_axis_m(end)]);

% Plot MTI Magnitude Variance
subplot(2,1,2);
plot(range_axis_m, variance_per_bin_mti);
hold on;
% Indicate expected range if variance method was intended
min_bin_idx_viz = round(para.expected_dist_min_m / para.RangeBin) + 1; 
max_bin_idx_viz = round(para.expected_dist_max_m / para.RangeBin) + 1;
min_bin_idx_viz = max(1, min_bin_idx_viz); 
max_bin_idx_viz = min(num_range_bins_mti, max_bin_idx_viz); 
xline(range_axis_m(min_bin_idx_viz), 'k--', 'Label', 'Min Expected Dist');
xline(range_axis_m(max_bin_idx_viz), 'k--', 'Label', 'Max Expected Dist');
title('MTI Magnitude Variance vs. Range (Ant 1)');
xlabel('Distance (m)'); ylabel('Variance'); grid on;
xlim([range_axis_m(1) range_axis_m(end)]);
hold off; 
drawnow;
%% --- Robust range-bin selection---------------------------------------

% (1) power summed over receivers
powRX     = squeeze( sum( abs(MTI_output) , 1) );    % (bins × time)
varBins   = var( powRX.' , 0 , 1 );                  % variance per bin

% (2) restrict to feasible human zone
searchMask          = range_axis_m >= para.expected_dist_min_m& range_axis_m <= para.expected_dist_max_m;
varBins(~searchMask)= 0;

% (3) power-floor mask to kill pure noise bins
meanPow   = mean(powRX,2);                           % avg power per bin
powThresh = 20*log10(max(meanPow)) - 40;             % −25 dB rule-of-thumb
powerMask = ( 20*log10(meanPow) > powThresh ).';
varBins(~powerMask)=0;

% (4) pick best-variance bin
[~, idxSorted]     = sort(varBins,'descend');
selected_range_bin_idx = idxSorted(1);

fprintf('[RANGE] picked bin %3d  (%.2f m)  var %.2e\n', ...
         selected_range_bin_idx, range_axis_m(selected_range_bin_idx), ...
         varBins(selected_range_bin_idx));

% --- show top-3 contenders for sanity check
figure('name','Variance check'); 
plot(range_axis_m, varBins); grid on; hold on;
stem(range_axis_m(idxSorted(1:3)), varBins(idxSorted(1:3)),'r','filled');
xlabel('Distance (m)'); ylabel('Variance'); legend('Var','top-3');
%% Processing the phase of the recieved signal
complex_signal_selected_bin = squeeze(Range_FFT_output(1, selected_range_bin_idx, :));
processed_phase = f_PhaseProcessing_Paper(complex_signal_selected_bin, para);
disp('Phase Extraction and Preprocessing complete.');

% --- Plot: VISUAL INSPECTION of Processed Phase ---
time_axis_phase = (0:length(processed_phase)-1) / para.fps; % Create time axis
figure('Name', ['Participant ' num2str(ParticipantNum) ' - Processed Phase Signal (Bin ' num2str(selected_range_bin_idx) ')']);
plot(time_axis_phase, processed_phase);
title(['Processed Phase Signal Used for Estimation (Bin ' num2str(selected_range_bin_idx) ')']);
xlabel('Time (s)'); ylabel('Phase (rad)'); grid on; xlim([time_axis_phase(1) time_axis_phase(end)]);
mean_ref_hr_approx = mean(Ref_HR_data_raw,'omitnan');
if ~isnan(mean_ref_hr_approx)
    expected_period_s = 60 / mean_ref_hr_approx;
    y_lims = ylim;
    text(time_axis_phase(1)+5, y_lims(2)*0.9, sprintf('Approx Ref HR Period: %.2f s', expected_period_s), 'Color', 'red');
end
drawnow;

%% 4-A.  Base-band & physiological band-passes  --------------------------
fs = para.fps;

% --- remove DC + drift below 0.05 Hz ---
[bHP,aHP]   = butter(2, 0.05/(fs/2), 'high');
phase_hp    = filtfilt(bHP,aHP, processed_phase);

resp_band  =  para.br_bpm_limits / 60;   % e.g. [10 35] → [0.167 0.583] Hz
[bBR,aBR]   = butter(6, resp_band/(fs/2), 'bandpass');
phase_br    = filtfilt(bBR,aBR, phase_hp);

heart_band =  para.hr_bpm_limits / 60;   % e.g. [40 160] → [0.67 2.67] Hz
% keep only 0.9–2.3 Hz (≈ 54–138 bpm) to dump respiration
heart_band(1) = max(heart_band(1), 0.90);
heart_band(2) = min(heart_band(2), 2.30);
[bHR,aHR] = butter(6, heart_band/(fs/2), 'bandpass');   
phase_hr    = filtfilt(bHR,aHR, phase_hp);

%% 4-B.  STFT on CLEAN bands  -------------------------------------------

parsHR = struct('stft_win',6,'overlap',0.80,'nfft_fac',8);
[HR_bpm_t,  tHR] = f_STFTRidge(phase_hr, fs, heart_band, parsHR);

parsBR = struct('stft_win',6,'overlap',0.90,'nfft_fac',8);
[BR_bpm_t,  tBR] = f_STFTRidge(phase_br, fs, resp_band,  parsBR);

% reference alignment
if ~isempty(Ref_HR_data_raw), RefHR_t = Ref_HR_data_raw(round(tHR)); end
if ~isempty(Ref_BR_data),     RefBR_t = Ref_BR_data(round(tBR));     end

figure('Name',['P' num2str(ParticipantNum) ' – Clean-phase HR / BR']);
subplot(2,1,1);
plot(tHR, HR_bpm_t,'b.-'); hold on;
if exist('RefHR_t','var'), plot(tHR,RefHR_t,'r--'); end
ylim(para.hr_bpm_limits+[-10 10]); grid on;
ylabel('Heart-rate (bpm)'); legend('Radar','Reference');

subplot(2,1,2);
plot(tBR, BR_bpm_t,'b.-'); hold on;
if exist('RefBR_t','var'), plot(tBR,RefBR_t,'r--'); end
ylim(para.br_bpm_limits+[-5 5]); grid on;
xlabel('Time (s)'); ylabel('Breathing-rate (bpm)'); legend('Radar','Reference');

fprintf('STFT on clean bands HR %.1f bpm | BR %.1f bpm\n', mean(HR_bpm_t), mean(BR_bpm_t));
%% ------------------------------------------------------------------------
%  6.  SUMMARY ‒ mean radar vs. mean reference (clean-band STFT)
% -------------------------------------------------------------------------
disp('Plotting summary …');

% ---------- compute means (omit any NaNs) -------------------------------
meanRefHR = mean(RefHR_t ,'omitnan');
meanRefBR = mean(RefBR_t ,'omitnan');
meanRadHR = mean(HR_bpm_t,'omitnan');
meanRadBR = mean(BR_bpm_t,'omitnan');

figure('Name',['P' num2str(ParticipantNum) ' – Mean HR / BR comparison'],...
       'Position',[100 100 700 300]);   % wide layout

%% ---------------- HEART-RATE summary ------------------------------------
subplot(1,2,1)
bar(1 , meanRefHR,'FaceColor',[.8 .2 .2],'DisplayName','Reference'); hold on
bar(2 , meanRadHR,'FaceColor',[.2 .2 .8],'DisplayName','Radar (clean-band STFT)');
ylabel('Heart-rate (bpm)')
title('Mean heart-rate')
xticks([1 2]); xticklabels({'Ref','Radar'})
ylim([para.hr_bpm_limits(1)-10 , para.hr_bpm_limits(2)+10])
grid on; legend location south

%% ---------------- BREATH-RATE summary -----------------------------------
subplot(1,2,2)
bar(1 , meanRefBR,'FaceColor',[.8 .2 .2],'DisplayName','Reference'); hold on
bar(2 , meanRadBR,'FaceColor',[.2 .2 .8],'DisplayName','Radar (clean-band STFT)');
ylabel('Breathing-rate (bpm)')
title('Mean breathing-rate')
xticks([1 2]); xticklabels({'Ref','Radar'})
ylim([para.br_bpm_limits(1)-5 , para.br_bpm_limits(2)+5])
grid on; legend location south

%% ------------------------------------------------------------------------
%  6-bis.  CDF of absolute error  (Radar vs. Reference)
% -------------------------------------------------------------------------
% --- keep only samples where BOTH radar & reference are available -------
goodHR = ~isnan(HR_bpm_t) & ~isnan(RefHR_t);
goodBR = ~isnan(BR_bpm_t) & ~isnan(RefBR_t);

absErrHR = abs(HR_bpm_t(goodHR) - RefHR_t(goodHR));
absErrBR = abs(BR_bpm_t(goodBR) - RefBR_t(goodBR));

figure('Name',['P' num2str(ParticipantNum) ' – CDF of absolute error'], ...
       'Position',[100 450 700 300]);

% ---------------- HEART-RATE CDF ----------------------------------------
subplot(1,2,1)
[fHR,xHR] = ecdf(absErrHR);
stairs(xHR, fHR,'LineWidth',1.6,'Color',[.1 .4 .8]); hold on
xline(5 ,'--k','5 bpm');              % helper reference line
xlabel('|error|  (bpm)')
ylabel('F(|error| ≤ x)')
title('Heart-rate error CDF')
grid on; ylim([0 1]); xlim([0 max(8,max(absErrHR)+1)])

% ---------------- BREATH-RATE CDF ---------------------------------------
subplot(1,2,2)
[fBR,xBR] = ecdf(absErrBR);
stairs(xBR, fBR,'LineWidth',1.6,'Color',[.1 .6 .1]); hold on
xline(2 ,'--k','2 bpm');              % helper reference line
xlabel('|error|  (bpm)')
ylabel('F(|error| ≤ x)')
title('Breathing-rate error CDF')
grid on; ylim([0 1]); xlim([0 max(6,max(absErrBR)+1)])

%% ------------------------------------------------------------------------
%  6-ter.  CDF of the measured rates themselves  (Radar vs Reference)
% -------------------------------------------------------------------------
% -------- build clean vectors where both signals exist -------------------
goodHR = ~isnan(HR_bpm_t) & ~isnan(RefHR_t);
goodBR = ~isnan(BR_bpm_t) & ~isnan(RefBR_t);

radarHR = HR_bpm_t(goodHR);   refHR = RefHR_t(goodHR);
radarBR = BR_bpm_t(goodBR);   refBR = RefBR_t(goodBR);

% -------- empirical CDFs -------------------------------------------------
[f_rHR ,x_rHR ] = ecdf(radarHR);
[f_RHR ,x_RHR ] = ecdf(refHR);
[f_rBR ,x_rBR ] = ecdf(radarBR);
[f_RBR ,x_RBR ] = ecdf(refBR);

figure('Name',['P' num2str(ParticipantNum) ' – CDF : Radar vs Reference'], ...
       'Position',[120 120 720 320]);

% ---------------- HEART-RATE --------------------------------------------
subplot(1,2,1)
stairs(x_RHR, f_RHR,'k-' ,'LineWidth',1.8,'DisplayName','Reference'); hold on
stairs(x_rHR, f_rHR,'b--','LineWidth',1.8,'DisplayName','Radar');
xlabel('Heart-rate (bpm)')
ylabel('F(X ≤ x)')
title('Heart-rate CDF')
grid on; legend('Location','southeast')
ylim([0 1])

% ---------------- BREATH-RATE -------------------------------------------
subplot(1,2,2)
stairs(x_RBR, f_RBR,'k-' ,'LineWidth',1.8,'DisplayName','Reference'); hold on
stairs(x_rBR, f_rBR,'g--','LineWidth',1.8,'DisplayName','Radar');
xlabel('Breathing-rate (bpm)')
ylabel('F(X ≤ x)')
title('Breathing-rate CDF')
grid on; legend('Location','southeast')
ylim([0 1])
%% ========================================================================
%  BATCH EVALUATION of the 50 participants
% ========================================================================
para = f_Parameter();                          

% ----------- pre-allocate result arrays ---------------------------------
Nsub       = 5;                                % total subjects
MAE_HR     = nan(Nsub,1);
RMSE_HR    = nan(Nsub,1);
RHO_HR     = nan(Nsub,1);
MAE_BR     = nan(Nsub,1);
RMSE_BR    = nan(Nsub,1);
RHO_BR     = nan(Nsub,1);

allErrHR   = [];                                % pooled errors (for CDF)
allErrBR   = [];

fprintf('\n  subj |  MAE_HR   RMSE_HR   ρ_HR  |  MAE_BR   RMSE_BR   ρ_BR\n');
fprintf(' ------+--------------------------------+-------------------------\n');

for p = 1:Nsub
    try
        [radHR, refHR, radBR, refBR] = runOneSubject(p, para, true);

        % ---------- make sure length & NaN handling is consistent -------
        goodHR     = ~isnan(radHR) & ~isnan(refHR);
        goodBR     = ~isnan(radBR) & ~isnan(refBR);
        if ~any(goodHR) || ~any(goodBR)
            warning('P%02d: no valid samples - skipped',p);  continue
        end

        errHR      = radHR(goodHR) - refHR(goodHR);
        errBR      = radBR(goodBR) - refBR(goodBR);

        MAE_HR(p)  = mean(abs(errHR));
        RMSE_HR(p) = sqrt(mean(errHR.^2));
        RHO_HR(p)  = corr(radHR(goodHR).', refHR(goodHR).','rows','complete');

        MAE_BR(p)  = mean(abs(errBR));
        RMSE_BR(p) = sqrt(mean(errBR.^2));
        RHO_BR(p)  = corr(radBR(goodBR).', refBR(goodBR).','rows','complete');

        allErrHR   = [allErrHR  ; abs(errHR(:))];
        allErrBR   = [allErrBR  ; abs(errBR(:))];

        fprintf('  %3d  | %7.2f  %7.2f  %+.2f | %7.2f  %7.2f  %+.2f\n', ...
                p, MAE_HR(p), RMSE_HR(p), RHO_HR(p), ...
                   MAE_BR(p), RMSE_BR(p), RHO_BR(p));
    catch ME
        warning('P%02d : %s', p, ME.message);
    end
end

% ----------- overall means (ignoring skipped subjects) ------------------
mMAE_HR  = mean(MAE_HR ,'omitnan');
mRMSE_HR = mean(RMSE_HR,'omitnan');
mRHO_HR  = mean(RHO_HR ,'omitnan');

mMAE_BR  = mean(MAE_BR ,'omitnan');
mRMSE_BR = mean(RMSE_BR,'omitnan');
mRHO_BR  = mean(RHO_BR ,'omitnan');

fprintf('\n>>>  Global-mean  HR  |  MAE %.2f  RMSE %.2f  ρ %.2f\n', ...
         mMAE_HR, mRMSE_HR, mRHO_HR);
fprintf(  '>>>  Global-mean  BR  |  MAE %.2f  RMSE %.2f  ρ %.2f\n', ...
         mMAE_BR, mRMSE_BR, mRHO_BR);

%% ----------------------------------------------------------------------
%  PLOT 1 ‒ Box-and-whisker of MAE across subjects
% -----------------------------------------------------------------------
figure('Name','Radar vs Reference - MAE per participant', ...
       'Position',[80 80 560 350]);
boxchart([ones(Nsub,1); 2*ones(Nsub,1)], ...
         [MAE_HR ; MAE_BR],'GroupByColor',[ones(Nsub,1);2*ones(Nsub,1)])
set(gca,'XTick',[1 2],'XTickLabel',{'Heart-rate','Breathing-rate'})
ylabel('MAE (bpm)')
title('Distribution of absolute error across participants')
grid on

%% ----------------------------------------------------------------------
%  PLOT 2 ‒ Global CDF of absolute error
% -----------------------------------------------------------------------
figure('Name','Global CDF of |error|','Position',[680 80 540 350]);
[fH,xH]=ecdf(allErrHR);   stairs(xH,fH,'b','LineWidth',1.7,'DisplayName','Heart-rate'); hold on
[fB,xB]=ecdf(allErrBR);   stairs(xB,fB,'g','LineWidth',1.7,'DisplayName','Breathing-rate');
xline(5,'--k','5 bpm');   xline(2,'--k','2 bpm');
xlabel('|error|  (bpm)');  ylabel('F(|error| ≤ x)')
title('Empirical CDF of absolute error  (all windows, all subjects)')
grid on; legend location southeast
ylim([0 1])

%% ----------------------------------------------------------------------
%  TABLE ‒ ranked by heart-rate MAE  (top-10 best)
% -----------------------------------------------------------------------
% [~,rankHR] = sort(MAE_HR,'ascend');
% best10     = rankHR(1:10);
% T = table(best10(:), MAE_HR(best10), RMSE_HR(best10), RHO_HR(best10), ...
%           MAE_BR(best10), RMSE_BR(best10), RHO_BR(best10), ...
%      'VariableNames',{'Subj','MAE_HR','RMSE_HR','rho_HR','MAE_BR','RMSE_BR','rho_BR'});
% disp('Top-10 participants on heart-rate MAE:');  disp(T)

%% 5-A.  Wavelet decomposition (MODWT) – band-limited phase  -------------
disp('Wavelet MODWT decomposition…');
[HR_bpm_wav, BR_bpm_wav, H_sig_wav, B_sig_wav, ...
 HR_bpm_t,   BR_bpm_t,   t_HR,      t_BR]      = ...
         f_VitalSigns_WaveletRobust(phase_hp, para);
fprintf('Wavelet global HR %.1f bpm | BR %.1f bpm\n', HR_bpm_wav, BR_bpm_wav);

% ---------- spike-robust smoothing (zero-lag) --------------------------- %%% MOD
HR_bpm_t = sgolayfilt( medfilt1(HR_bpm_t,3,'omitnan','truncate') , 3 , 9 );
BR_bpm_t = sgolayfilt( medfilt1(BR_bpm_t,3,'omitnan','truncate') , 3 , 9 );

% (optional sanity plot of the two sub-bands)
figure('Name',['P' num2str(ParticipantNum) ' – Wavelet bands']);
plot((0:numel(H_sig_wav)-1)/para.fps, H_sig_wav,'r'); hold on;
plot((0:numel(B_sig_wav)-1)/para.fps, B_sig_wav,'b');
legend('Heart band','Breath band'); xlabel('Time (s)'); grid on;

% ---------- 5-B.  Time-varying plots  -----------------------------------

RefHR_w = Ref_HR_data_raw( round( t_HR ) );   
RefBR_w = Ref_BR_data    ( round( t_BR ) );

figure('Name',['Participant ' num2str(ParticipantNum) ' – Wavelet + STFT HR/BR']);
subplot(2,1,1);
plot(t_HR,HR_bpm_t,'b.-','DisplayName','Radar HR'); hold on
if exist('RefHR_w','var'), plot(t_HR,RefHR_w,'r--','DisplayName','Reference'); end
ylabel('Heart-rate (bpm)'); grid on
ylim(para.hr_bpm_limits + [-10 10]); legend

subplot(2,1,2);
plot(t_BR,BR_bpm_t,'b.-','DisplayName','Radar BR'); hold on
if exist('RefBR_w','var'), plot(t_BR,RefBR_w,'r--','DisplayName','Reference'); end
ylabel('Breathing-rate (bpm)'); xlabel('Time (s)'); grid on
ylim(para.br_bpm_limits + [-5 5]); legend
drawnow

% ---------- 5-D.  Optional: show wavelet bands for sanity ---------------
t_vec = (0:numel(phase_hp)-1)/para.fps;
figure('Name',['Participant ' num2str(ParticipantNum) ' – Wavelet bands']);
subplot(3,1,1), plot(t_vec,phase_hp,'k'); title('HP-cleaned phase'); grid on
subplot(3,1,2), plot(t_vec,H_sig_wav,'r'); title('Heart-band (MODWT)'); grid on
subplot(3,1,3), plot(t_vec,B_sig_wav,'b'); title('Breath-band (MODWT)'); grid on


%% ------------------------------------------------------------------------
%  6-ter.  CDF of the measured rates themselves  (Radar vs Reference)
% -------------------------------------------------------------------------
% -------- build clean vectors where both signals exist -------------------
goodHR = ~isnan(HR_bpm_t) & ~isnan(RefHR_w);
goodBR = ~isnan(BR_bpm_t) & ~isnan(RefBR_w);

radarHR = HR_bpm_t(goodHR);   refHR = RefHR_w(goodHR);
radarBR = BR_bpm_t(goodBR);   refBR = RefBR_w(goodBR);

% -------- empirical CDFs -------------------------------------------------
[f_rHR ,x_rHR ] = ecdf(radarHR);
[f_RHR ,x_RHR ] = ecdf(refHR);
[f_rBR ,x_rBR ] = ecdf(radarBR);
[f_RBR ,x_RBR ] = ecdf(refBR);

figure('Name',['P' num2str(ParticipantNum) ' – CDF : Radar vs Reference'], ...
       'Position',[120 120 720 320]);

% ---------------- HEART-RATE --------------------------------------------
subplot(1,2,1)
stairs(x_RHR, f_RHR,'k-' ,'LineWidth',1.8,'DisplayName','Reference'); hold on
stairs(x_rHR, f_rHR,'b--','LineWidth',1.8,'DisplayName','Radar');
xlabel('Heart-rate (bpm)')
ylabel('F(X ≤ x)')
title('Heart-rate CDF')
grid on; legend('Location','southeast')
ylim([0 1])

% ---------------- BREATH-RATE -------------------------------------------
subplot(1,2,2)
stairs(x_RBR, f_RBR,'k-' ,'LineWidth',1.8,'DisplayName','Reference'); hold on
stairs(x_rBR, f_rBR,'g--','LineWidth',1.8,'DisplayName','Radar');
xlabel('Breathing-rate (bpm)')
ylabel('F(X ≤ x)')
title('Breathing-rate CDF')
grid on; legend('Location','southeast')
ylim([0 1])
%% --- Function Definitions ---
function para = f_Parameter()
    %% Parameters Setting
    para.AntNum = 4;
    para.freqslope = 40.845001220703125e6*1e6; 
    para.samplerate = 3000e3;
    para.bw =3.746303561822511e9;
    para.chirploops= 2;
    para.adcsample = 256;
    para.startfreq =60.25e9;
    para.c = 3e8;
    para.lambda = para.c/(para.bw/2 + para.startfreq); 
    para.rangeresol = para.c/(2*para.bw);
    para.rangemax = (para.samplerate*para.c)/(2*para.freqslope);
    para.tc = para.bw/para.freqslope; 
    para.FFTSize = 2^10; 
    para.RangeBin = (para.c * para.samplerate) / (2*para.freqslope*para.FFTSize) ; 
    para.fps = 20; 
    para.datalength = para.fps * 60 * 5; 
end

function OutputData = f_RangeFFT(InputData, para)
    OutputData = zeros(para.AntNum, para.FFTSize, para.datalength);
    for i1 = 1:para.AntNum 
        for i2 = 1:para.datalength 
            OutputData(i1,:,i2) = fft(squeeze(InputData(i1,:,i2)), para.FFTSize);
        end
    end
end

function OutputData = f_MTI(InputData, para)
    Clutter = zeros(para.AntNum, para.FFTSize, para.datalength); 
    OutputData = zeros(para.AntNum, para.FFTSize/2, para.datalength); 
    alpha = 0.01; 
    for i1 = 1:para.AntNum 
        for i2 = 2:para.datalength 
            Clutter(i1,:,i2) = alpha .* InputData(i1,:,i2) + (1-alpha) .* Clutter(i1,:,i2-1);
            OutputData(i1,:,i2) = InputData(i1,1:para.FFTSize/2,i2) - Clutter(i1,1:para.FFTSize/2,i2);
        end
    end
end

function [Total_Movement_Sum,Total_MovementFlag] = f_MovementMeasure(MTI)
        MTI_Sum_Rx_1 = squeeze(sum(abs(MTI(1,:,:)),2)); 
        Movement_Rx_1 = abs(MTI_Sum_Rx_1(2:end) - MTI_Sum_Rx_1(1:end-1));
        Movement_Rx_1 = [0; Movement_Rx_1(:)]; 
        
        MTI_Sum_Rx_2 = squeeze(sum(abs(MTI(2,:,:)),2));
        Movement_Rx_2 = abs(MTI_Sum_Rx_2(2:end) - MTI_Sum_Rx_2(1:end-1));
        Movement_Rx_2 = [0; Movement_Rx_2(:)];
        
        MTI_Sum_Rx_3 = squeeze(sum(abs(MTI(3,:,:)),2));
        Movement_Rx_3 = abs(MTI_Sum_Rx_3(2:end) - MTI_Sum_Rx_3(1:end-1));
        Movement_Rx_3 = [0; Movement_Rx_3(:)];
        
        MTI_Sum_Rx_4 = squeeze(sum(abs(MTI(4,:,:)),2));
        Movement_Rx_4 = abs(MTI_Sum_Rx_4(2:end) - MTI_Sum_Rx_4(1:end-1));
        Movement_Rx_4 = [0; Movement_Rx_4(:)];
        
        Total_Movement_Sum =  Movement_Rx_1 + Movement_Rx_2 + Movement_Rx_3 +  Movement_Rx_4;
        
        if length(Total_Movement_Sum) >= 100
            Total_Movement_Sum(1:100) = 0;
        else
            Total_Movement_Sum(:) = 0; 
        end
        Total_MovementFlag = Total_Movement_Sum >= 850000;
end

function processed_phase = f_PhaseProcessing_Paper(complex_signal_selected_bin, para)
    phase_raw = angle(complex_signal_selected_bin);
    phase_unwrapped = unwrap(phase_raw);
    current_phase_signal = detrend(phase_unwrapped, 'linear'); 
    processed_phase = current_phase_signal; 
    for p = 1:para.sliding_avg_passes
        processed_phase = movmean(processed_phase, para.sliding_avg_window_size);
    end
end

% ========================================================================
% Ridge-tracking STFT for a *single* band-limited signal
% ========================================================================
function [bpm_out, tvec] = f_STFTRidge(x, fs, band_Hz, p)
% x          : real signal (row or col)
% fs         : sample rate
% band_Hz    : [fmin fmax] in Hz of the band you expect a single dominant tone
% p          : struct with fields
%                • stft_win  (sec)
%                • overlap   (0-1)
%                • nfft_fac  (power-of-2 multiple, eg 4)

win   = round(p.stft_win * fs);
olap  = round(win * p.overlap);
nfft  = 2^nextpow2(win * p.nfft_fac);
[S,F,tvec] = spectrogram(x, hann(win), olap, nfft, fs);
P = abs(S).^2;

fmin = band_Hz(1); fmax = band_Hz(2);
keep = F>=fmin & F<=fmax;

bpm_out = NaN(1,numel(tvec));
last_bpm = NaN;

for k = 1:numel(tvec)
    slice = P(keep,k);               % power in band
    Fk    = F(keep);                 % same size

    % simple ridge: choose global max in band, optionally near last
    if ~isnan(last_bpm)
        idxNear = abs(Fk*60 - last_bpm) <= 8;  % ±8 bpm guard
        if any(idxNear), slice(~idxNear) = 0; end
    end

    [~,ix] = max(slice);
    if ~isempty(ix) && slice(ix) > 0
        last_bpm     = Fk(ix)*60;
        bpm_out(k)   = last_bpm;
    end
end

% bridge gaps, light smoothing
bpm_out = fillmissing(bpm_out,'pchip');
bpm_out = movmedian(bpm_out,3);
end

function [radHR, refHR, radBR, refBR] = runOneSubject(ParticipantNum, para, dbg)
% RUNONESUBJECT  – run *exactly* the pipeline in vital_sign_ahmed.m
%
% Usage
%   [radHR, refHR, radBR, refBR] = runOneSubject(ID , PARA , DBG)
%
% Inputs
%   ID     : integer 1…50 (which participant)
%   DBG    : true / false   show all diagnostic figures (default = true)
%
% Outputs (row-vectors, common length T)
%   radHR  : radar-derived heart-rate          [bpm]
%   refHR  : reference heart-rate              [bpm]
%   radBR  : radar-derived breathing-rate      [bpm]
%   refBR  : reference breathing-rate          [bpm]

% -------------------------------------------------------------------------
if nargin<2 || isempty(para), para = f_Parameter(); end
if nargin<3, dbg = true; end
assert(isscalar(ParticipantNum)&&ParticipantNum>=1&&ParticipantNum<=50,...
       'runOneSubject:BadID','ID must be 1…50');

if ~dbg                                   % hide figures in batch mode
    oldVis   = get(0,'DefaultFigureVisible');
    set(0,'DefaultFigureVisible','off');
    cobj = onCleanup(@()set(0,'DefaultFigureVisible',oldVis)); 
end

datastart = 21;
dataend   = 260;

%% Set Parameters
datastart = 21;
dataend = 260;
% Define expected heart rate and breathing rate physiological limits (BPM)
para.hr_bpm_limits = [40 160]; 
para.br_bpm_limits = [15 35];   

% Sliding average filter window size
para.sliding_avg_window_size = 4; 
para.sliding_avg_passes = 1; 

% Bandpass filter order
para.bpf_order = 4;

para.expected_dist_min_m = 0.35; 
para.expected_dist_max_m = 1.5; 

% Option to handle movement artifacts
para.handle_movement_artifacts = true; % Set to true to NaN estimates during movement

%% Read CSV Data
disp('Reading CSV data...');
script_path = fileparts(mfilename('fullpath'));
rawDir = fullfile(script_path,'../FMCW Radar/Rawdata/');
refDir = fullfile(script_path,'../Nihon Kohden/Heart Rate & Breathing Rate/');
if isfolder(rawDir), addpath(rawDir); end      % suppress “folder not found”
if isfolder(refDir), addpath(refDir); end

fid = fopen(['Rawdata_' num2str(ParticipantNum) '.csv'],'r');
assert(fid~=-1,'Could not open Rawdata_%d.csv',ParticipantNum);
tstream = textscan(fid,'%d','Delimiter',','); fclose(fid);
tstream = tstream{1};                               

streamsize = numel(tstream)/para.AntNum;
Rawdata = zeros(para.AntNum,streamsize);
for i = 1:para.AntNum
    Rawdata(i,:) = tstream((i-1)*streamsize+1:i*streamsize);
end
%disp('Rawdata loaded.');

Ref_BR_data     = readmatrix(['Ref_Breath_' num2str(ParticipantNum) '.csv']);
Ref_HR_data_raw = readmatrix(['Ref_Heart_'  num2str(ParticipantNum) '.csv']);

%% Reshape Radar RawData & Calculate Duration
RawData_reshaped = reshape(Rawdata, para.AntNum,...
                           para.adcsample*para.chirploops,...
                           para.datalength);
para.total_duration_sec = para.datalength / para.fps;

%% 1. Range FFT
Range_FFT_output = f_RangeFFT(RawData_reshaped, para);

%% 2. MTI
MTI_output = f_MTI(Range_FFT_output, para);            

%% 3. Movement measure (only for the plot – identical)
[Total_Movement_Sum,Total_MovementFlag] = f_MovementMeasure(MTI_output);
if dbg
    figure; plot(Total_Movement_Sum(datastart*para.fps:dataend*para.fps),'b');
    title('Movement measure'); grid on;
end

%% Robust range-bin selection (variance rule)
range_axis_m = (0:size(MTI_output,2)-1)*para.RangeBin;
powRX   = squeeze(sum(abs(MTI_output),1));
varBins = var(powRX.',0,1);
searchMask          = range_axis_m>=para.expected_dist_min_m & range_axis_m<=para.expected_dist_max_m;
varBins(~searchMask)= 0;
meanPow   = mean(powRX,2);
powThresh = 20*log10(max(meanPow))-40;
powerMask = (20*log10(meanPow)>powThresh).';
varBins(~powerMask)=0;
[~,idxSorted] = sort(varBins,'descend');
selected_range_bin_idx = idxSorted(1);

%% Phase processing
fs = para.fps;  
complex_signal_selected_bin = squeeze(Range_FFT_output(1,...
                                        selected_range_bin_idx,:));
processed_phase = f_PhaseProcessing_Paper(complex_signal_selected_bin,para);

[bHP,aHP] = butter(6,0.05/(para.fps/2),'high');
phase_hp  = filtfilt(bHP,aHP,processed_phase);

resp_band  = para.br_bpm_limits / 60;
[bBR,aBR]  = butter(6,resp_band /(para.fps/2),'bandpass');
heart_band =  para.hr_bpm_limits / 60;   % e.g. [40 160] → [0.67 2.67] Hz
% keep only 0.9–2.3 Hz (≈ 54–138 bpm) to dump respiration
heart_band(1) = max(heart_band(1), 0.90);
heart_band(2) = min(heart_band(2), 2.30);
[bHR,aHR] = butter(6, heart_band/(fs/2), 'bandpass');   
phase_br   = filtfilt(bBR,aBR,phase_hp);
phase_hr   = filtfilt(bHR,aHR,phase_hp);

%% STFT on clean bands
parsHR = struct('stft_win',6,'overlap',0.80,'nfft_fac',8);
parsBR = struct('stft_win',6,'overlap',0.90,'nfft_fac',8);

[HR_bpm_t,tHR] = f_STFTRidge(phase_hr,fs,heart_band,parsHR);
[BR_bpm_t,tBR] = f_STFTRidge(phase_br,fs,resp_band ,parsBR);

if ~isempty(Ref_HR_data_raw), RefHR_t = Ref_HR_data_raw(round(tHR)); else RefHR_t = nan(size(tHR)); end
if ~isempty(Ref_BR_data)    , RefBR_t = Ref_BR_data    (round(tBR)); else RefBR_t = nan(size(tBR)); end

%% Package outputs (row & equal length)
radHR = HR_bpm_t(:)';  refHR = RefHR_t(:)';
radBR = BR_bpm_t(:)';  refBR = RefBR_t(:)';

%% Optional debug plot
if dbg
    figure('Name',sprintf('P%d – clean-band HR / BR',ParticipantNum));
    subplot(2,1,1);
    plot(tHR,radHR,'b.-'); hold on; plot(tHR,refHR,'r--');
    ylabel('Heart-rate (bpm)'); legend Radar Reference; grid on;
    subplot(2,1,2);
    plot(tBR,radBR,'b.-'); hold on; plot(tBR,refBR,'r--');
    ylabel('Breathing-rate (bpm)'); xlabel('Time (s)');
    legend Radar Reference; grid on;
end
end
function [HR_bpm_peak , BR_bpm_peak , H_sig , B_sig , ...
          HR_bpm_t     , BR_bpm_t     , t_HR , t_BR] = ...
          f_VitalSigns_WaveletRobust( phase_hp , para )
% f_VitalSigns_WaveletRobust  –  MODWT-based separation of HR / BR
%
%   [HR_peak,BR_peak,Hsig,Bsig,HR_t,BR_t,tHR,tBR] = ...
%         f_VitalSigns_WaveletRobust( phase_hp , para )
%
%   INPUT
%       phase_hp    : 1×N   DC-/drift-removed phase (rad)   (see §4-A)
%       para        : struct  (needs fields  fps , hr_bpm_limits , br_bpm_limits)
%
%   OUTPUT
%       HR_bpm_peak : scalar   global-peak HR over the whole record  [bpm]
%       BR_bpm_peak : scalar   global-peak BR over the whole record  [bpm]
%       H_sig       : 1×N      reconstructed heart-band waveform     (rad)
%       B_sig       : 1×N      reconstructed breath-band waveform    (rad)
%       HR_bpm_t    : 1×T      time-varying HR (STFT ridge)          [bpm]
%       BR_bpm_t    : 1×T      time-varying BR (STFT ridge)          [bpm]
%       t_HR , t_BR : 1×T      time instants (s) of the STFT points
%
%   ------------------------------------------------------------------

% ---------- constants -------------------------------------------------
fs       = para.fps;

wname    = 'sym8';                    % << was 'db4'

Jmax     = floor(log2(numel(phase_hp))) - 1;
wt       = modwt(phase_hp, wname, Jmax);        % size Jmax×N
f_vec    = fs ./ 2.^((1:Jmax)+1);               % pseudo-centre freq

% ---------- choose levels automatically -------------------------------
hr_lim_hz = para.hr_bpm_limits/60;
br_lim_hz = para.br_bpm_limits/60;

bw        = f_vec./2;                            % half-bandwidth ≈ fs/2^(j+2)
overlapHR = min(f_vec+bw,hr_lim_hz(2)) - max(f_vec-bw,hr_lim_hz(1));
overlapBR = min(f_vec+bw,br_lim_hz(2)) - max(f_vec-bw,br_lim_hz(1));
idx_HR    = find(overlapHR > 0.5*bw);            % ≥50 % overlap
idx_BR    = find(overlapBR > 0.5*bw);            % ≥50 % overlap

% --- power screening keeps only energetic levels -----------------------
pow_lvl = mean(abs(wt).^2,2);
idx_HR  = idx_HR(pow_lvl(idx_HR) >= 0.40*max(pow_lvl(idx_HR)));   % 40 %
idx_BR  = idx_BR(pow_lvl(idx_BR) >= 0.25*max(pow_lvl(idx_BR)));   % 25 %

% ---------- reconstruct ------------------------------------------------
keepHR          = false(Jmax+1,1); keepHR(idx_HR) = true; keepHR(end)=true;
keepBR          = false(Jmax+1,1); keepBR(idx_BR) = true; keepBR(end)=true;
H_sig = imodwt(wt(keepHR,:),wname);
B_sig = imodwt(wt(keepBR,:),wname);

[bH,aH] = butter(4, hr_lim_hz/(fs/2),'bandpass');
[bB,aB] = butter(4, br_lim_hz/(fs/2),'bandpass');
H_sig   = filtfilt(bH,aH,H_sig);
B_sig   = filtfilt(bB,aB,B_sig);

% ---------- global-peak PSDs with harmonic shield  ---------------------
[pxB,fB] = periodogram(B_sig,[],[],fs);          % BR first
[~,kB]   = max(pxB);  BR_bpm_peak = fB(kB)*60;

[pxH,fH] = periodogram(H_sig,[],[],fs);
[~,idx]  = sort(pxH,'descend');                  % peaks sorted high→low
for m = 1:numel(idx)                            
    cand = fH(idx(m))*60;
    if abs(cand-2*BR_bpm_peak) > 5 && ...
       abs(cand-3*BR_bpm_peak) > 5 && ...
       abs(cand-4*BR_bpm_peak) > 5
        HR_bpm_peak = cand;     break            
    end
end
% (fallback – if all harmonics rejected)
if ~exist('HR_bpm_peak','var'), HR_bpm_peak = fH(idx(1))*60; end

% ---------- STFT-ridge for time-varying estimates ---------------------
parsHR.stft_win = 10;  parsHR.overlap = 0.80; parsHR.nfft_fac = 8;
parsBR.stft_win = 10;  parsBR.overlap = 0.90; parsBR.nfft_fac = 8;

[HR_bpm_t , t_HR] = f_STFTRidge( H_sig , fs , hr_lim_hz , parsHR );
[BR_bpm_t , t_BR] = f_STFTRidge( B_sig , fs , br_lim_hz , parsBR );

if nargout==0
    t = (0:numel(phase_hp)-1)/fs;
    figure('name','MODWT bands / STFT');
    subplot(3,1,1), plot(t,phase_hp,'k'); title('Input phase (HP)'); grid on
    subplot(3,1,2), plot(t,H_sig,'r',t,B_sig,'b'); grid on
    legend('Heart-band','Breath-band');
    subplot(3,1,3), plot(t_HR,HR_bpm_t,'r.-',t_BR,BR_bpm_t,'b.-');
    ylabel('bpm'); legend('HR','BR'); grid on
end
end

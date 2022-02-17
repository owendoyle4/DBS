% Original Author: Mike Cohen
% Modified by Owen Doyle
% This script is adapted from Analyzing Neural Time Series Data by Mike
% Cohen to handle this specific ECoG and DBS dataset instead of EEG data.

%Code is designed for one patient, one condition data.
pat = 'BI025';          %Patient ID
condition = 'contra';   %Condition 

% load NDATA
cd(strcat('/Users/owen/Box/DBS_data/sdata_files/',pat));
load(strcat('NDATA_',pat,'_lfp_15.mat'));

%Pick electrode
elec = 4;

baseline = (NDATA.(condition).cond_spec_baseline(:,elec))';
trial_data = (mean(squeeze(NDATA.(condition).data(elec,:,:)),2))'; %average of trials

%concatenate the baseline to the average of trials
length_data = 9300;
times = -(length(baseline)-1):300;

figure;
plot(baseline)
figure;
plot(trial_data)
%% Parameters and convolution
% wavelet parameters
min_freq = 2;
max_freq = 128; 
num_frex = 30;  %~frequency resolution

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/NDATA.sr:1;                                            
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters for DATA (use next-power-of-2)
n_wavelet     = length(time);                   
n_data        = length_data;                       
n_convolution = n_wavelet+n_data-1;             
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; % ^ # of cycles -> ^ width of gausian -> ^ width of wavelet [assuming constant freq]
                   %low cycles -> high temporal precesion convolution

                   
% initialize output time-frequency data FOR AVERAGE
tf_average_data = zeros(length(frequencies),length_data);

num_trials = length(squeeze(NDATA.(condition).data(1,1,:)));

for trial = 1: num_trials
    
    trial_data = squeeze(NDATA.(condition).data(elec,:,trial));
    data = [baseline, trial_data];
    
    % FFT of data (note: this doesn't change on frequency iteration)
    fft_data = fft(data,n_conv_pow2);  

    % initialize output time-frequency data
    tf_data = zeros(length(frequencies),length(data));
    
    for fi=1:length(frequencies)

        % create wavelet and get its FFT
        wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
        fft_wavelet = fft(wavelet,n_conv_pow2);

        % run convolution
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
        convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

        % put power data into time-frequency matrix
        tf_data(fi,:) = abs(convolution_result_fft).^2;
    end
    
    tf_average_data = tf_average_data + tf_data;
end

tf_average_data = tf_average_data/num_trials;

%% plot results for DATA
ytickskip = 2:4:num_frex; % This will be explained in the text.

%linear color scale - poor visualization due to 1/f
% figure
% imagesc(times,[],tf_data)
% set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal');
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')
% title('Color limit of 0 to 25 (linear)')
% colorbar

%log color scale
figure
subplot(211)
imagesc(times,[],log10(tf_average_data))
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Color limit of (log_1_0 units) [full baseline + trial average]')
colorbar

subplot(212)
imagesc(times,[],log10(tf_average_data))
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal', 'xlim', [-100 length(trial_data)]);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Color limit of (log_1_0 units) [last 100ms of baseline + trial average]')
colorbar

%% Figure 18.3 - dB Converted

% define baseline period
%baselinetime = [ -length(baseline) 0 ]; % in ms

% convert baseline window time to indices
baselineidx = [1 length(baseline)];

% dB-correct
baseline_power = mean(tf_average_data(:,baselineidx(1):baselineidx(2)),2);
dbconverted = 10*log10( bsxfun(@rdivide,tf_average_data,baseline_power));
% FYI: the following lines of code are equivalent to the previous line:
% dbconverted = 10*( bsxfun(@minus,log10(tf_average_data),log10(baseline_power)));
% dbconverted = 10*log10( tf_average_data ./ repmat(baseline_power,1,length(data)) );
% dbconverted = 10*( log10(tf_average_data) - log10(repmat(baseline_power,1,length(data))) );

dbconverted_short = dbconverted(:,9001:9300);
times_short = 1:300;
figure
contourf(times_short,frequencies,dbconverted_short,40,'linecolor','none')
set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log', 'xlim', [0 length(trial_data)]); %, 'clim', [-40 40]);
title(strcat("dB Converted [last 100ms of baseline + trial average] ", pat, " ",NDATA.chanlocs(elec)));
xlabel("Time (ms)");
ylabel("Frequency");
colorbar;

%% Figure 18.4

time2plot = 150; % in ms, 1-300 is the trial, everything < 0 is baseline

[~,timeidx] = min(abs(times-time2plot));

% plot frequencies
figure
subplot(211)
plot(frequencies,tf_average_data(:,timeidx))
title([ 'Power spectrum at ' num2str(timeidx) ' ms' ])
ylabel('Raw power (\muV^2)')
xlabel('Frequency (Hz)')
 
subplot(212)
plot(frequencies,dbconverted(:,timeidx))
ylabel('Baseline-normalized power (dB)')
xlabel('Frequency (Hz)')

%% Figure 18.6
% real data: percent change vs. baseline division
figure
baseline_power = mean(tf_average_data(:,baselineidx(1):baselineidx(2)),2);
pctchange = 100 * (tf_average_data-repmat(baseline_power,1,length(data)))./ repmat(baseline_power,1,length(data)); %****
subplot(221)
baselinediv = tf_average_data ./ repmat(baseline_power,1,length(data)); %****
plot(dbconverted(1:5:end),baselinediv(1:5:end),'.') % don't need all the datapoints to make a point
xlabel('DB'), ylabel('Baseline division')

% dB vs. baseline division
subplot(222)
plot(pctchange(1:5:end),baselinediv(1:5:end),'.')
xlabel('Percent change'), ylabel('Baseline division')

% Z-transform vs. percent change
subplot(223)
baseline_power = tf_average_data(:,baselineidx(1):baselineidx(2));
baselineZ = (tf_average_data-repmat(mean(baseline_power,2),1,size(tf_average_data,2))) ./ repmat(std(baseline_power,[],2),1,size(tf_average_data,2)); %****
plot(baselineZ(1:5:end),pctchange(1:5:end),'.')
xlabel('Z-transform'), ylabel('Percent change')

% Z-transform vs. dB
subplot(224)
plot(baselineZ(1:5:end),dbconverted(1:5:end),'.')
xlabel('Z-transform'), ylabel('DB')


%% Figure 18.7

dbconverted_short = dbconverted(:,9001:9300);
pctchange_short = pctchange(:,9001:9300);
baselinediv_short = baselinediv(:,9001:9300);
baselineZ_short = baselineZ(:,9001:9300);

times_short = 1:300;

figure

% plot dB-converted power
subplot(221)
imagesc(times_short,[],dbconverted_short)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal'); %,'xlim',[-100 300]) %,'clim',[-10 10])
title('dB change from baseline')
colorbar

% plot percent-change
subplot(222)
imagesc(times_short,[],pctchange_short)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal'); %,'xlim',[-100 300]) %,'clim',[-500 500])
title('Percent change from baseline')
colorbar

% divide by baseline
subplot(223)
imagesc(times_short,[],baselinediv_short)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal'); %,'xlim',[-100 300]) %,'clim',[-7.5 7.5])
title('Divide by baseline')
colorbar

% z-transform
subplot(224)
imagesc(times_short,[],baselineZ_short)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal'); %,'xlim',[-100 300]) %,'clim',[-3.5 3.5])
title('Z-transform')
colorbar


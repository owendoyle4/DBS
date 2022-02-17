function [] = wholeEMGpreprocessing(pat_ids, lpf_freq)

    if pat_ids == "all"
        patient_ids = ["BI020","BI025","BI029","BI031","BI033","BI035","BI037","BI038","BI039","BI040","BI041","BI099"];
    else
        patient_ids = pat_ids;
    end

    for pat_idx = 1:length(patient_ids)
        pat = patient_ids(pat_idx);
        disp(pat)
        % EMG is already downsampled to 1000Hz with a 500Hz LPF
        sr = 1000;
        cd('/Users/owen/Box/DBS_data/functions/preprocessing');
        [emgs] = getEMGs(pat);

        figure;
        plot(emgs);

        % 60Hz notch filter
        cd('/Users/owen/Box/DBS_data/functions/');
        filt_deg = 3;
        line_noise_freq=60;
        LNFwindow=10;
        tic; emgs = remove_line_noise(emgs,line_noise_freq,sr,LNFwindow); toc % A larger window makes for a much smaller "notch filter" effect.

        figure;
        plot(emgs);

        % 50 Hz HPF (4th degree filter)
        for num = 1:size(emgs,2)
            emgs(:,num) = HPF(emgs(:,num), sr, 50, 4); %(data, sr, freq, filter order)
        end

        emgs_50hz = emgs;

        % Full wave rectified (absolute value)
        for num = 1:size(emgs,2)
            emgs(:,num) = abs(emgs(:,num));
        end

        % Initialize different smoothing techniques
        emgs_lpf = emgs;
        emgs_rms = emgs;
        emgs_peak= emgs;

        % 5 or 10 Hz LPF both EMGs (4th degree filter)
        for num = 1:size(emgs,2)
            emgs_lpf(:,num) = LPF(emgs(:,num), sr, lpf_freq, 4);
        end

        % RMS Envelope both EMGs (only upper envelope)
        for num = 1:size(emgs,2)
            [emgs_rms(:,num),~] = envelope(emgs(:,num), 50, 'rms');
        end

        % Peak Evelope both EMGs (only upper envelope)
        for num = 1:size(emgs,2)
            [emgs_peak(:,num),~] = envelope(emgs(:,num), 50, 'peak');
        end

        %% plot contra (1) or ipsi (2) for all smoothing methods
        cd(strcat('/Users/owen/Box/DBS_data/sdata_files/',pat));
        move_temp = load(strcat('move_', pat, '_temp.mat'));
        
        figure;
        % plot(emgs_50hz); hold on
        % plot(emgs_peak(:,1), 'g', 'Linewidth', 1);hold on
        % plot(emgs_rms(:,1), 'b', 'Linewidth', 1); hold on
        plot(emgs_lpf(:,1), 'r', 'Linewidth', 1);
        legend;
        
        title('Final LPF vs. RMS Envelope vs. Peak Envelope');
        

        %%
        figure;
        plot(emgs_lpf(:,1), 'r', 'Linewidth', 1);
        figure;
        plot(emgs_lpf(:,2), 'r', 'Linewidth', 1);

        %%
        %update the EMGs used in the encoding model
        cd('/Users/owen/Box/DBS_data/functions/preprocessing');
        replaceEMG(pat, emgs_lpf, strcat(string(lpf_freq),"hz"))
        replaceEMG(pat, emgs_rms, 'rms')
        replaceEMG(pat, emgs_peak,'peak')
    end
end
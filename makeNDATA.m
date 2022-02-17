%Creates and saves NDATA. I would comment out the last section on saving if
%you are tweaking/adding features. This way everyone that shares the Box
%isn't spammed.
function [NDATA] = makeNDATA(pat, window)

%initialize NDATA struct
NDATA = struct;

%include patient ID
NDATA.patID = pat;

%load patient files
cd(strcat('/Users/owen/Box/DBS_data/sdata_files/',pat)) %Cd = current directory 
load(strcat('SDATA_',pat,'.mat'));

addpath /Users/owen/Box/DBS_data/functions;

%sample rate
NDATA.sr = SDATA.info.sr;

%channel labels
NDATA.chanlocs = SDATA.info.chan_labels(1:14);

addpath /Users/owen/Box/DBS_data/metadata/metadata_10_24_20/;
load('metadata.mat');
NDATA.chanlocs(1:length(metadata.(pat).ecog_label)) = metadata.(pat).ecog_label;

%Baseline duration, 4.5 second or 4500 ms
bl_dur = 4499; 

%Levels of the structure
conditions = ["contra", "ipsi", "bimanual"];
blocks = ["b1","b2"];
hands = ["contra", "ipsi"];

for condition = conditions
    %Create variable cond to index within SDATA.event files 
    if condition == "bimanual"
            cond = "bilat";
        else
            cond = condition;
    end

    % Select rows for desired condition to obtain ready, set, go, stop times
    cond_index = contains(SDATA.events.event_names(:,1), cond);
    
    % Get ready times to be the end of the baseline for this condition
    bl_stops = SDATA.events.event_times(cond_index,1);

    % initialize condition baseline
    condition_baseline = zeros(1,14);

    for block = blocks
        if strcmp("b1",block)
            b = 1;
        elseif strcmp("b2", block)
            b = 2;
        end

        bl_stop = bl_stops(b);
        bl_start = bl_stop-bl_dur;

        % getting baseline neural data for block baseline
        NDATA.(condition).(block).bl_start = bl_start;
        NDATA.(condition).(block).bl_stop = bl_stop;
        NDATA.(condition).(block).block_baseline_data = SDATA.data(bl_start:bl_stop,1:14);
        NDATA.(condition).(block).block_baseline = mean(SDATA.data(bl_start:bl_stop,1:14),1);

        % add 1/2 of block_baseline to condition baseline (two blocks so
        % they each contribute 1/2 to the condition-average baseline).
        condition_baseline = condition_baseline + 0.5 * mean(SDATA.data(bl_start:bl_stop,1:14),1);

        for hand = hands
            % exclude contra hand for ipsi condition and vice versa
            if (hand == condition) || (condition == 'bimanual') 
                [chan_data, move_times, window] = EMGPeakTimes(pat, condition, b, hand, window);
                NDATA.(condition).(block).(hand).peak_times = move_times;
                NDATA.(condition).(block).(hand).times = -window:window;
                NDATA.(condition).(block).(hand).data = move_times;

                %Use SegAndAvg to create the epochs of neural data. Function automatically segments over 1st dimension which is time (desired).
                [AVG, SEGS, REJ] = SegAndAvg(SDATA.data(:,1:14), move_times, -window:window-1);
                NDATA.(condition).(block).(hand).data = SEGS;
                NDATA.(condition).(block).(hand).data = permute(SEGS, [2,1,3]);
            end   
        end      
    end
    
    if condition ~= 'bimanual'
        % combine epoch data from both blocks
        e1 = NDATA.(condition).b1.(condition).data;
        e2 = NDATA.(condition).b2.(condition).data;
        NDATA.(condition).data = cat(3,e1,e2);
        
        % getting epoch dimension of data matrix (#of epochs) to make indexing easier
        NDATA.(condition).epoch = size(NDATA.(condition).data, 3);
        
    else
        % combine epoch data from both blocks for both hands
        ec1 = NDATA.(condition).b1.contra.data;
        ec2 = NDATA.(condition).b2.contra.data;
        ei1 = NDATA.(condition).b1.ipsi.data;
        ei2 = NDATA.(condition).b2.ipsi.data;
        NDATA.(condition).contra_data = cat(3,ec1,ec2);
        NDATA.(condition).ipsi_data = cat(3,ei1,ei2);
        
        %getting epoch dimension of data matrix (#of epochs) to make indexing easier
        NDATA.(condition).contra_epoch = size(NDATA.(condition).contra_data, 3);
        NDATA.(condition).ipsi_epoch = size(NDATA.(condition).ipsi_data, 3);
    end
          
    %combine base line data from both blocks
    bl1 = NDATA.(condition).b1.block_baseline_data;
    bl2 = NDATA.(condition).b2.block_baseline_data;
    NDATA.(condition).cond_spec_baseline = cat(1, bl1, bl2);
    
    %adding condition-average baseline to structure
    NDATA.(condition).condition_baseline = condition_baseline;

end

%getting time dimension of data matrix to make indexing easier
NDATA.pnts = window*2;
NDATA.times = [1:NDATA.pnts];

%Saving the struct as a .mat file in patient's SDATA_files
% cd(strcat('/Users/owen/Box/DBS_data/sdata_files/',pat));
% file_name = sprintf('NDATA_%s',pat);
% save(file_name, 'NDATA', '-v7.3');
end





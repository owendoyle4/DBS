patient_id = 'BI032';
c_emg = 15; %15 or 17
i_emg = 16; %16 or 18
% change directory to yours %
%cd(append('C:\Users\natal\Box\DBS_data\sdata_files\',patient_id))
cd(strcat('/Users/owen/Box/DBS_data/sdata_files/',patient_id))
% % % % % % % % % % % % % % % 

load(append('SDATA_',patient_id))
load('data_10hz')
load('events_labeled')
move_temp = load(strcat('move_', patient_id, '_temp.mat'));
move_temp = move_temp.(strcat("move_", patient_id, "_temp"));

% cd('/Users/owen/Box/DBS_data/sdata_files/')
% 
% load('move_t_dict_12-14')
% lpf = load(append('data_lpf'))
% all(lpf.data==data)


%%
num_blocks = length(SDATA.events.event_names);
figure;

for block = 1:num_blocks
    %index 
    start = SDATA.events.event_times(block,1); %ready
    go = SDATA.events.event_times(block,3); %go
    stop = SDATA.events.event_times(block,4); %stop
    
    name = SDATA.events.event_names(block);
    
    %choose hand chan based on block name
    if strcmp(name,'Hand_contra') || strcmp(name,'HO_contra') %contra
        hand_idx = c_emg:i_emg; %contra hand 
        title_name = (append(patient_id,' ',name));
        
    elseif strcmp(name,'Hand_ipsi') || strcmp(name,'HO_ipsi') %ipsi
        hand_idx = c_emg:i_emg; %ipsi hand 
        title_name = (append(patient_id,' ',name));
       
    elseif strcmp(name,'Hand_bilat') || strcmp(name,'HO_bilat') %bimanual
        hand_idx = c_emg:i_emg; %contra and ipsi
        title_name = (append(patient_id,' ',name));
    end
    c_color = 'b';
    i_color = 'r';
 
    
    %subplot for each block
    rows = (num_blocks/3);
    subplot(num_blocks/3, 3, block)
    
    plot(start:stop+4000, data(start:stop+4000, hand_idx(1)), 'LineWidth',1, 'Color', c_color); hold on
    plot(start:stop+4000, data(start:stop+4000, hand_idx(2)), 'LineWidth',1, 'Color', i_color); hold on
    
    xline(go, 'k', {' SDATA go'});
    xline(stop, 'k', {' SDATA stop'});
    
    if move_temp(block, 2) ~= 0
        xline(move_temp(block, 2), 'b', {' move_t go'});
        xline(move_temp(block, 3), 'b', {' move_t stop'});
    end
    
    title(title_name)
end
%legend
% add legend
lgd = legend('show');
lgd.Position(1) = 0;
lgd.Position(2) = 0.9;
lgd.String = ({'contra','ipsi'});

%%
block = 1;
start = SDATA.events.event_times(block,1); %ready
go = SDATA.events.event_times(block,3); %go
stop = SDATA.events.event_times(block,4); %stop
figure;
plot(start:stop+2000, data(start:stop+2000, hand_idx(1)), 'LineWidth',1, 'Color', c_color); hold on
plot(start:stop+2000, data(start:stop+2000, hand_idx(2)), 'LineWidth',1, 'Color', i_color); hold on
%% Correlation Time series
section_dur = 500; %ms
num_blocks = length(SDATA.events.event_names);
data_sections = struct();
figure;
sgtitle(patient_id);

disp(patient_id)
for block = 1:num_blocks
%     %only look at bimanual blocks
%     if move_temp(block, 1) == 2
    b_name = 'b'+string(block);

    start  = move_temp(block, 2);
    stop   = move_temp(block, 3);
    data_b = data(start:stop, c_emg:i_emg);
    data_b(:,1) = [zeros(500,1); data_b(1:length(data_b)-500,1)];
    

    S = size(data_b);
    section_durs = section_dur * ones(1,ceil(S(1)/section_dur));
    section_durs(end) = 1+mod(S(1)-1,section_dur);
    data_sections.(b_name).data = mat2cell(data_b,section_durs,S(2));

    num_sections = length(data_sections.(b_name).data);
    corrs = ones(1, num_sections);

    for section = 1:num_sections
        cell = data_sections.(b_name).data(section);
        cell_data = cell{1,1};
        corrs(1, section) = corr(cell_data(:,1),cell_data(:,2));
    end

    data_sections.(b_name).corrs = corrs;

    subplot(num_blocks/3, 3, block)
%     axis([0,30, -1, 1])
    plot(corrs);
    yline(nanmean(corrs));

    if block == 3 || block == 6 || block == 9
        disp(b_name)
        disp(nanstd(corrs))
    end
    ylim([ -1, 1]);
    title(b_name);

%     end
end

disp('---------')

figure;
for block = 1:num_blocks
%     %only look at bimanual blocks
%     if move_temp(block, 1) == 2
    b_name = 'b'+string(block);

    start  = move_temp(block, 2);
    stop   = move_temp(block, 3);
    data_b = data(start:stop, c_emg:i_emg);
    data_b(:,1) = [zeros(180,1); data_b(1:length(data_b)-180,1)];
    subplot(num_blocks/3, 3, block)
    plot(data_b)
end

%plot bimanual signals
figure;
plot_idx = 1;
sgtitle(patient_id);

for block = 1:num_blocks
    %only look at bimanual blocks
    if move_temp(block, 1) == 2
        start  = move_temp(block, 2);
        stop   = move_temp(block, 3);
        
        subplot(num_bi,1, plot_idx);
        plot(start:stop, data(start:stop, hand_idx(1)), 'LineWidth',1, 'Color', c_color); hold on
        plot(start:stop, data(start:stop, hand_idx(2)), 'LineWidth',1, 'Color', i_color); hold on
        
        plot_idx = plot_idx + 1;
    end
end
%% Compare FDI unimanual correlations to flex unimanual correlations

%concerning period:
% start = 122065;  %SDATA.events.event_times(2,3);
% stop = 127721; %SDATA.events.event_times(2,4);

%normal block
start = SDATA.events.event_times(5,3);
stop =  SDATA.events.event_times(5,4);
 
% start = 127500; 
% stop = SDATA.events.event_times(2,4);

block = 2;

%FDI, ipsi block 1
contra_fdi = SDATA.data(start:stop, 15);
ipsi_fdi   = SDATA.data(start:stop, 16);

%Flex, ipsi block 1
contra_flex = SDATA.data(start:stop, 17);
ipsi_flex   = SDATA.data(start:stop, 18);

disp(SDATA.events.event_names(block,1));
disp("contra:" + string(corr(contra_fdi, contra_flex)));
disp("ipsi:" + string(corr(ipsi_fdi, ipsi_flex)));

%% make move_times
%tag the new times that you want to change in the figure and update
%move_time with those new times
move_times = zeros(num_blocks, 3);

%% Save move_times (MAKE SURE ITS NOT ALL ZEROS STILL)
%save old version of events_labeled
% old_events_labeled = events_labeled;
% filename = strcat('move_', patient_id,'_temp');
% save old matrix as old_events_labeled for each patient folder
cd(append('/Users/owen/Box/DBS_data/sdata_files/', patient_id))
filename = append("move_", patient_id, "_temp");
move_BI032_temp = move_times;
save(filename, filename,'-v7.3')


%%
% save new matrix as events_labeled for each patient folder
% save('events_labeled','events_labeled','-v7.3')


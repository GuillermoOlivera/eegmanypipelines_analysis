base_path = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated\';
match = [ base_path '*.out.mat' ];
files  = ls(match);
len = size(ls(match), 1);

f1 = figure('visible','on');
f1.Position = [20 20 1250 650];
hax1 = axes;

for m=1:len
    file = files(m,:);
    load(fullfile(base_path, file));

    times = loaded_raw_data_from_eeglab.times;
    data = loaded_raw_data_from_eeglab.data;
    hold off;
    plot(hax1, times, data(1,:))
    hold on;
    plot(hax1, times, data(10,:))
    plot(hax1, times, data(20,:))
    plot(hax1, times, data(30,:))
    plot(hax1, times, data(40,:))
    plot(hax1, times, data(50,:))
    plot(hax1, times, data(60,:))
    plot(hax1, times, data(70,:))
    legend({file})

    idx_str = num2str(m);
    if numel(idx_str) == 1
        idx_str = ['0' idx_str];
    end
    filename = ['G:\_EEGManyPipelines\EMP_data\eeg_brainvision\images\all_sessions' idx_str '.png'];
    saveas(f1, filename)
end


%%
f = figure('visible','off');
% for r=1:length(all_segments_erp_summary.manmade_new_mean_electrodes)
for r=1:2
    plot(all_segments_erp_summary.manmade_new_mean_electrodes{r});
    idx_str = num2str(r);
    if numel(idx_str) == 1
        idx_str = ['0' idx_str];
    end
    legend({[ 'Electrode' idx_str ]})
    filename = ['G:\_EEGManyPipelines\EMP_data\eeg_brainvision\images\manmade_new_electrode_' idx_str '.png'];
    saveas(f, filename)
end

system('magick  -delay 200 -loop 0 images/*.png animation.gif');
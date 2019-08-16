%% play with data froim the Ziv lab
clear
clc

%% load the data
load('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\C6M4_Day3_A_am.mat')
load('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\results_1.mat')
info = readtable('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\frameLog_1.csv');
spikeTrials = spikeTrials{1}.spikeTrials;
fs = 20;
dt = 1/fs;
behave_trials = 2:6;

%%
pnl = panel();
pnl.pack(2,5);
for trial = 2:6
    pos = my_mvmt{trial};
    trial_ti = [info.begFrame(trial) info.endFrame(trial)];
    nSpikes = [];
    pnl(1,trial-1).select();
    title("trial "+trial)
    xlabel('#frame')
    ylabel('position')
    hold on
    plot(pos.position,'k');
    for cell = 1:size(spikeTrials,2)
        spikes = spikeTrials{trial,cell};
        nSpikes(cell) = size(spikes,2);
        if ~isempty(spikes)
            IX = spikes(1,:);
%             spikes(2,:)
%             spikes(3,:)
%             pause
%             plot(IX, pos.position(IX),'.','MarkerSize',10);
        end
    end
    pnl(2,trial-1).select();
    histogram(nSpikes)
    title("trial "+trial)
    xlabel('#spikes')
    ylabel('#cells')
end

%%
figure
hold on
% plot(pos.raw_x,pos.raw_y,'.')
% axis equal
plot(pos.position,'.','MarkerSize',10)
IX = find(pos.goodframes);
plot(IX,pos.position(IX),'.r','MarkerSize',10)
legend({'all frames';'good frames'})




%%









%%

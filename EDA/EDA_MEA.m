%% 23_10_25 Nate Stanko MEA Presentation
%% load/clean pre/post H data
h5disp('10dayOld_postHist.h5') %visualizes the large h5 file;
poDa = h5read('10dayOld_postHist.h5','/data_store/data0000/spikes'); %loads file
preDa = h5read('10dayOld_preHist.h5','/data_store/data0000/spikes'); %loads file


%collect relevant columns within each stucture
po_chan = poDa.channel;
po_spT = double(poDa.frameno);
pre_chan = preDa.channel;
pre_spT = double(preDa.frameno);


% normalize peak occurence values
po_spT = po_spT - min(po_spT); 
pre_spT = pre_spT - min(pre_spT);


% convert data to seconds using the 20kHz freq
po_spT = po_spT/20000;
pre_spT = pre_spT/20000;


% find total length of recording in minutes
po_totTime = max(po_spT)/60;
pre_totTime = max(pre_spT)/60;


%% Generate trace figure (with mean lines) for both pre and post conditions
figure
po_spTr = histcounts(po_spT,0:1:600);
pre_spTr = histcounts(pre_spT,0:1:600);
plot(po_spTr,'m', LineWidth=1.05) % post condition
po_mn = mean(po_spTr);
pre_mn = mean(pre_spTr);
hold on
plot(pre_spTr,'k', LineWidth=1.05) % pre condition
title('Overall Activity in the Retina before (black) and after (purple) the addition of Histamine', FontSize=20)
ylabel('Firing Rate (spikes/sec)', FontSize=15)
xlabel('Time (s)', FontSize=15)
hold on
line([0 600],[1111 1111],Linewidth=1.5,color="b",LineStyle = "--")
line([0 600],[1963 1963],Linewidth=1.5,color="b",LineStyle = '--')

%% Generate boxchart figure for both pre and post conditions
figure
subplot(1,2,1)
boxchart(pre_spTr)
ylim([500,3000])
ylabel('Firing Rate (activity/sec)', FontSize=15)
title('Activity in the Retina Pre Histamine', FontSize=20)
subplot(1,2,2)
boxchart(po_spTr)
ylim([500,3000])
ylabel('Firing Rate (activity/sec)', FontSize=15)
title('Activity in the Retina Post Histamine', FontSize=20)

%% Generate swarmchart figure for both pre and post conditions
figure
swarmchart(ones(1,600),pre_spTr,'k*')
hold on
swarmchart(ones(1,600),po_spTr,'m*')
ylabel('Firing Rate (activity/sec)', FontSize=15)
title('Overall Activity in the Retina before (black) and after (purple) the addition of Histamine', FontSize=20)


%% load/clean different MEA data for movie

h5disp("data.raw.h5") %loads the large h5 file
data = h5read("data.raw.h5",'/data_store/data0000/spikes'); %loads the relevant portion of the h5 file
channels = data.channel;

spikeTimes = double(data.frameno); %selects a portion of the relevant structure, ie at what time each peak occurs
spikeTimes = spikeTimes - min(spikeTimes); %normalizes peak occurance but not nessecarily to seconds
%data is aquired at 20kHz, want to convert to per second
spikeTimes = spikeTimes/20000; %converts data to seconds at which peaks occur
totTime = max(spikeTimes)/60; % will show the number of minutes that passed during data collection

uniqueChannels = unique(channels);
numChannels = size(uniqueChannels,1);

for i = 1:numChannels
    currentChannel = uniqueChannels(i);
    inds = channels==currentChannel;
    spikeTimeChannel = spikeTimes(inds);
    spikeRateChannel(i,:) = histcounts(spikeTimeChannel,0:1:1800);
end

%% load MEA mapping data for movie
map = h5read("data.raw.h5",'/data_store/data0000/settings/mapping');
[commonChannels,ia,ib] = intersect(map.channel,uniqueChannels);

%making sure metadata is only common channels
map.channel = map.channel(ia);
map.electrode = map.electrode (ia);
map.x = map.x (ia);
map.y = map.y (ia);

numChannels = size(commonChannels,1);
xuniq = size(unique(map.x),1);
yuniq = size(unique(map.y),1);
xDim = unique(map.x);
yDim = unique(map.y);

%% Make a movie

movieVar = nan(yuniq,xuniq,1800);

for i = 1:numChannels
    currentX = map.x(i);
    currentY = map.y(i);
    xInd = find(currentX == xDim);
    yInd = find(currentY == yDim);
    movieVar(yInd,xInd,:) = spikeRateChannel(i,:);
end

figure
for i = 1:1800
    imagesc(movieVar(:,:,i))
    M(i) = getframe;
end

movie(M,1,30);
%% end

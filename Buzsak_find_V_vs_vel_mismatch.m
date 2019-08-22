%%
clear
clc
load('C:\Tamir\work\Courses\MBL\project\figures\buzsaki\opt_UVW_pos_R=1__run2\results.mat');

%% 
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

%% plot with brush functionality
x = interp1(pos.TimeStamps, pos.OneDLocation, t ,'linear','extrap');
xy = interp1(pos.TimeStamps, pos.TwoDLocation, t ,'linear','extrap');
hlines = {};
figure

subplot(121)
hlines{1} = plot(vel,V,'.');
xlim([0 1])
xlabel('speed')
ylabel('V')
title('Latent variable vs. speed')

subplot(122)
hlines{2} = plot(xy(:,1),xy(:,2),'.');
xlabel('position X')
ylabel('position Y')
title('2D position')

% subplot(212)
% hlines{3} = plot(t,x,'.');
% xlabel('Time')
% ylabel('Linearized position')

setappdata(gcf, 'hlines', hlines);
% suptitle('Mismatch between latent variable V and speed peojected to behavior')

% set brush callback
bO = brush(gcf);
set(bO, 'ActionPostCallback', @brush_CB);

%%
subplot(121)
% vel_band = [0.0335 0.1612];
vel_band = [0 0.1612];
xline(vel_band(1));
xline(vel_band(2));
IX = find(vel>vel_band(1) & vel<vel_band(2));
m = median(V(IX));
yline(m);
% % % hax = gca;
% % % hline = findall(hax(ii_ax),'type','Line');
% % 
% % %% choose upper part
% % hlines{1}.BrushData = (vel>vel_band(1) & vel<vel_band(2) & V>m);
% % 
% % %% choose lower part
% % hline.BrushData = (vel>vel_band(1) & vel<vel_band(2) & V<m);

%% callback function for brush selection
function brush_CB(varargin)
    hline = findall(gca,'type','Line');
    tf = hline.BrushData;
    hax = findall(gcf,'type','axes');
    for ii_ax = 1:length(hax)
        hline = findall(hax(ii_ax),'type','Line');
        hline.BrushData = tf;
    end

%     % get selected data (with brush) - from V vs. vel
%     hfigdata = getappdata(gcf);
%     hlines = hfigdata.hlines{1};
%     tf = hlines.BrushData;
%     % set selected data (with brush) - for spikes voltage plots
%     hlines = hfigdata.hlines{2};
%     hlines.BrushData = tf;
%     hlines = hfigdata.hlines{3};
%     hlines.BrushData = tf;
end
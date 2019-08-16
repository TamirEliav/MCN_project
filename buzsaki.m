%%
clear
clc

%% load data
load('C:\Tamir\work\Courses\MBL\project\Datasets\buzsaki\i01_maze15_MS.001\i01_maze15_MS.001_BehavElectrData.mat');

%% create time vector
dt = 1 ./ Par.SamplingFrequency;
t = [1:Par.nTimebins] .* dt;

%% create activity matrix
% X = zeros(length(Clu.totClu), length(t));
% X(sub2ind(size(X),Spike.totclu, Spike.res)) = 1;
% sigma_sec = 1;
% sigma = round(sigma_sec/dt);
% hsize = sigma*5+1;
% ker = fspecial('gaussian',[1 hsize],sigma);
% X2 = imfilter(X,ker);

%% discretize neural data in bins
bin_size = 0.05;
bin_edges = 0:bin_size:t(end);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))./2;
Y = zeros(length(Clu.totClu), length(bin_centers));
for clu=1:length(Clu.totClu)
    IX = find(Spike.totclu==clu);
    Y(clu,:) = histcounts(t(Spike.res(IX)),bin_edges) ./ bin_size ;
end
sigma_sec = 0.1;
sigma = round(sigma_sec/bin_size);
hsize = sigma*5+1;
ker = fspecial('gaussian',[1 hsize],sigma);
Y2 = imfilter(Y,ker);

%% discretize behavioral data in bins
X = zeros(2,length(bin_centers)); % 2 covariates
X(1,:) = interp1(t, Track.X', bin_centers,'linear','extrap');
X(2,:) = interp1(t, Track.Y, bin_centers,'linear','extrap');

%% plot behavior
figure
plot(Track.X, Track.Y, '.');

%% PCA
IX = find(~Clu.isIntern);
[coeff,score,latent,tsquared,explained,mu] = pca(Y(IX,:)');
% [coeff,score,latent,tsquared,explained,mu] = pca(Y');

ncolor = 100;
cmap = jet(ncolor);

figure
c = interp1(linspace(min(X(1,:)),max(X(1,:)),ncolor), cmap, X(1,:));
scatter3(score(:,1),score(:,2),score(:,3),5,c);
h=colorbar;
h.Label.String = 'X pos';
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

figure
c = interp1(linspace(min(X(2,:)),max(X(2,:)),ncolor), cmap, X(2,:));
scatter3(score(:,1),score(:,2),score(:,3),5,c)
h=colorbar;
h.Label.String = 'Y pos';
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')


%% t-SNE
tic
[Y, Lose] = tsne(activity', 'NumDimensions',2, 'NumPCAComponents',0, 'Perplexity',500);
toc


%%


%%
figure
hold on
plot(t(Spike.res)./60, Spike.totclu,'.')
interneuron_IX = find(Clu.isIntern);
IX = boolean(Clu.isIntern(Spike.totclu));
plot(t(Spike.res(IX))./60, Spike.totclu(IX),'.')

%%
figure
hold on
for clu = 1:length(Clu.totClu)
    cla
    plot(Track.X, Track.Y, '.', 'Color',0.5*[1 1 1]);
    IX = find(Spike.totclu==clu);
    plot(Spike.X(IX),Spike.Y(IX),'.r')
    title(num2str(clu))
    pause
end


%%













%%

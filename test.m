%%
clear
clc

%% create real tuning (pos)
n1 = 10;
n2 = 10;
k1 = 20;
k2 = 20;
M = zeros(n1,n2,k1,k2);
Sigma=10;
n_to_k_factor=k1/n1;
x = 1:k1;
y = 1:k2;
for ii = 1:n1    
    for jj = 1:n2
        x0=(ii-0.5)*n_to_k_factor;
        y0=(jj-0.5)*n_to_k_factor;
        M(ii,jj,:,:) = exp(-((x-x0).^2+(y'-y0).^2)/(2*Sigma^2));
    end
end
M2 = reshape(M,[n1*n2 k1*k2]);
M2 = M2.*20; % set max FR

%% add tuning to other params (HD)
HD_nbins = 100;
HD_FR = zeros(size(M2,1),HD_nbins);
HD_bins = linspace(-pi,pi,HD_nbins);
preferred_angles = linspace(-pi,pi,size(M2,1));
for cell=1:size(HD_FR,1)
    HD_FR(cell,:) = circ_vmpdf(HD_bins, preferred_angles(cell));
end

%% plot HD tuning
figure
hold on
imagesc(HD_FR);
axis tight
axis equal


%% plot tuning maps for all cells
figure
pnl = panel();
pnl.pack(n1,n2);
pnl.margin=10
pnl.de.margin=4
for ii = 1:n1
    for jj = 1:n2
        pnl(ii,jj).select();
        imagesc(squeeze(M(ii,jj,:,:)))
        set(gca,'xtick',[],'ytick',[]);
    end
end

%% create behavior
T = 9000;
rng(0);
sdf = 5.*randn(2,T);
sdf(:,1) = [k1 k2]./2;
xy = cumsum(sdf,2);
xlimits = [0 k1+1];
ylimits = [0 k2+1];
for ii_t = 1:T
    if xy(1,ii_t) < xlimits(1)
        xy(1,ii_t:end) = 2*xlimits(1) - xy(1,ii_t:end);
    elseif xy(1,ii_t) > xlimits(2)
        xy(1,ii_t:end) = 2*xlimits(2) - xy(1,ii_t:end);
    end
    if xy(2,ii_t) < ylimits(1)
        xy(2,ii_t:end) = 2*ylimits(1) - xy(2,ii_t:end);
    elseif xy(2,ii_t) > ylimits(2)
        xy(2,ii_t:end) = 2*ylimits(2) - xy(2,ii_t:end);
    end
end

% smooth the behavior
xy = [smooth(xy(1,:)) smooth(xy(2,:))]';

% remove a band of behavior
% IX = find(xy(1,:) > 15 & xy(1,:) < 30);
% IX = find(sqrt(sum((xy-50).^2)) < 20);
% xy(:,IX(1:1.5:end))=[];

% calculate other behavior variables than position
vel= [0 sqrt(sum(diff(xy')'.^2))];
HD = [0 angle(diff(xy(1,:))+i.*diff(xy(2,:)))];
HD_IX = discretize(HD,HD_bins);

xy_binned = round(xy);
xy_binned(xy_binned<=0) = 1;
xy_binned(1,xy_binned(1,:)>=k1+1) = k1;
xy_binned(2,xy_binned(2,:)>=k2+1) = k2;

%% plot behavior
figure
subplot(1,2,1)
plot(xy(1,:), xy(2,:),'.')
axis equal
axis tight
subplot(1,2,2)
plot(xy_binned(1,:), xy_binned(2,:),'.')
axis equal
axis tight

%% plot nice behavior figure
figure
hold on
IX = 1:7000;
% plot(smooth(xy(1,IX)), smooth(xy(2,IX)),'-k')
plot(xy(1,IX), xy(2,IX),'-k')
axis equal
axis tight
xlabel('X position')
ylabel('Y position')

%% plot behavior with HD colors
HD_complex = exp(i.*HD);
pos_HD_bins = 1:5:100;
bins = interp1(pos_HD_bins, 1:length(pos_HD_bins), xy_binned,'nearest','extrap');
HD_mean_per_pos = angle(accumarray(bins', HD_complex', [], @mean));
HD_std_per_pos = angle(accumarray(bins', HD_complex', [], @std));
% HD_mean_per_pos = angle(accumarray(xy_binned', HD_complex', [], @mean));
% HD_std_per_pos = angle(accumarray(xy_binned', HD_complex', [], @std));
% HD_mean_per_pos = angle(accumarray(bins', HD', [], @circ_mean));
% HD_std_per_pos = angle(accumarray(bins', HD', [], @circ_std));

cmap = jet(HD_nbins);
figure
subplot(121)
hold on
scatter(xy(1,:),xy(2,:),5,cmap(HD_IX,:))
colormap hsv
axis equal
axis tight
subplot(122)
hold on
imagesc(HD_mean_per_pos)
colormap hsv
axis equal
axis tight
% subplot(223)
% hold on
% imagesc(HD_std_per_pos)
% axis equal
% axis tight
% colorbar


%% create activity matrix from position FR map
pos_IX = sub2ind([k1 k2], xy_binned(1,:), xy_binned(2,:));
FR = M2(:,pos_IX);

%% add low-dim activity
rng(0);
N = size(FR,1);
R = 2;
T = size(FR,2);
U = randn(N,R);
V = randn(R,T).*5;
h = fspecial('gaussian',[1 1000],250)*10;
V2 = imfilter(V,h);
FR = FR + U*V2;
FR(FR<0) = 0; % make sure no negative FR
figure
hold on
% plot(V')
plot(V2')

%% add activity dependent on other variables (HD/speed)
% FR = FR + HD_FR(:,HD_IX);

%% generate poisson spikes from FR
spikes = poissrnd(FR);

%% plot activity of single cell
figure
plot(FR(2,:))

%% validatre we get nice tuning curves
% time_spent = histcounts(
cell_num = 46;
A=accumarray(xy_binned(1,:)', FR(cell_num,:)',[k1 1],@mean,nan);
plot(A)


%% t-SNE
tic
[Y, Lose] = tsne(FR', 'NumDimensions',2, 'NumPCAComponents',0, 'Perplexity',500);
toc

%% plot tSNE results with pos colors
cmap = jet(k1);
figure
hold on
scatter(Y(:,1),Y(:,2),5,cmap(xy_binned(2,:),:))

%% plot tSNE results with HD colors
cmap = jet(HD_nbins);
figure
hold on
scatter(Y(:,1),Y(:,2),5,cmap(HD_IX,:))

%% t-SNE
tic
[Y2, Lose] = tsne(FR', 'NumDimensions',3, 'NumPCAComponents',0, 'Perplexity',300);
toc

%%
cmap = jet(k1);
figure
hold on
% scatter3(Y2(:,1),Y2(:,2),Y2(:,3),5,cmap(xy_binned(2,:),:))
scatter3(Y2(:,1),Y2(:,2),Y2(:,3),5)

%% PCA
cmap = jet(k1);
[coeff,score,latent,tsquared,explained,mu] = pca(FR');
figure
scatter3(score(:,1),score(:,2),score(:,3),5,cmap(xy_binned(1,:),:));
figure
scatter3(score(:,1),score(:,2),score(:,3),5,cmap(xy_binned(2,:),:));

%% PCA
IX = find(sqrt(sum((xy-50).^2)) < 50);
% [coeff,score,latent,tsquared,explained,mu] = pca(activity(:,IX)');
figure
scatter3(score(IX,1),score(IX,2),score(IX,3),5,cmap(xy_binned(1,IX),:));
figure
scatter3(score(IX,1),score(IX,2),score(IX,3),5,cmap(xy_binned(2,IX),:))

%% plot PCA coeff
figure
for ii = 1:n1
    for jj = 1:n2
        ii_cell = sub2ind([n1 n2],ii,jj); 
        subplot(n1,n2,ii_cell)
        imagesc(reshape((coeff(:,ii_cell)),[n1 n2]));
    end
end

%% convert position from 2D to sparse representation per bin
pos = zeros(k1*k2,T);
IX = sub2ind(size(pos), pos_IX, 1:length(pos_IX));
pos(IX) = 1;

%% try our model
Y = spikes;
% X = xy(1,:);
X = pos;
N = size(spikes,1);
T = size(spikes,2);
P = size(X,1);
R = 2;
U = randn(N,R);
V = randn(R,T);
W = randn(N,P);
nIter = 10;
tic
fval=[];
options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
U(:)=0;
V(:)=0;
for iter=1:nIter
    disp(iter)
    [W fval(1,iter)] = fminunc(@(W)optW(Y,U,V,W,X),W,options);
%     [U fval(2,iter)] = fminunc(@(U)optU(Y,U,V,W,X),U,options);
%     [V fval(3,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options);
end
toc


%% try simple PCA model
N = size(FR,1);
T = size(FR,2);
r = 2;
P = 1;
% X = xy(1,:);
Y = FR;
Y = Y - mean(Y,2);
U = randn(N,r);
V = randn(r,T);
% W = randn(n,k);
nIter = 10;
tic
loss=[];
for iter=1:nIter
    disp(iter)
    U = Y/V;
    V = U\Y;
%     [U fval(1,iter)] = fminunc(@(U)(norm(Y-U*V,'fro')), U);
%     [V fval(2,iter)] = fminunc(@(V)(norm(Y-U*V,'fro')), V);
    loss(iter) = norm(Y-U*V, 'fro');
%     figure
%     plot(V(1,:),V(2,:),'.')
end
toc
figure
plot(V(1,:),V(2,:),'.')

%% try only to explain X
Y = FR;
% Y = Y - mean(Y,2);
N = size(Y,1);
T = size(Y,2);
X = xy(1:2,:);
P = size(X,1);
% W = randn(n,k);
% r = 2;
% U = randn(n,r);
% V = randn(r,T);
nIter = 10;
tic
loss=[];
for iter=1:nIter
    disp(iter)
    U = Y/V;
    V = U\Y;
%     [U fval(1,iter)] = fminunc(@(U)(norm(Y-U*V,'fro')), U);
%     [V fval(2,iter)] = fminunc(@(V)(norm(Y-U*V,'fro')), V);
    loss(iter) = norm(Y-U*V, 'fro');
    figure
    plot(V(1,:),V(2,:),'.')
end
toc

%%
figure
for cell=1:100
%     imagesc(exp(reshape(W(cell,:)',[k1 k2])))
    imagesc(reshape(W(cell,:)',[k1 k2]));
    colorbar
    pause
end

%% optimization functions
function [f,g] = optW(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = sum(Yhat.*Y-exp(Yhat),'all');
    g = (Y-exp(Yhat))*X';
end
function [f,g] = optU(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = sum(Yhat.*Y-exp(Yhat),'all');
    g = (Y-exp(Yhat))*V';
end
function [f,g] = optV(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = sum(Yhat.*Y-exp(Yhat),'all');
    g = U'*(Y-exp(Yhat));
end



%%




%%
%%
clear
clc

%% TODO:
% add orthogonality constrain on UV
% add basline term b
% add regularization for smoothness in W
% add regularization for smoothness in V
% change V to be sparse and change the regularization on V to be sparse

%% create real tuning (pos)
n1 = 10;
n2 = 10;
k1 = 5;
k2 = 5;
M = zeros(n1,n2,k1,k2);
Sigma=1.5;
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
max_FR = 20;
M2 = M2.*max_FR; % set max FR

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
T = 500;
rng(0);
sdf = 2.*randn(2,T);
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
% IX = 1:5000;
IX = 1:T;
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
FR = log(M2(:,pos_IX));

%% add low-dim activity
rng(0);
N = size(FR,1);
% R = 2;
R = 1;
T = size(FR,2);
uv_pattern = 'rand_smooth';
% uv_pattern = 'linear_monotonic_EI';
switch uv_pattern
    case 'rand_smooth'
        U2 = rand(N,R);
        V2 = randn(R,T).*1;
        h = fspecial('gaussian',[1 1000],50);
        V2 = imfilter(V2,h,'symmetric');
%         U2 = U2 - mean(U2);
%         V2 = V2 - mean(V2);
    case 'linear_monotonic_EI'
        U2 = zeros(N,R);
        U2(1:(N/2),1) = 1;
        U2((N/2+1):end,2) = 1;
        V2 = [linspace(0,1,T);linspace(0,-1,T)];
        V2 = V2.*0.5;
end
figure
hold on
plot(V2')
%%
% FR(:)=0;
FR = FR + U2*V2;

%% add activity dependent on other variables (HD/speed)
% FR = FR + HD_FR(:,HD_IX);

%% generate poisson spikes from FR
rng(0);
spikes = poissrnd(exp(FR));

%% plot activity of single cell
figure
plot(spikes(1,:))

%% validate we get nice tuning curves
% time_spent = histcounts(
cell_num = 20;
A=accumarray(xy_binned(1,:)', spikes(cell_num,:)',[k1 1],@mean,nan);
figure
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
[coeff,score,latent,tsquared,explained,mu] = pca(spikes');
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
pos = [pos ; ones(1,T)];

%% try our model
Y = spikes;
X = pos;
N = size(spikes,1);
T = size(spikes,2);
P = size(X,1);
R = 1;
U = rand(N,R);
V = randn(R,T);
W = rand(N,P).*sqrt(log(20));
W(:,end) = 0; % zero base line
% U = U - mean(U(:));
% V = V - mean(V(:));
nIter = 10;
tic
options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
% options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% options.CheckGradients = true;
% U(:)=0;
% V(:)=0;
% W(:)=0;
% V = V2;
% U = U2;
tic
opttime=nan(4,nIter);
fval=nan(4,nIter);
for iter=1:nIter
    disp(iter);
%     A = sparse(-eye(numel(W)));
%     b = sparse(10*ones(numel(W),1));
%     lb = -10*ones(numel(W),1);
%     [W fval(1,iter)] = fmincon(@(W)optW(Y,U,V,W,X),W,A,b,[],[],[],[],[],options);
%     [W fval(1,iter)] = fmincon(@(W)optW(Y,U,V,W,X),W,[],[],[],[],lb,[],[],options);
%     [W fval(1,iter)] = fminunc(@(W)optW(Y,U,V,W,X),W,options); opttime(1,iter)=toc;
%     [U fval(2,iter)] = fminunc(@(U)optU(Y,U,V,W,X),U,options); opttime(2,iter)=toc;
    
    UW = UW_join(U,W);
    [UW fval(4,iter)] = fminunc(@(UW)optUW(Y,UW,V,X),UW,options); opttime(4,iter)=toc;
    [U,W] = UW_split(UW,R,P);
    
    [V fval(3,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options); opttime(3,iter)=toc;
end
toc
b = W(:,end);
W(:,end) = [];

%% plot optimization progress
figure
subplot(121)
plot(opttime(:)./60,'o-');
xlabel('steps')
ylabel('minutes')
subplot(122)
plot(fval(:),'o-');
xlabel('steps')
ylabel('loss')

%%
figure
hold on
plot(V2./max(V2),'k','linewidth',2)
plot(V./max(abs(V)),'r')
plot(-V./max(abs(V)),'b')


%%
[U3,S3,V3]=svd(U*V);
diag(S3)
figure
subplot(211)
% plot(-V3(:,1))
plot(V')
subplot(212)
plot(V2')

%% check UV vs. UV2
figure
subplot(221)
plot(V2,V,'.')
axis equal
xlabel('real')
ylabel('estimated')
title('V')
subplot(222)
plot(U2,U,'.')
% axis equal
xlabel('real')
ylabel('estimated')
title('U')
subplot(212)
hold on
plot(zscore(V))
plot(-zscore(V))
plot(zscore(V2),'k','linewidth',2)
title('V vs V2 (zscore)')
legend({'estimated','-estimated','real'})

%% plot W as image per cell
figure
for cell=1:size(W,1)
    imagesc(exp(reshape(W(cell,:)',[k1 k2])))
%     imagesc(reshape(W(cell,:)',[k1 k2]));
    title("cell"+cell)
    colorbar
%     pause
    pause(0.05)
end

%% plot baseline term
figure
hold on
imagesc(1:k1,1:k2,reshape(b,[n1 n2]))
plot(xy(1,:),xy(2,:),'.k')
colorbar
xlabel('neuron_X')
ylabel('neuron_Y')
legend('behavior')
title('baseline term per neuron by position')

%%
figure
subplot(411)
imagesc(U2*V2)
colorbar
title('real')

subplot(412)
imagesc(U*V)
title('estimated')
colorbar

subplot(413)
imagesc(zscore(U*V,1,2))
title('estimated (zscored)')
colorbar

subplot(414)
imagesc(exp(U2*V2))
title('spiking activity')
colorbar

%%
figure
histogram(max(exp(W)))
h=xline(max_FR);
h.Color = 'r';
h.LineWidth = 2;
xlabel('max FR');
ylabel('counts');
legend({'estimated';'real'})
h=gca;
% h.XScale = 'log';
title('max W vs. max real pos FR')

%% optimization functions
function [f,g] = optW(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -(Y-exp(Yhat))*X';
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(W,'fro')^2;
    g = g + 2*lambda*W;
end
function [f,g] = optU(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -(Y-exp(Yhat))*V';
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    g = g + 2*lambda*U;
end
function [f,g] = optV(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -U'*(Y-exp(Yhat));
    % add smoothness regularization
    Vpad = padarray(V,[0 1],'replicate','both');
    lambda = 1e3;
    f = f + lambda*sum(diff(Vpad,1,2).^2,'all');
    g = g + 2*lambda.*( 2*V -Vpad(:,1:end-2) -Vpad(:,3:end) );
    % add L2 regularization for V
    lambda = 1;
    f = f + lambda*norm(V,'fro')^2;
    g = g + 2*lambda*V;
end
function [f,g] = optUW(Y,UW,V,X)
    R = size(V,1);
    P = size(X,1);
    [U,W] = UW_split(UW,R,P);
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    gW = -(Y-exp(Yhat))*X';
    gU = -(Y-exp(Yhat))*V';
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gU = gU + 2*lambda*U;
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gW = gW + 2*lambda*W;
    g = UW_join(gU,gW);
end
function UW = UW_join(U,W)
    UW = [U W];
end
function [U,W] = UW_split(UW,R,P)
    U = UW(:,1:R);
    W = UW(:,(R+1):end);
end

%%




%%
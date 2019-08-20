%% play with data froim the Ziv lab
clear
clc

%% TODO
% try only WX, and looks for the residuals
% opt UWV together

%% load the data
load('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\C6M4_Day3_A_am.mat')
load('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\results_1.mat')
info = readtable('C:\Tamir\work\Courses\MBL\project\Datasets\Ziv\frameLog_1.csv');
spikeTrials = spikeTrials{1}.spikeTrials;
fs = 20;
dt = 1/fs;
behave_trials = 2:6;
N = length(cellRegistered);

%% explore data structure
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
%         cla
%         plot(pos.position,'k');
        spikes = spikeTrials{trial,cell};
        nSpikes(cell) = size(spikes,2);
%         if ~isempty(spikes)
%             IX = spikes(1,:);
% %             plot(IX, pos.position(IX),'.','MarkerSize',10);
% %             pause
%         end
    end
    pnl(2,trial-1).select();
    histogram(nSpikes)
    title("trial "+trial)
    xlabel('#spikes')
    ylabel('#cells')
end


%% start preparing data 
trial = 2;
pos = my_mvmt{trial};
T = length(pos.position);

%% create smooth velocty profile
smooth_p = 1e-2;
t = 1:T;
pp = csaps(1:T,pos.position,smooth_p);
pp1 = fnder(pp,1);
pos_csaps = fnval(pp,t);
vel_csaps = fnval(pp1,t).*fs;
dir0 = zeros(1,T);
dir1 = zeros(1,T);
dir2 = zeros(1,T);
vel_thr = 0.075;
dir0(abs(vel_csaps)< vel_thr) = 1;
dir1(    vel_csaps > vel_thr) = 1;
dir2(    vel_csaps <-vel_thr) = 1;
figure
subplot(121)
hold on
plot(t,pos_csaps)
plot(t,vel_csaps);
% fnplt(fnder(pp,0))
% fnplt(fnder(pp,1))
subplot(122)
hold on
plot(pos_csaps(logical(dir0)), vel_csaps(logical(dir0)),'.k');
plot(pos_csaps(logical(dir1)), vel_csaps(logical(dir1)),'.b');
plot(pos_csaps(logical(dir2)), vel_csaps(logical(dir2)),'.r');
legend({'dir0';'dir1';'dir2'})

%% prepare data!
nBinPos = 10; % number of spatial positions
X = zeros(nBinPos,T);
[~,EDGES,BIN] = histcounts(pos.position,nBinPos);
IX = sub2ind(size(X), BIN', 1:T);
X(IX) = 1;
X = [X;dir0;dir1;dir2]; % add directions (0 is no movement)
X = [X; ones(1,T)]; % add baseline term

P = size(X,1);

Y = zeros(N,T);
for cell = 1:N
    spikes = spikeTrials{trial,cell};
    if ~isempty(spikes)
        IX = spikes(1,:);
        Y(cell,IX) = 1;
    end
end

%% plot X/Y
figure
subplot(211)
imagesc(Y)
title('Y')
subplot(212)
imagesc(X)
title('X')

%% run the model !!!
R = 1;
U = rand(N,R);
V = randn(R,T);
W = rand(N,P);
W(:,end) = 0; % zero base line
% U = U - mean(U(:));
% V = V - mean(V(:));
nIter = 100;
tic
options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
% options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% options.CheckGradients = true;
% U(:)=0;
% V(:)=0;
% W(:)=0;
% V = V2;
% U = U2;
opttime2=nan(2,nIter);
fval2=nan(2,nIter);
for iter=1:nIter
    disp(iter);
    UW = UW_join(U,W);
    [UW fval2(1,iter)] = fminunc(@(UW)optUW(Y,UW,V,X),UW,options); opttime2(1,iter)=toc;
    [U,W] = UW_split(UW,R,P);
    [V fval2(2,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options); opttime2(2,iter)=toc;
end
tic
opttime=nan(2,nIter);
fval=nan(2,nIter);
UVW = UVW_join(U,V,W);
for iter=1:nIter
    disp(iter);
%     UW = UW_join(U,W);
%     [UW fval(1,iter)] = fminunc(@(UW)optUW(Y,UW,V,X),UW,options); opttime(1,iter)=toc;
%     [U,W] = UW_split(UW,R,P);
%     [V fval(2,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options); opttime(2,iter)=toc;

    [UVW fval(1,iter)] = fminunc(@(UVW)optUVW(Y,UVW,X,R,P,N,T),UVW,options); opttime(1,iter)=toc;
end
[U,V,W] = UVW_split(UVW,R,P,N,T);
toc

%% seperate b from W
% b = W(:,end);
% W(:,end) = [];

%% plot optimization progress
figure
subplot(121)
plot(opttime'./60,'o-');
xlabel('steps')
ylabel('minutes')
ylimits = get(gca,'ylim');
ylimits(1) = 0;
set(gca,'ylim', ylimits);
subplot(122)
plot(fval','o-');
xlabel('steps')
ylabel('loss')
% legend({'UW';'V'})

%% plot V vs. vel
figure
subplot(511)
hold on
plot(V)
plot(vel_csaps)
legend({'V';'velocity'})
axis tight
title('V vs. velocity')
xlabel('time')

subplot(512)
imagesc(U*V)
title('UV')
xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(513)
imagesc(W*X)
title('WX')
xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(514)
imagesc(fs*exp(U*V+W*X))
title('estimate FR')
xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(515)
imagesc(Y)
title('Y')
xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

%% plot W
figure
imagesc(W)
colorbar


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
    % add smoothness regularization on V
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
function [f,g] = optUVW(Y,UVW,X,R,P,N,T)
    [U,V,W] = UVW_split(UVW,R,P,N,T);
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    gW = -(Y-exp(Yhat))*X';
    gU = -(Y-exp(Yhat))*V';
    gV = -U'*(Y-exp(Yhat));
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gU = gU + 2*lambda*U;
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gW = gW + 2*lambda*W;
    % add smoothness regularization on V
    Vpad = padarray(V,[0 1],'replicate','both');
    lambda = 1e3;
    f = f + lambda*sum(diff(Vpad,1,2).^2,'all');
    gV = gV + 2*lambda.*( 2*V -Vpad(:,1:end-2) -Vpad(:,3:end) );
    % add L2 regularization for V
    lambda = 1;
    f = f + lambda*norm(V,'fro')^2;
    gV = gV + 2*lambda*V;
    
    g = UVW_join(gU,gV,gW);
end
function UVW = UVW_join(U,V,W)
    UVW = [U(:);V(:);W(:)];
end
function [U,V,W] = UVW_split(UVW,R,P,N,T)
    U = UVW(              1:(N*R)  );
    V = UVW( (N*R)     + (1:(R*T)) );
    W = UVW( (N*R+R*T) + (1:(N*P)) );
    U = reshape(U,[N R]);
    V = reshape(V,[R T]);
    W = reshape(W,[N P]);
end




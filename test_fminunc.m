%% play with fminunc
clear
clc

%% create data
% rng(0);
n = 3;
m = 1000;
r = 1;
% r = min(n,m);
U0 = randn(n,r);
V0 = randn(r,m);
W0 = [U0(:);V0(:)];
Y = randn(n,m);
switch n
    case 2
        Y = Y.*[4 1]';
        Y = [1 1;-1 1]'*Y;
    case 3
        Y = Y.*[10 2 .1]';
        Y = [1 1 0;-1 1 0; 0 0 1]'*Y;
%         Y = [1 1 1;0 1 0; 0 0 1]'*Y;
end
[coeff,score,latent,tsquared,explained,mu] = pca(Y');

switch n
%     case 3
%         figure
%         hold on
%         plot3(Y(1,:),Y(2,:),Y(3,:),'.')
    case {2,3}
        figure
        hold on
        plot(Y(1,:),Y(2,:),'.')
        refline(coeff(1,1)/coeff(1,2),0)
        refline(coeff(2,1)/coeff(2,2),0)
        axis equal
end

%%
% tic
% [W,fval,exitflag,output,grad,hessian] = fminunc(@(W) pca_opt(W,n,m,r,Y), W0);
% toc
tic
W = fminunc(@(W) pca_opt(W,n,m,r,Y), W0);
toc
[U,V] = reshapeW(W,n,m,r);

%%
nIter = 10;
U = U0;
V = V0;
tic
fval=[];
for ii_iter=1:nIter
%     [U fval(1,ii_iter), exitflag,output,grad,hessian] = fminunc(@(U)(norm(Y-U*V,'fro')), U);
    [U fval(1,ii_iter)] = fminunc(@(U)(norm(Y-U*V,'fro')), U);
    [V fval(2,ii_iter)] = fminunc(@(V)(norm(Y-U*V,'fro')), V);    
end
toc


%%
nIter = 10;
U = U0;
V = V0;
tic
fval=[];
loss=[];
options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
for iter=1:nIter
%     U
%     U = Y*V'*inv(V*V');
%     V = inv(U'*U)*U'*Y;
    [U fval(1,iter)] = fminunc(@(U)pcalossU(Y,U,V),U,options);
    [V fval(2,iter)] = fminunc(@(V)pcalossV(Y,U,V),V,options);
    loss(iter) = norm(Y-U*V, 'fro');
end
toc

%%
[f,g] = pcalossU(Y,U,V);
% [f,g] = pcalossV(Y,U,V);

%%
function [f,g] = pcalossU(Y,U,V)    
    M = Y-U*V;
    f = norm(M,'fro')^2;
    g = -2*Y*V'+2*U*V*V';
end

%%
function [f,g] = pcalossV(Y,U,V)    
    M = Y-U*V;
    f = norm(M,'fro')^2;
    g = -2*U'*Y + 2*U'*U*V;
end

%%
% pca_opt(W0,n,m,r,Y);
function cost = pca_opt(W,n,m,r,Y)
    [U,V] = reshapeW(W,n,m,r);
    cost = norm(U*V-Y,'fro');
end

function [U,V] = reshapeW(W,n,m,r)
    Usize = n*r;
    Vsize = m*r;
    U = W(1:Usize);
    V = W((Usize+1):end);
    U = reshape(U,n,r);
    V = reshape(V,r,m);
end




%%

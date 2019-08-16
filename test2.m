%%
clear
clc

%% using simple
% Y = [1:10] + [1:3]';
Y = [1:10] .* [1:3]';
X = 1:10;
W = Y/X;

%% using fminunc
W = fminunc(@(W)(norm(Y-W*X,'fro')), W);

%% using fminunc with gradient provided
W = fminunc(@(W)(norm(Y-W*X,'fro')), W);

%%
figure
subplot(221)
imagesc(Y)
subplot(222)
imagesc(W*X)
subplot(223)
hold on
plot(Y)
plot(W*X,'--')




%% play with GLM
mdl=fitglm(Y(1,:)', X','Distribution','poisson');
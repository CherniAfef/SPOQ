%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPOQ: Smooth lp-Over-lq ratio
%%% This code implements the SPOQ regularization function presented in
%%% "SPOQ $\ell_p$-Over-$\ell_q$ Regularization for Sparse Signal: Recovery 
%%% applied to Mass Spectrometry"
%%% Afef Cherni, IEEE member, 
%%% Emilie Chouzenoux, IEEE Member,
%%% Laurent Duval, IEEE Member,
%%% Jean-Christophe Pesquet, IEEE Fellow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Afef Cherni 04/03/2020


clc
clear all
close all
addpath(genpath('Data'));
addpath(genpath('Tools'));

%% Initialization
% load xtrue (sparse signal)
xtrue = load('x'); 
% load K (measurement matrix)
K = load('K');  
% build y = K*x
y = K*xtrue;
% add the gaussian noise with a standard deviation sigma
sigma = 0.1*max(y)/100;
noise = load('noise');
y = y + sigma*noise;
% choose SPOQ parameters
N = length(xtrue);
xi = 1.1*sqrt(N)*sigma;
eta = 2E-6;
alpha = 7E-7;
beta = 3E-2;
p = 0.25;
q = 2;
nbiter=2000;

%% SPOQ test
disp('_____________________________________________________')
disp(['You have choosen to test Trust-Region algorithm with l', num2str(p), '/l', num2str(q)]);
disp('_____________________________________________________')
tic;
[xrec,fcost,Bwhile,time,mysnr]=FB_PPXALpLq(K,y,p,q,2,alpha,beta,eta,xi,nbiter,xtrue);
tf= toc;
text=['Reconstruction en ', num2str(length(time)), ' iterations'];
disp(text);
text=['SNR = ', num2str(-10*log10(sum((xtrue-xrec).^2)/sum(xtrue.^2)))];
disp(text);
text=['Temps de reconstruction = ', num2str(sum(time)), 's.'];
disp(text);
disp('_____________________________________________________')

%% Results
figure();
plot(xtrue, 'r-o'); hold on; plot(xrec, 'b--'); hold off;
legend("Original signal", "Estimated signal");
title("Reconstruction results")
figure();
plot([0;cumsum(time(:))],mysnr,'-k','linewidth',2);
title("Convergence of the algorithm to recover dataset A");
xlabel("Time in second");
ylabel("SNR (dB)");
legend('TR-VMFB')
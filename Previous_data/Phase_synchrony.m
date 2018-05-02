%% load your data
clear; clc
load('FC_LoopGainModulatingWeightRaw0.4Weight1.mat')
close all
%% plot all channels
plot(data_struct)
data=data_struct;
t=0:0.5:750;
t(end)=[]
%% preprocessing- detrending
figure
Mean=repmat(mean(data),[1500 1]);
data=data-Mean;
plot(t,data)
%% simulated data
% t=0:0.001:2;
% y=sin(2*pi*10*t);
% data=repmat(y,[20,1,1])+0.1*randn(20,length(y));
%data=0.1*randn(20,length(y));
%% power-spectral density
srate=2000;
[pxx,f] = pwelch(data,500,300,500,srate);
plot(f,pxx)
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
%% band-pass filter and compute Hilbert transform
% D=data';
% srate=2000;
% locutoff=6;
% hicutoff=10;
% [smoothdata] = eegfilt(D,srate,locutoff,hicutoff);
%% compute the phase synchrony
Complex=hilbert(data)';
Complex=Complex./abs(Complex);
phase_synchrony=abs(mean(Complex));
phase_synchrony_mean=mean(phase_synchrony);


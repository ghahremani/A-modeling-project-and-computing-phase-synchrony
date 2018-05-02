%% load your data
clear; clc; close all
count_two=0;phase_synchrony_mean=[]; phase_synchrony_sd=[];
for gain={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.8', '1','1.5','2.0'}
    count_two=count_two+1;
    count_one=0;
    for G={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5','0.6','0.8','1.0'}
        count_one=count_one+1;
        load(['FC_LoopOPRaw' char(G) 'Weight' char(gain) '.mat'])
        data=data_struct;
        %% plot all channels
%         plot(data_struct)
%         
%         t=0:0.5:750;
%         t(end)=[];
        %% preprocessing- detrending
        Mean=repmat(mean(data_struct),[1500 1]);
        data=data-Mean;
        %plot(t,data)
        %% simulated data
        % t=0:0.001:2;
        % y=sin(2*pi*10*t);
        % data=repmat(y,[20,1,1])+0.1*randn(20,length(y));
        %data=0.1*randn(20,length(y));
        %% power-spectral density
%         srate=2000;
%         [pxx,f] = pwelch(data,500,300,500,srate);
%         plot(f,pxx)
%         xlabel('Frequency (Hz)')
%         ylabel('PSD (dB/Hz)')
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
        phase_synchrony_mean(count_one,count_two)=mean(phase_synchrony);
        phase_synchrony_sd(count_one,count_two)=sqrt(var(phase_synchrony));
    end
end

weight=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2.0];
plot(weight,phase_synchrony_mean(5,:))

figure
errorbar(weight,phase_synchrony_mean(5,:),phase_synchrony_sd(5,:))

%% ---------------------- save
save('Phase_synchrony_metrics_OP.mat','phase_synchrony_mean')
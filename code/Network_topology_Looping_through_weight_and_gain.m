% Network Topology Heat maps as function of global gain and connection
% weights
clc; close all
addpath('/Users/RChenLab/Documents/TVB_Distribution/demo_scripts/Github/ThalCorProject/BCT/BCT/2017_01_15_BCT')
count_two=0;Integ_within=[];
for gain={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.8', '1','1.5','2.0'}
    count_two=count_two+1;
    count_one=0;
    for G={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5','0.6','0.8','1.0'}
        count_one=count_one+1;
        Data=load(['FC_LoopOP' char(G) 'Weight' char(gain) '.mat'])
        net=Data.data_struct;
%         [nRows,nCols] = size(net);
%         net(1:(nRows+1):nRows*nCols) = 0;
        SC_Data=load('Struct_data.mat')
        struct=SC_Data.data_struct;
        % if you want to get rid of negative weights
        % optional
        % net(net<0)=0;
        %net=abs(net);
        % net=threshold_proportional(net, 0.90); % see Hwang, tried 0.1 to 0.15 and averaged the result
        %% ---------- Community structure through Louvain's algorithm
        nNodes=size(net,1);
        ci = zeros(nNodes,1);
        gamma = 1;
        tau = 0.1; 
        nReps = 10;

        for x = 1:500
            [ci_temp(:,x),q_temp(x,1)] = community_louvain(net,gamma,1:1:nNodes,'negative_asym'); 
        end

        %estimate a 'consensus' partition (tau and nReps can be altered to change threshold - see https://sites.google.com/site/bctnet/)
        D = agreement(ci_temp);
        ci = consensus_und(D,tau,nReps);
        Module(count_one,count_two)=length(unique(ci));
        q = nanmean(q_temp);

        %% ------------------------------- Participation index (BA)  First do the community structure estimation before this step
        BA = zeros(nNodes,1);
        BA = participation_coef_sign(net,ci);
        Integ_BA(count_one,count_two)=mean(BA)
        Thal_PC(count_one,count_two)=mean(BA([42,90]));
        Others=BA;
        Others([43, 91])=[];
        other_PC(count_one,count_two)=mean(Others);
        % subplot(2,4,6)
        % info_plot=BA
        % hist(info_plot)
        % thal.deg=info_plot([41, 42, 43, 89, 90, 91])
        % min(thal.deg)
        % max(thal.deg)
        % title('participation_coef')
        % 
        % hold on
        % y=0:0.001:20; % How much is long
        % x=ones(size(y))*min(thal.deg);
        % plot(x, y, 'r') % but is not large enough
        % % Thalamic nodes do not seem to be that high degree
        % 
        % hold on
        % y=0:0.001:20; % How much is long
        % x=ones(size(y))*max(thal.deg);
        % plot(x, y, 'r') % but is not large enough
        %% ---------------------------- Within-Module degree
        Z=module_degree_zscore(net,ci);
        Integ_within(count_one,count_two)=mean(Z);
        Thal_within(count_one,count_two)=mean(Z([42,90]));
        Others=Z;
        Others([43, 91])=[];
        other_within(count_one,count_two)=mean(Others);
        %% -------------------------- Integration: inverse modularity Q^-1
        Integ_q(count_one,count_two)=q^-1;

        %% --------------------------- Global efficiency
%         L=net.^(-1);  % is this correct? Shine do not explain how they convert their weight matrix into a distance matrix
%         [D,B] = distance_wei(L);
%         [lambda,efficiency,ecc,radius,diameter] = charpath(D);
        % efficiency wei it takes the weighted connection matrix local
        GE(count_one,count_two)=efficiency_wei(net);
        %% --------------------------- Clustering coefficient
        [C_pos,C_neg,Ctot_pos(count_one,count_two),Ctot_neg] = clustering_coef_wu_sign(net);  % (default) Onnela et al. formula, used in original clustering_coef_wu.m. Computed separately for positive & negative weights

        %% Plot grid_communities
        % [X,Y,INDSORT] = grid_communities(ci); % call function
        % figure
        % imagesc(net(INDSORT,INDSORT));           % plot ordered adjacency matrix
        % hold on;                                 % hold on to overlay community visualization
        % plot(X,Y,'r','linewidth',2);   


        %% --------------------- compute correlation between structural and functional data
        CORR(count_one,count_two)=corr(reshape(net,[96*96,1]),reshape(struct,[96*96,1]))
        
        %% ----------------------- Newman's modularity

        % gamma is the resolution parameter; gamma=1 classic modularity
        gamma=1;
        [Ci Q(count_one,count_two)] = modularity_und(net,gamma);

    end
end

%% This explores SC data and calculates which nodes fall into the same degree match
%% Degree calculation and find degree matched nodes (falling in Mean+-2SD)
%figure
% [is,os,str] = strengths_dir(struct)
% 
% subplot(2,4,2)
% info_plot=str
% hist(info_plot)
% thal.deg=info_plot([42, 43, 44, 90, 91, 92])
% MEAN=mean(thal.deg);
% SD=sqrt(var(thal.deg))
% min(thal.deg)
% max(thal.deg)
% title('strength measure')
% 
% hold on
% y=0:0.001:20; % How much is long
% x=ones(size(y))*(MEAN+2*SD);
% plot(x, y, 'r') % but is not large enough
% % Thalamic nodes do not seem to be that high degree
% 
% hold on
% y=0:0.001:20; % How much is long
% x=ones(size(y))*(MEAN-2*SD);
% plot(x, y, 'r') % but is not large enough
% [a,deg_matched_indeces]=find((MEAN+2*SD)>str & str>(MEAN-2*SD));
% IND=round(rand(1,6)*length(a))
% Nodes_rand=deg_matched_indeces(IND);
% matched=str(deg_matched_indeces(IND));
%% Plot as a function of global coupling
weight=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2];
gain=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1];
figure
subplot(1,5,1)
imagesc(weight,gain,Integ_q)
title('inverse modularity')
subplot(1,5,2)
imagesc(weight,gain,Integ_BA)
title('participation coefficient')
% subplot(1,5,3)
% Q=Q.^(-1);
% imagesc(weight,gain,Q)
% title('Newman modularity')
subplot(1,5,4)
imagesc(weight,gain,Ctot_pos)
title('clustering coefficient')
subplot(1,5,5)
imagesc(weight,gain,CORR)
title('SC FC correlation')
subplot(1,5,3)
imagesc(weight,gain,GE)
title('GE')
%% Line plot at g=0.4
figure
subplot(1,5,1)
plot(weight,Integ_q(5,:))
title('inverse modularity')
axis([0 2 1.5 6])
subplot(1,5,2)
plot(weight,Integ_BA(5,:))
title('participation coefficient')
axis([0 2 0.15 0.4])
subplot(1,5,4)
plot(weight,Ctot_pos(5,:))
title('clustering coefficient')
axis([0 2 0.05 0.35])
subplot(1,5,5)
plot(weight,CORR(5,:))
title('SC FC correlation')
axis([0 2 0.34 0.43])
subplot(1,5,3)
plot(weight,GE(5,:))
title('GE')



%% comparing thalamic nodes to other nodes in terms of PC and within-module degree
figure
subplot(1,2,1)
plot(weight,Thal_PC(5,:),'b',weight,other_PC(5,:),'r')
subplot(1,2,2)
plot(weight,Thal_within(5,:),'b',weight,other_within(5,:),'r')


%% Degree-matched random nodes
% clear;clc
% SC_Data=load('data_struct.mat')
% struct=SC_Data.data_struct;
% [is,os,str] = strengths_dir(struct)
% 
% subplot(2,4,2)
% info_plot=str
% hist(info_plot)
% thal.deg=info_plot([41, 42, 43, 89, 90, 91])
% MEAN=mean(thal.deg);
% SD=sqrt(var(thal.deg))
% min(thal.deg)
% max(thal.deg)
% title('strength measure')
% 
% hold on
% y=0:0.001:20; % How much is long
% x=ones(size(y))*(MEAN+2*SD);
% plot(x, y, 'r') % but is not large enough
% % Thalamic nodes do not seem to be that high degree
% 
% hold on
% y=0:0.001:20; % How much is long
% x=ones(size(y))*(MEAN-2*SD);
% plot(x, y, 'r') % but is not large enough
% XX=str; 
% [a,deg_matched_indeces]=find((MEAN+2*SD)>XX & XX>(MEAN-2*SD));
% TH=[41, 42, 43, 89, 90, 91];
% deg_matched_indeces
% for i=1:6
% deg_matched_indeces(find(deg_matched_indeces==TH(i)))=[];
% end
% deg_matched_indeces
% Nodes_rand=[]
% matche=[]
% for count=1:10
% IND=randi(length(deg_matched_indeces),1,6);
% Nodes_rand(count,:)=deg_matched_indeces(IND);
% matched(count,:)=str(deg_matched_indeces(IND));
% end
% 

%% ----------------- Save the data here
save('FC_metrics_OP.mat','Integ_BA','Integ_q','GE','CORR')
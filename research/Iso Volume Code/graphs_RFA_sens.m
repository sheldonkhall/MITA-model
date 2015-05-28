%% graphs_RFA_sens
%
% This script takes the full_results structure and experiment matrices and
% produces a set of visualisations of the results. Needs to be run after
% lesion_volume_local_dir.

close all
clear all

% set working directory
fp = '/data/GO-SMART/MITA-model/RFA-sensitivity-cool/';

% set folders to take results from in order {experiment-1,experiment-2,etc}
folders = {'20150217-170854-experiment-1',...
    '20150218-172328-experiment-2',...
    '20150219-111815-experiment-3',...
    '20150219-210431-experiment-4',...
    '20150222-155751-experiment-5',...
    '20150223-120141-experiment-6'};

% set parameters to display

% select_param = [true true true true true true true true true true true true];
select_param = [false false true true true false true false false true false true]; % select params to display cd ------------------------------
% select_param = [true false true true true true true false false false false false]; % temp -----------------------------------------

% set data file

% dfile = '/runlesion_area.mat';
dfile = '/runlesion_area.cd.mat';
% dfile = '/runlesion_area.temp.mat';

% cycle through folders and load experiments and results
store_results = cell(length(folders),1);
experiments = cell(1,length(folders));
for i=1:length(folders)
    
    eval(['load ' fp folders{i} dfile]) % change data -----------------------------
    store_results{i} = full_results;
    eval(['load ' fp 'experiment-' num2str(i) '.mat'])
    experiments{i} = A;
    
end

% for each time compute sensitivities and store
mn = zeros(size(A,2),length(j));
mvt = zeros(length(j),1);
vvt = mvt;
sd=mn;
for i=1:length(j)
    
    responses = cell(1,length(folders));
    
    % use results from every experiment for mean and sd
    for k=1:length(folders)
        responses{1,k} = [store_results{k}(i,:).area];
    end
    
    % collect main effects
    [mn(:,i) sd(:,i)] = Process_Results(experiments,responses);
    
    % collect response means
    mvt(i) = mean([responses{:}]);
    vvt(i) = std([responses{:}]);
  
    % collect all responses
    
    
end

% define parameter names
param_names = {'$\rho c$','$\rho_{vap} c_{vap}$','$k_0$','$\Delta k$',...
    '$\sigma_0$','$\Delta \sigma$','$\omega_b$','$C$',...
    '$\Delta T$','$\xi$','$\sigma_{vap}$','$G_0$'};

%% morris method at final time point
figure
plot(mn(:,end-1),sd(:,end-1),'or',...
    [0 max(abs(mn(:,end-1)))], [0 max(abs(mn(:,end-1)))],'b-',...
    [0 -max(abs(mn(:,end-1)))], [0 max(abs(mn(:,end-1)))],'b-','LineWidth',4)
% text(mn(:,end)+5e-7,sd(:,end),cellfun(@num2str,num2cell(1:size(A,2)),'UniformOutput',false))
text(mn(select_param,end-1)+2e-6,sd(select_param,end-1),param_names(select_param),'Interpreter','LaTex','FontSize',20)
xlabel('mean','Interpreter','LaTex','FontSize',20)
ylabel('standard deviation','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

% set(gcf,'PaperUnits', 'Inches', 'PaperPosition', [0 0 10 10])
print -deps figure1

%% mn and sd vs time

% select times
times=[1:length(j)-1]; %----------------------------
% times=[1:100]; %---------------------------
% times=[1:13]; %---------------------------

figure
subplot(3,1,1)
% plot(j,cumsum(abs(mn)'))
plot(j(times),abs(mn(select_param,times)'),'LineWidth',4)
axis([0 j(times(end)) 0 max(max(abs(mn(:,times)))) ])
grid on
ylabel('mean','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
subplot(3,1,2)
plot(j(times),sd(select_param,times)','LineWidth',4)
% legend(cellfun(@num2str,num2cell(1:size(A,2)),'UniformOutput',false))
legend({param_names{select_param}},'Interpreter','LaTex','FontSize',20,'Location','NorthEast')
axis([0 j(times(end)) 0 max(max(abs(mn(:,times)))) ])
grid on
ylabel('standard deviation','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

%% plot mean lesion area vs time
subplot(3,1,3)
plot(j(times),mvt(times),j(times),(mvt(times)+vvt(times)),'r',j(times),(mvt(times)-vvt(times)),'r','LineWidth',4)
grid on
xlabel('time (s)','Interpreter','LaTex','FontSize',20)
ylabel('mean ablation area (m$^2$)','Interpreter','LaTex','FontSize',20)
legend({'mean','$\pm$ standard deviation'},'Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
axis([0 j(times(end)) 0 max(max(abs(mvt(times)+vvt(times))))*1.5])

set(gcf,'PaperUnits', 'Inches', 'PaperPosition', [0 0 10 10])
print -depsc figure2

%% plot integrated morris method
figure
plot(trapz(j,abs(mn')),trapz(j,sd'),'or',...
    [0 max(trapz(j,abs(mn')))], [0 max(trapz(j,abs(mn')))],'b-',...
    [0 -max(trapz(j,abs(mn')))], [0 max(trapz(j,abs(mn')))],'b-','LineWidth',4)

text(trapz(j,abs(mn'))+0.0001,trapz(j,sd'),param_names,'Interpreter','LaTex','FontSize',20)
xlabel('mean','Interpreter','LaTex','FontSize',20)
ylabel('standard deviation','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

print -deps figure3

% ==============================================================================
% This is a top-level routine for generating all figures in the manuscript.
% Mechanistic modeling of continuous lyophilization.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Model Equations', 'Events','Exporting Graphics','Plotting', ...
    'Validation Data','Simulations','Calculations','Saved Data');

% Figure selection
Fig5 = 'off';  % model validation
Fig6 = 'off';  % complete simulation
Fig7 = 'off';  % spatiotemporal data for primary drying
Fig8 = 'off';  % spatiotemporal data for secondary drying
Fig9R1 = 'off';  % optimizing VISF
Fig9R2 = 'off';  % optimizing VISF
Fig9R3 = 'on';  % optimizing VISF
Fig9R4 = 'off';  % optimizing VISF
Fig10 = 'off';  % conventional primary drying
Fig11 = 'off';  % conventional secondary drying
Fig12 = 'off';  % condenser failure
FigA1 = 'off';  % for appendix A
FigA2 = 'off';  % for appendix A
FigA3 = 'off';  % for appendix A
data = 'off';
All = 'off';  % run all simulations
None = 'off';  % run nothing
saveplot = 0;  % set to 1 to save plots

% Check mode of simulation
switch All
case 'on'
    Fig5 = 'on'; Fig6 = 'on'; Fig7 = 'on'; Fig8 = 'on'; Fig9 = 'on'; Fig10 = 'on'; Fig11 = 'on'; Fig12 = 'on';
    FigA1 = 'on'; FigA2 = 'on'; FigA3 = 'on';
end

switch None
case 'on'
    Fig5 = 'off'; Fig6 = 'off'; Fig7 = 'off'; Fig8 = 'off'; Fig9 = 'off'; Fig10 = 'off'; Fig11 = 'off'; Fig12 = 'off';
    FigA1 = 'off'; FigA2 = 'off'; FigA3 = 'off';
end


%% Figure 5: Model validation
switch Fig5
case 'on'

% Continuous Freezing
load('Data_Con_Freeze.mat')
ip = get_inputdata;
ip = overwrite_inputdata(ip,'expcon1_freezing');
ip = input_processing(ip);

sol1 = Sim_Freezing(ip);
time1 = sol1.t;
Temp1 = sol1.T;
Tg = sol1.Tg;

fig = figure; 
tiledlayout(1,3,"TileSpacing","loose","Padding","compact");
nexttile(1)
plot(time1/3600,Temp1,'linewidth',2,'Color',[0 0.5 1, .9],'DisplayName','Model (Case 1)'); hold on; 
plot(Data(1:end,1),Data(1:end,2),'s','MarkerSize',4,'color','k','MarkerFaceColor',[0 0.2 1],'DisplayName','Experiment (Case 1)')
%plot(time1/3600,Tg,':','MarkerSize',6,'color',[0 0 0],'linewidth',1.5,'MarkerFaceColor',[0 0.5 1],'DisplayName','Gas temperature')
h = legend('location','best'); h.ItemTokenSize(1) = 10;
ylabel({'Product temperature (K)'}); xlabel('Time (h)')
graphics_setup('1by3')
text(.02,.1,'b','Units','normalized','FontSize', 11,'fontweight', 'bold');


% Continuous Primary Drying
% Case 2a
load('Data_Con_1stDrying_263K.mat');
ip = get_inputdata;
ip = overwrite_inputdata(ip,'expcon1_primdry');
ip.T02 = Data(1,2);
Data(:,1) = Data(:,1)-Data(1,1);
ip = input_processing(ip);
sol1 = Sim_1stDrying(ip);
time1 = sol1.t;
Temp1 = sol1.T;
Tp1 = Temp1(:,end); 

nexttile(2)
plot(time1/3600,Tp1,'linewidth',2.5,'Color',[0 0.8 0.17, .9],'DisplayName','Model (Case 2a)'); hold on; 
plot(Data(:,1),Data(:,2),'s','color','k','MarkerFaceColor',[0 0.62 0.17],'markersize',4,'DisplayName','Experiment (Case 2a)'); hold on
ylabel('Product temperature (K)'); xlabel('Time (h)');
h = legend('location','best'); h.ItemTokenSize(1) = 10;


% Case 2b
load('Data_Con_1stDrying_313K.mat');
ip = get_inputdata;
ip = overwrite_inputdata(ip,'expcon1_primdry');
% ip.Tb2 = @(x) min(263+(1/60)*x,313);
ip.Tb2 = 313;
ip.T02 = Data(1,2);
Data(:,1) = Data(:,1)-Data(1,1);

ip = input_processing(ip);
sol1 = Sim_1stDrying(ip);
time1 = sol1.t;
Temp1 = sol1.T;
S = sol1.S; 
Tp1 = Temp1(:,end); 
Tb1 = sol1.Tb;

plot(time1/3600,Tp1,'linewidth',2.5,'Color',[1 0.5 0.8, 1],'DisplayName','Model (Case 2b)'); hold on; 
plot(Data(:,1),Data(:,2),'^','color','k','MarkerFaceColor',[1 0.2 0.8],'markersize',4,'DisplayName','Experiment (Case 2b)');
ylabel('Product temperature (K)'); xlabel('Time (h)');
ylim([220 260])
h = legend('location','best'); h.ItemTokenSize(1) = 10;
graphics_setup('1by3');
set(gcf, 'units', 'centimeters', 'Position',  [10, 6, 8.1, 5.3]);
text(.02,.1,'c','Units','normalized','FontSize', 11,'fontweight', 'bold');


% Continuous Secondary Drying
% Case 3a
load('Data_Con_2ndDrying_1.mat');
ip = get_inputdata;
ip = overwrite_inputdata(ip,'expcon1_secdry');
ip = input_processing(ip);
sol = Sim_2ndDrying(ip);
time = sol.t;
cw = sol.cw; cw_avg = mean(cw,2);

nexttile(3)
plot(time/3600,cw_avg,'linewidth',2,'Color',[0 0.8 0.17, .9],'DisplayName','Model (Case 3a)'); hold on; 
errorbar(Data(:,1),Data(:,2),Data(:,3),'s','color','k','MarkerFaceColor',[0 0.62 0.17],'CapSize',8,'linewidth',.5,'markersize',4,'DisplayName','Experiment (Case 3a)'); hold on
ylabel('Average concentration (kg water/kg solid)'); xlabel('Time (h)');
h = legend('location','best'); h.ItemTokenSize(1) = 10;
graphics_setup('1by3')


% Case 3b
load('Data_Con_2ndDrying_2.mat');
ip = get_inputdata;
ip = overwrite_inputdata(ip,'expcon1_secdry');
ip.cw0 = 0.075;
ip = input_processing(ip);
sol = Sim_2ndDrying(ip);
time = sol.t;
cw = sol.cw; cw_avg = mean(sol.cw,2);
Temp = sol.T; Tp = Temp(:,end); Tb = sol.Tb;

plot(time/3600,cw_avg,'linewidth',2,'Color',[1 0.5 0.8, 1],'DisplayName','Model (Case 3b)'); hold on; 
errorbar(Data(:,1),Data(:,2),Data(:,3),'^','color','k','MarkerFaceColor',[1 0.2 0.8],'CapSize',8,'linewidth',.5,'markersize',4,'DisplayName','Experiment (Case 3b)'); hold on
ylabel({'Average concentration' ; '(kg water/kg solid)'}); xlabel('Time (h)');
h = legend('location','best'); h.ItemTokenSize(1) = 10;
graphics_setup('1by3')
text(.02,.1,'d','Units','normalized','FontSize', 11,'fontweight', 'bold');

if saveplot == 1; export_figures(fig,'Validation'); end

end

%% Figure 6: Complete simulation
switch Fig6
case 'on'

% Parameters
ip0 = get_inputdata;
% ip0.Vl = 3e-6;  % modify any inputs here
ip = input_processing(ip0);

% Simulation and obtain solutions
% sol 1 = freezing, sol 2 = primary drying, sol 3 = secondary drying
tic; [sol1, sol2, sol3] = Sim_Lyo(ip); toc;  

% Plotting
fig = figure; 
plot_all2(sol1,sol2,sol3,ip)

if saveplot == 1; export_figures(fig,'Complete_lyo'); end

end


%% Figure 7: Spatiotemporal data for primary drying
switch Fig7
case 'on'

% Parameters
tic
ip0 = get_inputdata; 
ip0.H2 = 0.02;  
ip0.hb2 = 30;
ip0.nz2 = 100;
ip = input_processing(ip0);
sol2 = Sim_1stDrying(ip);

% Plot
fig = figure;
plot3D_T_primdrying(ip,sol2)

if saveplot == 1; export_figures(fig,'Lyo3D_13'); end

end


%% Figure 8: Spatiotemporal data for secondary drying
switch Fig8
case 'on'

% Parameters
tic
ip0 = get_inputdata; 
ip0.H3 = 0.02;  
ip0.hb3 = 30;
ip0.nz3 = 100;
ip0.cw0 = linspace(.05,.2,ip0.nz3)';
ip = input_processing(ip0);
sol3 = Sim_2ndDrying(ip);

% Plot
fig = figure;
plot3D_cw_secdrying(ip,sol3)
% axis equal

if saveplot == 1; export_figures(fig,'Lyo3D_2'); end

end


%% Figure 9: Optimizing VISF
switch Fig9R1
case 'on'

% Plotting
fig = tiledlayout(1,2,"TileSpacing","loose","Padding","compact");


% Parameters
ip0 = get_inputdata;
P = [1e4; 5e3; 1e3; 1e2];
tP = 60;
Ptot{1} = [1e5, P(1), P(1); 0, tP, 1000]'; 
Ptot{2} = [1e5, P(2), P(2); 0, tP, 1000]'; 
Ptot{3} = [1e5, P(3), P(3); 0, tP, 1000]'; 
Ptot{4} = [1e5, P(4), P(4); 0, tP, 1000]'; 

% Plotting
pressure = {'10^4 Pa','5000 Pa','1000 Pa','100 Pa'};
linspec = {'-','-.',':','-','-diamond','-*'};
msize = [6;4;4.5;4.5;4.5;4];
lw = [1.5;1.5;1.5;.9;.5;1];
color = {[0,0,139]/256,[0,0,255]/256,[65,105,225]/256,[0,191,255]/256,'m'};
np = 4;

for i = 1:np
    ip0 = get_inputdata;
    ip0 = overwrite_inputdata(ip0,'StoVISF_freezing');
    ip0.Ptot = Ptot{i};
    ip = input_processing(ip0);
    sol = Sim_Freezing_VISF(ip);
    time = sol.t;
    Temp = sol.T;
    Tg = sol.Tg;

    nexttile(1)
    % iend = find(time>2.2,1);
    plot(time/3600,Temp,linspec{i},'Color',color{i},'MarkerFaceColor',color{i},'linewidth',lw(i),'markersize',msize(i)); hold on 
    ylabel({'Product temperature (K)'}); xlabel('Time (h)')
    graphics_setup('1w')
    xlim([0.1 .8])
    ylim([255 280])
    % text(.84,.5,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
    if i == np
        h = legend(pressure,'location','best');
        h.ItemTokenSize(1) = 15;
    end

end

P = linspace(100,1e4,50);
np = length(P);
mloss = zeros(np,1);
tn = zeros(np,1);

for i = 1:np
    Ptot = [1e5, P(i), P(i); 0, tP, 1000]'; 
    ip0 = get_inputdata;
    ip0 = overwrite_inputdata(ip0,'StoVISF_freezing');
    ip0.Ptot = Ptot ;
    ip = input_processing(ip0);
    sol = Sim_Freezing_VISF(ip);
    mloss(i) = sol.mloss;
    tn(i) = sol.tn;

end

nexttile(2)
colororder([0.5,0,0.5;.5,0,1])

yyaxis left
semilogx(P,tn,':','color',[0.5,0,0.5],'LineWidth',1.5)
ylabel('Time of first nucleation (s)'); 

yyaxis right
semilogx(P,mloss,'-','color',[0.5,0,1],'LineWidth',1.5)
ylabel('Mass loss (kg)','Rotation',270); 
h = legend('Time','Mass','Location','west');
h.ItemTokenSize(1) = 15;

xlabel('VISF Pressure (Pa)')
graphics_setup('1w')
% text(.84,.5,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');

if saveplot == 1; export_figures(fig,'VISF1'); end

end


%% Figure 9: Optimizing VISF
switch Fig9R2
case 'on'

% Plotting
fig = tiledlayout(1,3,"TileSpacing","loose","Padding","compact");

% Parameters
ip0 = get_inputdata;
P = [1e4; 5e3; 1e2];
tP = 60;
Ptot{1} = [1e5, P(1), P(1); 0, tP, 1000]'; 
Ptot{2} = [1e5, P(2), P(2); 0, tP, 1000]'; 
Ptot{3} = [1e5, P(3), P(3); 0, tP, 1000]'; 

% Plotting
pressure = {'5000 Pa','1000 Pa'};
linspec = {'-','-.',':','-','-diamond','-*'};
msize = [6;4;4.5;4.5;4.5;4];
lw = [1.5;1.5;1.5;.9;.5;1];
color = {'b',[0 0.62 0.17],'r','m',[0.9290 0.5940 0.1250]};
np = 3;
nmc = 1e3;  % change to 1e3 when creating plots
tn = zeros(np,nmc)';
Tn = zeros(np,nmc)';
mn = zeros(np,nmc)';

% rng(1)
for i = 1:np
    ip0 = get_inputdata;
    if i == 1
        ip0 = overwrite_inputdata(ip0,'stochastic_freezing');
    else
        ip0 = overwrite_inputdata(ip0,'StoVISF_freezing2');
        ip0.Ptot = Ptot{i};
    end

    ip = input_processing(ip0);
    tref = unique([(0:20:ip.tpost1)';ip.tpost1]);
    % Tref = zeros(length(tref),nmc);
    Temp = cell(nmc,1);
    time = cell(nmc,1);
    
    parfor j = 1:nmc
        rng(j,'twister')
        if i == 1
            sol = Sim_Freezing_Sto2(ip);
        else
            sol = Sim_Freezing_StoVISF2(ip);
        end
        
        time{j} = sol.t;
        Temp{j} = sol.T;
        Tg = sol.Tg;
        % 
        % Tref(:,j) = interp1(time,Temp,tref);
        Tn(j,i) = sol.Tn;
        tn(j,i) = sol.tn;
        mn(j,i) = sol.mn;
   
    end

    tref = unique(sort([tref;vertcat(tn(:))]));
    Tref = zeros(length(tref),nmc);

    parfor j = 1:nmc
        Tref(:,j) = interp1(time{j},Temp{j},tref);
    end

    [Tmean, Tmin, Tmax] = cal_CI(Tref',.95);
    nexttile(i)
    plot_T_MeanCI(tref'/3600,Tmean,Tmin-1e-3,Tmax+1e-3,{'Color', color{i},'Linewidth',1.5},{color{i}, 'FaceAlpha',.3,'LineStyle','none'});
    ylabel('Product temperature (K)')
    lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
    if i == 1
        text(.04,.1,'Uncontrolled nucleation','Units','normalized','FontSize', 9);
    else
        text(.04,.1,['VISF at ' num2str(P(i)) ' Pa'],'Units','normalized','FontSize', 9);
    end
    xlim([0.1 .8])
    ylim([255 280])
    graphics_setup('1w')
    Tsd = std(Tn);
    tsd = std(tn);


end


if saveplot == 1; export_figures(fig,'VISF2'); end

end


%% Figure 9: Optimizing VISF
switch Fig9R3
case 'on'

load('VISF_mass_n_30.mat')
load('VISF_temp_n_30.mat')
load('VISF_time_n_30.mat')
load('VISF_pressure_30.mat')
fig = tiledlayout(1,2,"TileSpacing","loose","Padding","compact");

for i = 1:width(tn)
    [tmean(i), tmin(i), tmax(i)] = cal_CI(tn(:,i),.95);
end

for i = 1:width(Tn)
    [Tmean(i), Tmin(i), Tmax(i)] = cal_CI(Tn(:,i),.95);
end

for i = 1:width(mn)
    [mmean(i), mmin(i), mmax(i)] = cal_CI(mn(:,i),.95);
end

nexttile(1)
plot_tP_MeanCI_log(P,tmean,tmin,tmax,{'Color', [.5 0 .5],'Linewidth',1.5},{'m', 'FaceAlpha',.2,'LineStyle','none'});
%plot_tP_MeanCI_log([P,1e5],[tmean(2:end),tmean(1)],[tmin(2:end),tmin(1)],[tmax(2:end),tmax(1)],{'Color', [.5 0 .5],'Linewidth',1.5},{'m', 'FaceAlpha',.2,'LineStyle','none'});
hold on;
yneg = tmean(1) - tmin(1);
ypos = tmax(1) - tmean(1);
% plot([P(1);P(end)],[tmean(1);tmean(1)],':r','linewidth',1.5)
% errorbar(1e4,tmean(1),yneg,ypos,'o','Color', [.5 0 .5],'MarkerFaceColor',[.5 0 .5])
% xlim([1e2 1.2e4])
lg = legend({'Mean','95% CI'},'location','west'); lg.ItemTokenSize(1) = 15;
graphics_setup('1w')


nexttile(2)
ip0 = get_inputdata;
ip0 = overwrite_inputdata(ip0,'VISF');
ip = input_processing(ip0);  

sol = Sim_Freezing_VISF_exp(ip);
time = sol.t; Temp = sol.T; Tg = sol.Tg;
plot(time,Temp,'linewidth',2,'Color',[0 0.5 1, .9]); hold on;
load('Data_VISF.mat')
plot(Data.T(:,2),Data.T(:,1),'sk','MarkerSize',6,'MarkerFaceColor',[0 0.2 1])
ylabel('Product temperature (K)'); xlabel('Time (s)')
h = legend('Model','Experiment','Location','best');
h.ItemTokenSize(1) = 15;
graphics_setup('1w')

% nexttile(2)
% plot_tP_MeanCI_log(P,Tmean(2:end),Tmin(2:end),Tmax(2:end),{'Color', [.5 0 .5],'Linewidth',1.5},{'m', 'FaceAlpha',.2,'LineStyle','none'});
% lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
% graphics_setup('1w')
% ylabel('Nucleation temperature (K)'); 
% 
% nexttile(3)
% plot_tP_MeanCI_log(P,mmean(2:end),mmin(2:end),mmax(2:end),{'Color', [.5 0 .5],'Linewidth',1.5},{'m', 'FaceAlpha',.2,'LineStyle','none'});
% lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
% graphics_setup('1w')
% ylabel('Mass of first nucleus (K)'); 

% nexttile(2)
% colororder([210/255,105/255,30/255;139/255,69/255,19/255])
% yyaxis left
% semilogx(P,Tmean,'-','color',[210/255,105/255,30/255],'LineWidth',1.5)
% ylabel('Mean nucleation temperature (K)'); 
% yyaxis right
% semilogx(P,mmean,':','color',[139/255,69/255,19/255],'LineWidth',1.5)
% ylabel('Mean mass of first nucleus (kg)','Rotation',270); 
% xlabel('VISF Pressure (Pa)')
% h = legend('Temperature','Mass','Location','east');
% h.ItemTokenSize(1) = 15;
% graphics_setup('1w')

if saveplot == 1; export_figures(fig,'VISF3'); end

end


%% Figure 9: Optimizing VISF
switch Fig9R4
case 'on'

ip0 = get_inputdata;
tP = 60;

P = logspace(2,5,30);
np = length(P);
for i = 1:np
    Ptot{i} = [1e5, P(i), P(i); 0, tP, 1000]'; 
end

nmc = 1e3;
tn = zeros(np,nmc)';
Tn = zeros(np,nmc)';
mn = zeros(np,nmc)';

% rng(1)
for i = 1:np
    ip0 = get_inputdata;

    ip0 = overwrite_inputdata(ip0,'StoVISF_freezing2');
    ip0.Ptot = Ptot{i};

    ip = input_processing(ip0);

    parfor j = 1:nmc
        rng(j,'twister')
        sol = Sim_Freezing_StoVISF2(ip);
        Tn(j,i) = sol.Tn;
        tn(j,i) = sol.tn;
        mn(j,i) = sol.mn;

    end

end

end


%% Figure 10: Conventional primary drying
switch Fig10
case 'on'

ip = get_inputdata;
ip = overwrite_inputdata(ip,'expbatch1_primdry');
ip = input_processing(ip);

% Case 4a
load('Data_Batch_1stDrying_258K_Temp.mat');
load('Data_Batch_1stDrying_258K_Itf.mat');
sol1 = Sim_1stDrying(ip);
time1 = sol1.t;
Temp1 = sol1.T;
S = sol1.S; 
Tp1 = Temp1(:,end);

fig = figure;
tiledlayout(1,2); 
nexttile(1); plot(time1/3600,Tp1,'linewidth',2,'Color',[0 0.8 0.17, .9],'DisplayName','Model (Case 4a)'); hold on; 
plot(Data_Val1_Temp_Case258(:,1),Data_Val1_Temp_Case258(:,2),'s','color','k','MarkerFaceColor',[0 0.62 0.17],'markersize',5,'DisplayName','Experiment (Case 4a)');
ylabel('Product temperature (K)'); xlabel('Time (h)'); graphics_setup('1by2.5'); xlim([0 12]) 
h = legend('location','best'); h.ItemTokenSize(1) = 10;
nexttile(2); plot(time1/3600,S*100,'linewidth',2,'Color',[0 0.8 0.17, 0.9],'DisplayName','Model (Case 4a)'); hold on; 
plot(Data_Val1_ITF_Case258(:,1),Data_Val1_ITF_Case258(:,2),'s','color','k','MarkerFaceColor',[0 0.62 0.17],'markersize',5,'DisplayName','Reference (Case 4a)'); 
ylabel('Sublimation front position (cm)'); xlabel('Time (h)'); graphics_setup('1by2.5'); xlim([0 12])
h = legend('location','best'); h.ItemTokenSize(1) = 10;


ip = get_inputdata;
ip = overwrite_inputdata(ip,'expbatch1_primdry');
ip.Tb2 = @(x) min(228.15+(.25/60)*x,268.15);
ip = input_processing(ip);

% Case 4b
load('Data_Batch_1stDrying_268K_Temp.mat');
load('Data_Batch_1stDrying_268K_Itf.mat');
sol1 = Sim_1stDrying(ip);
time1 = sol1.t;
Temp1 = sol1.T;
S = sol1.S; 
Tp1 = Temp1(:,end); Tb1 = sol1.Tb;

nexttile(1); plot(time1/3600,Tp1,'linewidth',2,'Color',[1 0.5 0.8, 1],'DisplayName','Model (Case 4b)'); hold on; 
plot(Data_Val1_Temp_Case268(:,1),Data_Val1_Temp_Case268(:,2),'^','color','k','MarkerFaceColor',[1 0.2 0.8],'markersize',4,'DisplayName','Experiment (Case 4b)');
ylabel('Product temperature (K)'); xlabel('Time (h)'); graphics_setup('1by2.5'); 
text(.03,.93,'(A)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');
h = legend('location','best'); h.ItemTokenSize(1) = 10;
nexttile(2); plot(time1/3600,S*100,'linewidth',2,'Color',[1 0.5 0.8, 1],'DisplayName','Model (Case 4b)'); hold on; 
plot(Data_Val1_ITF_Case268(:,1),Data_Val1_ITF_Case268(:,2),'^','color','k','MarkerFaceColor',[1 0.2 0.8],'markersize',4,'DisplayName','Reference (Case 4b)'); 
ylabel('Sublimation front position (cm)'); xlabel('Time (h)'); graphics_setup('1by2.5');
h = legend('location','southeast'); h.ItemTokenSize(1) = 10;
text(.03,.93,'(B)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');
% h = legend('position',[0,0,.8,.1],'Orientation','horizontal'); h.ItemTokenSize(1) = 10;

if saveplot == 1; export_figures(fig,'BatchPrimary'); end

end

%% Figure 11: Conventional secondary drying
switch Fig11
case 'on'

ip = get_inputdata;
ip = overwrite_inputdata(ip,'expbatch2_secdry');
ip = input_processing(ip);

% Case 5
cs_exp = load('Data_Batch_2ndDrying_Conc.mat').cs;
T_exp = load('Data_Batch_2ndDrying_Temp.mat').T;
sol2 = Sim_2ndDrying(ip);
time2 = sol2.t;
cw = sol2.cw; cw_avg = mean(sol2.cw,2);
Temp2 = sol2.T; Tp2 = Temp2(:,end); Tb2 = sol2.Tb;

fig = figure; 
tiledlayout(1,2);
nexttile; plot(time2/3600,Tp2,'linewidth',2,'Color',[0 0.5 1, .9],'DisplayName','Model (Case 5)'); hold on; 
plot(T_exp(:,1),T_exp(:,2),'s','MarkerSize',6,'color','k','MarkerFaceColor',[0 0.2 1],'DisplayName','Experiment (Case 5)')
h = legend('location','best'); h.ItemTokenSize(1) = 10;
ylabel({'Product temperature (K)'}); xlabel('Time (hours)')
text(.03,.93,'(A)','Units','normalized','FontSize', 9 ,'fontweight', 'bold' );
graphics_setup('1by2.5')
nexttile; plot(time2/3600,cw_avg,'linewidth',2,'Color',[0 0.5 1, .9],'DisplayName','Model (Case 5)'); hold on; 
errorbar(cs_exp(:,1),cs_exp(:,2),cs_exp(:,3),'s','color','k','MarkerFaceColor',[0 0.2 1],'CapSize',8,'linewidth',.5,'DisplayName','Experiment (Case 5)'); hold on
ylabel({'Average concentration' ; '(kg water/kg solid)'}); xlabel('Time (h)');
h = legend('location','best'); h.ItemTokenSize(1) = 10;
text(.03,.93,'(B)','Units','normalized','FontSize', 9 ,'fontweight', 'bold' );
graphics_setup('1by2.5')

if saveplot == 1; export_figures(fig,'BatchSecondary'); end

end


%% Figure 12: Condenser failure
switch Fig12
case 'on'

% Normal operation
ip = get_inputdata;
ip = input_processing(ip);
sol1 = Sim_1stDrying(ip);
time1 = sol1.t;
Temp1 = sol1.T;
Tp1 = mean(Temp1,2); 


fig = figure;
tiledlayout(1,3,"TileSpacing","loose","Padding","compact")
nexttile(3); plot(time1/3600,Tp1,'linewidth',2); hold on
ylabel({'Product temperature (K)'}); xlabel('Time (h)')
nexttile(2); plot(time1/3600,sol1.S*100,'linewidth',2); hold on
ylabel({'Sublimation front position (cm)'}); xlabel('Time (h)')
nexttile(1); plot(time1/3600,sol1.P,'linewidth',2); hold on
ylabel({'Pressure (Pa)'}); xlabel('Time (h)')
graphics_setup('1by3')

% Choked
ip = get_inputdata;
ip = input_processing(ip);
sol1 = Sim_1stDrying_Choked(ip);
time1 = sol1.t;
Temp1 = sol1.T;
Tp1 = mean(Temp1,2); Tb1 = sol1.Tb;
Ti = Temp1(:,1);

nexttile(1); plot(time1/3600,sol1.P,':','linewidth',2); hold on
ylabel({'Pressure (Pa)'}); xlabel('Time (h)')
ylim([0 25]); xticks(0:2:10); text(.83,.1,'a','Units','normalized','FontSize', 12,'fontweight', 'bold');
graphics_setup('1by3s')
nexttile(2); plot(time1/3600,sol1.S*100,':','linewidth',2); hold on
ylabel({'Sublimation front position (cm)'}); xlabel('Time (h)'); xticks(0:2:10);
text(.83,.1,'b','Units','normalized','FontSize', 12,'fontweight', 'bold');
h = legend({'Normal operation','Condenser failure'},'location','southoutside','Orientation','horizontal');
h.ItemTokenSize(1) = 15;
graphics_setup('1by3s')
nexttile(3); plot(time1/3600,Tp1,':','linewidth',2); hold on
ylabel({'Product temperature (K)'}); xlabel('Time (h)'); xticks(0:2:10);
text(.83,.1,'c','Units','normalized','FontSize', 12,'fontweight', 'bold');
graphics_setup('1by3s')
Nw = cal_Nw(sol1.T(:,1),sol1.S,sol1.P,ip);
j = ip.nvial*Nw*ip.Ac;

if saveplot == 1; export_figures(fig,'Condenser'); end

end


%% Figure A.1: Biot number
switch FigA1
case 'on'

outputs = cal_T_rt_cylinder(0,8,0.012,2.25,2108,917,3000,300,230,50,100,50,0);
T = outputs.T; Tavg = outputs.Tavg; Tlump = outputs.Tlump; r = outputs.r; 
t = outputs.t; err = outputs.err; Bi = outputs.Bi; dT = outputs.dT; Q = outputs.Q;
Qlump = outputs.Qlump;

fig = figure;
tiledlayout(1,2,"TileSpacing","loose","Padding","compact");

nexttile; plot(t,Tavg,'-b','linewidth',2); hold on; plot(t,Tlump,'-or','linewidth',1,'markersize',3,'MarkerIndices',1:10:length(t))
h = legend('1D','Lump','location','best'); h.ItemTokenSize(1) = 15;
xlabel('Time (s)'); ylabel('Average temperature (K)')
text(.85,.35,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
graphics_setup('1by2')

nexttile; plot(t,Q,'-b','linewidth',2); hold on; plot(t,Qlump,'-or','linewidth',1,'markersize',3,'MarkerIndices',1:10:length(t))
h = legend('1D','Lump','location','best'); h.ItemTokenSize(1) = 15;
xlabel('Time (s)'); ylabel('Total heat removed (W/m)')
text(.85,.35,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');
graphics_setup('1by2')

if saveplot == 1; export_figures(fig,'Biot_base'); end

end

%% Figure A.2: Biot number
switch FigA2
case 'on'

outputs = cal_T_rt_cylinder(0,65,0.012,2.25,2108,917,1000,300,230,50,100,50,0);
T = outputs.T; Tavg = outputs.Tavg; Tlump = outputs.Tlump; r = outputs.r; 
t = outputs.t; err = outputs.err; Bi = outputs.Bi; dT = outputs.dT; Q = outputs.Q;
Qlump = outputs.Qlump;

fig = figure;
tiledlayout(1,2,"TileSpacing","loose","Padding","compact");

nexttile; plot(t,Tavg,'-b','linewidth',2); hold on; plot(t,Tlump,'-or','linewidth',1,'markersize',3,'MarkerIndices',1:10:length(t))
h = legend('1D','Lump','location','best'); h.ItemTokenSize(1) = 15;
xlabel('Time (s)'); ylabel('Average temperature (K)')
text(.85,.35,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
graphics_setup('1by2')

nexttile; plot(t,Q,'-b','linewidth',2); hold on; plot(t,Qlump,'-or','linewidth',1,'markersize',3,'MarkerIndices',1:10:length(t))
h = legend('1D','Lump','location','best'); h.ItemTokenSize(1) = 15;
xlabel('Time (s)'); ylabel('Total heat removed (W/m)')
text(.85,.35,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');
graphics_setup('1by2')

if saveplot == 1; export_figures(fig,'Biot_max'); end

end

%% Figure A.3: Biot number
switch FigA3
case 'on'

Bi_all = [.01; .043; .35; 1];

fig = tiledlayout(1,2,"TileSpacing","loose","Padding","compact");
for i = 1:length(Bi_all)
    outputs = cal_T_rt_cylinder(1,50,1e-2,2,2000,917,1000,300,230,100,200,50,Bi_all(i));
    T = outputs.T; Tavg = outputs.Tavg; Tlump = outputs.Tlump; r = outputs.r; 
    t = outputs.t; err = outputs.err; Bi = outputs.Bi; dT = outputs.dT; Q = outputs.Q;
    Qlump = outputs.Qlump; eQ = outputs.eQ;
    
    nexttile(1)
    plot(t,err,'linewidth',1,'Displayname',['Bi = ' num2str(Bi)]); hold on
    h = legend('location','best'); h.ItemTokenSize(1) = 15;
    xlabel('Time (s)'); ylabel('Error (K)')
    text(.85,.35,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
    graphics_setup('1by2')

    nexttile(2)
    plot(t,abs(Q-Qlump),'linewidth',1,'Displayname',['Bi = ' num2str(Bi)]); hold on
    h2 = legend('location','best'); h2.ItemTokenSize(1) = 15;
    xlabel('Time (s)'); ylabel('Error (W/m)')
    text(.85,.35,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');
    graphics_setup('1by2')

end

if saveplot == 1; export_figures(fig,'Biot_all'); end

end

%% Collect Data
switch data
case 'on'
% Parameters
ip0 = get_inputdata;
tP = 60;

P = logspace(2,5,30);
np = length(P);
for i = 1:np
    Ptot{i} = [1e5, P(i), P(i); 0, tP, 1000]'; 
end

nmc = 1e3;
tn = zeros(np,nmc)';
Tn = zeros(np,nmc)';
mn = zeros(np,nmc)';

% rng(1)
for i = 1:np
    ip0 = get_inputdata;

    ip0 = overwrite_inputdata(ip0,'StoVISF_freezing2');
    ip0.Ptot = Ptot{i};

    ip = input_processing(ip0);

    parfor j = 1:nmc
        rng(j,'twister')
        sol = Sim_Freezing_StoVISF2(ip);
        Tn(j,i) = sol.Tn;
        tn(j,i) = sol.tn;
        mn(j,i) = sol.mn;

    end

end

end
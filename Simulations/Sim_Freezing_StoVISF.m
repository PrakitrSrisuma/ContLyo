function outputs = Sim_Freezing_StoVISF(ip)

%% Cooling, Pre-conditioning
% Correct the gas temperature profile
if length(ip.Tg) > 1
    if ip.Tg(end,2) < ip.tpost1
        Tg_add = [ip.Tg(end,1), ip.tpost1];
        ip.Tg = [ip.Tg; Tg_add];
    end
end

% Parameters for ODE solver
tf = ip.tpre1;  % final time
dt = ip.dt1;  % data collection frequency from the ODE solver
tspan = unique([(0:dt:tf)';tf]);  % define the time span
T0 = ip.T01;  % initial conditions

% Solve the ODEs
opts_ode = odeset('Event', @(t,y) event_cooling_complete(t,y,ip),'RelTol',ip.tol,'AbsTol',ip.tol);
[t1,y1] = ode15s (@(t1,y1) ODE_FreezingCoolPre(t1,y1,ip), tspan, T0, opts_ode);
t = t1;
T = y1;
mv = [];


%% VISF
tv = []; Tv = []; mv = [];

% Parameters for ODE solver
y0 = [ip.mw; y1(end)];

% Parameters
prob_nuc = rand;  % probability for nucleation
tf = ip.tpost1;
dt = 1;  % data collection frequency from the ODE solver

% Simulation
prob = 0;
ip.Ptot(:,2) = ip.Ptot(:,2) + t(end);
tspan = [t(end) t(end)+dt];

while prob<prob_nuc
opts_ode2 = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
[t1v,y1v] = ode15s (@(t,y) ODE_FreezingVISF(t,y,ip), tspan, y0, opts_ode2);
Tend = y1v(end,2);
mend = y1v(end,1);

    if Tend < ip.Tf
        J = ip.bn*(ip.Tf-Tend)^ip.kn;
        prob = prob + J*ip.Vl*dt;
    end
    tv = [tv;t1v(end)];
    Tv = [Tv;Tend];
    mv = [mv;mend];
    y0 = [mend;Tend];
    tspan = [tv(end); tv(end)+dt];
    
end
outputs.Tn = Tend;
outputs.tn = tv(end);
outputs.mv = mv;
outputs.tv = tv;
outputs.Pv = cal_P(tv,ip.Ptot);
outputs.Tv = Tv;
t = [t;tv];
T = [T;Tv];

% Check if the temperature is above nucleation temperature
if T(end) > ip.Tnuc
    warning('VISF temperature is higher than the target nucleation temperature.')
end

mloss = mv(1)-mv(end);

m_ice = zeros(length(T),1);

% Update mass after VISF
ip.mws_new = ip.mws - mloss;  % total mass after evaporation via VISF
ip.ms_new = ip.ms;  % mass of solute after VISF (assumed non-volatile)
ip.mw_new = ip.mws_new-ip.ms_new;  % mass of water after evaporation via VISF
ip.m0 = ip.mw_new;  % total mass after evaporation via VISF


%% Nucleation
% Find the equilibrium temperature
K1 = 1;
K2 = -ip.Tf-ip.Tnuc-ip.dHfus*ip.mw/(ip.Cpws*ip.mws);
K3 = ip.dHfus*ip.mw*ip.Tf/(ip.Cpws*ip.mws) - ip.ms*(ip.Kf/ip.Ms)*ip.dHfus/(ip.Cpws*ip.mws) + ip.Tf*ip.Tnuc;
Teq_ini = 0.5*(-K2-sqrt(K2^2-4*K1*K3));  % new equilibrium temperature
mi = ip.mw - ip.ms*(ip.Kf/ip.Ms)/(ip.Tf-Teq_ini);  % mass of ice after first nucleation
 
% Collect data
tf = ip.tpost1;
t = [t;t(end)];
T = [T;Teq_ini];
m_ice = [m_ice;mi];
tspan = unique([(t(end):dt:tf)';tf]);


%% Solidification
% Solve the ODEs
opts_ode3 = odeset('Event', @(t,y) event_freezing_complete(t,y,ip),'RelTol',ip.tol,'AbsTol',ip.tol);
% [t2,y2] = ode15s (@(t,y) ODE_FreezingNucl(t,y,ip), tspan, mi, opts_ode3);
[t2,y2] = ode15s (@(t,y) ODE_FreezingNucl(t,y,ip), tspan, mi, opts_ode3);

outputs_tmp = cal_freezing_interface(y2,ip);
outputs.r = outputs_tmp.r; outputs.l = outputs_tmp.l; outputs.t_rl = t2; outputs.H = outputs_tmp.H;

% Collect data
Teq = ip.Tf - ((ip.Kf/ip.Ms)*(ip.ms_new./(ip.m0-y2)));
ip.mi_fin = y2(end);
t = [t;t2];
T = [T;Teq];
m_ice = [m_ice;y2];


%% Final Cooling
y0 = T(end);
tspan = unique([(t(end):dt:tf)';tf]);

% Solve the ODEs
opts_ode4 = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
[t3,y3] = ode15s (@(t,y) ODE_FreezingCoolPost(t,y,ip), tspan, y0, opts_ode4);  
t = [t;t3];
T = [T;y3];
m_ice = [m_ice; m_ice(end)*ones(length(t3),1)];


%% Export
P_profile = ip.Ptot;
P_profile = [[P_profile(1,1),0]; P_profile];
P_profile = [P_profile; [P_profile(1,1),t(end)]];

outputs.Tg = cal_Tg(t,ip.Tg);
outputs.Tw = cal_T(t,ip.Tc1);
outputs.t = t;
outputs.T = T;
outputs.S = ip.S0*ones(length(outputs.t),1);
outputs.cw = ip.cw0*ones(length(outputs.t),1);
outputs.mi = m_ice;
outputs.P  = cal_P(t,P_profile);

return
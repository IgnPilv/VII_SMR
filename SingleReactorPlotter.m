% ----------------------------------------------------------------------- %
% SingleReactorPlotter.m
% ----------------------------------------------------------------------- %
% SingleReactorPlotter.m
% MultipleReactorOptimizer.m
% calcz.m
% ODEsolver.m
% ----------------------------------------------------------------------- %
close all
clear all
clc
global P0 T0 Ft0 rhoc U a Rgas alpha rlist glbp Wspan2 rhob Ac sb hi
% -- DECLARATIONS OF INPUT PARAMETERS
% Wspan = linspace(0,20000,50000); % Range for the weight of the catalyst kg (xaxis of produced graphs)
Wspan = (0:5:70000); % Range for the weight of the catalyst kg (xaxis of produced graphs)

rlist = [];

% Simulation options [0/1]
% r1,r2,r3,heating,reactionheat,dP, inletflow scale
% glbp = [1 0 0 1 1 1 1];
glbp = [1 1 1 1 1 1 1];


P0 = 20*10^5; % Pa inlet pressure
T0 = 580+273.15; % K inlet temperature
Ta0 = 1200+273.15; % K inlet coolant temperature
% a-ch4 b-h2o c-co d-h2 e-co2
% Fa0 = glbp(7)*1292.0; % kmol/h inlet ch4 flow rate
% Fb0 = glbp(7)*5350.4; % kmol/h inlet h2o flow rate 3226
% Fc0 = glbp(7)*149.60; % kmol/h inlet co flow rate
% Fd0 = glbp(7)*367.20; % kmol/h inlet h2 flow rate
% Fe0 = glbp(7)*1; % kmol/h inlet co2 flow rate

Fa0 = glbp(7)*1274.0; % kmol/h inlet ch4 flow rate
Fb0 = glbp(7)*5352.5; % kmol/h inlet h2o flow rate 3226
Fc0 = glbp(7)*147.51; % kmol/h inlet co flow rate
Fd0 = glbp(7)*362.07; % kmol/h inlet h2 flow rate
Fe0 = glbp(7)*0; % kmol/h inlet co2 flow rate

Ft0 = Fa0 + Fb0 + Fc0 + Fd0 + Fe0; % kmol/h inlet total flow rate
% z = 1; % ideal gas for now
Rgas = 8.3144598; % kJ/kmol
sb = 5.670367*10^-8; % W?m?2?K?4
% --Tube Data
Dt = 0.0762; % m tube diameter
% Dt = 0.020; % m tube diameter

Nt = 1200; % number of tubes (used for pressure drop)
Ac = 3.1415*Dt*Dt/4*Nt; % m² tube x-section area (used for pressuredrop)
a = 4/Dt; % 1/m area to volume ratio

% --Furnace data
U = 100; % W/(m² K) overall heat transfer coefficient


% --Catalyst data
e = 0.605; % - voidage
rhoc = 2355.2; % kg/m³ particle density
rhob = 1362.0; % kg/m³ bed density
e = 0.605
Dp = 0.0174131; % m particle diameter

% % -- PRESSURE DROP CALCULATIONS
m = (Fa0*16.04+Fb0*18+Fc0*28.01+Fd0*2+Fe0*44)/3600; % kg/s Feed mass flow (used for pressure drop)
nu = 0.000027915; % kg/(m s) viscosity) Average value retrieved FROM THERMOWORKBENCH (used for pressure drop)
rho0 = P0*(m/(Ft0*1000/3600))/(T0*Rgas); % kg/m³ (used for pressure drop)
% 3.2533 kg/m³
G = m/(Ac); % kg/(m² s) mass flux (total mass / total area) (used for pressure drop)
Beta0 = (G/(rho0*Dp)) * ((1-e)/(e*e*e)) * (((150*(1-e)*nu)/Dp)+(1.75*G)); % (used for pressure drop)
alpha = (2*Beta0)/(Ac*rhob*P0); % (used for pressure drop)
u = m/rho0/Ac; %m/;s

% k = 0.1422; %W/m K
% Re = G*Dt/nu; % Reynolds
% Cpcorr = 700.743; %kJ/kg K
% Pr = Cpcorr*nu/k; %Prandtl
% Nu = 0.023*Re^0.8*Pr^0.4;
% hi = Nu*k/Dt

% % -- INPUT PARAMETERS FOR ODE SOLVER
i0 = [Fa0; Fb0; Fc0; Fd0; Fe0; P0; T0; Ta0; 0; 0; T0;]; % Initial values for the dependent variables.
% % -- CALL TO THE SOLVER FUNCTION
[w0,y0]=ode45(@ODEsolver,Wspan,i0); % Range of solved values for a range of catalyst mass stored in this array
% % -- INTERPOLATION OF OUTPUT DATA
% disp(y0(:,1));
Wf = interp1q(1-y0(:,1)./Fa0,w0,0.90) % kg needed mass of catalyst for 4750kmol/h of H2
% Wf = interp1q(y0(:,1)./Fa0,w0,50.0) % kg needed mass of catalyst for 4750kmol/h of H2
Tf = interp1(w0,y0(:,7),Wf)-273.15; % C outlet temperature
Taf = interp1(w0,y0(:,8),Wf)-273.15; % C outlet coolant temperature
Twf = interp1(w0,y0(:,11),Wf)-273.15; % C outlet coolant temperature
Pf = interp1(w0,y0(:,6),Wf); % Pa outlet pressure
Lf = Wf/(Ac*(1-e)*rhoc); % m Length of reactor
conv = 1-interp1(w0,y0(:,1),Wf)./Fa0; % - conversation of methane
Faf = interp1(w0,y0(:,1),Wf); % kmol/h outlet CH4
Fbf = interp1(w0,y0(:,2),Wf); % kmol/h outlet H2O
Fcf = interp1(w0,y0(:,3),Wf); % kmol/h outlet CO
Fdf = interp1(w0,y0(:,4),Wf); % kmol/h outlet H2
Fef = interp1(w0,y0(:,5),Wf); % kmol/h outlet CO2
Vf = Wf /(rhoc*(1-e)); % m³ total volume of reactor
y0(:,9)=y0(:,9)./3600*10^-6;
y0(:,10)=y0(:,10)./3600*10^-6;
Qf = interp1(w0,y0(:,9),Wf) % total heat generated J
Qff = interp1(w0,y0(:,10),Wf) % total heat generated J
dPbar = (P0-Pf)/100000; % bar pressure drop
Ftf = interp1(w0,y0(:,1),Wf)+interp1(w0,y0(:,2),Wf)+interp1(w0,y0(:,3),Wf)+interp1(w0,y0(:,4),Wf)+interp1(w0,y0(:,5),Wf);
% -- TABULAR RESULTS
fprintf(['INPUT:\nNumber of tubes: %i. Tube diameter: %1.3f m. Inlet temperature: %1.3f C. Inletcoolant temperature : %1.3f C. Inlet pressure : %1.3f Pa.\n'],Nt,Dt,T0-273,Ta0-273,P0 );


fprintf(['CH4 feed: %1.3f kmol/h. H2O feed: %1.3f kmol/h. CO feed:: %1.3f kmol/h.\n'],Fa0,Fb0,Fc0);
fprintf(['H2 feed: %1.3f kmol/h. CO2 feed: %1.3f kmol/h.\n'],Fd0,Fe0);
fprintf(['OUTPUT:\nCatalyst mass : %1.3f kg. Outlet Temperature : %1.3f C. Coolant Outlet Temperature : %1.3f C. \nOutlet Pressure : %1.3f Pa. Length : %1.3f m. Conversion : %1.3f. Volume : %1.3f m³.\nCH4: %1.3f. H2O: %1.3f. CO: %1.3f. H2: %1.3f. CO2: %1.3f. kmol/h\n'],Wf,Tf,Taf,Pf,Lf,conv,Vf,Faf,Fbf,Fcf,Fdf,Fef);
% T = table(Nt,T0-273,P0,Ta0-273,Fb0,Fp0,Wf,Tf,Pf,Taf,Lf,conv,V,dPbar); % Generate a table of results
% T.Properties.VariableNames = {'Nt' 'T0' 'P0' 'Ta0' 'Fb0' 'Fp0' 'W' 'Tf' 'Pf' 'Taf' 'L' 'X' 'Vf' 'dPf'}

% L0 = w0./(Nt*Ac*(1-e)*rhoc); % Length of the reactor
% -- PLOTTER

test = [Wspan',w0,y0,1-(y0(:,1))./Fa0];

fig1 = figure;
graph1 = plot(w0,(y0(:,1)),w0,(y0(:,2)),w0,(y0(:,3)),w0,(y0(:,4)),w0,(y0(:,5)));
legend('CH4 kmol/h','H2O kmol/h','CO kmol/h','H2 kmol/h','CO2 kmol/h')
xlabel('W, kg')
% xlabel('L, m')
ylabel('Flow, kmol/h')
title('Reactor flow rate Profiles')

fig2 = figure;
hold all;
dummy1 =(y0(:,1))+(y0(:,2))+(y0(:,3))+(y0(:,4))+(y0(:,5));
graph2 = plot(w0,(y0(:,1))./dummy1,w0,(y0(:,2))./dummy1,w0,(y0(:,3))./dummy1,w0,(y0(:,4))./dummy1,w0,(y0(:,5))./dummy1);
legend('CH4','H2O','CO','H2','CO2')
xlabel('W, kg')
ylabel('Mole Fraction')
title('Reactor flow rate Profiles')

fig3 = figure;
hold all;
graph3 = plot(w0,(y0(:,10)),w0,(y0(:,9)),w0,-(y0(:,9))+(y0(:,10)));
legend('Furnace','Reaction','Sensible')
xlabel('W, kg')
ylabel('Q, MW')
title('Heat consumption')

axelf =(y0(:,1));
fig4 = figure;
hold all;
graph4 = plot([Lf Lf], [0 conv], [0 Lf],[conv conv],w0/(Ac*(1-e)*rhoc),1-axelf./Fa0,'*');
xlabel('Length of reactor tube, m')
ylabel('Conversion')
title('Methane Conversion')

fig5 = figure;
hold all;
graph5 = plot([Lf Lf], [0 Taf],w0/(Ac*(1-e)*rhoc),(y0(:,7))-273,'*',w0/(Ac*(1-e)*rhoc),(y0(:,8))-273,'*',w0/(Ac*(1-e)*rhoc),(y0(:,11))-273,'*');
legend('','Reactor','Furnace','Wall')
xlabel('Length of reactor tube, m')
ylabel('Temperature, C')
title('Temperature Profile')

fig6 = figure;
hold all;
graph6 = plot(w0,(y0(:,6))/100000);
legend('Reactor')
xlabel('W, kg')
ylabel('Pressure, bar')
title('Pressure Profile')

set(fig1,'units','normalized','outerposition',[0 0.5 .333 .5])
set(fig2,'units','normalized','outerposition',[0 0 .333 .5])
set(fig3,'units','normalized','outerposition',[0.333 0 .333 .5])
set(fig4,'units','normalized','outerposition',[0.333 0.5 .333 .5])
set(fig5,'units','normalized','outerposition',[0.666 0.5 .333 .5])
set(fig6,'units','normalized','outerposition',[0.666 0 .333 .5])
rlist=rlist';
% % --Display reaction rates
% % fprintf([repmat('%7.4e    ',1,size(rlist,2)),'\n'],rlist')
% % L0 = w0./(Ac*(1-e)*rhoc) % Length of the reactor
% disp(1-(y0(1:10,1)./Fa0));
avgQ = Qff/(Nt*Dt*3.14*Lf)*1000 %kW/m2
% Twf = interp1(w0,y0(:,11),Wf)-273.15 % C outlet coolant temperature
acpF = Qff/1;
Atubes = Nt*Dt*3.14*Lf
Acp15 = Nt*Dt*1.5*Lf

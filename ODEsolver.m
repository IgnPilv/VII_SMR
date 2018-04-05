% ----------------------------------------------------------------------- %
% ODEsolver.m
% Cumene production packed bed reactor ODE solver
% ----------------------------------------------------------------------- %
%
% Malcolm Bell, Sophie Fitchett, Adam Machon, Jack Palmer, Ignas Pilvinis
% S1433709 S1405595 S1426444 S1439172 S1452897
% 07/11/2017
% ----------------------------------------------------------------------- %
% SingleReactorPlotter.m
% MultipleReactorOptimizer.m
% calcz.m
% ODEsolver.m
% ----------------------------------------------------------------------- %
function dYfuncvecdW = ODEsolver(W,Yfuncvec);
Fa = Yfuncvec(1); % kmol/h ch4 flow rate
Fb = Yfuncvec(2); % kmol/h h2o flow rate
Fc = Yfuncvec(3); % kmol/h co flow rate
Fd = Yfuncvec(4); % kmol/h h2 flow rate
Fe = Yfuncvec(5); % kmol/h co2 flow rate
P = Yfuncvec(6); % Pa
T = Yfuncvec(7); % K Temperature
Ta = Yfuncvec(8); % K Coolant Temperature
Qr = Yfuncvec(9); % W Heat generated
Qf = Yfuncvec(10); % W Heat generated
Tw = Yfuncvec(11); % W Heat generated
global P0 T0 Ft0 rhoc U a Rgas alpha rlist glbp Wspan2 rhob Ac sb hi
% Explicit equations
Ft = Fa + Fb + Fc + Fd + Fe; % kmol/h total molar flowrate
% Heat of reactions
%Shotmate equation parameters for CH4 H2O CO H2 CO2
coeffA = [-0.703029,30.092,25.56759,33.066178,24.99735];
coeffB = [108.4773,6.832514,6.09613,-11.363417,55.18696];
coeffC = [-42.52157,6.793435,4.054656,11.432816,-33.69137];
coeffD = [5.862788,-2.53448,-2.671301,-2.772874,7.948387];
coeffE = [0.678565,0.082139,0.131021,-0.158558,-0.136638];
coeffF = [-76.84376,-250.881,-118.0089,-9.980797,-403.6075];
coeffG = [158.7163,223.3967,227.3665,172.707974,228.2431];
coeffH = [-74.8731,-241.8264,-110.5271,0,-393.5224];

t1 = ((T)/1000); % 1/K 827C or 1100K temperature for mean capacity

Cpbank = coeffA + coeffB.*t1 + coeffC.*t1^2 + coeffD.*t1^3 + coeffE./t1^2; % J/(mol K) heat capacities
Cpa = Cpbank(1); % J/(mol K) heat capacity
Cpb = Cpbank(2); % J/(mol K) heat capacity
Cpc = Cpbank(3); % J/(mol K) heat capacity
Cpd = Cpbank(4); % J/(mol K) heat capacity
Cpe = Cpbank(5); % J/(mol K) heat capacity
dCp1std = 46.4252; % J/mol K basis h2o
dCp2std = 3.22087; % J/mol K basis h2o
dCp3std = 49.6461; % J/mol K basis ch4
dCp1 = Cpc + 3*Cpd - Cpb - Cpa; % J/mol K basis h2o
dCp2 = Cpd + Cpe - Cpc - Cpb; % J/mol K basis h2o
dCp3 = 4*Cpd + Cpe - 2*Cpb - Cpa; % J/mol K basis ch4
Hrx1 = 205700; % J/mol at 25 basis h2o
Hrx2 = -41550; % J/mol at 25 basis h2o
Hrx3 = 164150; % J/mol at 25 basis ch4
dHrx1 = Hrx1 + (dCp1+dCp1std)/2*(T-(298.15)); % J/mol Heat of reaction (reference 25degC)
dHrx2 = Hrx2 + (dCp2+dCp2std)/2*(T-(298.15)); % J/mol Heat of reaction (reference 25degC)
dHrx3 = Hrx3 + (dCp3+dCp3std)/2*(T-(298.15)); % J/mol Heat of reaction (reference 25degC)

% Entropy
sbank = coeffA.*log(t1)+coeffB.*t1+(coeffC.*t1.^2)./2+(coeffD.*t1.^3)./3-coeffE./(2.*t1.^2)+coeffG;
sTa = sbank(1); % J/mol entropy
sTb = sbank(2); % J/mol entropy
sTc = sbank(3); % J/mol entropy
sTd = sbank(4); % J/mol entropy
sTe = sbank(5); % J/mol entropy

% Enthalpy
h298bank = [-74400.000 -241830.000 -110530.000 0.000 -393910.000]; % J/mol enthalpy at standard conditions
hdiffbank = (coeffA.*t1+coeffB.*t1.^2./2+coeffC.*t1.^3./3+coeffD*t1.^4./4-coeffE./t1+coeffF-coeffH)*1000;
hbank = h298bank + hdiffbank; % J/mol enthalpy at standard conditions
hTa = hbank(1); % J/mol enthalpy
hTb = hbank(2); % J/mol enthalpy
hTc = hbank(3); % J/mol enthalpy
hTd = hbank(4); % J/mol enthalpy
hTe = hbank(5); % J/mol enthalpy

Gbank = hbank - T.*sbank; % J/mol Gibbs free
gTa = Gbank(1); % J/mol gibbs free energy
gTb = Gbank(2); % J/mol gibbs free energy
gTc = Gbank(3); % J/mol gibbs free energy
gTd = Gbank(4); % J/mol gibbs free energy
gTe = Gbank(5); % J/mol gibbs free energy

dG1 = 3*Gbank(4)+Gbank(3)-Gbank(2)-Gbank(1);
dG2 = Gbank(5)+Gbank(4)-Gbank(2)-Gbank(3);
dG3 = 4*Gbank(4)+Gbank(5)-2*Gbank(2)-Gbank(1);

% Reaction equilibrium constants
K1 = exp(-dG1/(Rgas.*T));
K2 = exp(-dG2/(Rgas.*T));
K3 = exp(-dG3/(Rgas.*T));

% Reaction rate constants
k1 = glbp(1)*4.22*10^15 * exp(-240100/(Rgas.*T)); % rate coefficients [kmol bar0.5 / kgcat h]
k2 = glbp(2)*1955000 * exp(-67130/(Rgas.*T)); % rate coefficients [kmol / kgcat h]
k3 = glbp(3)*1.02*10^15 * exp(-243900/(Rgas.*T)); % rate coefficients [kmol bar0.5 / kgcat h]



% K1 = exp(30.42-27106./T) % reaction rate constant [bar2]
% K2 = exp(-3.798+4160./T); % reaction rate constant -
% K3 = exp(34.218-31266./T) % reaction rate constant [bar2]
% K3 = K1.*K2;

% Adsorbtion constants
Ka = 6.65*10^-4*exp(38280/(Rgas.*T)); % ch4 bar-1
Kb = 6.12*10^-9*exp(82900/(Rgas.*T)); % h2o - 
Kc = 8.23*10^-5*exp(70650/(Rgas.*T)); % co bar-1
Kd = 1.77*10^5*exp(-88680/(Rgas.*T)); % h2 bar-1
% Ke = 1; % co2 -

% Stoichiometry
% Cto = P0./(Rsi*T0*1000*z); % mol/L initial feed concentration
pch4 = (Fa / Ft) * (P/100000); %/P0) * (T0/T); % bar ch4 concentration bar
ph2o = (Fb / Ft) * (P/100000); %/P0) * (T0/T); % bar h2o concentration bar
pco = (Fc / Ft) * (P/100000); %/P0) * (T0/T); % bar co concentration bar
ph2 = (Fd / Ft) * (P/100000); %/P0) * (T0/T); % bar h2 concentration bar
pco2 = (Fe / Ft) * (P/100000); %/P0) * (T0/T); % bar co2 concentration bar

% Rate laws
den = (1 + Kc*pco + Kd*ph2 + Ka*pch4 + Kb*ph2o/ph2); % denominator term
r1 = ((k1/(ph2^2.5))*(pch4*ph2o-((ph2^3)*pco/K1)))/((den)^2); % first reaction [kmol/kgcat/h]
r2 = ((k2/ph2)*(pco*ph2o-(ph2*pco2/K2)))/((den)^2); % second reaction [kmol/kgcat/h]
r3 = ((k3/(ph2^3.5))*(pch4*(ph2o^2)-((ph2^4)*pco2/K3)))/((den)^2); % third reaction [kmol/kgcat/h]


Wspan2 = [Wspan2 W];
rlist = [rlist [W r1 r2 r3]'];

% Differential equations / Mole Balance
dFadW = -r1-r3; %[kmol/kgcat/h]
dFbdW = -r1-r2-2*r3;
dFcdW = r1-r2;
dFddW = 3*r1+r2+4*r3;
dFedW = r2+r3;
dPdW = glbp(6)*(-alpha/2)*(P0/(P/P0))*(T/T0)*(Ft/Ft0);


% dTdW = ((0.75*U*3600/rhob * a * (Tw - T)) + glbp(5)*(-1*r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);
% dTadW = 0;
% dQrdW = -1*(-r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3);
% dQhdW = (sb*3600/rhob * a * (Ta.^4 - Tw.^4));
% dTwdW =  ((sb*3600/rhob * a * (Ta.^4 - Tw.^4)) - (0.75*U*3600/rhob * a * (Tw - T)))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);


dTdW = ((U*3600/rhob * a * (Tw - T)) + glbp(5)*(-1*r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);
dTadW = 0;
dQrdW = -1*(-r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3);
dQhdW = (U*3600/rhob * a * (Tw - T));
dTwdW = ((sb*3600*a*(Ta.^4-Tw.^4)/(rhob))-(U*3600/rhob * a * (Tw - T)))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);;

% dTdW = ((glb-p(4)*U*3600/rhoc * a * (Ta - T)) + glbp(5)*(-1*r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);
% dTdW = ((-1*glbp(4)*U*3600*a*(T-Ta)/(rhob)) + glbp(5)*(-1*r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3))/(Fa*1000*Cpa + Fb*1000*Cpb + Fc*1000*Cpc + Fd*1000*Cpd + Fe*1000*Cpe);
% dTadW = 0;
% dQrdW = (-r1*1000 * dHrx1 - r2*1000 * dHrx2 - r3*1000 * dHrx3);
% dQhdW = (U*3600/rhob * a * (Ta - T));
dYfuncvecdW = [dFadW; dFbdW; dFcdW; dFddW; dFedW; dPdW; dTdW;dTadW;dQrdW;dQhdW;dTwdW];
% v02: All scenarios

clear; close all;
addpath(genpath('../../Data'))
addpath(genpath('../../../Data'))
%% init
RCP = 2;

if RCP == 1
    emissions_filename = 'SSP119_emissions_202305042205.csv';
    magicc_filename    = 'SSP119_magicc_202305042205.csv';

    e_indx = 164;
    OHC_indx = 3;
    C_indx = 4;
    iE_indx = 5;
    LND_indx = 6;
    OCN_indx = 7;
    Fall_indx = 8;
    FCO2_indx = 13;
    TAS_indx = 18;
elseif RCP == 2
    emissions_filename = 'SSP126_emissions_202305072210.csv';
    magicc_filename    = 'SSP126_magicc_202305072210.csv';

    e_indx = 164;
    OHC_indx = 3;
    C_indx = 4;
    iE_indx = 5;
    LND_indx = 6;
    OCN_indx = 7;
    Fall_indx = 9;
    FCO2_indx = 14;
    TAS_indx = 8;
elseif RCP == 3
    emissions_filename = 'SSP245_emissions_202305072211.csv';
    magicc_filename    = 'SSP245_magicc_202305072211.csv';

    e_indx = 164;
    OHC_indx = 3;
    C_indx = 6;
    iE_indx = 5;
    LND_indx = 16;
    OCN_indx = 17;
    Fall_indx = 7;
    FCO2_indx = 12;
    TAS_indx = 18;
elseif RCP == 4
    emissions_filename = 'SSP370_emissions_202305072212.csv';
    magicc_filename    = 'SSP370_magicc_202305072212.csv';

    e_indx = 164;
    OHC_indx = 15;
    C_indx = 6;
    iE_indx = 5;
    LND_indx = 16;
    OCN_indx = 17;
    Fall_indx = 7;
    FCO2_indx = 12;
    TAS_indx = 18;
elseif RCP == 5
    emissions_filename = 'SSP434-over_emissions_202305072212.csv';
    magicc_filename    = 'SSP434-over_magicc_202305072212.csv';

    e_indx = 164;
    OHC_indx = 15;
    C_indx = 6;
    iE_indx = 5;
    LND_indx = 16;
    OCN_indx = 17;
    Fall_indx = 7;
    FCO2_indx = 12;
    TAS_indx = 18;
elseif RCP == 6
    emissions_filename = 'SSP585_emissions_202305072213.csv';
    magicc_filename    = 'SSP585_magicc_202305072213.csv';

    e_indx = 164;
    OHC_indx = 15;
    C_indx = 6;
    iE_indx = 5;
    LND_indx = 16;
    OCN_indx = 17;
    Fall_indx = 7;
    FCO2_indx = 12;
    TAS_indx = 18;
end
%% Load emissions
e_dat = csvread(emissions_filename,0,10);

t_e = e_dat(1,:)';
E_magicc = e_dat(e_indx,:)'/1e3*12/44;

%% Load magicc output
m_dat = csvread(magicc_filename,0,12);
t_m = m_dat(1,:)';

OHC_magicc = m_dat(OHC_indx,:)';
C_magicc = m_dat(C_indx,:)';
iE_magicc = m_dat(iE_indx,:)';
LND_magicc = m_dat(LND_indx,:)';
OCN_magicc = m_dat(OCN_indx,:)';
Fall_magicc = m_dat(Fall_indx,:)';
FCO2_magicc = m_dat(FCO2_indx,:)';
TAS_magicc = m_dat(TAS_indx,:)';

FnonCO2_magicc = Fall_magicc-FCO2_magicc;

%% plot emissions
figure;
plot(t_e,E_magicc,'b-','LineWidth',1.5), hold on
plot(t_m,iE_magicc,'r.-','LineWidth',1.5), hold on
title('Emissions');
legend('Input','Inverse Emissions','Location','Best');

%% plot sinks
implied_sink = nan(size(OCN_magicc));
implied_sink(end-length(C_magicc)+2:end) = E_magicc(end-length(C_magicc)+2:end) - diff(C_magicc);
figure;
plot(t_m,LND_magicc+OCN_magicc,'b-','LineWidth',1.5), hold on
plot(t_m,implied_sink,'r-.','LineWidth',1.5), hold on
title('Combined sinks');
legend('MAGICC','Implied (E-GATM)','Location','Best');


%% plot rest
figure;
subplot(2,2,1);
plot(t_m,C_magicc,'b-','LineWidth',1.5), hold on
title('Concentrations');

subplot(2,2,2);
plot(t_m,LND_magicc,'b-','LineWidth',1.5), hold on
plot(t_m,OCN_magicc,'r-.','LineWidth',1.5), hold on
title('Sinks');
legend('LAND','OCEAN','Location','Best');

subplot(2,2,3);
plot(t_m,FCO2_magicc,'b-','LineWidth',1.5), hold on
plot(t_m,FnonCO2_magicc,'r-.','LineWidth',1.5), hold on
title('Forcings');
legend('CO2','Non-CO2','Location','Best');

subplot(2,2,4);
plot(t_m,TAS_magicc,'b-','LineWidth',1.5), hold on
title('Surface temperature');

%% save
if RCP == 1
    save('SSP119_output_v2');
elseif RCP == 2
    save('SSP126_output_v2');
elseif RCP == 3
    save('SSP245_output_v2');
elseif RCP == 4
    save('SSP370_output_v2');
elseif RCP == 5
    save('SSP434_output_v2');
elseif RCP == 6
    save('SSP585_output_v2');
end



%% Compare with GCB 2023
dat = csvread('GCB_model_global_2023.csv',1);

disp('Analysing GCB2023 data...')
t        = dat(:,1);
E_FF     = dat(:,2);
E_LUC    = dat(:,3);
G_ATM    = dat(:,4);
S_OCEAN  = dat(:,5);
S_LAND   = dat(:,6);
%S_CEMENT = dat(:,7);
B_IM     = dat(:,7);
GDP      = dat(:,8);
DGDP     = dat(:,9);
GGDP     = dat(:,10);
SOI      = dat(:,11);
C        = dat(:,12);
E        = dat(:,13);

%% plot
LND_magicc2 = LND_magicc(t_m<2023) + E_LUC(t>1994);
tmp = E_magicc(t_e<2023);
tmpt = t_e(t_e<2023);
E_magicc2   = E_magicc(tmpt>1958) + E_LUC;
iE_magicc2 = iE_magicc(t_m<2023) + E_LUC(t>1994);

E_LU_m = e_dat(49,:)'/1e3*12/44;
E_FF_m = e_dat(50,:)'/1e3*12/44;
E_to_m = E_LU_m(t_e>2014)+E_FF_m(t_e>2014);

LND_magicc3 = LND_magicc(t_m>2014) + E_LU_m(t_e>2014);

figure;
subplot(2,2,1);
plot(t,C,'k-','LineWidth',2), hold on
plot(t_m,C_magicc,'b-.','LineWidth',1.5), hold on
title('Concentrations');
legend('GCB','MAGICC','Location','Best')

subplot(2,2,2);
plot(t,E,'k-','LineWidth',2), hold on
plot(t_e,E_magicc,'b-.','LineWidth',1.5), hold on
plot(t,E_magicc2,'r--','LineWidth',1.5), hold on
plot(t_m,iE_magicc,'g--','LineWidth',1.5), hold on
plot(1995:2022,iE_magicc2,'y--','LineWidth',1.5), hold on
plot(2015:2100,E_to_m,'c-','LineWidth',1.5), hold on
title('Emissions');
legend('GCB','MAGICC','MAGICC+LUC','MAGICCinv','MAGICCinv+LUC','Location','Best');

subplot(2,2,3);
plot(t,S_LAND,'k-','LineWidth',2), hold on
plot(t_m,LND_magicc,'b-.','LineWidth',1.5), hold on
plot(1995:2022,LND_magicc2,'r--','LineWidth',1.5), hold on
plot(2015:2100,LND_magicc3,'g--','LineWidth',1.5), hold on
title('Land sink');
legend('GCB','MAGICC','MAGICC+LUC','MAGICC+MAGICC-LUC','Location','Best');

subplot(2,2,4);
plot(t,S_OCEAN,'k-','LineWidth',2), hold on
plot(t_m,OCN_magicc,'b-.','LineWidth',1.5), hold on
title('Ocean sink');
legend('GCB','MAGICC','Location','Best')



%% plot
LND_magicc2 = LND_magicc(t_m<2023) + E_LUC(t>1994);
tmp = E_magicc(t_e<2023);
tmpt = t_e(t_e<2023);
E_magicc2   = E_magicc(tmpt>1958) + E_LUC;
iE_magicc2 = iE_magicc(t_m<2023) + E_LUC(t>1994);

E_LU_m = e_dat(49,:)'/1e3*12/44;
E_FF_m = e_dat(50,:)'/1e3*12/44;
E_to_m = E_LU_m(t_e>2014)+E_FF_m(t_e>2014);

LND_magicc3 = LND_magicc(t_m>2014) + E_LU_m(t_e>2014);

figure;
subplot(2,2,1);
plot(t,C,'k-','LineWidth',2), hold on
plot(t_m,C_magicc,'b-.','LineWidth',1.5), hold on
title('Concentrations');
legend('GCB','MAGICC','Location','Best')

subplot(2,2,2);
plot(t,E,'k-','LineWidth',2), hold on
plot(t_e,E_magicc,'b-.','LineWidth',1.5), hold on
%plot(t,E_magicc2,'r--','LineWidth',1.5), hold on
plot(t_m,iE_magicc,'g--','LineWidth',1.5), hold on
%plot(1995:2022,iE_magicc2,'y--','LineWidth',1.5), hold on
%plot(2015:2100,E_to_m,'c-','LineWidth',1.5), hold on
title('Emissions');
legend('GCB','MAGICC','MAGICC+LUC','MAGICCinv','MAGICCinv+LUC','Location','Best');

subplot(2,2,3);
plot(t,S_LAND,'k-','LineWidth',2), hold on
plot(t_m,LND_magicc,'b-.','LineWidth',1.5), hold on
plot(1995:2022,LND_magicc2,'r--','LineWidth',1.5), hold on
%plot(2015:2100,LND_magicc3,'g--','LineWidth',1.5), hold on
title('Land sink');
legend('GCB','MAGICC','MAGICC+LUC','MAGICC+MAGICC-LUC','Location','Best');

subplot(2,2,4);
plot(t,S_OCEAN,'k-','LineWidth',2), hold on
plot(t_m,OCN_magicc,'b-.','LineWidth',1.5), hold on
title('Ocean sink');
legend('GCB','MAGICC','Location','Best')


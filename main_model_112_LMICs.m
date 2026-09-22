
clear global all
clearvars
global vacc_switch_time  B wm wi1 wi2 u um beta d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv v1_old v2_old v3_old sc1_old sc2_old sc3_old sc2n_old sc3n_old wv_old reintro wA; %immig deaths 


%%  ---- model input data 
load('demographic_data_1980_2060.mat');                      % country-specific demographic data - birth and death rates, population and age-distribution
load('fixed_parameters_for_all_countries_1000_sim.mat');        % fixed model parameters: R0, wm, hosp, wv, d3
load("vaccine_coverage.mat");                                              % country-specific vaccine coverage                                           
load("current_vacc_status_1000_sim.mat")                                   % country-specific vaccine response rate, time of vaccination from 2004 to 2024

%% --- Projected vaccine impact under different strategies - run one strategy at a time ----------------
load("all_keep_current_vacc_status_1000_sim.mat")   % Country-specific vaccine response rates and vaccination timing starting in 2025 under each country's current vaccination strategy.
%load("all_uses_6_10_14_sch_1000_sim.mat")          % Country-specific vaccine response rates and vaccination timing starting in 2025, assuming all countries use a 3-dose schedule at 6, 10, and 14 weeks.
%load("all_switch_to_6_10_sch_1000_sim.mat")        % Country-specific vaccine response rates and vaccination timing starting in 2025, assuming all countries use a 2-dose schedule at 6 and 10 weeks.

% Output the country-specific number of vaccination doses administered.
% For each country, values will remain zero until the year the vaccine is introduced.
annual_dose1_sum = zeros(112,57);     % Annual sum of infants receiving the first dose
annual_dose2_sum = zeros(112,57);     % Annual sum of infants receiving the second dose
annual_dose3_sum = zeros(112,57);     % Annual sum of infants receiving the third dose


% Output country-specific population estimates from the model for 2004–2060
annual_pop_U1 = zeros(112,57);                  % Annual population of children under 1 year
annual_pop_U5 = zeros(112,57);                  % Annual population of children under 5 years
annual_pop_overall = zeros(112,57);             % Annual total population


% - -  results age stratified from 2004 to 2060 - each years is having 20 age groups
% 1     2   3   4   5   6   7        8      9        10      11     12       13       14     15      16      17      18      19     20
%<1Y	1Y	2Y	3Y	4Y	5Y	6_10Y	11_15Y	16_20Y	21_25Y	26_30Y	31_35Y	36_40Y	41_45Y	46_50Y	51_55Y	56_60Y	61_65Y	66_70Y	>70Y
% Store simulation outputs for moderate-to-severe and non-severe cases
MS_stratified = zeros(1000,1140);       % 1,000 simulations of moderate-to-severe cases
NS_stratified = zeros(1000,1140);       % 1,000 simulations of non-severe cases
MS_stratified_mean = zeros(1,1140);     % Mean across 1,000 simulations of moderate-to-severe cases
NS_stratified_mean = zeros(1,1140);     % Mean across 1,000 simulations of non-severe cases
 
 
% Vaccination strategy switching year:
% the year a country without a vaccine introduces one, or a country with an existing vaccination strategy switches to a different schedule.
vacc_switch_time = 540;   % Vaccination strategy switch time, fixed at December 2024 (540 months from January 1980 to December 2024)

 
% --- Total simulation length from January 1980 to December 2060 ---
t0 = round(12*24);              % Burn-in period: January 1980 to December 2003
simulation_part = 684;          % Actual simulation period: January 2004 to December 2060
tmax = t0 + simulation_part;    % Total simulation length: January 1980 to December 2060  Monthly time step (total = 972 months)




%%%% for the global simulations %%%%%%%%%%%%%
%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];                                                      %age distribution 
al=length(age);                                                                     %length age distribution
u=[1/1.0*ones(1,24) 1/12*ones(1,4) 1/(12*5)*ones(1,13) 1/(12*25)];                  %interval between age groups 
datesim=(datenum([1980 1 1 0 0 0]):30.44:datenum([2060 12 31 0 0 0]))';             %length of simulations from 1980 to 2060
datepop=[(1980:1:2060)' 12*ones(81,1) ones(81,1) zeros(81,3)];                      %annaul population 


% Loop over all 112 low- and middle-income countries (LMICs)
for icountry =1:112

    % Loop over 1,000 simulation runs
    for isim =1:1000
        



%country-specific demographics set up
N=pop_pyramid_1980(icountry,:);                                 %country-specific population at 1980 
N=(N/exp(.08));
B=interp1(datenum([1980 1 1 0 0 0; datepop+1]),[birth_rate_1980_2060(1,icountry); birth_rate_1980_2060(:,icountry)]/1000,datesim);    %annual country-specific crude birth rate
B=[B zeros(tmax,al-1)];
cdr=interp1(datenum([1980 1 1 0 0 0; datepop+1]),[death_rate_1980_2060(1,icountry); death_rate_1980_2060(:,icountry)]/1000,datesim);  %annual country-specific crude death rate



% - - - country-specific immigration and rate of aging  for each country- - - - - -
cdr=log(1+cdr)/12;
immig=pop_set(icountry,1)*ones(length(cdr),1); 
immig=log(1+immig)/12;
um=pop_set(icountry,2)*ones(1,al);
um=log(1+um-immig-pop_set(icountry,3))/12;  



%%% INFECTION PARAMETERS  ---  fixed %%%
dur=1/4.3;                          %duration of infectiousness
d1=1/dur;                           %rate of recovery from primary infection (per week)
d2=2*d1;                            %rate of recovery from subsequent infection (per week)
rr1=0.62;                           %relative risk of second infection 
rr2=0.35;                           %relative risk of third infection 
ri2=0.5;                             %relative infectiousness of secondary infection
wi1=1/2.999;                         %rate of waning immunity following primary infection
wi2=1/2.999;                         %rate of waning immunity following 2nd infection
wA=0;                                %rate of waning of immunity to symptomatic infections (S2->S1)
 

ptrans=R0_para(isim,icountry)/dur;     %probability of transmission given contact
wm=wm_para(isim);                      %rate of waning maternal immunity (distribution gamma(6,1))
c1=100*ones(al);                        %Homogeneous mixing 
beta=(ptrans/100)*c1;                   %transmission rate matrix
h=hosp_para(isim);                      %proportion of severe diarrhea cases hospitalized
hosp1=0.13*h*ones(1,al);                %proportion of primary infections with severe diarrhea who are hospitalized -Velazquez et al
hosp2=0.03*h*ones(1,al);                %proportion of secondary infections with severe diarrhea who are hospitalized--Velazquez et al
hosp3=d3_para(isim)*h*ones(1,al);        %proportion of subsequent infections with severe diarrhea who are hospitalized
delta1=0.41*h*ones(1,al);               %proportion of primary infections that are symptomatic
delta2=0.35*h*ones(1,al);                %proportion of secondary infections that are symptomatic
delta3=0.21*h*ones(1,al);                %rate of detection of subsequent infection
ri3=0.1;                                 %relative infectiousness of asymptomatic infection
reintro=0;                               %Allows for possible constant low background risk of infection



%%% VACCINATION PARAMETERS

v1_old=zeros(tmax,al); v2_old=zeros(tmax,al); v3_old=zeros(tmax,al); %initialize vaccination rate across all ages

%% --- Current vaccination strategy ---
% Countries without a vaccination strategy have no vaccination outcomes.
avacc_old = [base_age_dose1_admin(icountry) base_age_dose2_admin(icountry) base_age_dose3_admin(icountry)]; % Country-specific month of vaccination administration for each dose
sc1_old = sc1_base(isim,icountry);      % Country-specific probability of responding to the first dose
sc2_old = sc2_base(isim,icountry);       % Country-specific probability of responding to the second dose
sc3_old = sc3_base(isim,icountry);       % Country-specific probability of responding to the third dose
sc2n_old = sc2n_base(isim,icountry);     % Country-specific probability of responding to the second dose given no response to the first dose
sc3n_old = sc3n_base(isim,icountry);     % Country-specific probability of responding to the third dose given no response to either the first or second dose

%  - - - - current vaccination coverage - countries without vaccine have coverage coevrage of 0------
if avacc_old(1)>0
v1_old(:,avacc_old(1))=base_vac_dose1_scalar(icountry)*vacc_cov_baseline(:,icountry);   %country-specific dose 1 coverage
end
if avacc_old(2)>0
v2_old(:,avacc_old(2))=base_vac_dose2_scalar(icountry)*vacc_cov_baseline(:,icountry);   %country-specific dose 2 coverage
end
if avacc_old(3)>0
v3_old(:,avacc_old(3))=base_vac_dose3_scalar(icountry)*vacc_cov_baseline(:,icountry);    %country-specific dose 3 coverage
end

wv_old=1/wv_para(isim); %vaccine-induced immunity duration - fixed across all countries

 
%% --- projected vaccination strategy from Jan 2025 ---

v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages

avacc=[dose1_month(icountry) dose2_month(icountry) dose3_month(icountry)]; 
sc1=sc1_par(isim,icountry);                 % Country-specific probability of responding to the first dose 
sc2=sc2_par(isim,icountry);                 % Country-specific probability of responding to the second dose
sc3=sc3_par(isim,icountry);                 % Country-specific probability of responding to the third dose
sc2n=sc2n_par(isim,icountry);               % Country-specific probability of responding to the second dose given no response to the first dose
sc3n=sc3n_par(isim,icountry);               % Country-specific probability of responding to the third dose given no response to either the first or second dose


%% ---projected vaccination coverage January 2025 using fixed coverage of either rotavirus or DPT3 coverage as of 2024. --
if avacc(1)>0
v1(:,avacc(1))=dose1_frac(icountry)*[0*vacc_cov(1:vacc_switch_time,icountry); vacc_cov(vacc_switch_time+1:end,icountry)];       %country-specific dose 1 coverage
end
if avacc(2)>0
v2(:,avacc(2))=dose2_frac(icountry)*[0*vacc_cov(1:vacc_switch_time,icountry); vacc_cov(vacc_switch_time+1:end,icountry)];       %country-specific dose 2 coverage
end
if avacc(3)>0
v3(:,avacc(3))=dose3_frac(icountry)*[0*vacc_cov(1:vacc_switch_time,icountry); vacc_cov(vacc_switch_time+1:end,icountry)];       %country-specific dose 3 coverage
end


wv=1/wv_para(isim); %vaccine-induced immunity duration - fixed across all countries


%Initialize vector to keep track of the number of people in each state
St0=zeros(33*al,1);
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-[ones(1,al-11) zeros(1,10)]]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-11) zeros(1,10)]; %Infectious_1 (primary) 
St0(3*al+1:4*al,1)=zeros(1,al); %Recovered_1
St0(4*al+1:5*al,1)=zeros(1,al); %Susceptible_1
St0(5*al+1:6*al,1)=zeros(1,al); %Infectious_2 (2nd time)
St0(6*al+1:7*al,1)=zeros(1,al); %Recovered_2
St0(7*al+1:8*al,1)=zeros(1,al); %Susceptible-Resistant
St0(8*al+1:9*al,1)=zeros(1,al); %Asymptomatic Infectious_3 (subsequent)
St0(9*al+1:10*al,1)=zeros(1,al); %Temp Resistant
St0(10*al+1:11*al,1)=zeros(1,al); %Maternal immunity
St0(11*al+1:12*al,1)=zeros(1,al); %
St0(12*al+1:13*al,1)=zeros(1,al); %
St0(13*al+1:14*al,1)=zeros(1,al); %
St0(14*al+1:15*al,1)=zeros(1,al); %
St0(15*al+1:16*al,1)=zeros(1,al); %SV0
St0(16*al+1:17*al,1)=zeros(1,al); %IV0
St0(17*al+1:18*al,1)=zeros(1,al); %RV0
St0(18*al+1:19*al,1)=zeros(1,al); %SV1
St0(19*al+1:20*al,1)=zeros(1,al); %IV1
St0(20*al+1:21*al,1)=zeros(1,al); %RV1
St0(21*al+1:22*al,1)=zeros(1,al); %SV2
St0(22*al+1:23*al,1)=zeros(1,al); %IV2
St0(23*al+1:24*al,1)=zeros(1,al); %RV2
St0(24*al+1:25*al,1)=zeros(1,al); %MV1
St0(25*al+1:26*al,1)=zeros(1,al); %MV2
St0(26*al+1:27*al,1)=zeros(1,al); %MV3
St0(27*al+1:28*al,1)=zeros(1,al); %MV4
St0(28*al+1:29*al,1)=zeros(1,al); %MV5
St0(29*al+1:30*al,1)=zeros(1,al); %MV6
St0(30*al+1:31*al,1)=zeros(1,al);  %V1
St0(31*al+1:32*al,1)=zeros(1,al);  %V2
St0(32*al+1:33*al,1)=zeros(1,al);  %V3

clear St lambda H %clear outcome variables which may be in memory

options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time St]=ode45('ode_equations',1:tmax,St0,options);



%Discard the burn-in period
time(1:t0,:)=[]; %delete output from from burn-in period
St(1:t0,:)=[]; 

lambda=zeros(tmax-t0,al);
for t=1:size(St,1)
     lambda(t,:)=((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); 
     
end

% model-estimated moderate-to-severe and non-severe cases 
Cu=zeros(tmax-t0,al); Hu=zeros(tmax-t0,al); Cv=zeros(tmax-t0,al); Hv=zeros(tmax-t0,al); Incid=zeros(tmax-t0,al); Prev=zeros(tmax-t0,al);
for i=1:al
    Cu(:,i)=max(0,delta1(i)*St(:,al+i).*lambda(:,i)+delta2(i)*rr1*St(:,4*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,7*al+i).*lambda(:,i));              %unvaccinated non-severe cases
    Cv(:,i)=max(0,delta1(i)*St(:,15*al+i).*lambda(:,i)+delta2(i)*rr1*St(:,18*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,21*al+i).*lambda(:,i));          %vaccinated non-severe cases
    Hu(:,i)=max(0,hosp1(:,i).*St(:,al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,4*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,7*al+i).*lambda(:,i));       %unvaccinated moderate-to-severe cases
    Hv(:,i)=max(0,hosp1(:,i).*St(:,15*al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,18*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,21*al+i).*lambda(:,i));  %vaccinated moderate-to-severe cases
end


H=Hu+Hv;        %combined age stratified moderate-to-severe
C=Cu+Cv;        %combined age stratified non-severe cases


%number of doses adminstered
annual_dose1_sum(icountry,:)=sum(reshape(sum(St(:,630+base_age_dose1_admin(icountry):42:1386),2),12,57));     %sum vaccination age group at 2month
annual_dose2_sum(icountry,:)=sum(reshape(sum(St(:,630+base_age_dose2_admin(icountry):42:1386),2),12,57));     %sum vaccination age group at 3month
annual_dose3_sum(icountry,:)=sum(reshape(sum(St(:,630+base_age_dose3_admin(icountry):42:1386),2),12,57));     %sum vaccination age group at 4month

%% results from population  2004 to 2060
annual_pop_U1(icountry,:)=mean(reshape(sum(St(:,[1:12,43:54,85:96,127:138,169:180,211:222,253:264,295:306,337:348,379:390,421:432,463:474,505:516,547:558,589:600,631:642,673:684,715:726,757:768,799:810,841:852,883:894,925:936,967:978,1009:1020,1051:1062,1093:1104,1135:1146,1177:1188,1219:1230,1261:1272,1303:1314,1345:1356]),2),12,57));
annual_pop_overall(icountry,:)=mean(reshape(sum(St(:,[1:42, 43:84, 85:126, 127:168, 169:210, 211:252, 253:294, 295:336, 337:378, 379:420, 421:462, 463:504, 505:546, 547:588, 589:630, 631:672, 673:714, 715:756, 757:798, 799:840, 841:882, 883:924, 925:966, 967:1008, 1009:1050, 1051:1092, 1093:1134, 1135:1176, 1177:1218, 1219:1260, 1261:1302, 1303:1344, 1345:1386]),2),12,57));
annual_pop_U5(icountry,:)=mean(reshape(sum(St(:,[1:27, 43:69, 85:111, 127:153, 169:195, 211:237, 253:279, 295:321, 337:363, 379:405, 421:447, 463:489, 505:531, 547:573, 589:615, 631:657, 673:699, 715:741, 757:783, 799:825, 841:867, 883:909, 925:951, 967:993, 1009:1035, 1051:1077, 1093:1119, 1135:1161, 1177:1203, 1219:1245, 1261:1287, 1303:1329, 1345:1371]),2),12,57));





%%% *********** cases (moderate-to-severe and non severe cases)  from 2004 to 2060 ************************

mod_to_severe_2004_2060=([sum(H(:,1:12),2) sum(H(:,13:24),2) sum(H(:,25),2) sum(H(:,26),2) sum(H(:,27),2) sum(H(:,28),2) sum(H(:,29),2) sum(H(:,30),2) sum(H(:,31),2) ...
sum(H(:,32),2) sum(H(:,33),2) sum(H(:,34),2) sum(H(:,35),2) sum(H(:,36),2) sum(H(:,37),2) sum(H(:,38),2) sum(H(:,39),2) sum(H(:,40),2) sum(H(:,41),2) sum(H(:,42),2)]);

non_severe_2004_2060=([sum(C(:,1:12),2) sum(C(:,13:24),2) sum(C(:,25),2) sum(C(:,26),2) sum(C(:,27),2) sum(C(:,28),2) sum(C(:,29),2) sum(C(:,30),2) sum(C(:,31),2) ...
sum(C(:,32),2) sum(C(:,33),2) sum(C(:,34),2) sum(C(:,35),2) sum(C(:,36),2) sum(C(:,37),2) sum(C(:,38),2) sum(C(:,39),2) sum(C(:,40),2) sum(C(:,41),2) sum(C(:,42),2)]);



%  - -  - combining the output for age stratefied moderate-to-severe and non-severe cases   - - - - - - -
%age stratified moderate-to-severe cases between 2004 and 2060
[nrows,ncols]=size(mod_to_severe_2004_2060(:,:));
rowscount=12;
temp_MS=reshape(mod_to_severe_2004_2060(:,:),rowscount,1,[]);
temp_MS=sum(temp_MS,1,'omitnan')/h;
temp_MS=reshape(temp_MS,[],ncols);
tempval_MS=temp_MS';
MS_stratified(isim,:)=tempval_MS(:)';  
MS_stratified_mean(1,:)=mean(MS_stratified,1);                  %mean of 1000 simulations


%age stratified non-severe cases between 2004 and 2060
[mrows,mcols]=size(non_severe_2004_2060(:,:));
rowscount=12;
temp_NS=reshape(non_severe_2004_2060(:,:),rowscount,1,[]);
temp_NS=sum(temp_NS,1,'omitnan')/h;
temp_NS=reshape(temp_NS,[],ncols);
tempval_NS=temp_NS';
NS_stratified(isim,:)=tempval_NS(:)';
NS_stratified_mean(1,:)=mean(NS_stratified,1);                  %mean of 1000 simulations


end

save(sprintf('results_curr_vac_novac_zero_from_2004_%02d',icountry),'MS_stratified','NS_stratified','MS_stratified_mean','NS_stratified_mean','annual_dose1_sum','annual_dose2_sum','annual_dose3_sum','annual_pop_U1','annual_pop_U5','annual_pop_overall');


end

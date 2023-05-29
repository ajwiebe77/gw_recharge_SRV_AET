% rechargeTotal30.m
%
% Andrew J. Wiebe, 11 Jan 2023; modifies rechargeTotal29.m, 29 Jul 2022; modifies rechargeTotal29.m, 07 Dec 2020, modifies rechargeTotal26.m, 24 Nov 2020; modifies rechargeTotal24.m, 23 Nov 2020; modifies rechargeTotal23.m, 5 Nov 2020; 21 Oct 2020; modifies: rechargeTotal22.m, 4 Sep 2020, 21 Jul 2020; modifies rechargeTotal21test.m, 17 Jul 2020; modifies rechargeTotal21.m, 9 Jul 2020
%
% Code for GNU Octave (Eaton et al., 2018).
%
% Overall Objective: Read in random rainfall time series, generate random AET/P estimates,
%                    and calculate the recharge for the Alder Creek watershed 
%                    for each realization based on a stochastic annual vadose zone water budget.
%
% Method:
% Use annual streamflow and BFI values derived from the Water Survey of Canada gauge within the watershed (WSC, 2019).
% Assume uniform snowfall across the watershed and use the annual values from the Environment Canada 
%      weather station at Roseville (Government of Canada, 2019)
% Calculate the ratio of AET / P from the Budyko curve, with or without random perturbation within 
%      observed MOPEX envelope (Budyko, 1961; Gentine et al., 2012)
% 
% Note: 11 Jan 2023: replaced all but last three uses of nrealz with maxrealz; set nrealz to convergence number of iterations; 
%       now has the calculations stop when convergence is reached. Need to recalculate the precipitation histogram based on the 
%       number of realizations employed
% 
% This script calculates: 
%            1) the Thiessen polygon, weighted average total precipitation from three virtual stations
%               (after adding Roseville snowfall to the annual rain sums at the virtual stations),
%            2) an estimate of the AET/P ratio based on the Budyko curve and observed variability in
%                47 MOPEX watersheds (Duan et al., 2006) where 0.787 <= PET/P <= 0887, and
%            3) an estimate of total recharge to the watershed over n = 46 years.
%
% References:
% Barlow, P.M., Cunningham, W.L., Zhai, T., and Gray, M., 2015. U.S. Geological Survey Groundwater Toolbox, 
%    a graphical and mapping interface for analysis of hydrologic data (version 1.0) – User guide for 
%    estimation of base flow, runoff, and groundwater recharge from streamflow data. U.S. Geological Survey
%    Techniques and Methods, book 3, chap. B10, 27 p. https://doi.org/10.3133/tm3B10.
% Budyko, M., 1961. The heat balance of the Earth’s surface. Natl. Weather Serv., U.S. Dep. of Commer.,
%    Washington, D.C., USA.
% Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%    high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
% Gentine, P., D’Odorico, P., Lintner, B.R., Sivandran, G., Salvucci, G., 2012. Interdependence of climate,
%    soil, and vegetation as constrained by the Budyko curve. Geophys. Res. Lett. 39, L19404.
%    https://doi.org/10.1029/2012GL053492.
% Government of Canada, 2019. Historical data: Rainfall, snowfall, and air temperature data for the Roseville,
%    ON, weather station data [computer files]. 
%    http://climate.weather.gc.ca/historical_data/search_historic_data_e.html.
% Matrix Solutions Inc. (Matrix), S.S. Papadopulos and Associates Inc. (SSPA), 2014b. Region of Waterloo
%    Tier Three Water Budget and Local Area Risk Assessment. Final Report, Sep. 2014. 
%    Prepared for: Region of Waterloo. 
%    https://www.sourcewater.ca/en/source-protection-areas/resources/Documents/Grand/RMOW-September-2014-WQRA_chpt-1-10.pdf.
% Neff, B.P., Day, S.M., Piggott, A.R., Fuller, L.M., 2005. Base flow in the Great Lakes basin.
%    Date Posted: 29 Nov 2005. U.S. Geological Survey Sci. Inv. Rep. 2005-5217.
%    https://pubs.usgs.gov/sir/2005/5217/pdf/SIR2005-5217.pdf.
% Roberson, W., 2017. remove rows with all zeros. MATLAB Answers.
%    https://www.mathworks.com/matlabcentral/answers/40390-remove-rows-with-all-zeros.
% Water Survey of Canada (WSC), 2019. Daily volumetric flow rates for the New Dundee gauging 
%    station (02GA030) for 1965 to 2018. https://wateroffice.ec.gc.ca/mainmenu/historical_data_index_e.html.
%

tic; % 11 Jan 2023

clear all; clc;

pkg load statistics;

window = 500; % window size for checking percentage change in metrics with a set of realizations
converge = 0.1; % convergence threshold in percent (0.1% = 0.001)
convergence = 0;

aet_option = 3; % 1 -- budyko curve directly; 2 -- budykoRNDnormal (normally distributed about Budyko curve); 3 -- empirical (MOPEX5)

%%% FORMAT: statsRnd = [sum1973, sum1974, ..., sum2018]
statsRnd = dlmread("C:\\Octave640working\\topic3journal\\Rnd_rainfall27\\stats_Summary115373compiled.txt", " ,/:\t", 5, 1); % the zeroth column is all zeros
statsRnd = statsRnd(any(statsRnd,2) != 0,:); % remove zero timestamp rows (Roberson, 2017)

%%% assign contains a list of rows in statsRnd that relate to each of the five virtual rain gauges
assign = dlmread("C:\\Octave640working\\topic3journal\\assignments_VirtualStations27v3compiled.txt", " ,/:\t", 5, 1); % * the zeroth column is all zeros
assign = assign(any(assign,2) != 0,:); % remove zero timestamp rows (Roberson, 2017)

%%% use ETo calculations with average wind speed = 1.6 m/s, from the SOWC stations
% rville_ETo = dlmread("E:Octave-5.1.0.0\\topic3journal\\Roseville_ETo_ws1_6ms.txt", " ,/:\t", 3, 0); % format: [year ETo_mm]
load Roseville_ETo_ws1_6ms.mat; % format: [year ETo_mm]
pet_annual = rville_ETo([24:end],2);  % select only the data for years 1973 to 2018
pet_avg = sum(pet_annual) / length(pet_annual);

load Roseville_rain_snow_totPrecip1973-2018infilled.mat; % format: [day month year rainfall_mm snowfall_mm totalPrecipitation_mm]
for yr = 1973:2018
   idx2 = rville_precip(:,3) == yr;
   rville_snow(yr - 1972,1) = sum(rville_precip(idx2,5));
end

% load sigma_bracket5.txt;
% load sigma_bracket8.txt; % fewer brackets; each one has at least 100 points
load sigma_bracket10.txt;

bfi = 0.56; % based on analysis with PART in USGS Groundwater Toolbox (Barlow et al., 2015)
bfi_results = dlmread("C:\\Octave640working\\topic3journal\\Alder_PART_BFI_results.txt", " ,/:\t", 3, 0);
bfi_annual = bfi_results(:,2);
Q_tot_avg = 236; % revised the 216 mm from Wiebe (2020) Appendix Q
Qscalingfactor = Q_tot_avg / 140.5;
Q_sw_fraction = (1- bfi) * Q_tot_avg; 
Q_tot_annual_results = dlmread("C:\\Octave640working\\topic3journal\\WSC_Qannual.txt", " ,/:\t", 2, 0);
Q_tot_annual = Q_tot_annual_results(:,2);

n = 46;

maxrealz = length(assign);

%%% Exclude AET from the saturated zone via a correction factor derived from applying PET to areas mapped as bog/swamp/wetland/marsh/open water
vzAET_factor = 0.93; % 0.89; % 0.93; % 0.91

recharge = zeros(maxrealz,8); %[row1 row2 row3 row4 row5 p_ws_avg, aet_p, recharge]

percent_avgrch_change = zeros(maxrealz,1);

convStats = zeros(maxrealz, 4); % format = [avg_recharge_avg_i, percent_diff_wrt_avg_i, percent_diff_mass_loadings_change(i,1), factor_change(i,1)];

q_K26_well = 6756 * 365; % pumping volume in m3

for i = 1:maxrealz
   if(convergence == 0)
	   p_list = zeros(n, 5);

	   %%% Read each row; pick the first five valid time series of the six
	   temp = assign(i,:)';
	   temp = temp(any(temp,2) != 0,:); % remove zero timestamp rows (Roberson, 2017)
	   
	   p_list(:,1) = statsRnd(temp(1,1), :)' + rville_snow; % v26: removed the PET/P ratios from statsRnd
	   % p_list(:,2) = statsRnd(temp(2,1), [3:48])' + rville_snow;
	   p_list(:,3) = statsRnd(temp(3,1), :)' + rville_snow;
	   % p_list(:,4) = statsRnd(temp(4,1), [3:48])' + rville_snow;
	   p_list(:,5) = statsRnd(temp(5,1), :)' + rville_snow;
	   
	   p_ws_avg_annual = calcWSavgPrecip357annual(p_list(:,[1,3,5])); % 3 Thiessen polygons; will now be a list of 46 amounts
	   p_ws_avg = sum(p_ws_avg_annual) / n;
	   
	   if aet_option == 1
	   
		  aet_p = vzAET_factor * budyko(pet_avg, p_ws_avg);
		  aet_p_annual = ones(n,1) * aet_p;
		  
	   elseif aet_option == 2
		  
		  for j = 1:n
			 %%% normal scatter about the Budyko Curve
			 aet_p_annual(j,1) = vzAET_factor * budykoRNDnormalMOPEX2(pet_annual(j,1), p_ws_avg_annual(j,1));
			 
		  end
		  
		  % filename = sprintf("Rnd_aet_p_30\\aet_p_annual_aetopt2_%d.txt", i);
		  % save(filename, "aet_p_annual");
		  
		  aet_p = sum(aet_p_annual) / n;

	   elseif aet_option == 3
		  
		  for j = 1:n
			 aet_p_annual(j,1) = vzAET_factor * budykoRNDnormalMOPEX5(pet_annual(j,1), p_ws_avg_annual(j,1), sigma_bracket(:,[1:4]));
		  end
		  
		  % filename = sprintf("Rnd_aet_p_30\\aet_p_annual_aetopt3_%d.txt", i);
		  % save(filename, "aet_p_annual");
		  
		  aet_p = sum(aet_p_annual) / n;
		  
	   end
	   
	   recharge(i,[1:5]) = temp([1:5],1)';
	   recharge(i,6) = p_ws_avg;
	   recharge(i,7) = aet_p;

	   recharge_star_annual = p_ws_avg_annual - (p_ws_avg_annual .* aet_p_annual) - Qscalingfactor * (1 - bfi_annual) .* Q_tot_annual;
	   recharge(i,8) = sum(recharge_star_annual);
	   
	   if i >= 2 * window  % quantify change in metrics (avg change in value with additional 'window' realizations)
		  
		  recharge_avg_list_i = recharge([1:i],8) / n;
		  recharge_avg_list_prev = recharge([1:(i-window)],8) / n;

		  quan0025_avg_recharge_i = quantile(recharge_avg_list_i,[0.025]);
		  quan0975_avg_recharge_i = quantile(recharge_avg_list_i,[0.975]);
		  quan0025_avg_recharge_prev = quantile(recharge_avg_list_prev,[0.025]);
		  quan0975_avg_recharge_prev = quantile(recharge_avg_list_prev,[0.975]);
		  
		  avg_recharge_avg_i = sum(recharge_avg_list_i) / i;
		  avg_recharge_avg_prev = sum(recharge_avg_list_prev) / (i-window);
		  
		  %metric = change when additional 'window' averages considered
		  percent_avgrch_change(i,1) = 100*(avg_recharge_avg_prev - avg_recharge_avg_i) / avg_recharge_avg_prev;
		  
		  % metric = cap zone area change
		  area1_K26 = q_K26_well / (quan0025_avg_recharge_i / 1000); % area in m2
		  area2_K26 = q_K26_well / (quan0975_avg_recharge_i / 1000); % area in m2
		  areaK26avg = q_K26_well / (avg_recharge_avg_i / 1000); % area in m2
		  areaDiff_K26 = max(abs(area1_K26 - areaK26avg), abs(areaK26avg - area2_K26)) / 1e6;   % difference in sq. km
		  areaK26_avg = areaK26avg / 1e6;
		  percent_diff_wrt_avg_i = 100 * areaDiff_K26 / areaK26_avg;
		  
		  area1_K26 = q_K26_well / (quan0025_avg_recharge_prev / 1000); % area in m2
		  area2_K26 = q_K26_well / (quan0975_avg_recharge_prev / 1000); % area in m2
		  areaK26avg = q_K26_well / (avg_recharge_avg_prev / 1000); % area in m2
		  areaDiff_K26 = max(abs(area1_K26 - areaK26avg), abs(areaK26avg - area2_K26)) / 1e6;   % difference in sq. km
		  areaK26_avg = areaK26avg / 1e6;
		  percent_diff_wrt_avg_prev = 100 * areaDiff_K26 / areaK26_avg;
		  
		  percent_capzoneareachange_change(i,1) = 100*(percent_diff_wrt_avg_prev - percent_diff_wrt_avg_i)/percent_diff_wrt_avg_prev;
		  
		  % metric = mass loadings change	  
		  percent_diff_mass_loadings_i = 100 * max(abs(quan0975_avg_recharge_i - avg_recharge_avg_i), abs(quan0025_avg_recharge_i - avg_recharge_avg_i)) / avg_recharge_avg_i;
		  percent_diff_mass_loadings_prev = 100 * max(abs(quan0975_avg_recharge_prev - avg_recharge_avg_prev), abs(quan0025_avg_recharge_prev - avg_recharge_avg_prev)) / avg_recharge_avg_prev;
		  percent_diff_mass_loadings_change(i,1) = 100 * (percent_diff_mass_loadings_prev - percent_diff_mass_loadings_i)/percent_diff_mass_loadings_prev;
		  
		  % metric = max cumulative recharge / min cumulative recharge
		  factor_i = quantile(recharge([1:i],8),[0.975]) / quantile(recharge([1:i],8),[0.025]);
		  factor_prev = quantile(recharge([1:(i-window)],8),[0.975]) / quantile(recharge([1:(i-window)],8),[0.025]);
		  factor_change(i,1) = 100*(factor_prev - factor_i)/factor_prev;	  
		  
		  convStats(i,:) = [avg_recharge_avg_i, percent_diff_wrt_avg_i, percent_diff_mass_loadings_i, factor_i];
		  
	   end
	   
	   if or(mod(i,500) == 0, i == maxrealz)
		  i
	   end
	   
	   % if and(i >= 2 * window, i > 5000, not(convergence)) % use this regularly
	   % if and(i >= 2 * window, not(convergence)) % test with this (only ~ 1000 realizations)
	   % if i >= 2 * window % test with this (only ~ 1000 realizations)
	   if i >= 4 * window
		  %%% determine whether average change in each metric is less than convergence criterion
		  conv1 = sum(abs(percent_avgrch_change([i-window+1]:i,1)))/window < converge;
		  conv2 = sum(abs(percent_capzoneareachange_change([i-window+1]:i,1)))/window < converge; % 
		  conv3 = sum(abs(percent_diff_mass_loadings_change([i-window+1]:i,1)))/window < converge; %
		  conv4 = sum(abs(factor_change([i-window+1]:i,1)))/window < converge;
		  
		  if and(conv1,conv2,conv3,conv4)
			 printf("All metrics < %f%% over %d iterations after %d iterations.\n",converge,window,i);
			 convergence = i;
			 nrealz = convergence;
		  end
		  
		  % if and(or(mod(i,500) == 0, i == maxrealz),conv1,conv2,conv3,conv4)
			 % printf("All metrics < %f%% over %d iterations at %d iterations.\n",converge,window,i);
		  % end
		  
	   end
	end
end

% recharge = recharge(any(recharge,2) != 0,:); % remove zero timestamp rows (Roberson, 2017)
recharge = recharge(1:nrealz,:); % trim the matrix; ignore zeros after convergence
convStats = convStats(1:nrealz,:); % trim the matrix; ignore zeros after convergence

%%% --- RESULTS ---------------------------------------------------------------------

vzAET_factor100 = 100 * vzAET_factor;

if aet_option == 1
   filename = sprintf("rechargeTotal30_357_results_no_scatter_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   diary(filename);
elseif aet_option == 2
   filename = sprintf("rechargeTotal30_357_results_scatter_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   diary(filename);
elseif aet_option == 3
   filename = sprintf("rechargeTotal30_357_results_mopex_empirical_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   diary(filename);
end

printf("\n");
if aet_option == 3 
   printf("AET/P estimates via budykoRNDnormalMOPEX5\n");
end

maximum_possible_realizations = maxrealz
nrealz
aet_option
window
convergePercent = converge
convcheck = [conv1 conv2 conv3 conv4]
if and(conv1,conv2,conv3,conv4)
   printf("All metrics < %f%% after %d iterations.\n",converge, convergence);
   convergence = i;
end

bfi
Q_tot_avg

p_ws_avg_overall = sum(recharge(:,6)) / nrealz
stdev_p_ws_avg = std(recharge(:,6));

aet_p_avg = sum(recharge(:,7)) / nrealz
stdev_aet_p = std(recharge(:,7))

recharge_avg = recharge(:,8) / n;
avg_recharge_avg = sum(recharge_avg) / nrealz
stdev = std(recharge_avg)

min_avg_recharge = min(recharge_avg)
quan0025_avg_recharge = quantile(recharge_avg,[0.025])
quan0975_avg_recharge = quantile(recharge_avg,[0.975])
max_avg_recharge = max(recharge_avg)

%%% --- Mass Loadings ---

printf("\n");
worst_case_percent_diff_mass_loadings_wrt_avg = 100 * max(abs(quan0975_avg_recharge - avg_recharge_avg), abs(quan0025_avg_recharge - avg_recharge_avg)) / avg_recharge_avg
printf("\n");

%%% --- Well Capture Zone ---

area1_K26 = q_K26_well / (quan0025_avg_recharge / 1000); % area in m2
area2_K26 = q_K26_well / (quan0975_avg_recharge / 1000); % area in m2
areaK26avg = q_K26_well / (avg_recharge_avg / 1000); % area in m2
areaDiff_K26 = max(area1_K26 - areaK26avg, areaK26avg - area2_K26) / 1e6   % difference in sq. km

areaK26_avg = areaK26avg / 1e6

printf("\n");
percent_diff_wrt_avg = 100 * areaDiff_K26 / areaK26_avg
printf("\n");

%%% --- Cumulative Recharge ---

dif_cum_tot = quantile(recharge(:,8),[0.975]) - quantile(recharge(:,8),[0.025]); % difference in mm
dif_cum_tot_avg = dif_cum_tot / n  % average difference per year in mm
printf("\n");

factor_cum_quan0975_div_by_cum_quan0025 = quantile(recharge(:,8),[0.975]) / quantile(recharge(:,8),[0.025])
printf("\n");
diary off;

if aet_option == 1
   filename = sprintf("rechargeTotal30_357_budyko_no_scatter_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   save(filename, "recharge");
elseif aet_option == 2
   filename = sprintf("rechargeTotal30_357_budykoMOPEX2_scatter_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   save(filename, "recharge");
elseif aet_option == 3
   filename = sprintf("rechargeTotal30_357_budykoMOPEX5_mopex_empirical_aetvzfactor0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
   save(filename, "recharge");
end

filename = sprintf("convergenceStats30_vzaet_factor_0_%d_aetopt%d.txt", vzAET_factor100, aet_option);
save(filename, "convStats");

toc; % 11 Jan 2023
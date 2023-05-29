% gen_alderRain27.m
% 
% Andrew James Wiebe, 30 Nov 2020, 26 Nov 2020; modifies gen_alderRain26.m, 26 Nov 2020, 24 Nov 2020; modifies gen_alderRain25.m, 23 Nov 2020; gen_alderRain24.m, 20 Nov 2020; modifies gen_alderRain23.m, 19 Oct 2020; modifies gen_alderRain22.m, 5 Aug 2020, 21 Jul 2020; modifies gen_alderRain21test.m, 15 Jul 2020; modifies gen_alderRain20.m, 22 Jan 2020
%
% Code for GNU Octave (Eaton et al., 2018).
%
% Objective: Generate synthetic 46-yr time series based on a parametric method that uses
%            Markov transition probabilities for intervals between rainy days and 
%            a mixed exponential function for rainfall amounts.
% Overview: Use Markov probabilities for spacing and mixed exponential function to generate rainfall time series. 
%           Check sets of seven (six time-series plus one observed) time series for acceptable spatial correlation
%           (based on observed Pearson Product-Moment and Spearman Rank Correlation coefficients). 
%           Save the synthetic time series from sets of six that have acceptable spatial correlation, 
%           and a summary file containing a list of assignments of rainfall time series to five virtual stations.
%
% Method:
% Step 1: Generate six random time series based on Markov transition probabilities (Basinger et al., 2010) 
%         for rainy day spacing and a mixed exponential model (Li et al., 2013) for daily amounts.
% Step 2: Use the Iman and Conover (1982) method (Tarpanelli et al., 2012) (the "IC method") to 
%         modify correlation coefficients among the seven time series (the six random time series plus
%         the observed Roseville rainfall) to be similar to the observed Pearson and Spearman coefficients.
%         Reject any sets of time series with correlation coefficients falling outside the observed Spearman
%         Rank Correlation coefficient range.
% Step 3: Save the data for sets of five or six time series that have acceptable spatial correlation 
%         coefficients. A summary file called stats_SummaryX.txt (where X is the total number of valid time series 
%         generated) is saved; this file contains total precipitation sums for each year (1973 to 2018) 
%         for each time series (information for one time series per row). Finally, a list of all sets of 
%         valid time series is saved in "assignments_VirtualStations26vX.txt" (where X = version). Each row 
%         of this file lists a valid set of five or six time series. A zero is an empty placeholder if there are
%         only five valid time series in the set.
%
% Notes:
% 1. Time series that correlate similarly with Roseville compared to the SOWC field data
%    are adjusted in blocks of 3 (or 4 for the first set) years at a time, because the observed 
%    field dataset was 3 years in length.
% 2. The outer "for i" loop generates 5 or 6 random time series during each successful iteration.
% 3. If this file is used more than once to generate additional time series, the "c" and "assignRow"
%    and "version" parameters should be incremented to the next numbers after the respective numbers
%    already processed. Also, the user needs to splice together (separately) the stats_Summary*.txt files 
%    and the assignments_virtualStations*.txt files before running rechargeTotal27.m.
% 4. A folder called "Rnd_rainfall27" needs to be created in the folder containing this file prior 
%    to running this script.
% 5. Difference between this version and gen_alderRain20.m: Removed the constraint that the PET/P ratio at each 
%    station needed to fall within observed range of PET/P ratios +/-10% at Roseville. Those ratios are
%    constrained by the observed precipitation at Roseville, which is specific to Roseville.
% 6. The mixed exponential function parameters for the Roseville rainfall depths were updated for this version (lm_mixed7depth3.m).
% 7. The list variable is still needed to ensure that only the acceptably correlated time series numbers are recorded.
%
% References:
% Basinger, M., Montalto, F., Lall, U., 2010. A rainwater harvesting system reliability model based on 
%    nonparametric stochastic rainfall generator. J. Hydrol. 392, 105-118.
%    https://doi.org/10.1016/j.jhydrol.2010.07.039.
% Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%    high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
% Iman, R.L., Conover, W.J., 1982. A distribution-free approach to inducing rank correlation among input 
%    variables. Commun. Stat. Simul Comput. 11(3), 311-334. https://doi.org/10.1080/03610918208812265.
% Roberson, W., 2017. remove rows with all zeros. MATLAB Answers.
%    https://www.mathworks.com/matlabcentral/answers/40390-remove-rows-with-all-zeros.
% Tarpanelli, A., Franchini, M., Brocca, L., Camici, S., Melone, F., Moramarco, T., 2012. A simple approach
%    for stochastic generation of spatial rainfall patterns. J. Hydrol 472-473, 63-76.
%    https://doi.org/10.1016/j.jhydrol.2012.09.010.

clear all;

% apply the observed range of Spearman daily spatial correlation coefficients (among WS2, WS3, ..., WS7, Roseville)
max_correl = 0.825; % updated - omits missing data for WS5, WS2, and WS6 (rainCorrel14v6.m)
min_correl = 0.5196; 

% 3-year block final index values for 1976, 1979, ..., 2015, 2018 (first block is 4 years)
breaks = [1461, 2556, 3652, 4748, 5844, 6939, 8035, 9131, 10227, 11322, 12418, 13514, 14610, 15705, 16801]';

% rville_precip = dlmread("E:\\Octave-5.1.0.0\\topic3journal\\Roseville_rain_snow_totPrecip1973-2018infilled.txt", " ,/:\t", 5, 1);
% rville_precip = rville_precip(any(rville_precip,2) != 0,:); % remove zero timestamp rows (Roberson, 2017)
load Roseville_rain_snow_totPrecip1973-2018infilled.mat; % format: [day month year rainfall_mm snowfall_mm totalPrecipitation_mm]

load Roseville_ETo_ws1_6ms.mat; % format: [year ETo_mm]
rville_ETo = rville_ETo([24:end],:); % select only the data for years 1973 to 2018

load sowc_minus1_Rville_24hr_correl_pearson_spearman_rev5.mat; % observed Pearson and Spearman correlation (7x7) matrices for WS2,WS3,WS4,WS5,WS6,WS7,Roseville; omits missing data for WS5, WS2, and WS6 (rainCorrel14v5.m)

load markov_probabilities_Roseville_infill.mat pw pww pdd; % load overall rainfall probability (pw) and transition probabilities
pww(366,1) = pww(365,1); % Set Julian Day 366 values to Julian Day 365 values for consistency (too little data for leap years)
pdd(366,1) = pdd(365,1); % Set Julian Day 366 values to Julian Day 365 values for consistency (too little data for leap years)

assign = zeros(10000,6);

n = 10000;
ny = 46;

times = rville_precip(:,[1:3]);
maxdaily = 100; % maximum daily rainfall depth allowed at Roseville (observed maximum is 97 mm)

statsRnd = zeros(n, 46); %statsRnd = [sum1973, sum1974, ..., sum2018]

ndays = length(rville_precip);
load mixed7depth3_pbeta1beta2.txt;
p_fit_depths = p_fit';

version = 1; % 3; %2; 1; % ********** update by one each time this is run *********
c = 1; % 83431; % 18978; %1; % ********** update with next value before re-running ***************
assignRow = 1; %14744; % 3360; %1; % ********** update with next value before re-running ***************

for i = 1:n
	
	if mod(i,100) == 0
       i
    end	   
	
	% select = 1;
	
	rain_yr = zeros(6,ny);
	
	%%% --- Step 1: generate a set of six rainfall time series
	rnd_rain = zeros(ndays,7);
	
	for stn = 1:6
	   % rnd_rain(:,stn) = rndMixedExpRain(ndays, p_fit_gaps, p_fit_depths); % average rainfall is too high
	   rnd_rain(:,stn) = markovGaps_mixedExpRain(times, p_fit_depths, pw, pww, pdd, maxdaily);
	end
	
	%%% --- Step 2: Calculate spatial correlation coefficients and determine number of valid time series ---
	rnd_rain(:,7) = rville_precip(:,4);
	rnd_rain2 = zeros(length(rnd_rain), 7);
	
	%%% --- Compare in 3 or 4 year blocks using the indices in "breaks"
	minSpearSim = 1;
	maxSpearSim = 0;
	minPearSim = 1;
	maxPearSim = 0;
    for j = 1:size(breaks,1)
		if j == 1
			[rnd_rain2([1:breaks(j,1)],:), pearSim, spearSim] = methodIC5(rnd_rain([1:breaks(j,1)],:), Rpear);
			minSpearSim = min( min(spearSim(spearSim > 0), minSpearSim));
			minPearSim = min( min(pearSim(pearSim > 0), minPearSim));
			maxSpearSim = max( max(spearSim(spearSim < 0.99), maxSpearSim));
			maxPearSim = max( max(pearSim(pearSim < 0.99), maxPearSim));
		else
			[rnd_rain2([breaks(j-1,1)+1:breaks(j,1)],:), pearSim, spearSim] = methodIC5(rnd_rain([breaks(j-1,1)+1:breaks(j,1)],:), Rpear);
			minSpearSim = min( min(spearSim(spearSim > 0), minSpearSim));
			minPearSim = min( min(pearSim(pearSim > 0), minPearSim));
			maxSpearSim = max( max(spearSim(spearSim < 0.99), maxSpearSim));
			maxPearSim = max( max(pearSim(pearSim < 0.99), maxPearSim));
		end
	end
	
	flags = ones(7);
	for j = 1:6
		for k = j+1:7
			if not(and(spearSim(j,k) >= min_correl, spearSim(j,k) <= max_correl))
				flags(j,k) = 0;
				flags(k,j) = 0;
			end
		end
	end
	
	temp = sum(flags);
	
	%%% Check correlation coefficient
	if sum(temp(1,[1:6]) > 5) >= 5 % If there are at least five valid time series out of six
		printf("correlation OK\n");
	   
		list = zeros(1,6);
		% valid = zeros(length(rnd_rain2),6);
	   
		for j = 1:6
			if sum(flags(:,j)) >= 6 %if all correlation coefficients for this time series are within the observed range
				
				%%% --- Calculate annual rain sums ---
				for yr = 1973:2018
					%%sum the simulated precip for the year "yr"
					idx3 = rville_precip(:,3) == yr;

					rain_yr(j,yr - 1972) = sum(rnd_rain2(idx3,j)); % rain_yr values have averages between 900 and 940 mm -- too high!
					
					% P_tot = sum(rnd_rain2(idx3,j)) + sum(rville_precip(idx3,5)); % random rain (within correlation envelope) + Roseville snow
					% PET = rville_ETo(yr - 1972,2);
					% ratio(yr - 1972,1) = PET / P_tot;
				end

				% maxRatio = max(ratio);
				% minRatio = min(ratio);

				list(1,j) = 1; % all ratios accepted in this version if correlation among time series OK
			end
		end
	
		%%% --- Step 3: Save lists of valid time series ---
		for j = 1:6
			if list(1,j) == 1
				% statsRnd(c,1) = minRatio;
				% statsRnd(c,2) = maxRatio;
				statsRnd(c,[1:46]) = rain_yr(j,:);
			 
				filename = sprintf('Rnd_rainfall27\\rainseries%d.txt',c);
				rainseries = rnd_rain2(:,j);
				save(filename,'rainseries');
				
				assign(assignRow,j) = c;
				c = c + 1;
			end
		end		 
		assignRow
		assignRow = assignRow + 1;

	else
		if and(maxSpearSim > max_correl, minSpearSim < min_correl)
			printf("%g valid time series. Correlation NOT OK: correl both too high and too low.\n", sum(temp(1,[1:6]) > 5));
		elseif maxSpearSim > max_correl
			printf("%g valid time series. Correlation NOT OK: correl too high.\n", sum(temp(1,[1:6]) > 5));
		else
			printf("%g valid time series. Correlation NOT OK: correl too low.\n", sum(temp(1,[1:6]) > 5));
		end
	end
   
end

numAssign = c - 1;
filename = sprintf("Rnd_rainfall27\\stats_Summary%d.txt",numAssign);
save(filename, "statsRnd");

filename = sprintf("assignments_VirtualStations27v%d.txt", version);
save(filename, "assign");

%%% -------------------------------------------------------------
%%% TEST - check average annual rainfall for each row in statsRnd
% for yr = 1973:2018
   % idx2 = rville_precip(:,3) == yr;
   % rville_snow(yr - 1972,1) = sum(rville_precip(idx2,5));
% end

% check = zeros(assignRow - 1,1);
% for i =1:(assignRow - 1)
   % totp = statsRnd(i,:)' .+ rville_snow;
   % check(i,1) = sum(totp)/ny;
% end

% avg_of_all_Ptot_avgs = sum(check) / (assignRow - 1)

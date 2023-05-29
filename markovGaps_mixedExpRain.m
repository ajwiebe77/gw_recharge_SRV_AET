% markovGaps_mixedExpRain.m
% 
% Andrew James Wiebe, 30 Nov 2020; based on code from markovrainseq and using results generated via markov.m
%
% Code for GNU Octave (Eaton et al., 2018).
%
% Objective: Generate a sequence of rainy and non-rainy days based on the Markov chain 
%            method described by Basinger et al. (2010); generate rainfall amounts
%            based on a mixed exponential model (e.g., Li et al., 2013). 
%
% Notes:
% Uses the results of markov.m.
% Uses the rainfall generation code from gen_alderRain27test1 for generating rainfall depths.
% The script that calls this function should use the data from Julian Day 365 for Julian Day 366
% due to limited historical data for leap years.
%
% Input variables: 
%    times - a list of days (format = [day, month, year])
%    rainparams - a vector of three values describing a mixed exponential model for the
%       relative frequency distribution of the intervals between days with rainfall.
%    pw - the overall probability of rainfall, for consideration of the first day in the time series
%    pww - the probability that a rainy day ("wet" - "w") will follow a rainy day ("wet" - "w"), 
%          based on the Julian Day under consideration
%    pdd - the probability that a dry day ("d" - no rain) will follow a dry day ("d" - no rain),
%          based on the Julian Day under consideration
%
% Output variable:
%    rain - a vector of zero and nonzero rainfall amounts corresponding to the transition
%           probabilities between rainy and non-rainy days, and the mixed exponential
%           distribution of non-zero rainfall amounts described by rainparams.
% 
% Notes:
%    The mixed exponential model is based on Li et al. (2013) equations 5a and 5b:
%       f(x) = (p/beta1) * exp(-x/beta1) + ((1-p)/beta2)*exp(-x/beta2),
%    where p represents the mixing fraction, and beta1 and beta2 are scale parameters for two different
%    exponential distributions. Constrains: 0 <= p <= 1, beta1 > 0, beta2 > 0.
%    The equation
%       rt = -bt * ln(vt)
%    was used to simulate the interval or rainfall amounts. In this case, t = time (day index in "rain" vector),
%    bt = either beta1 or beta2 (chosen based on the mixing parameter p), and vt is a uniform random number
%    between 0 and 1.
%
%    The depth generation lines ("depth =") use the natural logarithm, i.e., "log" function in Octave.
%
% References:
%   Basinger, M., Montalto, F., Lall, U., 2010. A rainwater harvesting system reliability model based on 
%      nonparametric stochastic rainfall generator. J. Hydrol. 392, 105-118.
%      https://doi.org/10.1016/j.jhydrol.2010.07.039.
%   Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%      high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
%   Li, Z., Brissette, F., Chen, J., 2013. Finding the most appropriate precipitation probability 
%      distribution for stochastic weather generation and hydrological modelling in Nordic watersheds. 
%      Hydrol. Process. 27, 3718-3729. https://doi.org/10.1002/hyp.9499.
%
function rain = markovGaps_mixedExpRain(times, rainparams, pw, pww, pdd, maxdaily)

	ndays = length(times);
	rainday = 0;
	rain = zeros(ndays,1);

	p = rand();
	if p <= pw
		rainday = 1;

		%%% Assign rainfall amount
		vt = rand(); % uniform random number
		temp = rand(); % uniform random number
		if temp < rainparams(1,1)
		   bt = rainparams(1,2);
		else
		   bt = rainparams(1,3);
		end
		
		rain(1,1) = -bt*log(vt);
		
		if rain(1,1) > maxdaily
			while rain(1,1) > maxdaily
				vt = rand();
				rain(1,1) = -bt*log(vt);
			end
		end
	end

	prevrainday = rainday;
	
	for j = 2:ndays
		jd = getJulianDay(times(j,1),times(j,2),times(j,3));
		
		p = rand();
		if and(prevrainday == 1, p <= pww(jd,1))
			rainday = 1; % specify another rainy day
		elseif and(prevrainday == 0, p <= (1 - pdd(jd,1))) % P(dry,rainy) = (1 - P(dry,dry))
			rainday = 1; % specify a rainy day after a non-rainy day
		else
			rainday = 0;
		end
		
		if rainday == 1
			%%% Assign rainfall amount
			vt = rand(); % uniform random number
			temp = rand(); % uniform random number
			if temp < rainparams(1,1)
			   bt = rainparams(1,2);
			else
			   bt = rainparams(1,3);
			end
			
			rain(j,1) = -bt*log(vt);
			
			if rain(j,1) > maxdaily
				while rain(j,1) > maxdaily
					vt = rand();
					rain(j,1) = -bt*log(vt);
				end
			end
		end
		
		prevrainday = rainday;
	end

return;
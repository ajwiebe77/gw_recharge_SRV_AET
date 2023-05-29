% rndMixedExpRain.m
%
% Andrew J. Wiebe, 21 Jul 2020
% 
% Code for GNU Octave (Eaton et al., 2018).
%
% Objective: generate a daily rainfall timeseries of length equal to ndays based on
%            mixed exponential models for the probability (relative frequency) distributions
%            of both the intervals between rainy days and the daily rainfall amounts.
%
% Input variables: 
%    ndays is total number of days simulated
%    spacingparams is a vector of three values describing a mixed exponential model for the
%       relative frequency distribution of the intervals between days with rainfall.
%    rainparams is a vector of three values describing a mixed exponential model for the
%       relative frequency distribution of the intervals between days with rainfall.
%
% Output variable:
%    rain is a vector of zero and nonzero rainfall amounts corresponding to a sequence of ndays in length
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
%   Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%      high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
%   Li, Z., Brissette, F., Chen, J., 2013. Finding the most appropriate precipitation probability 
%      distribution for stochastic weather generation and hydrological modelling in Nordic watersheds. 
%      Hydrol. Process. 27, 3718-3729. https://doi.org/10.1002/hyp.9499.
%
function rain = rndMixedExpRain(ndays, spacingparams, rainparams)

rain=zeros(ndays,1); % initialize output vector with zeros

vt = rand(); % draw a uniformly distributed random number
temp = rand();
if temp < spacingparams(1,1)
   bt = spacingparams(1,2);
else
   bt = spacingparams(1,3);
end

intvl = ceil(-bt*log(vt));

rday=1+intvl; % choose first rainy day

vt = rand(); % draw a uniformly distributed random number
temp = rand();
if temp <= rainparams(1,1)
   bt = rainparams(1,2);
else
   bt = rainparams(1,3);
end

depth = -bt*log(vt); % try taking the ceiling; Li et al (2013) equation 5b; uses natural logarithm, or "log" function in Octave
rain(rday,1) = depth; % assign a rainfall amount to the first rainy day

while rday<=ndays

	vt = rand(); % draw a uniformly distributed random number
	temp = rand();
	if temp < spacingparams(1,1)
	   bt = spacingparams(1,2);
	else
	   bt = spacingparams(1,3);
	end

	intvl = max(1,ceil(-bt*log(vt)));
   
	rday=rday+intvl; % advance the index of the current rainy day
	
	vt = rand(); % draw a uniformly distributed random number
	temp = rand();
	if temp <= rainparams(1,1)
	   bt = rainparams(1,2);
	else
	   bt = rainparams(1,3);
	end

	depth = -bt*log(vt);
	rain(rday,1) = depth; % assign a rainfall amount

end

rain = rain(1:ndays,1);

return
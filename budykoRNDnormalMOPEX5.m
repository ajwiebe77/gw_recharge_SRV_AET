% budykoRNDnormalMOPEX5.m
%
% Andrew J. Wiebe, 02 Dec 2020; modifies budykoRNDnormalMOPEX2.m, 20 Nov 2020; modifies budykoRNDnormalMOPEX.m, 9 Jul 2020; modifies budykoRNDnormal, 11 Dec 2019
% 
% Code for GNU Octave (Eaton et al., 2018).
%
% Calculate a watershed's average ratio of AET / P (i.e., actual evapotranspiration (ET) divided by precipitation)
%    given pet (annual potential ET) and p (annual total precipitation, i.e., rain plus snow).
% Add some randomness using the 'normrnd' function, where the mean value is the Budyko curve estimate 
%     (Gentine et al., 2012) and the standard deviation is based on the scatter of 
%     annual MOPEX points about the Budyko curve for 43 watersheds where 0.787 <= PET/P <= 0.887
%     (calculated from the MOPEX dataset; Duan et al., 2006). The PET and P values used for this calculation
%     were the long term averages.
%
% Input parameters:
%    PET = potential evapotranspiration for a watershed
%    P = total precipitation for a watershed
%    mopex_bracket = list of bracket definitions with min(PET/P ratio), max(PET/P ratio), mean for bracket, standard deviation for bracket
%
% Output:
%     aet_p is an estimate of the ratio of AET (avg) to P (avg) based on the Budyko curve and
%     observed variability.
%
% References:
% Duan, Q., et al., 2006. Model Parameter Estimation Experiment (MOPEX): An overview of science strategy 
%    and major results from the second and third workshops. J. Hydrol. 320, 3-17.
%    https://doi.org/10.1016/j.jhydrol.2005.07.031.
% Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%    high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
% Gentine, P., D’Odorico, P., Lintner, B.R., Sivandran, G., Salvucci, G., 2012. Interdependence of climate,
%    soil, and vegetation as constrained by the Budyko curve. Geophys. Res. Lett. 39, L19404.
%    https://doi.org/10.1029/2012GL053492.
%

function aet_p = budykoRNDnormalMOPEX5(pet, p, mopex_bracket)

   phi = pet/p;
   % aet_p_curve = sqrt(phi * tanh(1/phi) * (1 - exp(-phi)));
   
   % sigma = 0.11; % use maximum observed standard deviation from 45 MOPEX watersheds where 0.787 <= PET/P <= 0.887
   
   done = 0;
   
   for i = 1:length(mopex_bracket)
      if and(done == 0, phi < mopex_bracket(i,2))
	     % use the mean and standard deviation for the bracket defined by the i'th row
		 aet_p = normrnd(mopex_bracket(i,3), mopex_bracket(i,4), 1, 1); % generate one random value based on the normal distribution
		 done = 1 ;
      end
   end
   % aet_p = normrnd(aet_p_curve, sigma, 1, 1); % generate one random value based on the normal distribution
end

 % 0.5 0.6 0.50807 0.052625 18
 % 0.6 0.7 0.54952 0.079736 185
 % 0.7 0.8 0.60311 0.075727 473
 % 0.8 0.9 0.65423 0.075305 674
 % 0.9 1 0.68956 0.078162 458
 % 1 1.1 0.70996 0.095555 252
 % 1.1 1.2 0.71589 0.074168 99
 % 1.2 1.3 0.6943 0.11449 29
 % 1.3 2 0.63108 0.10859 26
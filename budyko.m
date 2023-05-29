% budyko.m
%
% Andrew J. Wiebe, 18 Jan 2019
% 
% Code for GNU Octave (Eaton et al., 2018).
%
% Calculate a watershed's average ratio of AET / P (i.e., actual evapotranspiration (AET) divided by precipitation)
%    given pet (potential ET) and p (total precipitation, i.e., rain plus snow).
% Uses the Budyko curve equation by Gentine et al. (2012).

% Input parameters:
%    pet = the long-term average PET estimate (mm/unit watershed area)
%    p   = the long term average total precipitation for the watershed
%
% Output:
%     aet_p is the estimate of the ratio of AET (avg) to P (avg) based on the Budyko curve.
%
% Reference:
% Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%    high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
% Gentine, P., D’Odorico, P., Lintner, B.R., Sivandran, G., Salvucci, G., 2012. Interdependence of climate,
%    soil, and vegetation as constrained by the Budyko curve. Geophys. Res. Lett. 39, L19404.
%    https://doi.org/10.1029/2012GL053492.

function aet_p = budyko(pet, p)

   phi = pet/p;
   
   aet_p = sqrt(phi * tanh(1/phi) * (1 - exp(-phi)));

end
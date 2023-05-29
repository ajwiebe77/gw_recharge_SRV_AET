% AJW, 11 Jan 2019
% Use the Villarini et al (2008) function with c1 = 1 (see Villarini et al (2010))

% input parameters
% c = vector with two rows: c = [c2; c3]
% (If vector lenght was 3, c(1) = nugget (but fixed at 1, here))
% c(1) = c2 = decorrelation distance (distance at which spatial correlation coefficient reaches 1/e on the curve; i.e., where process decorrelates)
% c(2) = c3 = a shape factor

% x = distance (km)
% 
% References:
%   Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%      high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
%   Villarini, G., Mandapaka, P.V., Krajewski, W.F., Moore, R.J., 2008. Rainfall and sampling uncertainties: A rain gauge perspective.
%      J. Geophys. Res. 113, D11102. https://doi.org/10.1029/2007JD009214.
%   Villarini, G., Smith, J.A., Baeck, M.L., Sturdevant-Rees, P., Krajewski, W.F., 2010. Radar analyses of extreme rainfall and flooding 
%      in urban drainage basins. J. Hydrol. 381, 266-286. https://doi.org/10.1016/j.jhydrol.2009.11.048.

function r = Villarini (x, c, consts)

    r = 1 * exp (-1 * power(x / c(1),c(2))); % really c(1) = c2 and c(2) = c3
	
end
% AJW, 11 Jan 2019
% Use the Li et al (2013) mixed exponential equation.
% 
% Code for GNU Octave (Eaton et al., 2018).
%
% input parameters:
% x = rainfall depth (mm)
% c = vector with three rows: c = [p; beta1; beta2]
% c(1) = nugget (fixed, here)
% c(2) = decorrelation distance (distance at which spatial correlation coefficient reaches 1/e on the curve; i.e., where process decorrelates)
% c(3) = a shape factor
% 0 <= p <= 1; beta1 > 0; beta2 > 0
%
% References:
%   Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%      high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
%   Li, Z., Brissette, F., Chen, J., 2013. Finding the most appropriate precipitation probability distribution for stochastic weather
%      generation and hydrological modelling in Nordic watersheds. Hydrol. Process. 27, 3718-3729. https://doi.org/10.1002/hyp.9499.

function f = mixedEXP(x, c, consts)

	for i = 1:size(x,1)
		f(i,1) = (c(1)/c(2))*exp(-x(i)/c(2)) + ((1-c(1))/c(3))*exp(-x(i)/c(3));
	end
	
end

% calcWSavgPrecip357annual.m
%
% Andrew J. Wiebe, 18 Aug 2020
%
% Code for GNU Octave (Eaton et al., 2018).
%
% Objective: Calculate the average precipitation over the Alder Creek watershed from three 
%            virtual stations by weighting with their representative subcatchment areas.
%
% Input parameters:
%    totP3 = a matrix where each column represents one station and has n years of total precipitation
%            sums (random rainfall plus Roseville snowfall). Assume 3 stations (i.e., 3 columns in totP3).
%
% Output parameter:
%    p = list of average annual watershed total precip amounts, weighted by subwatershed area.
%
% References:
%   Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%      high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
%

function p = calcWSavgPrecip357annual(totP3)

   area = [25.85297449;... % Thiessen polygon area for station 1
           45.61051953;... % Thiessen polygon area for station 2
           6.596530245;... % Thiessen polygon area for station 3
   ];
	
   n = length(totP3);

   p = zeros(n,1); % annual average precip estimate, via weighted (by area) average
   
   for i = 1:n
      p(i,1) = area' * totP3(i,:)' / sum(area);
   end
   
end
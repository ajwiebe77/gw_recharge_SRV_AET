# gw_recharge_SRV_AET
GNU Octave code for a stochastic vadose zone water budget. A range of groundwater recharge estimates is generated using a Monte Carlo approach based on stochastic spatiotemporal rainfall variability and empirically-derived temporal evapotranspiration variability.

README.txt

Author: Andrew J. Wiebe, 29 May 2023

Licence: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)

** NOTE: **
** 28 Oct 2024 **
** This archive has been superceded by the updated version at https://doi.org/10.17632/stbtg4hnv9.1: **
** Wiebe, A., 2024. Scripts for calculating groundwater recharge via a stochastic vadose zone water budget, Mendeley Data, V1, doi: 10.17632/stbtg4hnv9.1. **
** Please refer to the updated version instead of this one. **


INTRODUCTION

This archive contains code related to the Open Source scientific computation software GNU Octave (https://octave.org/) for a stochastic vadose zone water budget. A range of groundwater recharge estimates is generated using a Monte Carlo approach based on stochastic spatiotemporal rainfall variability and empirically-derived temporal evapotranspiration variability. The approach is a modified version of the one outlined in Wiebe (2020). 

This project is useful for calculating the variability in recharge rates due to spatial/temporal variability in rainfall and annual variability in actual evapotranspiration. These components of groundwater recharge uncertainty are rarely considered.

QUICK GUIDE

What you need to do to run this code:

Download and install GNU Octave (https://octave.org/).

Download the .m files related to this repository.

Compile long-term regional precipitation (rainfall and snowfall) and temperature data for your watershed for at least one station.

Compile local rainfall data (several stations)for your watershed (several years of data are recommended). Calculate Pearson and Spearman Rank correlation coefficients between the stations.

Fit a distribution to the relative frequency histogram of the long-term rainfall data (e.g., exponential or mixed exponential). The example given here used a mixed exponential function (there also needs to be a function written to generate random rainfall amounts based on the distribution).

Calculate potential or reference evapotranspiration based on temperature (and additional data, if available). A tool such as the ETo Calculator (Raes, 2009) may be used.

Calculate the ratio of average annual potential evapotranspiration (PET) to average annual precipitation (P) for your watershed using long-term data (several decades). Select MOPEX watersheds (for example, the US MOPEX data: Duan et al., 2006;  https://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/) similar to your watershed's PET/P ratio and estimate the mean and standard deviation for normal distributions of AET/P within different PET/P ranges (e.g., 0.1 intervals).



METHODS - VADOSE ZONE WATER BUDGET

Equations (from Wiebe, 2020):

(1) Change in VZ storage = I - AETvz - R,

where I is infiltration, AETvz is vadose zone actual evapotranspiration, and R is recharge.

(2) R = P - AETvz - QoL - change in VZ storage,

where P is precipitation, and QoL is overland flow.

(3) R_tot = sum from yr=1 to yr=n {P_yr - AETvz_yr - Q_oL,yr},

where n is the number of years to consider for each Monte Carlo realization, R_tot is the total groundwater recharge over all years for the realization, and yr is the index of the year.

(4) R_tot = sum from yr=1 to yr=n {P_yr - P_yr(AETvz_yr / Pyr) - (1-BFI_yr)Qtot_yr},

where BFI is baseflow index (i.e., the fraction of annual total streamflow constituted by groundwater baseflow), and Qtot is annual total streamflow.


ASSUMPTIONS:

1) Vadose zone storage change over n years is negligible (0 mm)

2) n years captures PET variability

3) the Budyko curve may reasonably be used for the watershed (few anthropogenic influences)

4) Ratios of AETtot_yr / P_yr ratios may reasonably be represented by random values generated based on reference watersheds

5) a correction factor for AETvz AETtot may be estimated

6) snowfall is relatively uniform over the watershed and a spatial analysis of this component of total precipitation is not desired.



REFERENCES

Duan, Q., et al., 2006. Model Parameter Estimation Experiment (MOPEX): An overview of science strategy and major results from the second and third workshops. J. Hydrol. 320, 3-17. https://doi.org/10.1016/j.jhydrol.2005.07.031.

Raes, D., 2009. The ETo Calculator: Evapotranspiration from a reference surface. Reference Manual Version 3.1. Food and Agriculture Organization of the United Nations, Land and Water Division, Rome, IT. http://www.fao.org/land-water/databases-and-software/eto-calculator/en/. (Accessed 01.02.2018).

Wiebe, A.J., 2020. The influences of spatially variable rainfall and localized infiltration on groundwater recharge in a water management context, PhD dissertation, University ofWaterloo,Waterloo, ON, Canada, http://hdl.handle.net/10012/16476.



AUTHOR/CONTACT INFORMATION

This code was written by Andrew J. Wiebe:

https://www.researchgate.net/profile/Andrew-Wiebe-3


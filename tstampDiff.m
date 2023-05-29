% tstampDiff.m
% AJW, 20 Jul 2020; modifies tstampDiff.m, 19 Mar 2015
%
% A function that returns the difference between two timestamps (in hours), 
%    i.e., dd2-mh2-yyyy2 hh2:mm2 - dd1-mh1-yyyy1 hh1:mm1
% Returns a positive value if timestamp 2 > timestamp 1
% If one of the timestamps is unreasonable, return error code -9999
% Assumes a 24h clock
function diff = tstampDiff(dd1, mh1, yyyy1, hh1, mm1, dd2, mh2, yyyy2, hh2, mm2)
   
   % --------------------------------------------------------------
   % VARIABLES

   hours_2 = 0;
   hours_1 = 0;

   % --------------------------------------------------------------

   if is_leap_year(yyyy1)
      n_days_yr1 = 366; % leap year
      n_days_month1 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
   else
      n_days_yr1 = 365;
      n_days_month1 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
   end

   if is_leap_year(yyyy2)
      n_days_yr2 = 366; % leap year
      n_days_month2 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
   else
      n_days_yr2 = 365;
      n_days_month2 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
   end
   % --------------------------------------------------------------


   if (or(mh1 > 12, mh2 > 12, dd1 > n_days_month1(mh1), dd2 > n_days_month2(mh2), mm1 > 60, mm2 > 60 ))
      result = -9999;
      return;
   end


   %case 1 - only minutes are different

   if and(yyyy2 - yyyy1 == 0, mh2 - mh1 == 0, dd2 - dd1 == 0, hh2 - hh1 == 0)
      diff = (mm2 - mm1) / 60;
   elseif and(yyyy2 - yyyy1 == 0, mh2 - mh1 == 0, dd2 - dd1 == 0)

      %case 2 - hours differ
      diff = (hh2 + mm2/60) - (hh1 + mm1/60);

   elseif and(yyyy2 - yyyy1 == 0, mh2 - mh1 == 0)

      %case 3 - days differ

      diff = ((dd2 - 1) * 24 + hh2 + mm2/60) - ((dd1 - 1) * 24 + hh1 + mm1/60); % no. of hours into month minus no. of hours into month

   elseif yyyy2 - yyyy1 == 0

      %case 4 - months differ

      if mh2 > mh1

         hours_2 = (dd2 - 1) * 24 + hh2 + mm2/60;  %no. of hours into mh2 current month
      
         for i = (mh2 - 1):-1:(mh1 + 1)
            hours_2 = hours_2 + 24 * n_days_month2(i);  % add all hours in months between the two months mh2 and mh1
         endfor


         % Now, add the hours from timestamp1 to the end if its month
         diff = hours_2 + (24 * n_days_month1(mh1) - ((dd1 - 1) * 24 + hh1 + mm1/60));

      else
         % so, mh1 > mh2

         hours_1 = (dd1 - 1) * 24 + hh1 + mm1/60;  %no. of hours in mh1 current month
      
         for i = (mh1 - 1):-1:(mh2 + 1)
            hours_1 = hours_1 + 24 * n_days_month1(i);  % add all hours in months between the two months mh2 and mh1
         endfor

         % Now, add the hours from timestamp1 to the end if its month
         diff = hours_1 + (24 * n_days_month2(mh2) - ((dd2 - 1) * 24 + hh2 + mm2/60));
         diff = diff * - 1;
      end

   else

      % case 5 - years differ

      if yyyy2 > yyyy1

         hours_2 = (dd2 - 1) * 24 + hh2 + mm2/60;  %no. of hours in mh2 current month

         for i = (mh2 - 1):-1:1
            hours_2 = hours_2 + 24 * n_days_month2(i);  % add all hours in months in the current year
         endfor
      
         for year = (yyyy2 - 1):-1:(yyyy1 + 1)

            if and(mod(year, 4) == 0, not(or(mod(year, 100) == 0, mod(year, 400) == 0)))
               n_days_yr = 366; % leap year
               n_days_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
            else
               n_days_yr = 365;
               n_days_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            end


            for m = 12:-1:1
               hours_2 = hours_2 + 24 * n_days_month(m);  % add all hours in months in intervening years
            endfor
         endfor

         %
         for i = 12:-1:(mh1 + 1)
            hours_2 = hours_2 + 24 * n_days_month1(i);  % add all hours in months between (mh1 + 1) and the end if its year
         endfor

         % Now, add the hours from timestamp1 to the end if its month
         diff = hours_2 + (24 * n_days_month1(mh1) - ((dd1 - 1) * 24 + hh1 + mm1/60));

      else

         % So, yyyy1 > yyyy2

         hours_1 = (dd1 - 1) * 24 + hh1 + mm1/60;  %no. of hours in mh1 current month

         for i = (mh1 - 1):-1:1
            hours_1 = hours_1 + 24 * n_days_month1(i);  % add all hours in months in the current year
         endfor
      
         for year = (yyyy1 - 1):-1:(yyyy2 + 1)

            if and(mod(year, 4) == 0, not(or(mod(year, 100) == 0, mod(year, 400) == 0)))
               n_days_yr = 366; % leap year
               n_days_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; 
            else
               n_days_yr = 365;
               n_days_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            end

            for m = 12:-1:1
               hours_1 = hours_1 + 24 * n_days_month(m);  % add all hours in months in intervening years
            endfor
         end

         %
         for i = 12:-1:(mh2 + 1)
            hours_1 = hours_1 + 24 * n_days_month2(i);  % add all hours in months between (mh1 + 1) and the end if its year
         endfor

         % Now, add the hours from timestamp1 to the end if its month
         diff = hours_1 + (24 * n_days_month2(mh2) - ((dd2 - 1) * 24 + hh2 + mm2/60));
         diff = diff * - 1;
      end


   end

end
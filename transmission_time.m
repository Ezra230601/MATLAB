function [time_tx, dtS] = transmission_time(time_rx, range, Eph, SP3, err_tropo, err_iono, dtR, frequencies, obs_comb, lambda)

%https://ji-rss-aa.pw/link/XBQcUEKXtjiqLDe7
%?mu=0
% INPUT:
%   time_rx   = reception time
%   range     = observed range
%   sat       = satellite index
%   Eph       = ephemeris
%   SP3       = structure containing precise ephemeris data
%   sbas      = SBAS corrections
%   err_tropo = tropospheric delays
%   err_iono  = ionospheric delays
%   dtR       = receiver clock offset
%
% OUTPUT:
%   time_tx = transmission time
%   dtS     = satellite clock error
%
% DESCRIPTION:
%   Compute the signal transmission time.
%SBAS clock offsets

global glp

time_tx_RAW = time_rx - (range - err_tropo - err_iono) / glp.vlight + dtR;
if isempty(SP3)
    % GLONASS
    if strcmp( char( Eph(31) ),'R' )
        dtS = sat_clock_error_correction(time_tx_RAW, Eph) ;
        dtS = sat_clock_error_correction(time_tx_RAW - dtS, Eph);
    else
        dtrel = relativistic_error_correction(time_tx_RAW, Eph);
        %group delay correction term
        %--- this correction term is only for the benefit of "single-frequency"
        %    (L1 P(Y) or L2 P(Y)) users
        if nargin>9 && ~strcmp(obs_comb,'IONO_FREE')
            
            tgd = Eph(28);
            if (frequencies(1) == 2)
                gamma = (lambda(1,2)/lambda(1,1))^2;
                tgd = gamma*tgd;
            end
        else
            tgd=0;
        end
        
        dtS = sat_clock_error_correction(time_tx_RAW, Eph);
        dtS = dtS + dtrel - tgd ;
        dtS = sat_clock_error_correction(time_tx_RAW - dtS, Eph);
        dtS = dtS + dtrel - tgd ;
        
    end
else
    % SP3 data not available in this version
    % dtS = interpolate_SP3_clock(time_tx_RAW, SP3, sat);
    fprintf('SP3 data is not available in this version!\n');
end
time_tx = time_tx_RAW - dtS;

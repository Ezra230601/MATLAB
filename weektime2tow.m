function [sow] = weektime2tow(week, time)

% SYNTAX:
%   [sow] = weektime2tow(week, sow);
%
% INPUT:
%   week = GPS week
%   time = GPS time (continuous since 6-1-1980)
%
% OUTPUT:
%   sow  = GPS seconds-of-week
%
% DESCRIPTION:
%   Conversion from GPS time in continuous format (similar to datenum) to
%   GPS time in week, seconds-of-week.



sow = time - week*7*86400;

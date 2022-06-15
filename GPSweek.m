%*******************************************************
% function [GPS_wk, GPS_sec_wk] = GPSweek(Y,M,D,H,min,sec)
%
% DESCRIPTION:
%
%  This function finds GPS week and GPS second of the week based on the
%  input calendar date and time. 
%  
% ARGUMENTS:
%
%  Y - year (4 digits)
%  M - month
%  D - day
%  H - hours (in military time)
%  min - minutes
%  sec - seconds
%  
% OUTPUT:
%
%  GPS_wk - GPS week
%  GPS_sec_wk - GPS seconds of the week
%  
% CALLED BY:
%
%  loadInput
%
% FUNCTIONS CALLED:
%
%  None
%
% MODIFICATIONS:    
% 
%             ??  :  Kristine Larson - Original
%       02-05-02  :  Lisa Reeh
% 
% 
% Colorado Center for Astrodynamics Research
% Copyright 2002 University of Colorado, Boulder
%*******************************************************
function [GPS_wk, GPS_sec_wk] = GPSweek(Y,M,D,H,min,sec)
format long g;
if nargin==4
	min=0;
	sec=0;
end;
% UT = H+(min/60)+(sec/3600);
% if M > 2
% 	y=Y;
% 	m=M;
% else
%     y=Y-1;
% 	m=M+12;
% end;
% 
% JD = fix(365.25*y) + fix(30.6001*(m+1)) + D + (UT/24) + 1720981.5;
JD = juliandate([Y,M,D,H,min,sec]);
GPS_wk = fix((JD-2444244.5)/7);
% GPS_sec_wk = round( ( ((JD-2444244.5)/7)-GPS_wk)*7*24*3600);
GPS_sec_wk = round( ( ((JD-2444244.5)/7)-GPS_wk)*7*24*3600,2);
% Ensure that 1 < GPS_wk < 1024
GPS_wk = mod(fix((JD-2444244.5)/7), 1024);
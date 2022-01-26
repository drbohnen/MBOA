clear all; close all

%% Grid Input, Verbose and Plot Production Settings 
mapgrid=('mata_exam.grd');      % load a NetCDF GMT file
                                %
direc=('./');  % directory path to grid file,  include trailing "/". 
                        % For example ('/myfiles/grids/')   
                        % 
plotlevel=2;    % 1 makes map with seamounts (displayed & saved graphic)
                % 2 makes map plus 3-D view of each edifice (saved graphics)
                % 
verb=0;         % verb =1 for verbose, verb = 0 for relatively quiet 
                % 
writestats=1;   % set to one to write formated stats file 
                %
figsavetype=1;  % file format to save the graphics 
                % 0 = matlab fig file 
                % 1 = PNG 
                % 2 = postscript 

%% Inputs that affect the selection of basal contours 
contype=0;      % 0= use matlab's contourm
                % 1 = use GMT's grdcontour (usually faster). Must have GMT installed
                % and in the system path. Only tested for OSX, Linux 
                % 
intpflag=0;     % set to 1 to fill small NaN data gaps in the grid    
                %
clstep=40;      % Integer contour interval (set to 1/3 to 1/2 the
                % minimum height mound you want to detect 
                % 
clstart=10;     % rounds start level of contouring to the multple of X that
                % is closest to the minimum elevation point in the grid.  
                % Typically an integer between 1 m and clstep, but could 
                % be <1 m for sub-meter scale topography. 
                % 
ptsmin=15;      % Minimum number of points in the contour to have it 
                % considered as a base (10-15 works well for ship-based multibeam, 
                % min vlaue is 8).  
                %
maxdim=12;      % Maximum dimension of a basal contour in km. Defined as  
                % max distance between 2 points on the closed contour.
                %
elp_ratio=4;    % Maximum eliipical ratio to keep as possible basal contour.
                %
Efit=0.4;       % Threshold misfit to ellipical model. Smaller Efit=better fit,
                % Exclude contours with values larger than Efit thresfold  
                %
minHR=tand(2);  % min HR radius (average) for a base to be concidered as a
                % closed contour (typically 2-3 deg)
                % 
maxHR=tand(60); % max HR radius (average) for a base consideration of a 
                % closed contour (typically 60 deg)

                
%% These parameters are used by the base modification routine
Naz=180;            % Number of profiles to consider in a 1/2 circle 
                   % (90 gives 2 deg sectors, 180 gives 1 deg sectors)
                   %
Rstep=150;         % Range step for A/O and slope checks in meters 
                   % Commonly 1.5-3 time the grid resolution, 
                   % Minimum value = 1*gird spacing, Maximum Value = 10* grid spacing  
                   % 
AOcrit = 1;        % X-sec Area to perimeter ratio check, 
                   % AOcirit=0 no checking; AOcirit=1 standard checking; 
                   % AOcirit =2 include base of profile in calculation;
                   % 
MaxHdiff = 5;      % Fractional difference in the basal heigh forward verses 
                   % reverse profile directions 
                   % (typically 0.25-1 ; max = 5). Set to 5 initally to
                   % disregard 
                   % 
AdthresR= 5;       % Base can extend "AdthresR" times the initial radius of 
                   % the edifice along each azimuth. 
                   % (typically 0.25-2; maximum value is 5 )
                   % 
SlopePortion=0.25; % Fractional portion of the profile near the edge where slope is 
                   % checked along each azimuth (0.25 typical) 
                   % 
SlopeRatioThres=0; % End modification if the slope near the edge of the  
                   % profile becomes less than this amount of the 
                   % initial slope.  Set to 0 to disregard 

%% These parameters smooth the final basal picks  
FT=0;  % zero= mean,    one = median  two = maximum   three = minimum
FL=7;  % filter length for averaging or median filter (Use an odd number 3, 5 or 7) 


%% Run the program 
callscript_MBOA1_1   % scr


%% Comments 
% 
% Input file for MBOA 1.1 (Modified Basal Outlining Alg) 
% Primary Authors: Del Bohnenstiehl (NCSU) & Julia Howll (USC). 
% Contacts: drbohnen@ncsu.edu 
% 
% Notes: 
% Software Requirements: 
% MATLAB   Version 7.11          (R2010b)
% Mapping Toolbox  Version 3.2   (R2010b)
% Statistics Toolbox Version 7.4 (R2010b)
% Developed and tested using Mac OSX 10.6 and Ubunto v110 & 11.  
% 
% Hardware Requirements 
%  For working with large grids, 
% with significant relief we recommended OSX (which allocates swap/virtual 
% memory as needed) or Linux with a Swap Partition with 10's GB allocated. 
%%-------------------------------------------------------------------------
% estimate_radar_range: Calculated Radar range of the target with given
% beat frequency and parameters
% @parameters: beat freequency of each target 
% To Find: BSw fir 1 m resolution, TCh Chirp time for Max Range, define the
% frequency shifts
% Inputs: Fb = [0, 1.1e6, 13e6,24e6 ]
%%
function R = estimate_radar_range(Fb)
    c = 3e8;
    Rres = 1; %Resolution of the Radar
    % TODO : Find the Bsweep of chirp for 1 m resolution
    BSw = c/(2*Rres); %Bandwidth of the radar for 1 m resulution 
    RMax = 300;
    % TODO : Calculate the chirp time based on the Radar's Max Range
    % The sweep time can be computed based on the time needed for the 
    % signal to travel the maximum range. In general, for an FMCW radar 
    % system, the sweep time should be at least 5 to 6 times the round 
    % trip time. This example uses a factor of 5.5
    Tch = 5.5*2*RMax/c;    
    % TODO : define the frequency shifts and calculate the range for each
    % object
    R = c*Tch*Fb/(2*BSw);
    % Display the calculated range
    %disp(R);       
end
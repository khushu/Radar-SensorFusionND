%%-------------------------------------------------------------------------
% estimate_velocity_from_dopler_shift: function computes dopler velocity
% from DS frequency
% @parameters: dopler frequency of the target
%%
function velocity = estimate_velocity_from_dopler_shift(dopler_frequency)
    % Doppler Velocity Calculation
    c = 3*10^8;         %speed of light
    frequency = 77e9;   %frequency in Hz
    
    % TODO : Calculate the wavelength
    %Wavelength calculation of the radar wave
    lamda = c/frequency;

    % TODO : Define the doppler shifts in Hz using the information from above
    % TODO : Calculate the velocity of the targets  fd = 2*vr/lambda

    velocity = dopler_frequency * lamda / 2;    
end

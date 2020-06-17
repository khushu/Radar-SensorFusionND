% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Define the number of training cells and guard cells
% Start sliding the window one cell at a time across the complete FFT 1D array. Total window size should be: 2(T+G)+CUT
% For each step, sum the signal (noise) within all the leading or lagging training cells
% Average the sum to determine the noise threshold
% Using an appropriate offset value scale the threshold
% Now, measure the signal in the CUT, which is T+G+1 from the window starting point
% Compare the signal measured in 5 against the threshold measured in 4
% If the level of signal measured in CUT is smaller than the threshold measured, then assign 0 value to the signal within CUT.

% Close and delete all currently open figures
close all;
% Data_points
Ns = 1000;
% Generate random noise
s=abs(randn(Ns,1));
%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
tar = [100 ,200, 300, 700];
s(tar)=[8 9 4 11];
%plot the output
%plot(s);
% TODO: Apply CFAR to detect the targets by filtering the noise.
% 1. Define the following:
% 1a. Training Cells
% 1b. Guard Cells 
T = 25;
G = 2;
% Offset : Adding room above noise threshold for desired SNR 
offset=3;

% Vector to hold threshold values 
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];

% 2. Slide window across the signal length
for i = 1:(Ns-2*(G+T))

    CutIndex = i+G+T;
    LagTraningIndex = i:1:i+T-1;
    LeadTraningIndex = CutIndex+G+1:1:CutIndex+G+T;
    treshold_index = [ LagTraningIndex LeadTraningIndex ];    
    % 2. - 5. Determine the noise threshold by measuring it within the training cells
    threshold = sum(s(treshold_index))/(2*T);
    %Multiply the treshold with an offset
    threshold = threshold*offset;    
    % 6. Measuring the signal within the CUT
    signal = s(CutIndex);
    % 8. Filter the signal above the threshold
    if (s(CutIndex)<threshold)        
        signal = 0;
    end
    %Construct the signal and the treashold
    threshold_cfar = [threshold_cfar {threshold}];
    signal_cfar = [signal_cfar, {signal}];
end

% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')
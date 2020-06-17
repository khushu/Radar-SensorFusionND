clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s (chapter its 70 m/s 
% Velocity Resulution  = 3 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;
Rres = 1; %Resolution of the Radar
RMax = 200;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
% Innitial velocity = constant = TargetValocity
% Target range is reducing (since the radar is assumed to accelerate
% towards the Target vehicle
TargetRange = 100; 
TargetValocity = 50;

%% FMCW Waveform Generation

% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

B = c/(2*Rres); %Bandwidth of the radar for 1 m resulution
Tchirp = 5.5*2*RMax/c; % Chirp Time   
AlphaSlope = B/Tchirp; % Chirp Slope

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)     
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    % Range of the vehicle Innitial range + m 
    r_t(i) = TargetRange + TargetValocity * t(i);
    
    %Todo, to verify if it is - or + the target range (Closing in or out)
    
    %Delayed version of the signal would be the time travel by the signal 
    td(i) = 2 * r_t(i)/c;
    
    TTravelDelay = t(i) - td(i); 
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*((fc*t(i))+ (AlphaSlope*t(i)*t(i)/2)));
    Rx(i) = cos(2*pi*((fc*(TTravelDelay)) + AlphaSlope* (TTravelDelay*TTravelDelay)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT

% *%TODO* :
% Reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
% Range and Doppler FFT respectively.
Mix=reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Y = fft(Mix,Nr);

%Todo: Chirp time * Bandwidth (Time Bandwidth product, for estimation of length of the FFT) 
L = Tchirp * B;
%L = Nr; 

%TODO* :
% Take the absolute value of FFT output
% Taking half the samples and and Normalizing and taking absolute
signal_fft = abs(Y(1:(L/2)+1)/L);

%TODO: making the samples equal to the length of the signal 
f = B*(0:(L/2))/L;

%plotting the FFT and the range
figure ('Name','Range from First FFT');
subplot(2,1,1);
plot(f,signal_fft);
legend('FFT of the mixed signal'); 
R = (c*Tchirp*f)/(2*B);
subplot(2,1,2);
plot(R,signal_fft);
legend('Peak at the target range'); 

axis ([0 200 0 1]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure ('Name','Dopler from Second FFT'),surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.

% 1. Determine the number of Training cells for each dimension Tr and Td. 
% Similarly, pick the number of guard cells Gr and Gd.

Tr = 10;    %Training Range
Td = 3;    %Training Doppler
Gr = 2;     %Guard Range
Gd = 2;     %Guard Doppler

% *%TODO* :
% offset the threshold by SNR value in dB
Offset = 7;

% 2. Slide the Cell Under Test (CUT) across the complete cell matrix
% Select the grid that includes the training, guard and test cells. 
% Grid Size = (2Tr+2Gr+1)(2Td+2Gd+1).

%Total cell count for averaging
TotalCellCount = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);
 

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
%changed to just a vlaue 
TrainingSum = 0; 


%CFAR Signal value 
CFARSignal = zeros(Nr/2,Nd);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.


for j = 1:(Nr/2)-(2*(Gd+Td))
    %Starting the CUT Y from j + Gd + Td index
    CutYIndex = j + Td + Gd;
    
    for i = 1:(Nd-2*(Gr+Tr))        
        
        %Starting the CUT X from i + Gr + Tr index
        CutXIndex = i + Tr + Gr;
        
        % ***ttttttttt***| A: Indexes for Section A 
        % ***ttttttttt***| CutXIndex-Tr-Gr:1:CutXIndex+Gr+Tr X CutYIndex-Td-Gd:1:CutYIndex-Gd-1 
        % ***ttttttttt***| 
        % ***ttttttttt***| 
        % ***ttgggggtt***| B: Indexes for Section B 
        % ***ttggCggtt***| CutXIndex-Tr-Gr:1:CutXIndex-Gr-1 X CutYIndex-Gd:1:CutYIndex+Gd
        % ***ttgggggtt***| CutXIndex+Gr+1:1:CutXIndex+Tr+Gr X CutYIndex-Gd:1:CutYIndex+Gd
        % ***ttttttttt***|
        % ***ttttttttt***| C: Indexes for Section C 
        % ***ttttttttt***| CutXIndex-Tr-Gr:1:CutXIndex+Gr+Tr X CutYIndex+Gd+1:1:CutYIndex+Td+Gd
        % ***ttttttttt***| 
        
        %A + C portion :
        TrainingIndexX_AC   = CutXIndex-Tr-Gr:1:CutXIndex+Gr+Tr; 
        TrainingIndexY_A    = CutYIndex-Td-Gd:1:CutYIndex-Gd-1; 
        TrainingIndexY_C    = CutYIndex+Gd+1:1:CutYIndex+Td+Gd;
        
        %B portion 
        LagTraningIndexX_B  = CutXIndex-Tr-Gr:1:CutXIndex-Gr-1;
        LeadTraningIndexX_B = CutXIndex+Gr+1:1:CutXIndex+Gr+Tr;
        TrainingIndexY_B    = CutYIndex-Gd:1:CutYIndex+Gd;
        
        TrainingXIndexA = [ TrainingIndexX_AC ];
        TrainingYIndexA = [ TrainingIndexY_A ];
        
        TrainingXIndexB = [ LagTraningIndexX_B LeadTraningIndexX_B];
        TrainingYIndexB = [ TrainingIndexY_B ];
        
        TrainingXIndexC = [ TrainingIndexX_AC ];
        TrainingYIndexC = [ TrainingIndexY_C ];        
      
        
        %Calculating the sum of the training
        TrainingSum = sum(db2pow(RDM(TrainingYIndexA, TrainingXIndexA)),'all');
        TrainingSum = TrainingSum + sum(db2pow(RDM(TrainingYIndexB, TrainingXIndexB)),'all');
        TrainingSum = TrainingSum + sum(db2pow(RDM(TrainingYIndexC, TrainingXIndexC)),'all');        
      
       
        % Sum all the cells and devide by total number of cells 
        % 5. Measure and average the noise across all the training cells. 
        % This gives the threshold
        AverageNoise = TrainingSum/TotalCellCount;

        
        % 6. Add the offset (if in signal strength in dB) to the threshold to keep 
        % the false alarm to the minimum.    
        AverageNoise = pow2db(AverageNoise) + Offset;
        AverageNoise = db2pow(AverageNoise); 
        
        % 7. Determine the signal level at the Cell Under Test.

        % 8. If the CUT signal level is greater than the Threshold, assign a value 
        % of 1, else equate it to zero.


        % 9. Since the cell under test are not located at the edges, due to the 
        % training cells occupying the edges, we suppress the edges to zero. 
        % Any cell value that is neither 1 nor a 0, assign it a zero.
        
        % *%TODO* :
        % The process above will generate a thresholded block, which is smaller 
        %than the Range Doppler Map as the CUT cannot be located at the edges of
        %matrix. Hence,few cells will not be thresholded. To keep the map size same
        % set those values to 0. 
       
        CUTSignalValue = db2pow(RDM(CutYIndex, CutXIndex));
        if (CUTSignalValue <= AverageNoise)
            CUTSignalValue = 0;
        else 
            CUTSignalValue = 1;
        end         

        %Appending the CUT signal value
        CFARSignal(CutYIndex,CutXIndex) = CUTSignalValue; 

    end 
end


% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure ('Name','2D CFAR'),surf(doppler_axis,range_axis,CFARSignal);
colorbar;


 
 
## Radar-SensorfusionND
Radar target generation and detection:
- This is the 4th project as part of the Udacity's Sensor fusion nano degree program with Radar sensor. 

**Summary of the results and observations:**
- FMCW Waveform Design as per given specification 
- Simulation Loop to generate the mix signal after receiving the reflected signal
- Detection of target at almost the right location using Range FFT (1st FFT) 
- Doppler  FFT (2nd FFT) (was already implemented)
- Calculated 2D CFAR dynamic threshold for removing the noise 
- in 2D CFAR The object is detected at the correct location, 100m away without noise.
- Velocity in the doppler map is not accurate but is almost centered with correct value 
- Velocity doppler map is spread across a range of values with an accuracy of ~+-10 m/s centered around 50 m/s


**Implementation steps for the 2D CFAR process:**
- Chose values for Guard cells and Training cells
- Chose a value for the Offset
- Made the CFAR signal to all zero to suppress non thresholded values
- Designed loop to go through the 2D RDM signal to calculate the CFAR signal
  - For each iteration Calculate the indexes of the training cells
  - To calculate the dynamic CFAR threshold from the Range Doppler Map, convert the signal value from db2pow
  - Average the overall training indexes from RDM (sum of all RDM indexes/total number of cells) for the Cell Under Test (CUT)
  - Since the average is in logarithmic scale, add Offset to it and convert it back to db.
  - For each value of RDM Signal for the CUT, perform thresholding based on the calculated CFAR Threshold (0, if < CFAR threshold; else  1)
  - Append the signal to the CFAR signal at the CUT location
- Plot the result.



**Selection of Training cells:**

- If we increase the training doppler, then the spread across velocity increases from +-10 m/s to beyond
- If we reduce the training cells < 3, then the object splits with multiple peaks in doppler map, but one of them is very close to the actual velocity
- Overall, this parameter is quite sensitive, and may be limitation of the method to identify the velocity precisely
  - **Tr = 10, Td = 3**

**Selection of Guard cells:**
- Guard cell of around 1-4 is a good value
- If we reduce the guard cell to zero, there are many peaks all over the spectrum, and the algorithm has many false alarms
  - **Gr=2, Gd=1**


**Selection of Offset value:**
- Offset value of 9 gives a very good result, with minimal spread in the doppler velocity, < +/- 10m/s
- Offset of 10 and above results in no output
- Above 6 and below 9 is good
- If we reduce the offset < 6, there are many false alarms 
  - **Offset = 9**

**Steps taken to suppress the non-thresholded cells at the edges:**
- The signal is initialized to zero first and only the thresholded values are allowed to pass through

![summary of results](/results/summary.png)

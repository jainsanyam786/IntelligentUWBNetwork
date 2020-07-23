%Setting parameters for Root Raised Cosined Transmitter Fiter

Nsym = 6;           % Filter span in symbol durations
beta = 0.25;         % Roll-off factor
sampsPerSym = 4;     % Upsampling factor

dataToTransmit = ones(16,1);
dataToTransmit(3:3:16)=-1;

%Visualize the origibal data
figure(1);
stem(dataToTransmit);
title("Data To Transmit");
xlabel("Symbol index");ylabel("Amplitude");
xticks(0:1:18);
axis([0 18 -1.5 1.5]);

%upsampling data for transmission
Datalength = 64;
sig=zeros(Datalength,1);
sig(4:4:Datalength)=1;
sig(12:12:Datalength)=-1;

%assumed data rate
DataRate = 1000;

%sampling frequency
sampleFrequency = DataRate * sampsPerSym;

%intializing RaisedCosineTransmitFilter
txfilter = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym);


% Normalize to obtain maximum filter tap value of 1
b = coeffs(txfilter);
txfilter.Gain = 1/max(b.Numerator);



% Time vector sampled at symbol rate in milliseconds
tx = 1000 * (0: Datalength - 1) / DataRate;

% Visualize the impulse response
% fvtool(rctFilt, 'Analysis', 'impulse');


% Filter
yo = txfilter([sig ; zeros(Nsym/2,1)]);

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*DataRate);

% Correct for propagation delay by removing filter transients
yo = yo(fltDelay*sampleFrequency+1:end);

% Time vector sampled at sampling frequency in milliseconds
to = 1000 * (0: Datalength*sampsPerSym - 1) / sampleFrequency;

% Plot data
figure(2);
Sigplot = stem(tx, sig, 'kx'); hold on;
% Plot filtered data
FilterData = plot(to, yo, 'r--'); 
% Set axes and labels
axis([0 30 -1.7 1.7]);  xlabel('Time (ms)'); ylabel('Amplitude');
xticks(0:5:30); yticks(-1.5:0.25:1.5);
legend([Sigplot,FilterData],'Transmitted Data', 'RRCFilter Output', 'Location', 'southeast')
hold off;

%Setting parameters for Root Raised Cosined Filter Fiter


rxfilter  = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol', sampsPerSym, ...
  'DecimationFactor',       1);

rxfilter.Gain = 1/txfilter.Gain;

yr = rxfilter([yo; zeros(Nsym*sampsPerSym/2, 1)]);

% Correct for propagation delay by removing filter transients
yr = yr(fltDelay*sampleFrequency+1:end);

%Downsampling the received data to get transmitted data
receivedData = downsample(yr,4);

% Plot filtered data.
figure(3);
Sigplot = stem(tx, receivedData, 'kx'); hold on;
hold on;
plot(to, yr, 'b--'); hold off;
% Set axes and labels.
axis([0 30 -1.7 1.7]);  xlabel('Time (ms)'); ylabel('Amplitude');
legend('Received Data','Rcv Filter Output', 'Location', 'southeast');


%Downsampling the data from first Downsampling to get original data
receivedData = downsample(yr,4);
extractedOriginalData = downsample(receivedData,4,3);
figure(4);
stem(extractedOriginalData);
title("Extracted Original Data");
xlabel("Symbol index");ylabel("Amplitude");
xticks(0:1:18);
axis([0 18 -1.5 1.5]);
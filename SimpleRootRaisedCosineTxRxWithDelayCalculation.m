differArray = zeros(1,14);
delayArray = zeros(1,14);
j=1;
for i = linspace(-30,100,14)
    outPut = test(i);
    differArray(j) = outPut(1);
    delayArray(j) = outPut(2);
    j=j+1;
end

x = linspace(-30,100,14);
figure();
diff = plot(x,differArray,'m--');hold on;
xlabel('SNR (dB)');ylabel('Diff');
title("Diff b/w 1st and 2nd Max w.r.t to SNR");
legend(diff,'Diff b/w 1st and 2nd Max','Location', 'southeast')
hold off;

x = linspace(-30,100,14);
figure();
delay = plot(x,delayArray, 'k--'); 
title("Delay w.r.t to SNR");
axis([-40 100 0 12]);  xlabel('SNR (dB)'); ylabel('Delay (ms)');
xticks(-40:10:100); yticks(0:1:12);
grid on;
legend(delay,'Delay', 'Location', 'southeast')
hold off;

function varOutput = test(SNR)
%Setting parameters for Root Raised Cosined Transmitter Fiter
varOutput=zeros(1,2);

Nsym = 6;           % Filter span in symbol durations
beta = 0.25;         % Roll-off factor
sampsPerSym = 4;     % Upsampling factor

DataRate = 1000;
%sampling frequency
sampleFrequency = DataRate * sampsPerSym;

dataToTransmit = ones(16,1);
dataToTransmit(3:3:16)=-1;
dataLength=16;
tx = 1000 * (0: dataLength - 1) / DataRate;

%Visualize the origibal data
if (SNR == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    stem(tx,dataToTransmit);
    title("Data To Transmit");
    xlabel("Time (ms)");ylabel("Amplitude");
    axis([0 18 -2.0 2.0]);
end

%intializing RaisedCosineTransmitFilter
txfilter = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym);

%fvtool(txfilter, 'Analysis', 'impulse')

% Normalize to obtain maximum filter tap value of 1
b = coeffs(txfilter);
txfilter.Gain = 1/max(b.Numerator);
dataToTransmitmodified = [dataToTransmit;(zeros(6,1))];
modifiedDatalength = 22;
mtx = 1000 * (0: modifiedDatalength - 1) / DataRate;
% Visualize the impulse response
% fvtool(rctFilt, 'Analysis', 'impulse');

% Time vector sampled at sampling frequency in milliseconds
to = 1000 * (0: modifiedDatalength*sampsPerSym - 1) / sampleFrequency;
yo = txfilter(dataToTransmitmodified);
%AWGNoise
yp=awgn(yo,SNR,'measured');
% Plot data
sig = upsample(dataToTransmitmodified,4);

if (SNR == 100||SNR == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    Sigplot = stem(to,sig, 'mx'); hold on;
    % Plot filtered data
    FilterData = plot(to,yp, 'k--'); 
    % Set axes and labels
    title("Modified Data with RRC Output at Transmitter for SNR" + SNR);
    axis([0 23 -inf inf]);  xlabel('Time (ms)'); ylabel('Amplitude');
    grid on;
    legend([Sigplot,FilterData],'Data with delay', 'RRCFilter Output', 'Location', 'southeast')
    hold off;
end

%Setting parameters for Root Raised Cosined Filter Fiter
rxfilter  = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol', sampsPerSym, ...
  'DecimationFactor',       1);

rxfilter.Gain = 1/txfilter.Gain;

yr = rxfilter(yp);

%Downsampling the received data to get transmitted data
receivedData = downsample(yr,4);

% Plot filtered data.

if (SNR == 100||SNR == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    stem(mtx, receivedData, 'kx'); hold on;
    % hold on
    plot(to,yr, 'b--'); hold off;
    % Set axes and labels.
    title("Plot for SNR" + SNR);
    axis([0 25 -inf inf]);  xlabel('Time (ms)'); ylabel('Amplitude');
    grid on;
    legend('Received Data','Rcv Filter Output', 'Location', 'southeast');
end

if (SNR == 100||SNR == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    stem(receivedData);
    title("Extracted Data" + SNR);
    xlabel("Symbol index");ylabel("Amplitude");
    grid on;
    axis([0 25 -inf inf])
end    

% Correct for propagation delay by removing filter transients
fltDelay = Nsym / (DataRate);
yrCorreted = yr(fltDelay*sampleFrequency+1:end);
receivedDataCorrected = downsample(yrCorreted,4);
to = 1000 * (0: dataLength*sampsPerSym - 1) / sampleFrequency;

if (SNR == 100||SNR == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    plot(to,yrCorreted, 'r--'); hold on;
    stem(tx,receivedDataCorrected,'kx');
    % Set axes and labels.
    axis([0 18 -inf inf]);  xlabel('Time (ms)'); ylabel('Amplitude');
    grid on;
    title("Extracted Original Data with corrected RRC output with SNR " + SNR);
    xlabel("Time in (ms)");ylabel("Amplitude");
    legend('Corrected Rcv Filter Output','Extracted Original Data', 'Location', 'southeast');
end
 
%constructing original wave form with Transmitter output
fltDelay1 = Nsym /(2*DataRate);
orignalWaveForm1 = yo(fltDelay1*sampleFrequency+1:end);

orignalWaveForm=downsample(orignalWaveForm1,4);
outputWaveForm = downsample(yr,4);
[c,lags] = xcorr(orignalWaveForm,outputWaveForm);
c = normalize(c,'range');

if (SNR == 100||SNR == -20)
    figure();
    stem(lags,c);
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    title("Normalize Corrleation coffiecient with time delay for SNR " + SNR);
    ylabel('Corrleation coffiecient)');xlabel('Time delay');
    axis([-25 25 0 1.5]);
    grid on;
end

%Time at which coffiecient is maximum that gives delay.
[~,index1] = max(abs(c));
disp("Delay is " + abs(lags(index1)) + "ms" + " for SNR " + SNR);

d = sort(c,"descend");
diff = d(1) - d(2);

varOutput(1) = diff;
varOutput(2) = abs(lags(index1));
end



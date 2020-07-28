function delayComputationAndComparison()
differArray = zeros(3,11);
delayArray = zeros(3,11);

snrValues = -30:5:20;

%creating data for random sequence [1 1 -1 1 1 -1 ....]
dataToTransmitRN = ones(13,1);
dataToTransmitRN(3:3:13)=-1;


%creating data for Barker code sequence 
barker=comm.BarkerCode('Length',13,'SamplesPerFrame',13);
dataToTransmitBR = -1*(barker());


%creating data for PN code sequence
pnSequence = comm.PNSequence('Polynomial','x^4+x+1',... 
                'InitialConditions',[1 0 0 1],'SamplesPerFrame',13);
dataToTransmitPN = pnSequence();
dataToTransmitPN (dataToTransmitPN == 0) = -1;

%Visualize the origibal data

figure();
    ax = [0 18 -2.0 2.0];
    subplot(3,1,1);
    tx = 1000 * (0: length(dataToTransmitRN) - 1) / DataRate;
    stem(tx,dataToTransmitRN);
    title("Random [ 1 1 -1 1 -1 -1 ..] sequence (Symbols = 16)");
    xlabel("Time (ms)");ylabel("Amplitude");
    axis(ax);
    
    subplot(3,1,2);
    tx = 1000 * (0: length(dataToTransmitBR) - 1) / DataRate;
    stem(tx,dataToTransmitBR);
    title("Barker sequence (Symbols = 13)");
    xlabel("Time (ms)");ylabel("Amplitude");
    axis(ax);
    
    subplot(3,1,3);
    tx = 1000 * (0: length(dataToTransmitPN) - 1) / DataRate;
    stem(tx,dataToTransmitPN);
    title("PN sequence (Symbols = 16)");
    xlabel("Time (ms)");ylabel("Amplitude");
    axis(ax);
    



for  i  = 1:1:3
    for j = 1:1:length(snrvalues)
        if (i==1)
            outPut = computeDelay("Random [ 1 1 -1]",dataToTransmitRN,snrValues(j));
            differArray(i,j) = outPut(1);
            delayArray(i,j) = outPut(2);
            
        elseif (i==2)
            outPut = computeDelay("Barker Code",dataToTransmitBR,snrValues(j));
            differArray(i,j) = outPut(1);
            delayArray(i,j) = outPut(2);
            
        elseif (i==3)
            outPut = computeDelay("PN Code",dataToTransmitPN,snrValues(j));
            differArray(i,j) = outPut(1);
            delayArray(i,j) = outPut(2);
        end
    end
end
end


function varOutput = computeDelay(name,dTTrans,snr)
%Setting parameters for Root Raised Cosined Transmitter Fiter
varOutput=zeros(1,2);

Nsym = 6;           % Filter span in symbol durations
beta = 0.25;         % Roll-off factor
sampsPerSym = 4;     % Upsampling factor

DataRate = 1000;
%sampling frequency
sampleFrequency = DataRate * sampsPerSym;
dataToTransmit = dTTrans;
dataLength=length(dataToTransmit);
tx = 1000 * (0: dataLength - 1) / DataRate;

dataToTransmitmodified = [dataToTransmit;(zeros(6,1))];
modifiedDatalength = length(dataToTransmitmodified);

% Scaling the sequence, so that the mean power of a transmitted symbol
% is one
scalinFactor = 1/sqrt(mean(abs(dataToTransmitmodified).^2));
dataToTransmitmodified = dataToTransmitmodified*scalinFactor;

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

% Time vector sampled at symbol rate in milliseconds
mtx = 1000 * (0: modifiedDatalength - 1) / DataRate;
% Visualize the impulse response
% fvtool(rctFilt, 'Analysis', 'impulse');

% Time vector sampled at sampling frequency in milliseconds
to = 1000 * (0: modifiedDatalength*sampsPerSym - 1) / sampleFrequency;

yo = txfilter(dataToTransmitmodified);

%AWGNoise
yp=awgn(yo,snr,'measured');
% Plot data
sig = upsample(dataToTransmitmodified,4);

if (snr == 10||snr == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    Sigplot = stem(to,sig, 'mx'); hold on;
    % Plot filtered data
    FilterData = plot(to,yp, 'k--'); 
    % Set axes and labels
    title("Plot for " + name+ "sequence with SNR " + snr);
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

if (snr == 10||snr == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    stem(mtx, receivedData, 'kx'); hold on;
    % hold on
    plot(to,yr, 'b--'); hold off;
    % Set axes and labels.
    title("Received data Plot for " + name + " sequence with SNR " + snr);
    axis([0 25 -inf inf]);  xlabel('Time (ms)'); ylabel('Amplitude');
    grid on;
    legend('Received Data','Rcv Filter Output', 'Location', 'southeast');
end

if (snr == 10||snr == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    stem(receivedData);
    title("Extracted Data");
    title("Plot for SNR" + snr);
    xlabel("Symbol index");ylabel("Amplitude");
    grid on;
    axis([0 25 -inf inf])
end    

% Correct for propagation delay by removing filter transients
fltDelay = Nsym / (DataRate);
yrCorreted = yr(fltDelay*sampleFrequency+1:end);
receivedDataCorrected = downsample(yrCorreted,4);
to = 1000 * (0: dataLength*sampsPerSym - 1) / sampleFrequency;

if (snr == 100||snr == -20)
    figure();
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    plot(to,yrCorreted, 'r--'); hold on;
    stem(tx,receivedDataCorrected,'kx');
    % Set axes and labels.
    axis([0 18 -inf inf]);  xlabel('Time (ms)'); ylabel('Amplitude');
    grid on;
    title("Extracted Original Data with corrected RRC output for " + name + " with SNR" + snr);
    xlabel("Time (ms)");ylabel("Amplitude");
    legend('Corrected Rcv Filter Output','Extracted Original Data', 'Location', 'southeast');
end
%calculating delay using correlation function 

%constructing original wave form with Transmitter output
fltDelay1 = Nsym /(2*DataRate);
orignalWaveForm1 = yo(fltDelay1*sampleFrequency+1:end);
orignalWaveForm=downsample(orignalWaveForm1,4);
outputWaveForm = downsample(yr,4);
[c,lags] = xcorr(orignalWaveForm,outputWaveForm);
c = normalize(c,'range');

if (snr == 100||snr == -20)
    figure();
    stem(lags,c);
    ax=gca;
    ax.XTickMode = 'auto';
    ax.XTickMode = 'auto';
    title("Normalize Corrleation coffiecient with time delay for " + name + " sequence with SNR " + snr);
    ylabel('Corrleation coffiecient)');xlabel('Time delay');
    axis([-25 25 0 1.5]);
    grid on;
end
%Time at which coffiecient is maximum that gives delay.
[~,index1] = max(c);
disp("Delay is " + abs(lags(index1)) + "ms" + " for SNR " + snr);

d = sort(c,"descend");
diff = d(1) - d(2);
varOutput(1) = diff;
varOutput(2) = abs(lags(index1));
end
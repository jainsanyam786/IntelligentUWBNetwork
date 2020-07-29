function RMSComputationAndComparison()

global combinedDelayValue;
snrValues = -30:5:25;
interations = 1000;

combinedDelayValue = zeros(3,length(snrValues),interations);
combinedRMSE = zeros(3,length(snrValues));
combinedStdOfError = zeros(3,length(snrValues));
combinedmeanAbsError = zeros(3,length(snrValues));

DataRate = 1000;
%creating data for random sequence [1 1 -1 1 1 -1 ....]
dataToTransmitRN = ones(16,1);
dataToTransmitRN(3:3:16)=-1;

%creating data for Barker code sequence 
barker=comm.BarkerCode('Length',13,'SamplesPerFrame',13);
dataToTransmitBR = -1*(barker());

%creating data for PN code sequence
pnSequence = comm.PNSequence('Polynomial','x^4+x+1',... 
                'InitialConditions',[1 0 0 1],'SamplesPerFrame',15);
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
    title("PN sequence (Symbols = 15)");
    xlabel("Time (ms)");ylabel("Amplitude");
    axis(ax);
    

for  i  = 1:1:size(combinedDelayValue,1)
    for j = 1:1:size(combinedDelayValue,2)
        if (i==1)
            updatedelay(dataToTransmitRN,snrValues(j),i,j);
        elseif (i==2)
           updatedelay(dataToTransmitBR,snrValues(j),i,j);
        elseif (i==3)
           updatedelay(dataToTransmitPN,snrValues(j),i,j);
        end
    end
end

save("Delays","combinedDelayValue");
yhat = 6*ones(1,size(combinedDelayValue,3));
for i = 1:1:size(combinedDelayValue,1)
    for j = 1:1:size(combinedDelayValue,2)
        delay = combinedDelayValue(i,j,:);
        delay = delay(:)';
        MeanAbsError = mean(abs(delay- yhat));
        RMSE = sqrt(mean((delay- yhat).^2));
        sdOfError = std(delay- yhat);
        combinedRMSE(i,j) = RMSE;
        combinedStdOfError(i,j) = sdOfError;
        combinedmeanAbsError(i,j) = MeanAbsError;
    end
end

figure();
subplot(3,1,1)
plot(snrValues,combinedmeanAbsError(1,:),'m--*',snrValues,combinedmeanAbsError(2,:),"b-o",snrValues,combinedmeanAbsError(3,:),"k:sq");
xlabel('SNR (dB)');ylabel('Mean Absolute Error (ms)');
axis([-40 30 0 5]);
title("Mean absolute error for a Random, Barker and PN Sequence w.r.t to SNR ");
grid on;
legend("Random Sequence","Barker Sequence","PN Sequence",'Location', 'best',"Linewidth",1.5);

subplot(3,1,2)
plot(snrValues,combinedRMSE(1,:),'m--*',snrValues,combinedRMSE(2,:),"b-o",snrValues,combinedRMSE(3,:),"k:sq");
xlabel('SNR (dB)');ylabel('RMSE (ms)');
title("RMSE for a Random, Barker and PN Sequence w.r.t to SNR ");
axis([-40 30 0 5]);
xticks(-40:10:30);
grid on;
legend("Random Sequence","Barker Sequence","PN Sequence",'Location', 'best',"Linewidth",1.5);


subplot(3,1,3)
plot(snrValues,combinedStdOfError(1,:),'m--*',snrValues,combinedStdOfError(2,:),"b-o",snrValues,combinedStdOfError(3,:),"k:sq");
xlabel('SNR (dB)');ylabel('Standard Deviation (ms)');
title("Delay w.r.t to SNR ");
title("Standard Deviation for Error for a Random, Barker and PN Sequence w.r.t to SNR ");
axis([-40 30 0 5]);
xticks(-40:10:30);
grid on;
legend("Random Sequence","Barker Sequence","PN Sequence",'Location', 'best',"Linewidth",1.5);
end

function updatedelay(dTTrans,snr,i,j)
    %Setting parameters for Root Raised Cosined Transmitter Fiter
    global combinedDelayValue;

    Nsym = 6;           % Filter span in symbol durations
    beta = 0.25;         % Roll-off factor
    sampsPerSym = 4;     % Upsampling factor

    DataRate = 1000;
    %sampling frequency
    sampleFrequency = DataRate * sampsPerSym;
    dataToTransmit = dTTrans;
    dataToTransmitmodified = [dataToTransmit;(zeros(6,1))];

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

    % Normalize to obtain maximum filter tap value of 1
    b = coeffs(txfilter);
    txfilter.Gain = 1/max(b.Numerator);
    yo = txfilter(dataToTransmitmodified);

    %Setting parameters for Root Raised Cosined Filter Fiter
    rxfilter  = comm.RaisedCosineReceiveFilter(...
      'Shape',                  'Square root', ...
      'RolloffFactor',          beta, ...
      'FilterSpanInSymbols',    Nsym, ...
      'InputSamplesPerSymbol', sampsPerSym, ...
      'DecimationFactor',       1);

    rxfilter.Gain = 1/txfilter.Gain;
    for k = 1:1:size(combinedDelayValue,3)

        yp=awgn(yo,snr,'measured');
        yr = rxfilter(yp);
        %constructing original wave form with Transmitter output
        fltDelay1 = Nsym /(2*DataRate);
        orignalWaveForm1 = yo(fltDelay1*sampleFrequency+1:end);
        orignalWaveForm=downsample(orignalWaveForm1,4);
        outputWaveForm = downsample(yr,4);
        [c,lags] = xcorr(orignalWaveForm,outputWaveForm);
        c = normalize(c,'range');
        %Time at which coffiecient is maximum that gives delay.
        [~,index1] = max(c);
        combinedDelayValue(i,j,k) = abs(lags(index1));
    end
end
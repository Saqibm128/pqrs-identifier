function qrs = pt_qrs_detect(ecg,fs)
    %%initialize qrs
    qrs = [];
    %%initialize window size info
    WINDOW_SIZE_TIME = .150; %%empirically chosen window size for algorithm based on paper, 150ms
    WINDOW_SIZE = WINDOW_SIZE_TIME * fs;
    
    %%filter signal according to article
    firstFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',5,'CutoffFrequency2',11, ...
         'SampleRate',fs);
     %%from digital filter guide on MatLab website
    filteredData = firstFilt.filter(ecg);
    noiseData = ecg - filteredData;
    
    approxDiff = diff(filteredData) * fs; %% x(t + eps) - x(t) / eps
    %%eps = 1/fs
    %%this is a differentiation of the data
    squaredDiff = approxDiff .^ 2;
    %%squared differentiation
    
    allIndices = 1:(length(ecg) - WINDOW_SIZE);
    integratedSquares = zeros(length(allIndices),1);
    for i = 1:(length(ecg) - WINDOW_SIZE - 1)
        indicesToConsider = i:i+WINDOW_SIZE;
        integratedSquares(i) = trapz(squaredDiff(indicesToConsider));
    end

    
    %%  learning phase 1, get an SPK and NPK to start with
    [peaksOfIntegratedSquares, isLocs] = findpeaks(integratedSquares);
    [sortedPIS, indPIS] = sort(peaksOfIntegratedSquares);
    
    %%sort in ascending order
    noisePeak = mean(sortedPIS(1:3)); %get likely noise peaks (NPK)
    signalPeak = mean(sortedPIS(length(sortedPIS) - 3: length(sortedPIS))); %get likely signal peaks (SPK)
    
    %% learning phase 2, get an RR value, assume RR value is never less than 1/4 the fs (i.e. we dont have 240 bpm)
    

%     indicesToUse = []; %%use only indices above noise peak + signal peak / 2
%     for i = 1:length(derivIS)
%         if (integratedSquares(i) > noisePeak + .25 * (signalPeak - noisePeak))
%             indicesToUse = [indicesToUse i];
%         end
%     end
    secondFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',10,'CutoffFrequency2',15, ...
         'SampleRate',fs);
     %integratedSquares = secondFilt.filter(integratedSquares);
    dy = gradient(integratedSquares, 1);
    [peaksIS, index] = findpeaks(dy);
    qrs1 = [];
    for i = 1:length(index)
        if peaksIS(i) > mean(dy) + max(peaksIS) / 2
            qrs1 = [qrs1, (index(i) + WINDOW_SIZE/2) ./fs]; %%double difference funct use and lost window here
        end
    end
    j = 1;
    while j < length(qrs1)
        if (abs(qrs1(j) - qrs1(j + 1)) < .2) %%if qrs1 shows two lines
            qrs1(j + 1) = (qrs1(j) + qrs1(j + 1)) / 2;
            j = j + 1;
        else
            qrs = [qrs qrs1(j)];
            j = j + 1;
    end
    

%     for i = index
%         if integratedSquares(i) > noisePeak
%             qrs = [qrs i];
%         end
%     end
    

    
%     for j = 1:10
%         toRemove = [];
%         for i = 1:length(peaksOfIntegratedSquares)- 1
%             if (isLocs(i + 1) - isLocs(i) < WINDOW_SIZE) %%remove close peaks
%                 toRemove = [toRemove  i + 1];
%             end
%         end
%         adjust = 0;
%         for i = toRemove
%             isLocs(i - adjust) = [];
%             peaksOfIntegratedSquares(i - adjust) = [];
%             adjust = adjust + 1;
%         end
%     end
    

    
    
 %%stores distance between QRS R peaks
    
    
%     %%function to partition peaks data into signal and noise
%     function [noise, signal] = splitSN(peaks, mid)
%         noise = [];
%         signal = [];
%         for i = peaks
%             if (i > mid)
%                 signal = [signal i];
%             else
%                 noise = [noise i];
%             end
%         end
%     end
    
%     %%function which will be called for each window
%     function [threshold1, threshold2, SPK, NPK] = runningThreshold(peaksOfFilteredData, peaksOfNoiseData, oldSPK, oldNPK)
%         %%assume a good threshold is from the top peak
%         SPK = .125 * max(peaksOfFilteredData) + .875 * oldSPK;
%         NPK = .125 * max(peaksOfNoiseData) + .875 * oldNPK;
%         %%calculations from paper
%         threshold1 = NPK + .25 * (SPK - NPK);
%         threshold2 = .5 * threshold1;
%     end

end
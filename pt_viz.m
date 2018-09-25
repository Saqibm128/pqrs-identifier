function pt_viz(ecg,fs)
    
    % run your beat detection code
    qrs = pt_qrs_detect(ecg,fs);
    
    % plot the data
    figure;
    tt = (1:length(ecg))./fs;
    plot(tt,ecg);
    
    % plot the results
    hold on
    for k = 1:length(qrs)
        plot(qrs([k k]),[-500 500],'g-');
    end
    
end
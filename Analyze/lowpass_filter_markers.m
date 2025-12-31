function filtered_data = lowpass_filter_markers(data_3d, cutoff_freq, sampling_rate)
    % data_3d: [时间帧数, 标记点数, 3]
    % cutoff_freq: cutoff frequency(Hz)
    % sampling_rate: sampling frequancy(Hz)
    
    [~, num_markers, num_dims] = size(data_3d);
    filtered_data = zeros(size(data_3d));
    
    % Butterworth lowpass filter
    nyquist_freq = sampling_rate / 2;
    Wn = cutoff_freq / nyquist_freq;
    
    if Wn >= 1
        error('The cut-off frequency must be less than the Nyquist frequency(%.1f Hz)', nyquist_freq);
    end
    
    [b, a] = butter(4, Wn, 'low');  % 4-order Butterworth lowpass filter
    
    % Each point and Each dimension
    for marker = 1:num_markers
        for dim = 1:num_dims
            signal = squeeze(data_3d(:, marker, dim));
            
            filtered_signal = filtfilt(b, a, signal);
            
            filtered_data(:, marker, dim) = filtered_signal;
        end
    end
    
    fprintf('Low-pass filtering completed: Cut-off frequency %.1f Hz\n', cutoff_freq);
end
function [power f] = power_fft(x, sample_rate)
    data_length = length(x);
    % window = hamming(data_length);
    window = .5 - .5*cos(2*pi.*linspace(0,1,data_length));
    x_windowed = x.*window';
    x_detrend = detrend(x_windowed); % should be safe to detrend at this point
    y = fft(x_detrend);
    f = (0 : data_length - 1) * (sample_rate / data_length);     % frequency range
    power = abs(y).^2 / data_length;    % power of the DFT
end
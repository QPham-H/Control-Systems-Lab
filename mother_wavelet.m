function [ out ] = mother_wavelet( time )
out = zeros(1,length(time));
for i = 1:length(time);
    if(abs(time(i)) < 0.5)
        out(i) = sin(2*pi*(time(i) + 0.5));
    else
        out(i) = 0;
    end
end
end
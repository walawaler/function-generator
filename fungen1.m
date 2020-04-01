n = 700;
f = 450;
[sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f);
dc = 0;
waveType = 'Triangle';
a = 1;
%function fungen1(n, dc, t, a, f, w0, y, waveType)
w0 = 2*pi.*f;
t = linspace(0,2,w0.*pointsPerPeriod);
y = zeros(1,length(t));


switch waveType 
    case 'Sine'
        sinwave(n,dc,t,a,f,w0,y)
    case 'Square'
        sqr(n,dc,t,a,f,w0,y);
    case 'Triangle'
        Triangle(n,dc,t,a,f,w0,y);
    case 'Sawtooth'
        sawtooth(n,dc,t,a,f,w0,y);
    case 'Full Wave Rectified'
        rectified(n,dc,t,a,f,w0,y);
end


function sinwave(n,dc,t,a,f,w0,y)
    y = dc + a.*sin(w0.*t);
    plot(t,y); xlim([0 4./f]); ylim([min(y), max(y)]);
    [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f)
        sound(y,length(y)./max(t))
end
function sqr(n,dc,t,a,f,w0,y)
    for i = 1:2:n
        an = 4/((i.*pi));
        y =  y + an*sin(w0.*t.*i);
    end
        y = y.*a;
        y = y + dc;
        plot(t,y); xlim([0 4./f]); ylim([min(y), max(y)]);
        [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f);
        sound(y,length(y)./max(t))
end
function Triangle(n,dc,t,a,f,w0,y) 
    for i = 1:2:n
    an = 0.5 - (4./pi.^2);
        y =  y + (1./i.^2).*an*cos(w0.*i.*t);
    end
    y = y .* a;
    y = y + dc;
    plot(t,y); xlim([0 4./f]); ylim([min(y), max(y)]);
    [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f);
        sound(y,length(y)./max(t))
end 
function sawtooth(n,dc,t,a,f,w0,y) 
    for i = 1:n
    an = 0.5 - (1./pi);
        y =  y + (1./i).*an*sin(w0.*i.*t);   
    end
    y = y.*a;
    y = y + dc;
    plot(t,y); xlim([0 4./f]); ylim([min(y), max(y)]);
    [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f);
        sound(y,length(y)./max(t))
end
function rectified(n,dc,t,a,f,w0,y) 
    for i = 1:n
    an = 2./pi - 4./pi;
        y =  y + (1./(4.*i.^2-1)).*an*cos(w0.*i.*t);
    end
    y = y.*a;
    y = y + dc;
    plot(t,y); xlim([0 4./f]); ylim([min(y), max(y)]);
    [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(f);
        sound(y,length(y)./max(t))
end
function [sr, pointsPerPeriod, maxBufferSize, maxrate] = CalculateWaveformParameters(frequency)
            %This buffer size is to prevent the buffer size from getting
            %too big.  The smaller the wavelength, or in other words the
            %longer the wavelength, the larger the buffer.  With
            %frequencies like 0.0001 Hz, the buffer size starts getting
            %into the hundreds of megabytes size, so to prevent that we
            %limit it.  Observe that 8192 is a power of 2.  This is no
            %accident.
            
            maxBufferSize = 8192;
            
            %First things first let's get the hardware info and find the output sample rate.
            maxrate = 44100; %Most devices support this.  To do:  figure out way to find this number programatically future versions.
            
            %Next let's compute the number of points per period.
            %Lets make sure it's an integer power of 2.
            pointsPerPeriod = maxrate / frequency;
            pointsPerPeriod = pow2(floor(log(pointsPerPeriod) / log(2.0)));
            
            %Let's keep the buffer size from growing too big.
            %This is to control the memory consumption of the program with large wavelengths.
            if (pointsPerPeriod > maxBufferSize)
                pointsPerPeriod = maxBufferSize
            end
            
            %Also let's make sure our buffer isn't too small either.
            %Anything less than 2 samples per buffer is worthless.  It'll
            %just be a DC wave.  You will have to figure out a mechanism to
            %make sure those 2 samples land on the maxima of a waveform
            %however.
            if (pointsPerPeriod < 2)
                pointsPerPeriod = 2;
            end
            
            %Now let's compute the output sample rate
            sr = pointsPerPeriod * frequency;
            
        end

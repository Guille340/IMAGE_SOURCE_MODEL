
function R = surfaceReflection(f,theta,windspeed,met)

% f: frequency [Hz]
% theta: grazing angle of incidence [deg]
% windspeed: wind speed [m/s]

c = 1500; % speed of sound in water [m/s]
switch met
    case 'coherent' % Medwin & Clay (1998)
        h = sqrt(0.003 + 0.00512*windspeed);
        k = 2*pi*f/c;
        R = -exp(-2*(k*h*cos(theta*pi/180)).^2); % sea surface reflection factor    
    case 'incoherent' % Coates (1988)
        f = f*1e-3; % frequency [kHz]
        windspeed = windspeed * 1.94384;
        f2 = 378/windspeed^2;
        f1 = f2 * sqrt(10);
        E = 10*log10((1 + (f/f1).^2)./(1 + (f/f2).^2)) -(1 + (90 - windspeed)/60) * (theta/30)^2; % Beckmann-Spizzichino loss coefficient
        R = -10.^(E/20); % sea surface reflection factor       
    otherwise
        error('Input argument MET not recognised')
end
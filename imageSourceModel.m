%  tl = IMAGESOURCEMODEL(xr,zr,zs,f,hw,theta,Rb,Rs,varargin)
% 
%  DESCRIPTION
%  Calculates the complex transmission loss along a horizontal transect (XR,ZR) 
%  from a point source for a single frequency F using the Image Source method. 
%  The model accepts environmental parameters such as sound speed C, depth of 
%  water column HW, surface reflection factor RS, and bottom reflection factor 
%  RB. The transmission loss can be calculated for one or more reflection 
%  orders (higher order = higher accuracy).
%
%  INPUT VARIABLES
%  - xr: vector of receiver ranges [m]
%  - zr: receiver depth [m]
%  - zs: source depth [m]
%  - f: frequency [Hz]
%  - hw: water column depth [m]
%  - theta: vector of elevation angles (ref. normal). Ignored if RB and RS
%    are 1-element vectors [deg]
%  - Rb: vector of bottom reflection factors at THETA. Each value must be a 
%    complex number with absolute value between 0 and 1. RB can be a 1-element
%    vector if the reflection factor is constant with THETA.
%  - Rs: vector of surface reflection factors at THETA. Each value must be a 
%    complex number with absolute value between 0 and 1. RS can be a 1-element
%    vector if the reflection factor is constant with THETA.If information 
%    about the surface state is lacking, a value of Rs = -1 still provides 
%    accurate results.
%
%  VARIABLE INPUT ARGUMENTS
%  Specified in ('string',value) pairs within the function call. Type any
%  of the following strings followed by a suitable value:
%  - 'order': order of reflection group. The value must be an integer
%    higher or equal than 0. Each reflection group consists of four 
%    reflections (b... + bs... + sb... + bsb..., b = bottom, s = surface),
%    except group 0 that only has two (direct + ghost). The total number 
%    of simulated rays for a given order is nref = 4*order + 2.
%  - 'soundspeed': speed of sound in water [m/s]. If unknown, a value of 
%    1500 m/s will generally provide accurate results.
%
%  OUTPUT VARIABLES
%  - tl: transmission loss vector [-]. If more than one order is specified,
%    tl is an array where each column corresponds to each given order.
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  1) tl = imageSourceModel(f,hw,xr,zr,zs,theta,Rb,Rs)
%      � 'order'-> ordref = 100, 'soundspeed' -> c = 1500 m/s
%  2) tl = imageSourceModel(...,STR,VAL)
%      � ...,'order',ordref
%      � ...,'soundspeed',c
%
%  CONSIDERATIONS & LIMITATIONS
%  - Vincenty formula is slower than Haversine but accurate at every 
%    distance (1 mm error withing considered elliptic earth model)
%
%  REFERENCES
%  - Brekhovskikh & Lysanov (1990). "Propagation of sound in shallow
%    water", from Fundamentals of Ocean Acoustics, 2nd ed.
%  - Lepper et al. (2014). "Establishing the sensitivity of cetaceans
%    and seals to acoutic deterrent devices in Scotland", section 2.3.2,
%    p. 23-25, Introduction to Source Image Model (ImageTL).
%  - Robinson et al. (2011). "Measurement of underwater noise arising from
%    marine aggregate dredging operations", Appendix C, p. 103-104.
%
%  See also lloydMirror.m, surfaceReflection.m

%  VERSION 2.0 from 14 Mar 2021
%  - Added input argument THETA and updated RB and RS into vectors to
%    accept angle-dependent reflection factors. RB and RS should always
%    input as angle-dependent for more realistic results, since sound
%    arriving at the seabed at angles close to the normal experience large
%    attenuations.
%
%  VERSION 1.1 from 13 Mar 2021
%  - Simple update to accept a vector of receiver depths. Still a single 
%    receiver depth can be specified (back-compatibility with ver 1.1)
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  27 Feb 2019

function tl = imageSourceModel(xr,zr,zs,f,hw,theta,Rb,Rs,varargin)

if rem(nargin - 8,2) 
    error('Name/value pair arguments must come in pairs')
else
    ordref = 100; % bottom reflection order (no. boundary reflections = 4 * ordref + 1)
    c = 1500; % speed of sound in water [m/s]
    m = 1;
    while m <= nargin - 8
       if ~ischar(varargin{m}), error('Expected string parameter'); end
       switch varargin{m}               
           case 'order'
               ordref = varargin{m+1}; % if a vector, the pressure output is an array of the form [M length(ordref)]
           case 'soundspeed'
               c = varargin{m+1};
           otherwise
              error('Unknown parameter string');
       end
       m = m + 2;
    end
end

M = length(xr); % number of points in range vector;
xr = xr(:); % convert range vector into column
zr = zr(:); % convert receiver depth vector into column 
theta = theta(:)'; % convert vector of elevation angles into row
k = 2*pi*f/c; % wave number [rad/m]
nbot = max(ordref); % number of reflection groups
n = 1:nbot; % vector of number of reflection groups
Mb = round(250 * 1024^2 / ((nbot+1) * 8 * 12)); % number of range points to include in each processing block
Q = ceil(M/Mb); % number of processing blocks

hwb = waitbar(0,'Processing range-dependent acoustic pressure ...','Name','Image Source Model');
p = zeros(M,nbot+1);
for q = 1:Q
    
    % Progress Bar
    hwb = waitbar((q-1)/Q,hwb,'Processing range-dependent acoustic pressure ...');

    % Extract Range Block
    i1 = Mb*(q-1) + 1;
    if q < Q, i2 = Mb*q; else, i2 = M; end
    xb = xr(i1:i2);
    
    % Calculate Travelled Distances
    r1 = sqrt(xb.^2 + (zr - zs).^2);
    r2 = sqrt(xb.^2 + (zr + zs).^2);
    r1n = sqrt(xb.^2 + (2*n*hw - zr - zs).^2);
    r2n = sqrt(xb.^2 + (2*n*hw + zr - zs).^2);
    r3n = sqrt(xb.^2 + (2*n*hw - zr + zs).^2);
    r4n = sqrt(xb.^2 + (2*n*hw + zr + zs).^2);
    
    % Calculate Seabed Incidence Angles
    theta1n = atan(xb./(2*n*hw - zr - zs)) * 180/pi;
    theta2n = atan(xb./(2*n*hw + zr - zs)) * 180/pi;
    theta3n = atan(xb./(2*n*hw - zr + zs)) * 180/pi;
    theta4n = atan(xb./(2*n*hw + zr + zs)) * 180/pi;
    
    % Calculate Bottom Reflection Coefficients
    if numel(Rb) == 1
        Rb1n = Rb;
        Rb2n = Rb;
        Rb3n = Rb;
        Rb4n = Rb;
    else
        Rb1n = interp1(theta,Rb,theta1n);
        Rb2n = interp1(theta,Rb,theta2n);
        Rb3n = interp1(theta,Rb,theta3n);
        Rb4n = interp1(theta,Rb,theta4n);
    end
    
    % Calculate Surface Reflection Coefficients
    if numel(Rs) == 1
        Rs1n = Rs;
        Rs2n = Rs;
        Rs3n = Rs;
        Rs4n = Rs;
    else
        Rs1n = interp1(theta,Rs,theta1n);
        Rs2n = interp1(theta,Rs,theta2n);
        Rs3n = interp1(theta,Rs,theta3n);
        Rs4n = interp1(theta,Rs,theta4n);
    end
    
    % Calculate Sound Pressures
    pd = exp(1j*k*r1)./r1; % pressure from direct sound [Pa]
    pg = Rs * exp(1j*k*r2)./r2; % pressure from ghost reflection [Pa]#
    if nbot % if max reflection order is different than 0
        pr1 = Rs1n.^(n-1) .* exp(1j*k*r1n)./r1n .* Rb1n.^n; % pressure from reflections of order 2*n - 1 (bottom first)
        pr2 = Rs2n.^n .* exp(1j*k*r2n)./r2n .* Rb2n.^n; % pressure from reflections of order 2*n (surface first)
        pr3 = Rs3n.^n .* exp(1j*k*r3n)./r3n .* Rb3n.^n; % pressure from reflections of order 2*n (bottom first)
        pr4 = Rs4n.^(n+1) .* exp(1j*k*r4n)./r4n .* Rb4n.^n; % pressure from reflections of order 2*n + 1 (surface first)
    else
        pr1 = [];
        pr2 = [];
        pr3 = [];
        pr4 = [];
    end
    
    % Total Pressure with Range (All Reflection Orders)
    p(i1:i2,:) = cumsum([pd + pg , pr1 + pr2 + pr3 + pr4],2);
    
    % Clear variables
    clear r1 r2 r1n r2n r3n r4n
    clear pd pg pr1 pr2 pr3 pr4
end
hwb = waitbar(1,hwb,'Processing complete!');
pause(0.2)
close(hwb)

% Total Pressure with Range (Selected Reflection Orders)
p = p(:,ordref+1); % complex received sound pressure (source pressure = 1 Pa) [Pa]
tl = 1./p; % complex transmission loss (= source pressure / received pressure) [-]

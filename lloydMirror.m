%  tl = LLOYDMIRROR(r,zr,zs,f)
%
%  DESCRIPTION: calculates the pressure transmission loss at frequency F at the 
%  ranges in vector R, for a source and receiver at depths ZS and ZR.
%  
%  INPUT VARIABLES
%  - r: vector of ranges for TL [m]
%  - zr: receiver depth [m]
%  - zs: source depth [m]
%  - f: frequency [Hz]
%  - form: formula used for calculations. Two options:
%    ¬ 'exact': exact formulation (DEFAULT)
%    ¬ 'short': approximate formulation.
%
%  OUTPUT VARIABLES
%  - tl: transmission loss at frequency F in a semi-infinite water environment.
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - This is a very powerful method for simplifying transmission loss
%    calculations within a frequency band, but it has its limitations. After
%    an in-depth analysis of its performance I found out that the accuracy
%    of the approximation is closely related to the range resolution and
%    the characteristics of the averaging window. The main conclusions I 
%    arrived at for evenly-spaced range points (relevant for AcTUP and
%    most modelling scenarios) were:
%
%    1. The range resolution step dr must be dr <= rmin/(37.5*bpo) to get
%    an error at range of less than 0.3 dB, where rmin is the minimum range
%    that meets the error condition, and bpo is the number of bands per
%    octave (i.e. 3 indicates an equivalent frequency average over a
%    third-octave band).
%
%    2. An "overfit" approximation may occur when the relative averaging band
%    is too narrow. To avoid overfitting, the number of bands per octave must
%    be bpo < fc * rmin * 1e-5. A larger bpo is preferred to a small
%    value that may result in "oversmoothing" (see next point).
%
%    3. An "oversmoothed" approximation will occur when the relative averaging 
%    band is too wide. To avoid oversmoothing, the number of bands per octave 
%    must be bpo >= 3 to bring the standard error down to 0.5 dB or less.
%
%  See also tlavg_ex2.m, tlavg.m
%  

%  REVISION 1.2 [17 Feb 2021]
%  - Added option to choose between approximate (old) and exact formulation.
%
%  REVISION 1.1 [02 Feb 2021]
%  - Added function help
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  12 May 2020

function tl = lloydMirror(xr,zr,zs,f,varargin)
    
form = 'exact';
if nargin == 5
    form = varargin{1};
end

xr = xr(:);
f = f(:)';
c = 1500; % speed of sound in water [m/s]
k = 2*pi*f/c; % wave number [rad/m]

r1 = sqrt(xr.^2 + (zr - zs)^2); % direct travelled distance [m]
r2 = sqrt(xr.^2 + (zr + zs)^2); % reflected travelled distance [m]

switch form
    case 'exact'
        pd = exp(1j*k*r1)./r1; % pressure from direct sound [Pa]
        pg = -exp(1j*k*r2)./r2; % pressure from ghost reflection [Pa]
        p = pd + pg; % total pressure [Pa]
        tl = 1./p; % transmission loss [-]
    case 'short'
        tl = r./(2*sin(k*zr*zs./r)) .* exp(-1j*k.*(r1+r2)/2); % trans. loss [-]
end

    
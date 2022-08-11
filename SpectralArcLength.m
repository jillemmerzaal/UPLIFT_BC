function sparc = SpectralArcLength(mov,fsamp,params)
%
% SPECTRALARCLENGTH(acc,fsamp,parameters) computes the smoothness of a given
% movement profile using the spectral arc length method.
%
% ARGUMENTS
%  mov    ... Movement profile over time (vector)
%  fsamp  ... Sampling rate of the signal. 
%  params ... This contains the parameters to be used for spectral arc 
%             length computation (optional). This is a 1x2 column vector. 
%           - PARAMETER(1): The amplitude threshold to be used to choose
%           the cut-off frequency. The default value is chosen to be 0.05.
%           - PARAMETER(2): Maximum cut-off frequency for the spectral arc
%           length calcualtion. (DEFAULT VALUE = 10HZ) NOTE: 20Hz IS USED
%           TO REPRESENT THE MAXIMUM FREQUENCY COMPONENT OF A MOVEMENT.
%           THIS WAS CHOSEN TO COVER BOTH NORMAL AND ABNORMAL MOTOR
%           BEAHVIOUR. YOU CAN USE A VALUE LOWER THAN 20Hz IF YOU ARE AWARE
%           OF THE MAXIMUM FREQUENCY COMPONENT IN THE MOVEMENT OF INTEREST.
%           - PARAMETER(3): Zero padding index. This parameter controls the
%           resolution of the movement spectrum calculated from the speed
%           profile. (DEFAULT VALUE = 4). NOTE: IT IS NOT ADVISABLE TO USE
%           VALUES LOWER THAN 4.
%
% RETURN
%  sparc  ... This is smoothness of the given movement.
%
%
% For any queries about the method or the code, or if you come across any
% bugs in the code, feel free to contact me at siva82kb@gmail.com
% Sivakumar Balasubramanian. July 02, 2014.


% Revision history
%
% Updated by Philippe C. Dixon May 2016
% - This code computes the updated SPARC as described by Balasubraminian 2015
%   Checked vs output of python code here: https://github.com/siva82kb/SPARC


% Error checking / Set defaults 
%
if nargin == 0                          % Phil D test using python data
    mov = getMovement;
    fsamp = 100;
    params = [0.05, 10, 4];             % default parameter settings
elseif nargin == 1
    fsamp = 100; 
    params = [0.05, 10, 4];
elseif nargin == 2                     
    params = [0.05, 10, 4];
end



% Check if the input argument are of the appropriate dimensions
%
sz = size(mov);
if (sz(1) == 1) && (sz(2) > 1)
    mov = mov';
end


% Parameters
%
if length(params) ~= 3
    error('parameter is a vector with two elements.');
end



% Calculate the spectrum of the speed profile
%
%N = length(movement);
%Nfft = 2^(ceil(log2(N))+parameters(3))                       % number of zero pads
Nfft = 4096;                                                  % always 4096
f = 0:fsamp/Nfft:fsamp-fsamp/Nfft;

% Normalize spectrum with respect to the DC component.
%
Mf = abs(fft( mov, Nfft ));                                   % V(w)
Mf = Mf'/max(Mf);                                             % Vhat(w) = V(w)/V(0)

% Compute sal without any filtering (original approach)
%
freq = 0:(1/fsamp)*(1/Nfft):(1/fsamp)*((Nfft-1)/Nfft);
inxFc = find( freq(1:end) <= params(1), 1, 'last' );
dArcLengths = sqrt(  (1/(inxFc-1))^2 + diff(Mf).^2);
sal = -sum(dArcLengths);


% Indices to choose only the spectrum within the given cut off frequency Fc.
% NOTE: This is a low pass filtering operation to get rid of high frequency
% noise from affecting the next step (amplitude threshold based cut off for
% arc length calculation).
fc = params(2);
fc_inx = find(f <= fc);                                       % same as Python
Mf_sel = Mf(fc_inx);


% Calculates centroid of power spectrum
%
Mdf = sum( f(fc_inx).*Mf_sel) / sum(Mf_sel);

if nargin==0
    figure
    plot(f(fc_inx),Mf_sel)
    vline(Mdf,'k')
end

% Choose the amplitude threshold based cut off frequency.
% Index of the last point on the magnitude spectrum that is greater than
% or equal to the amplitude threshold.
amp_th = params(1);
inx = find(Mf_sel >=amp_th);
fc_inx = inx(1):inx(end);
Mf_sel = Mf_sel(fc_inx);


% Calculate the arc length and smoothness (updated version)
%
inxFc = fc_inx(end);     % this is Kc from Bala 2012
dArcLengths = sqrt( ( 1/(inxFc-1) )^2 + diff(Mf_sel).^2);

sparc = -sum(dArcLengths);

if nargin==0
    hline(params(1))
    x = f(fc_inx(end));
    y =params(1)+0.02;
    text(x,y,'Vbar')
    text(5,0.5,['SPARC = ',num2str(sparc)])
end

function movement = getMovement

% data from the python example to use as test

movement = [0.00673795,  0.00744286,  0.0082133 ,  0.00905444,  0.00997174,...
    0.010971  ,  0.01205832,  0.01324017,  0.01452331,  0.01591489,...
    0.01742237,  0.01905359,  0.02081669,  0.02272022,  0.02477302,...
    0.0269843 ,  0.02936358,  0.03192072,  0.03466586,  0.03760945,...
    0.0407622 ,  0.0441351 ,  0.04773932,  0.05158626,  0.05568748,...
    0.06005467,  0.06469962,  0.06963416,  0.07487015,  0.08041939,...
    0.08629359,  0.09250431,  0.09906293,  0.10598052,  0.11326784,...
    0.12093525,  0.12899263,  0.13744932,  0.14631404,  0.15559481,...
    0.16529889,  0.17543266,  0.1860016 ,  0.19701016,  0.20846169,...
    0.22035839,  0.23270121,  0.24548977,  0.2587223 ,  0.27239556,...
    0.2865048 ,  0.30104365,  0.31600413,  0.33137653,  0.34714942,...
    0.36330957,  0.37984196,  0.39672973,  0.41395417,  0.43149472,...
    0.44932896,  0.46743265,  0.48577972,  0.50434234,  0.52309091,...
    0.54199419,  0.56101928,  0.58013178,  0.59929579,  0.61847408,...
    0.63762815,  0.65671838,  0.67570411,  0.69454384,  0.71319529,...
    0.73161563,  0.74976159,  0.76758965,  0.78505618,  0.80211764,...
    0.81873075,  0.83485268,  0.8504412 ,  0.86545491,  0.87985338,...
    0.89359735,  0.9066489 ,  0.91897166,  0.9305309 ,  0.94129377,...
    0.95122942,  0.96030916,  0.96850658,  0.97579769,  0.98216103,...
    0.9875778 ,  0.99203191,  0.99551011,  0.998002  ,  0.99950012,...
    1.0       ,  0.99950012,  0.998002  ,  0.99551011,  0.99203191,...
    0.9875778 ,  0.98216103,  0.97579769,  0.96850658,  0.96030916,...
    0.95122942,  0.94129377,  0.9305309 ,  0.91897166,  0.9066489 ,...
    0.89359735,  0.87985338,  0.86545491,  0.8504412 ,  0.83485268,...
    0.81873075,  0.80211764,  0.78505618,  0.76758965,  0.74976159,...
    0.73161563,  0.71319529,  0.69454384,  0.67570411,  0.65671838,...
    0.63762815,  0.61847408,  0.59929579,  0.58013178,  0.56101928,...
    0.54199419,  0.52309091,  0.50434234,  0.48577972,  0.46743265,...
    0.44932896,  0.43149472,  0.41395417,  0.39672973,  0.37984196,...
    0.36330957,  0.34714942,  0.33137653,  0.31600413,  0.30104365,...
    0.2865048 ,  0.27239556,  0.2587223 ,  0.24548977,  0.23270121,...
    0.22035839,  0.20846169,  0.19701016,  0.1860016 ,  0.17543266,...
    0.16529889,  0.15559481,  0.14631404,  0.13744932,  0.12899263,...
    0.12093525,  0.11326784,  0.10598052,  0.09906293,  0.09250431,...
    0.08629359,  0.08041939,  0.07487015,  0.06963416,  0.06469962,...
    0.06005467,  0.05568748,  0.05158626,  0.04773932,  0.0441351 ,...
    0.0407622 ,  0.03760945,  0.03466586,  0.03192072,  0.02936358,...
    0.0269843 ,  0.02477302,  0.02272022,  0.02081669,  0.01905359,...
    0.01742237,  0.01591489,  0.01452331,  0.01324017,  0.01205832,...
    0.010971  ,  0.00997174,  0.00905444,  0.0082133 ,  0.00744286];




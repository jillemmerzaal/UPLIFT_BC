function [lyapExp,eLag, eDim] = DivergenceExponent(df, fs)

%%
% 
% DivergenceExponent calculates the lyapunov exponent
% 
% # Selects the relevant parameters signal and creates a time axis.
% # Calculates the power spectrum to 
% # Determines the dominant frequency of the signal.
% # Determines the range over which the Lyapunov Exponent is calculated (i.e. half a cycle)
% # Reconstructs the number of state spaces and timelag that fully capture the signal
% # Calculates the Lyapunov Exponent with the predetermined
%   parameters. 

%%
% 
%  written by:  JILL EMMERZAAL, Ph.D. Internal rehabilitation KU Leuven
%  Version:     v1.0
%  Date:        27/09/2022 


% 1) relevant parameters
signal = df;
N = length(signal);
k = (0:N-1);
dt = 1/fs;
f = k*(1/(N*dt));

% 2) Fourir analysis to calculate power
Y = fft(signal);
Power = abs(Y)/N;

% 3) determine dominant frequency of the signal
tf = islocalmax(Power(1:(N/2)), 'MaxNumExtrema',1);
T2 = f(tf);

% 4) calculate the range over which the lyapunov exponent is calculated
L1 = round(0.5*T2*fs);
eRange=[0, L1];

% 5) reconstruction of the state spaces
[~,eLag,eDim] = phaseSpaceReconstruction(signal);

%calculate the lyapunov exponent
lyapExp = lyapunovExponent(signal,fs, eLag, eDim, 'MinSeparation', ceil(T2), 'ExpansionRange',eRange);


end
function [lyapExp,eLag, eDim] = DivergenceExponent(df, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% relevant parameters
signal = df;
N = length(signal);
k = [0:N-1];
dt = 1/fs;
f = k*(1/(N*dt));

% Fourir analysis to calculate power
Y = fft(signal);
Power = abs(Y)/N;

% determine dominant frequency of the signal
tf = islocalmax(Power(1:(N/2)), 'MaxNumExtrema',1);
T2 = f(tf);

%calculate the range over which the lyapunov exponent is calculated
L1 = round(0.5*T2*fs); 
eRange=[0, L1];

% reconstruction of the state spaces
[~,eLag,eDim] = phaseSpaceReconstruction(signal);

%calculate the lyapunov exponent
lyapExp = lyapunovExponent(signal,fs, eLag, eDim, 'MinSeparation', ceil(T2), 'ExpansionRange',eRange);




end
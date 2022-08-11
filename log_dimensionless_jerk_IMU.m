function [ldlj] = log_dimensionless_jerk_IMU(accls,t, fs)
%log_dimensionless_jerk_IMU calculates the smoothness of movement directly
%                           from acceleration signals.

% Input parameters
%   -acc:       nx3 matrix from an acceleration signal gravity subtracted
%               and rotated to the earth coordinate system. and mean
%               subtracted (i.e. the mean of the signals should be around zero!)
%               this reduces the amount of overestimation of the peak
%               acceleration
%	-t:         1x2 vector of start and stop times

% Output variables 
%   - jerk:     The smoothness of the acceleration signal

% Script information 
%   written by: Jill Emmerzaal, PhD candidate Human movement biomechanics
%               research group KU Leuven
%   version:    v1.0
%   date:       20-04-2021

% Reference:
%   Melendez-Calderon, A., Shirota, C., & Balasubramanian, S. (2020). 
%   Estimating Movement Smoothness from Inertial Measurement Units.

%-------------------------------------------------------------------------%

%seperate start and stop time
t1 = t(1);
t2 = t(2);

%Sample time
dt = 1 / fs;
N = length(accls(t1:t2));

%Movement duration.
mdur = N * dt;

%mean subtracted acceleration
accms = accls(t1:t2,:) - (1/(t2-t1)) * trapz(accls(t1:t2,:));

%Gravity subtracted mean square ampitude
mamp = (norm(accms, 'fro').^2)/N;

%Derivative of the accelerometer signal
daccls = diff(accls(t1:t2,:));

% scaling factor peak absolute acceleration (Euclidische norm)
a_peak = norm(accms, 'fro');


% a_peak = max(abs(accms)); 

%Derivative of the accelerometer signal
daccls = [0,0,0; (diff(accms) * fs)]';

%Corrected jerk
mjerk = sum(vecnorm(daccls,2).^2)*dt;

% factors for the LDLJ calculation
f = [-log(mdur), log(mamp), -log(mjerk)];


% LDLJ calculation
ldlj = f(1) + f(2) + f(3); 

end



% def log_dimensionless_jerk_imu_factors(accls, gyros, grav, fs):
%     """
%     Returns the individual factors of the log dimensionless jerk metric
%     used for IMU data.
%     Parameters
%     ----------
%     accls : np.array
%             The array containing the accelerometer profile. This is a
%             multi-dimensional with the rows corresponding to the time samples
%             and the columns corresponding to the x, y, and z components.
%     gyros : np.array
%             The array containing the gyroscope profile. This is 
%             multi-dimensional with the rows corresponding to the time samples
%             and the columns corresponding to the x, y, and z components.
%     grav  : np.array
%             Gravity vector. This is a 3x1 array with the x, y, and z component
%             of gravity. 
%     fs    : float
%             The sampling frequency of the data.
%     Returns
%     -------
%     -ln(T) : float
%                Duration scaling factor.
%     +ln(A) : float
%                Amplitude scaling factor.
%     -ln(J)   : float
%                Jerk cost.
%     Notes
%     -----
%     Examples
%     --------
%     """
%     # Sample time
%     dt = 1. / fs
%     _N = len(accls)
% 
%     # Movement duration.
%     mdur = _N * dt
% 
%     # Gravity subtracted mean square ampitude
%     mamp = np.power(np.linalg.norm(accls), 2) / _N
%     if gyros is not None:
%       mamp = mamp - np.power(np.linalg.norm(grav), 2)
%     
%     # Derivative of the accelerometer signal
%     _daccls = np.vstack((np.zeros((1, 3)), np.diff(accls, axis=0) * fs)).T
%     
%     # Get corrected jerk if gyroscope data is available.
%     if gyros is not None:
%       _awcross = np.array([np.cross(_as, _ws)
%                            for _as, _ws in zip(accls, gyros)]).T
%     else:
%       _awcross = np.zeros(np.shape(_daccls))
%     
%     # Corrected jerk
%     _jsc = _daccls - _awcross
%     mjerk = np.sum(np.power(np.linalg.norm(_jsc, axis=0), 2)) * dt
% 
%     return - np.log(mdur), np.log(mamp), - np.log(mjerk)

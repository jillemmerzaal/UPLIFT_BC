function [acorr_norm, symm] = Symmetry(df)
%SYMMETRY will caclulate the regularity/symmetry of a cyclical movement
%pattern using the unbiassed autocorrelation. 
%
%   Input: signal 
%   Output1: normalised autocorraltion signal 
%   Output2: normalised magnitude of the first dominant peak after 
%   zerophase to denote the measure of regularity of the signal


acorr_df = xcov(df, 'unbiased');

[pos_pks, pos_locs] = findpeaks(acorr_df, 'MinPeakHeight', 0, 'MinPeakDistance',99);
nPeaks = numel(pos_locs);
TT = table(pos_locs,pos_pks);
[~,I] = max(TT.pos_pks(10:height(TT)-10));
I = I + 9;
zerophase_i = TT.pos_locs(I);

acorr_norm = acorr_df./acorr_df(zerophase_i);

symm_idx = TT.pos_locs(I + 1);
symm = acorr_norm(symm_idx); 



end
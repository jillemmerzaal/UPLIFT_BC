function [acorr_norm, symm] = Symmetry(df, plot_or_not, name)
%SYMMETRY will caclulate the regularity/symmetry of a cyclical movement
%pattern using the unbiassed autocorrelation. 
%
%   Input: signal 
%   Output1: normalised autocorraltion signal 
%   Output2: normalised magnitude of the first dominant peak after 
%   zerophase to denote the measure of regularity of the signal


acorr_df = xcov(df, 'unbiased');


[pos_pks, pos_locs] = findpeaks(acorr_df, 'MinPeakHeight', 0, 'MinPeakDistance',60);
%nPeaks = numel(pos_locs);
TT = table(pos_locs,pos_pks);
[~,I] = max(TT.pos_pks(5:height(TT)-5));
I = I + 4;
zerophase_i = TT.pos_locs(I);

acorr_norm = acorr_df./acorr_df(zerophase_i);


[~,I_symm] = max(TT.pos_pks(I+1:I+3));

symm_idx = TT.pos_locs(I + I_symm);
symm = acorr_norm(symm_idx); 

if plot_or_not
%     nexttile
%     plot(acorr_norm)
%     hold on
%     plot(pos_locs, acorr_norm(pos_locs), 'o')
%     plot([zerophase_i, symm_idx], [acorr_norm(zerophase_i) ,symm], 'b*')
%     title(name)
end



end
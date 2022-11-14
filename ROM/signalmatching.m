function [iStart,iStop, signal] = signalmatching(df,extrema, plot_or_not, subj_name, fileName)
%%
%
%  We find the peak magnitude, and based on that mannitude
%  we scale sine signal to match the height of the ROM
%  repetition
%  After that, we match the sine signal with the full ROM
%  signal and match the individual repetitions.
%


[~, peakMag] = peakfinder(df, [],[],extrema,false);
[~, peakZero] = peakfinder(df, [], [], extrema*-1, false);

iStart = zeros(size(peakMag));
iStop = zeros(size(peakMag));

temp = df;

for idx = 1:size(peakMag,1)
    if extrema <0
        amplidude = abs(peakMag(idx)) + mean(peakZero);
    elseif extrema > 0
        amplidude = abs(peakMag(idx)) - mean(peakZero);
    end

    shift =(max(peakZero));
    t=0:0.01:1;
    f=1;

    x = (sin(pi*f*t) * amplidude * extrema) + shift;

    [iStart(idx,1), iStop(idx,1)] = findsignal(temp, x,'TimeAlignment','dtw','Metric','absolute');

    temp(iStart(idx,1):iStop(idx,1)) = 0;
    signal(idx,:) = x;
end

if plot_or_not
    figure;
    t = tiledlayout(1,3, 'TileSpacing','Compact');

    ax1 = nexttile;
    plot(signal', 'o-', 'MarkerSize',2, 'MarkerFaceColor','auto')
    title('Signal')

    ax2 = nexttile([1 2]);
    plot(df); hold on
    for idx = 1:size(iStart,1)
        plot(iStart(idx):iStop(idx), df(iStart(idx):iStop(idx)), 'o', 'MarkerSize',2)
        title('Signal matches found')
    end

    linkaxes([ax1 ax2],'y')
    sgtitle(['Subject: ',subj_name, ' File: ', fileName])
end



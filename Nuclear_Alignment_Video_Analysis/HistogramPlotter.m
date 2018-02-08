function [FinalStd, MaxAngle] = HistogramPlotter(AngleVector, WrinkleAngle, plot_toggle)

% This makes a histogram with angles rotated so that the cell angles are
% centered around the wrinkle angle (wrinkle angle = 90 degress in plot)
HistogramVector=WrinkleNorm(AngleVector,WrinkleAngle);
HistogramVectorRad = HistogramVector*pi()/180;

if plot_toggle == 1
    figure(1);
    [t,r] = rose(HistogramVectorRad,30);
    r = r./max(r);
    polar(t,r,'-');
end

% This makes a histogram with angles rotated so that the cell angles are
% centered around the angle of max orientation (which should be close to
% the wrinkle angle)
[FinalStd, MaxAngle, HistogramMaxVector] = AngleSpread2(AngleVector);
HistogramVectorRad = HistogramMaxVector*pi()/180;
if plot_toggle == 1
    figure(2);
    [t,r] = rose(HistogramVectorRad,30);
    r = r./max(r);
    polar(t,r,'-');
end

% FinalStd tells you the overall degree of alignment (52 = completely
% random, as you get get lower than 52 you increase orientation)

% MaxAngle tells you the angle that the cells are most oriented to



function [FinalStd, MaxAngle, HistogramMaxVector]  = AngleSpread2(AngleVector)

%  Determines the minimum standard deviation of an angle vector by rotating
%  by 1 degree at a time for 180 degrees
%
%
AngleVector = AngleVector(~isnan(AngleVector));
angles = zeros(size(AngleVector,1),180);
%

StdVector = zeros(180,1);
MinStd = 1000000;
NonTruncStdVector = zeros(180,1);

for i = 1:180;
    AlignedAngleVector = AngleAlign(AngleVector);
    AlignedAngleVector = AlignedAngleVector - 90;
    if max(AlignedAngleVector) > 90
        disp('error max angle')
        continue
    end
    if min(AlignedAngleVector) < -90
        disp('error min angle')
    end
      
    NonTruncStdVector(i) = (sum(AlignedAngleVector.^2)/size(AlignedAngleVector,1))^0.5;
    StdVector(i) = 52/(1+543*NonTruncStdVector(i)^(-1.96));
    
    if StdVector(i) < MinStd;
        MaxAngleVector = AngleVector;
        MinStd = StdVector(i);
        MaxAngle = i;
    end;
    angles(:,i) = AlignedAngleVector;
    AngleVector = AngleVector - 1;
    
end;
%
%
ModMaxAngleVector = AngleAlign(MaxAngleVector);
%HistogramMaxVector = ModMaxAngleVector - 90;
HistogramMaxVector = ModMaxAngleVector;
FinalStd = std(HistogramMaxVector);
% 
%**************************************************************************

%{
AngleVector = AngleVectorCopy(~isnan(AngleVectorCopy));
angles = zeros(size(AngleVector,1),180);
f = figure
clear M
for i = 1:180;         
    AlignedAngleVector = AngleAlign(AngleVector);
    AlignedAngleVector = AlignedAngleVector - 90;
    angles(:,i) = AlignedAngleVector;
    AngleVector = AngleVector - 1;
    rose(AlignedAngleVector*(pi/180))
    M(i) = getframe(f);
    clf;
end;
%
%}
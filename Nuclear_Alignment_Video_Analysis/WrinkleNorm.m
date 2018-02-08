function HistogramVector = WrinkleNorm(AngleVector,WrinkleAngle);

AngleVector = AngleVector(~isnan(AngleVector));

x = size(AngleVector);

AngleVector = AngleVector-WrinkleAngle;

for y = 1:x;
    if AngleVector(y) < 0;
        AngleVectorMod(y) = AngleVector(y) + 360;
    else AngleVectorMod(y) = AngleVector(y);
    end;

    % Make all angles between 0 and 180
    if AngleVectorMod(y) > 180
        AngleVectorMod(y) = AngleVectorMod(y) - 180;
    end

    if AngleVectorMod(y) < 90;
        AngleVectorMod2(y) = 90 - AngleVectorMod(y);
    else AngleVectorMod2(y) = AngleVectorMod(y)*(-1) + 270;
    end;
end
    

HistogramVector=AngleVectorMod2;

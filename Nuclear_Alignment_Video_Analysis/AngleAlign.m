function AngleVectorMod2 = AngleAlign(AngleVector)

x = length(AngleVector);

AngleVectorMod = zeros(x(:,1),1);
AngleVectorMod2 = zeros(x(:,1),1);

    for y = 1:x;
        % Make all angles positive
        if AngleVector(y) < 0;
            AngleVectorMod(y) = AngleVector(y) + 360;
        else AngleVectorMod(y) = AngleVector(y);
        end
        
        % Make all angles between 0 and 180
        if AngleVectorMod(y) > 180
            AngleVectorMod(y) = AngleVectorMod(y) - 180;
        end
        
        % Centers the angles about 90 Degrees
        if AngleVectorMod(y) < 90;
            AngleVectorMod2(y) = 90 - AngleVectorMod(y);
        else AngleVectorMod2(y) = AngleVectorMod(y)*(-1) + 270;
        end
    end
    
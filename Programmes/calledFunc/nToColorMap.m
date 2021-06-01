function [col] = nToColorMap(n)
%% https://www.particleincell.com/2014/colormap/
    cmin = 0;
    cmax = 255;

    f = linspace(0,1,n);
    a=(1-f)/0.2;	%invert and group
    X=fix(a);	%this is the integer part
    Y=fix(cmax*(a-X)); %fractional part from 0 to 255
    
    col = {};
    for i=1:length(X)
        switch X(i)
            case 0
                r=cmax;g=Y(i);b=cmin;
            case 1
                r=cmax-Y(i);g=cmax;b=cmin;
            case 2
                r=cmin;g=cmax;b=Y(i);
            case 3
                r=cmin;g=cmax-Y(i);b=cmax;
            case 4
                r=Y(i);g=cmin;b=cmax;
            case 5
                r=cmax;g=cmin;b=cmax;
        end

        col{i} = [r g b]/cmax;
    end, clear i,
    
    col = col';
end
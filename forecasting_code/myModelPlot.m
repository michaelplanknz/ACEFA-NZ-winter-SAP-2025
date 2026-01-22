function myModelPlot(t, Y, clr)        

FaceAlpha = 0.5;

alpha1 = 0.5;
alpha2 = 0.25;

clr_band1 = alpha1*clr + 1-alpha1;
clr_band2 = alpha2*clr + 1-alpha2;

[x, y] = getFillArgs(t, Y(1, :), Y(5, :)  );
fill(x, y, clr_band2, 'LineStyle', 'none' , 'FaceAlpha', FaceAlpha )
hold on
[x, y] = getFillArgs(t, Y(2, :), Y(4, :)  );
fill(x, y, clr_band1, 'LineStyle', 'none' , 'FaceAlpha', FaceAlpha )
plot(t, Y(3, :), 'Color', clr) 



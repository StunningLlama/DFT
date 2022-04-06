%# Usage: ppm(fname,red,green,blue)
%#
%# Output- color ppm image in file "fname" (view with "xli fname").
%# Input- red, green, blue: 2d data of red, green, blue intensities
function ppm(fname,red,green,blue)
%# Enlarge image
for en=1:4
red=smooth(red); green=smooth(green); blue=smooth(blue);
end
pixmx=255;
height=size(red,1); width=size(red,2);
mx=max(max([red, green, blue]));
mn=min(min([red, green, blue]));
fid=fopen(fname,'w');
fprintf(fid,'P3\n');
fprintf(fid,'%d %d\n',width,height);
fprintf(fid,'%d\n',pixmx);
dat=([reshape(red',1,width*height); ...
reshape(green',1,width*height); ...
reshape(blue',1,width*height)] ...
-mn)/(mx-mn)*pixmx;
fprintf(fid,'%d ',dat);
fclose(fid);
end
function [flag] = is_inside2(area,x,y,z,p)
toll=0.0000001;
[area_new] = area_polygon2(p,x,y,z);
if area_new-area<toll
    flag=1;
else
    flag=0;
end

end
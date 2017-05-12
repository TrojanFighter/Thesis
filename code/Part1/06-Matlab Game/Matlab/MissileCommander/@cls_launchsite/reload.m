function ls = reload(ls)
global interceptors
if ls.hitcount < 3

amunition = cell(1,1);
for i = 1:10
    offset = ls.xoffset +.02 +.015*i; 
    amunition{i} = cls_interceptor(offset,.6);
    ls.interceptors(i) = i;
end

interceptors = [interceptors,amunition];
ls.maxcount = max(size(interceptors));
ls.launchcount = max(size(interceptors))-9;

end

end
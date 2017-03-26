function ls = launchinterceptor(ls,tx,ty)
global interceptors;


if ls.launchcount <= ls.maxcount && ls.hitcount < 3
    interceptors{ls.launchcount} = launch(interceptors{ls.launchcount},tx,ty);
    ls.launchcount = ls.launchcount+1;
end
% ls.launchcount = ls.launchcount+1;



end
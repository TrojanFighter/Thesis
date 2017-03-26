function highlight(ls,state)
switch state
    case 'on'
        set(ls.launchfacility,'FaceColor',[.6 .6 .6]);
    case 'off'
        set(ls.launchfacility,'FaceColor',[.5 .5 .5]);
end

end
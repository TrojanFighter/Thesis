function x = get(m,s)
switch s
    case 'Position'
        x = m.position;
    case 'stopped'
        x = m.stop;
    case 'hit'
        x = m.hit;
    otherwise
end

end
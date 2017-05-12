function delete(m)

for i=1:size(m.trail)
   if ~isnan(m.trail(i,3))
       try
       delete(m.trail(i,3));
       m.trail(i) = NaN;
       catch
%           warning('there was an error deleting missile trail');
       end
   end
end

try
    delete(m.graphic)
end

try
    delete(m.explosion(3))
end

end


function delete(abl)

for i = 1:max(size(abl.graphicshandles))
    delete(abl.graphicshandles(i));
end


end
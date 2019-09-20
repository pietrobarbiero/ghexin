idx = [];
for i = 1 : size(names)
    newidx = find(strcmp(names(i), genenames));
    idx = [idx;newidx];
end
A3 = cell(10,1);

parfor ix = 1:10
    for jx = 1:10
        A3{ix}(jx) = ix + jx;
    end
end


%%

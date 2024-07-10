for i = 1:39
    for j = 1:16
        f(j) = allRes(j).msdR(i);
        tab(i) = mean(f, 'all');
    end
end


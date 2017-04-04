function m = mfunc(Tind,Tmultiset)
    m = 1;
    % If the likelihood matches, they may be the same tree.
    llike = Tmultiset{Tind}.Lliketree;
    for ii = find(1:length(Tmultiset) ~= Tind)
        Ttemp = Tmultiset{Tind};
        if Ttemp.Lliketree == llike
             % Check Tree size
            if length(Ttemp.Allnodes) == length(Tmultiset{ii}.Allnodes)
                tm = treematch(Ttemp,Tmultiset{ii},[],[]);
                m = m + tm;
            end
        end
    end
end
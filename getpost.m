function [out,llikes,priors] = getpost(Tmultiset,X)
    K = length(Tmultiset);
    priors = zeros(K,1);
    llikes = zeros(K,1);
    for ii=1:K
        llikes(ii) = Tmultiset{ii}.Lliketree;
        % TODO: Make this more efficient with an indicator variable!
        priors(ii) = prior_eval(Tmultiset{ii},X);
    end
    post = llikes + priors;
    a = max(post);
    out = a + log(sum(exp(post - a)));
end
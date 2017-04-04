% T: an object of class T
% probs: a vector of the original proposal probabilities of the steps
% nbirths: number of births the Tree has available
% swappossible: [] or 1 indicating no/yes to a swap step on tree T
function [p_g,p_p,p_c,p_s] = propprobs(T,probs,nbirths,swappossible)
    if length(T.Allnodes) >= 5
        if nbirths > 0 % birth step possible
            if ~isempty(swappossible) % If swap is possible
                p_g = probs(1);
                p_p = probs(2);
                p_c = probs(3);
                p_s = probs(4);
            else % swap not possible
                prob_t = probs(1:3)./sum(probs(1:3));
                p_g = prob_t(1);
                p_p = prob_t(2);
                p_c = prob_t(3);
                p_s = 0;
            end
        else % If birth step is not possible
            if ~isempty(swappossible)  % swap possible
                prob_t = probs(2:4)./sum(probs(2:4));
                p_g = 0;
                p_p = prob_t(1);
                p_c = prob_t(2);
                p_s = prob_t(3);
            else % swap not possible
                prob_t = probs(2:3)./sum(probs(2:3));
                p_g = 0;
                p_p = prob_t(1);
                p_c = prob_t(2);
                p_s = 0;
            end
        end
    elseif length(T.Allnodes) > 1 % Only grow, death, and change steps available
        p_s = 0;
        if nbirths > 0 % grow step possible
            prob_t = probs(1:3)./sum(probs(1:3));
            p_g = prob_t(1);
            p_p = prob_t(2);
            p_c = prob_t(3);
        else % If grow step is not possible
            prob_t = probs(2:3)./sum(probs(2:3));
            p_g = 0;
            p_p = prob_t(1);
            p_c = prob_t(2);
        end
    else % root node
        p_g = 1;
        p_p = 0;
        p_c = 0;
        p_s = 0;
    end
end
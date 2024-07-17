
function ComputeKullbackLeiblerDivergence(P,Q)
    return sum([p*log(p/q) for (p,q) in zip(P,Q) if p>0 ])
end
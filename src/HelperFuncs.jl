function lpToAngle(lps, r0=3.8)
    return @. acos(exp(-r0/lps))
end
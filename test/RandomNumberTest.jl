using Plots, StatsBase

K =13.6
x= 0:0.01:π


d = RosenbluthChains.Cosine_BondAngle_Sampler(K, 30)


Data = SimData("$TestPath/tmp/", 1.0, 10, 8, 10, 1)

rand_num =  [begin RosenbluthChains.give_rand(d, Data);  Data.rand_val; end for _ in 1:500_000]


x=0:0.01:maximum(rand_num)
hist = fit(Histogram, rand_num, x)
StatsBase.normalize(hist, mode=:density)
hist =hist.weights ./ sum(hist.weights*(x[2]-x[1]))

x = (x[2:end].+x[1:end-1])/2.0
theory = RosenbluthChains.Cosine_BondAngle_func.(K, x)/(sum(RosenbluthChains.Cosine_BondAngle_func.(K, x))*(x[2]-x[1]))

if DOPLOTS
    fig= plot(x, theory, label="theory")
    plot!(x, hist, label="histogram")
    Plots.savefig(fig, "$TestPath/tmp/RandNum.png")
end

@test ComputeKullbackLeiblerDivergence(hist,theory) <0.1
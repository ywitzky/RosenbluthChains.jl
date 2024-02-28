using Plots

K =13.6
x= 0:0.01:π

d = RosenbluthChains.Cosine_BondAngle_Sampler(K, 30)

fig= plot(x, RosenbluthChains.Cosine_BondAngle_func.(K, x)/(sum(RosenbluthChains.Cosine_BondAngle_func.(K, x))*π/(length(x))), label="dist")
histogram!([RosenbluthChains.give_rand(d) for _ in 1:500_000], normalize=:true, label ="random number histogram")

Plots.savefig(fig, "./tmp/RandNum.png")
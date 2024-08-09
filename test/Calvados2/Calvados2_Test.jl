using Plots

#include("./HarmonicBonds_Test.jl")
#include("./Debye_Ashbaught_Hatch_Test.jl")

### data.txt was taken from https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/05_24_20/1_init_proteinsRgs.pkl and contains the data used to create figure 7 in the original calvados2 paper

NTests= 1#67

file = open("./Calvados2/Calvados2_05_24_20.csv", "r")
lines = readlines(file)
close(file)
NSeq=length(lines)

Sequences= Vector{String}(undef,NSeq)
Names= Vector{String}(undef,NSeq)
RGs= zeros(NSeq)
Temp= zeros(NSeq)
Salt = zeros(NSeq)
pH = zeros(NSeq)

for (i, line) in enumerate(lines)
    tmp = split(rstrip(line),",")
    Names[i] = tmp[1]
    Sequences[i] = strip(tmp[2])
    RGs[i] = parse(Float32,tmp[3])
    Temp[i] = parse(Float32,tmp[4])
    Salt[i] = parse(Float32,tmp[5])
    pH[i] = parse(Float32,tmp[6])
end 

N = length.(Sequences)
indices = sortperm(N) ### sort according to length
N = N[indices]
Names = Names[indices]
Sequences = Sequences[indices]
RGs = RGs[indices].*10.0 ### they measured in nm
Temp = Temp[indices]
Salt = Salt[indices]
pH = pH[indices]

### https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/05_24_20/residues.csv only for test purpose
OneToLambda_05_24_20 = Dict('D'=>0.05517, 'E'=>0.001712, 'H'=>0.4488, 'K'=>0.1785, 'F'=>0.8823, 'W'=>0.9926, 'Y'=>0.9318, 'A'=>0.2985, 'R'=>0.7307, 'N'=>0.4255, 'C'=>0.5571, 'Q'=>0.3618, 'G'=>0.6645, 'I'=>0.5362, 'L'=>0.6051, 'M'=>0.5480, 'P'=>0.3452, 'S'=>0.5616, 'T'=>0.3135, 'V'=>0.2805, '#'=>0.4625, '&'=>0.3713)

AaToWeight = Dict('A'=>71.078, 'R'=> 156.186, 'N'=> 114.103, 'D'=> 115.087, 'C'=> 103.143, 'Q'=> 128.129, 'E'=>129.114,  'G'=> 57.051, 'H'=>137.139, 'I'=>113.158, 'L'=>113.158, 'K'=>128.172, 'M'=>131.196, 'F'=>147.174, 'P'=>97.115, 'S'=> 87.077, 'T'=>101.104, 'W'=>186.210, 'Y'=> 163.173, 'V'=>99.131, 'X'=>75, '\t'=>0, '\n'=>0 , '#'=>167, '&'=>243, '@'=>241.15) #units of Da which are atomare weights  = 1.660539*10**-27kg, UNRES transfers X to G which then weights also 75Da


N_Batch = 40_000
println("Only the $(NTests) shortest proteins are compared to the results in the papers. Test whether we have less then 2.5% disagreement to paper")
Intervals = collect(5:500)
RGs_bl = zeros(length(Intervals))
RG_Errs = zeros(length(Intervals))


@testset "Calvados2 Model" begin
    for i in 1:NTests
        N_Trial = 48
        Data = SimData("../tmp/", 1.0,N[i], N_Trial, N_Batch, 1)

        Calvados2_Model = Calvados2(Sequences[i],Temp[i]; OneToLambda=OneToLambda_05_24_20, SaltConcentration=Salt[i], pH=pH[i])
        
        ### according to https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/code/analyse.py in def calcRG
        Masses = [AaToWeight[c] for c in Sequences[i]] 
        Masses[1] += 2.0
        Masses[end] += 16.0
 
        Result = RunSim(Data, Calvados2_Model, RG_Measurement(N_Batch,Masses ))

        ### numerical stability
        Result.Weights ./= maximum(Result.Weights)
        Total_Weights = sum(Result.Weights)

        ### thats how they compute it in Calvados2 papaer
        #https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/code/analyse.py
        Result.RGs = sqrt.(Result.RGs)
        Result.REEs = sqrt.(Result.REEs)

        RG_avg, RG_err = ComputeMeanError(Result.RGs, Result.Weights; NIntervals=100)

        for (ind, inter) in enumerate(Intervals)
            RGs_bl[ind], RG_Errs[ind] = ComputeMeanError(Result.RGs, Result.Weights; NIntervals=inter)
        end

        println("Name: $(Names[i])")
        println("RG: $RG_avg±$RG_err ,  Calvados2: $(RGs[i])\n")

        @test abs(RGs[i] - RG_avg)< (2.5*RG_err) ||  abs(RGs[i] - RG_avg) < RGs[i]*0.025    ###  unlikely to be outside of this margin + error in their measurements

        fig = Plots.plot(Intervals, RG_Errs, ylabel="estimated error")
        Plots.savefig(fig, "./tmp/Calvados2_Errorblock_$(Names[i])_$N_Trial.pdf")

        b = 10.0 .^ (-8:0.1:0)
        fig = Plots.histogram(Result.Weights; bin=b, xscale=:log10, yscale=:identity,  xlim=extrema(b), ylabel="Weights")
        Plots.savefig(fig, "./tmp/Calvados2_$(Names[i])_Weights_$N_Trial.pdf")

        ### @TODO trim the tests to smaller size
    end
end

N_Batch = 10_000
NTests = 1
@testset "PERM: Calvados2 Model" begin
    for i in 1:NTests
        N_Trial = 32
        Data = SimData("../tmp/", 1.0,N[i], N_Trial, N_Batch, 1)

        Calvados2_Model = Calvados2(Sequences[i],Temp[i]; OneToLambda=OneToLambda_05_24_20, SaltConcentration=Salt[i], pH=pH[i])
        
        ### according to https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/code/analyse.py in def calcRG
        Masses = [AaToWeight[c] for c in Sequences[i]] 
        Masses[1] += 2.0
        Masses[end] += 16.0

        val = 5.0
        ratio = 100.0
        PERM = PermData(N[i];K=4, cMax=val, cMin=val/ratio, PreAverage=10_000)
 
        Result = RunSim(Data, Calvados2_Model, RG_Measurement(N_Batch,Masses ), PERM)

        println("RGs: $(Result.RGs[end-10:end])")

        println("Weight: $(Result.Weights[end-10:end])")
        println("max: $(maximum(Result.Weights))")

        ### numerical stability
        Result.Weights ./= maximum(Result.Weights)
        Total_Weights = sum(Result.Weights)


        ### do same computation as in calvados2 paper
        Result.RGs = sqrt.(Result.RGs)

        RG_avg, RG_err = ComputeMeanError(Result.RGs, Result.Weights; NIntervals=100)

        for (ind, inter) in enumerate(Intervals)
            RGs_bl[ind], RG_Errs[ind] = ComputeMeanError(Result.RGs, Result.Weights; NIntervals=inter)
        end

        println("Name: $(Names[i])")
        println("RG: $RG_avg±$RG_err ,  Calvados2: $(RGs[i])\n")

        println("unique $(PERM.Total_Grown) vs. $(N_Batch)")
        println("prune $(PERM.TotalPrune) vs. $(PERM.TotalEnrich)")
        println("\n\n\n\n")

        @test abs(RGs[i] - RG_avg)< (2.5*RG_err) ||  abs(RGs[i] - RG_avg) < RGs[i]*0.025    ###  unlikely to be outside of this margin + error in their measurements

        fig = Plots.plot(Intervals, RG_Errs, ylabel="estimated error")
        Plots.savefig(fig, "./tmp/Calvados2_Errorblock_$(Names[i])_$N_Trial.pdf")

        b = 10.0 .^ (-8:0.1:0)
        fig = Plots.histogram(Result.Weights; bin=b, xscale=:log10, yscale=:identity,  xlim=extrema(b), ylabel="Weights")
        Plots.savefig(fig, "./tmp/PERM_Calvados2_$(Names[i])_Weights_$N_Trial.pdf")

        ### @TODO trim the tests to smaller size
    end
end
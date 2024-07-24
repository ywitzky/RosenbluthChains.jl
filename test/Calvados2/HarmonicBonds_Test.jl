using StatsBase, RosenbluthChains, Plots

# Function to calculate harmonic energy based on distance
function harmonic_energy(distance::Float64, k::Float64, r0::Float64)
    return 0.5 * k * (distance - r0)^2
end

function getTheoryhist(bins::Vector{R}, k::R, r0::R;  kbT=1.0) where {R<:Real}
    energies = [harmonic_energy(dist, k, r0) for dist in bins]
    Z = sum(exp.(-energies / kbT))
    theory_hist =  exp.(-energies / kbT) / Z 
    theory_hist ./= sum(theory_hist)
    return theory_hist./(bins[2]-bins[1])
end


# emulate the creation of 100_000 bond distances and their weights
N = 100_000
Weights = zeros(N)
distances = zeros(N)
radii = zeros(N)

NBeads = 4
NTrial = 8
PolymerPerBatch=1
NBatches = 1
kbT = 1

My_Simulation_Data = SimData("../tmp/", 1.0, NBeads, NTrial, PolymerPerBatch, NBatches)

bins = collect(0.0:0.1:10.0)
x = (bins[2:end].+bins[1:end-1])./2.0

@testset "HarmonicBondLength_Slow " begin
    for (k, r0, Δr) in zip([10.0,20.0,30.0], [2.0,5.0, 8.0 ], [2.0, 3.0, 4.0])

        Test_Polymer_Model = SimulationParameters(  RosenbluthChains.HarmonicBondLength_Slow(k, r0, 2.0)
        ,RandBondAngle(), RandTorsion(), IdealChain())

        New_Result = RunSim(My_Simulation_Data, Test_Polymer_Model, NoMeasurement())


        My_Simulation_Data.id = 3
        ### do some set up
        My_Simulation_Data.tmp1 .= My_Simulation_Data.xyz[3].-My_Simulation_Data.xyz[2]
        ×( My_Simulation_Data.xyz[2]-My_Simulation_Data.xyz[1], My_Simulation_Data.tmp1,My_Simulation_Data.crossproduct)
        My_Simulation_Data.current .= My_Simulation_Data.tmp1
        for i in 1:N ### we overwrite bead 4 all the time
            My_Simulation_Data.RosenbluthWeight = 1 ### reset the weight
            My_Simulation_Data.LogRosenbluthWeight=0.0
            My_Simulation_Data.id = 4

            ### emulate the step which selects the radius
            RosenbluthChains.SetTrialRadius(My_Simulation_Data, Test_Polymer_Model)
            RosenbluthChains.SetTrialBondAngle(My_Simulation_Data, Test_Polymer_Model)
            RosenbluthChains.SetTrialTorsionAngle(My_Simulation_Data, Test_Polymer_Model)
            RosenbluthChains.ComputeTrialPositions(My_Simulation_Data,Test_Polymer_Model)
            RosenbluthChains.ChooseTrialPosition(My_Simulation_Data, Test_Polymer_Model)

            My_Simulation_Data.xyz[4] .= My_Simulation_Data.trial_positions[My_Simulation_Data.tid]
            Weights[i] = My_Simulation_Data.RosenbluthWeight
            radii[i] = My_Simulation_Data.trial_radius[My_Simulation_Data.tid]

            distances[i] = RosenbluthChains.norm(My_Simulation_Data.xyz[4]-My_Simulation_Data.xyz[3]) ### in case you wouldnt trust the implementation which builds on top of .trial_radius
        end

        radii_hist = fit(Histogram, radii, weights(Weights), bins)
        StatsBase.normalize(radii_hist, mode=:density)
        radii_hist =radii_hist.weights ./ sum(radii_hist.weights*(bins[2]-bins[1]))


        distances_hist = fit(Histogram, distances, weights(Weights), bins)
        StatsBase.normalize(distances_hist, mode=:density)
        distances_hist =distances_hist.weights ./ sum(distances_hist.weights*(bins[2]-bins[1]))
        
        theory_hist = getTheoryhist(x, k, r0)

        fig = Plots.bar(x, radii_hist, label="radii")
        Plots.plot!(x, distances_hist, label="distances")
        Plots.plot!(x, theory_hist, label="theory")
        Plots.savefig(fig, "/uni-mainz.de/homes/ywitzky/Code_Projects/rosenbluthchains/test/tmp/HarmonicBondLength_Slow_$(k)_$(r0).pdf")

        @test ComputeKullbackLeiblerDivergence(radii_hist,theory_hist) <0.1
        @test ComputeKullbackLeiblerDivergence(distances_hist,theory_hist) <0.1
    end
end

### only test the random number generator since computeTrialPositions evaluates in the same way in the slow variant.
@testset "HarmonicBondLength " begin
    xgrids = collect(0.01:0.01:10) ### interpolation grid used to determin cumulative distribution functions
    for (k, r0) in zip([10.0,20.0,30.0], [2.0,5.0, 8.0 ])

        BondPotential =  HarmonicBondLength(k, r0, xgrids)

        ### get a random radius for each iteration
        for i in 1:N
            RosenbluthChains.SetTrialRadius(My_Simulation_Data,BondPotential)
            radii[i] = My_Simulation_Data.trial_radius[rand(1:NTrial)]
        end
        radii_hist = fit(Histogram, radii, bins)
        StatsBase.normalize(radii_hist, mode=:density)
        radii_hist =radii_hist.weights ./ sum(radii_hist.weights*(bins[2]-bins[1]))

        ### reuse distance to only test the random number generator
        distances.= BondPotential.inv_cdf.(rand(N))
        distances_hist = fit(Histogram, distances, weights(Weights), bins)
        StatsBase.normalize(distances_hist, mode=:density)
        distances_hist =distances_hist.weights ./ sum(distances_hist.weights*(bins[2]-bins[1]))

        theory_hist = getTheoryhist(x, k, r0)

        fig = Plots.bar(x, radii_hist, label="radii")
        Plots.plot!(x, theory_hist, label="theory")
        Plots.savefig(fig, "/uni-mainz.de/homes/ywitzky/Code_Projects/rosenbluthchains/test/tmp/HarmonicBondLength_$(k)_$(r0).pdf")

        @test ComputeKullbackLeiblerDivergence(radii_hist,theory_hist) <0.1
        @test ComputeKullbackLeiblerDivergence(distances_hist,theory_hist) <0.1
    end
end
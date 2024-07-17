using StatsBase

# Function to calculate harmonic energy based on distance
function harmonic_energy(distance::Float64, k::Float64, r0::Float64)
    return 0.5 * k * (distance - r0)^2
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

bins = collect(0, 0.1,10.0)

@testset "HarmonicBondLength_Slow " begin

    for (k, r0, Δr) in zip([10,20,30], [2.0,5.0, 8.0 ], [2.0, 3.0, 4.0])

        Test_Polymer_Model = SimulationParameters(  HarmonicBondLength_Slow(k, r0, 2.0)
        ,RandBondAngle(), RandTorsion(), IdealChain())

        New_Result = RunSim(My_Simulation_Data, New_Polymer_Model, NoMeasurement())



        My_Simulation_Data.id = 3
        ### do some set up
        My_Simulation_Data.tmp1 .= My_Simulation_Data.xyz[3].-My_Simulation_Data.xyz[2]
        crossprod( My_Simulation_Data.xyz[2]-My_Simulation_Data.xyz[1], My_Simulation_Data.tmp1,My_Simulation_Data.crossproduct)
        My_Simulation_Data.current .= My_Simulation_Data.tmp1
        for i in 1:N ### we overwrite bead 4 all the time
            My_Simulation_Data.RosenbluthWeight = 1 ### reset the weight
            My_Simulation_Data.LogRosenbluthWeight=0.0
            My_Simulation_Data.id = 4

            ### emulate the step which selects the radius
            RosenbluthChains.SetTrialRadius(My_Simulation_Data, New_Polymer_Model)
            RosenbluthChains.SetTrialBondAngle(My_Simulation_Data, New_Polymer_Model)
            RosenbluthChains.SetTrialTorsionAngle(My_Simulation_Data, New_Polymer_Model)
            RosenbluthChains.ComputeTrialPositions(My_Simulation_Data,New_Polymer_Model)
            RosenbluthChains.ChooseTrialPosition(My_Simulation_Data, New_Polymer_Model)

            My_Simulation_Data.xyz[4] .= My_Simulation_Data.trial_positions[My_Simulation_Data.tid]
            Weights[i] = My_Simulation_Data.RosenbluthWeight
            radii[i] = My_Simulation_Data.trial_radius[My_Simulation_Data.tid]
            ### data.xyz stores vector3 types which represent a normal 3d vector for which -,+,* are defined
            distances[i] = RosenbluthChains.norm(My_Simulation_Data.xyz[4]-My_Simulation_Data.xyz[3]) ### in case you wouldnt trust the implementation which builds on top of .trial_radius

            radii_hist = fit(Histogram, radii, bins)
            StatsBase.normalize(radii_hist, mode=:density)

            distances_hist = fit(Histogram, radii, bins)
            StatsBase.normalize(distances_hist, mode=:density)

            energies = [harmonic_energy(dist, k, r0) for dist in bins]
            Z = sum(exp.(-energies / kbT))
            theory_hist =  exp.(-energies / temperature) / Z 
            theory_hist ./= sum(theory_hist)


            @test ComputeKullbackLeiblerDivergence(radii_hist,theory_hist) <0.1
            @test ComputeKullbackLeiblerDivergence(distances_hist,theory_hist) <0.1

        end
    end
end
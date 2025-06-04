[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Rosenbluth Chains

Implementation of the basic Rosenbluth-Rosenbluth method used to generate and measure ensembles of short polymers (<200 Monomers) according to bond, angle and torsion potentials.

### Basic Usage with Predefined Potentials

```julia
BondLengthPotential = FixedBondParameters(r)
BondAnglePotential = RosenbluthChains.GaussianLp_Cosine_BondAngles(lp, lp*σ,N)
TorsionPotential = RandTorsion()
NonBondedPotential = LJ_Repulsion(ones(Int64,N),Dict(1 => r), r ) ### WCA with sigma=r
SimParam = SimulationParameters( BondLengthPotential, BondAnglePotential, TorsionPotential, NonBondedPotential)

FloatNum = 1.0 ### used to define float type in simulation
SimData = SimData("PathToData", FloatNum, PolymerLength, N_Trial_Moves,  BatchSize, NumberOfBatches)

Measurement =  RG_Measurement(PathForTmpData, BatchSize*NumberOfBatcher; T=Float64)
RunSim(SimulationData,SimParam, Measurement)
```

### Defining new Measurements
New measurements can be defined as a struct as a subtype of the abstract type *AbstractMeasurement* and extending the functions 

- InitMeasurement(data::SimData, param::SimulationParameters, Meas::New_Measurement) 
- MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Meas::New_Measurement)
- MeasureAfterBatch(data::SimData, param::SimulationParameters, Meas::New_Measurement)
- SaveMeasurements(data::SimData, param::SimulationParameters, Meas::New_Measurement)

The exact times at which these functions are called can be found in RunSim(...) and mainLoop(...) in *./src/RosenbluthMethod.jl*. Weights and positons of polymers will be deleted after each generation and have to be stored within the new struct if needed for analysis after batch generation. An example of an implementaion can be found  at *./src/Measurements.jl* for RG_Measurement() . 

### Defining new Potentials
New potentials have to be defined as a subtype of the apropriate abstract type *AbstractBondParam*, *AbstractBondAngleParam*, *AbstractTorsionAngleParam* or *AbstractSelfAvoidanceParameters* and define the either following function to include the energetic contribution

- GetTrialBoltzmannWeight(data::SimData, param::NewPot)

or implement the energetic contribution implicitly by drawing random bond lengths, angles, or torsion angles to according distributions and set the contribution to the weight to 0.

- SetTrialRadius(data::SimData,param::NewPot)
- SetTrialBondAngle(data::SimData,param::NewPot)
- SetTrialTorsionAngle(data::SimData,param::NewPot)

An example for both options can be found for *HarmonicBondLength* and *HarmonicBondLength_Slow* in *./src/Models/Calvados2/HarmonicBonds.jl*.

### Citation 
If you use RosenbluthChains, please cite the following papers: 

- Witzky, Yannick, Friederike Schmid, and Arash Nikoubashman. "From Heteropolymer Stiffness Distributions to Effective Homopolymers: A Conformational Analysis of Intrinsically Disordered Proteins." arXiv preprint arXiv:2504.11027 (2025).


### Caveats

- Test cases in runtest.jl are statistic in nature. Therefore it is often tested, whether the data lies within the 2.5σ intervall around the expected result such that the tests may fail. Generating more chains will increase the precision but increase the runtime of the tests.
- Type piracy has been used for Vector3 and Matrix3 defined by StaticArrays.
- The implementation of the **P**rune **E**nriched **R**osenbluth **M**ethod (PERM) has been started in PERM.jl but is not 100% successfully tested as of now. More info when running *./test/PERM_Test.jl*.
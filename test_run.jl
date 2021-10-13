# include the SSA for the nRNA and some other packages.
include("./DSSA-delayed-telegraph.jl"); using .SSAdelayedtelegraph;
using Plots, DelimitedFiles, LaTeXStrings, CSV;

# change the working dir
cd("/home/jamesholehouse/github/Detailed-nRNA-Model/")

# define some useful stuff
par_set = [2.11,0.609,0.0282]; # params should be float64.
τ = 150.0; # must be float64.
total_time = 1000000.0; # must be float64.
sp = 1.0; # must be float64.
sims = 1; # take a single trajectory. Must be integer.

# simulate
@time data=SSAdt(sims, par_set, τ, total_time, sp);

# take the individual trajectories of each species.
n_traj = data[2,1,:]; g_0 = data[1,1,:];
g_1 = 1 .- data[1,1,:];

writedlm("test_data/data1.csv",n_traj[1:end])


# plot first 1000 sample points
time = [sp*i for i in 1:1000];
plt = plot(time,n_traj[1:1000], label = L"$N$");
xlabel!("time/s"); ylabel!("Molecule #s");
gui(plt)

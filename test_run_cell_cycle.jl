# include the SSA for the nRNA and some other packages.
include("./Stoch-sims-cell-cycle.jl"); using .SSAcellcycle;
using Plots, DelimitedFiles, LaTeXStrings, Distributions, Parameters, StatsBase, LinearAlgebra, Colors;
Plots.theme(:dao) # plot theme

# change the working dir
cd("/home/jamesholehouse/github/Detailed-nRNA-Model/")

# define some useful stuff
par_set = [2.11,0.609,0.0282, 0.66, 0.05]; # params should be float64.
τ = 150.0; # must be float64.
N₀ = 10; # must be Int64 and ≧ 2.
total_time = 1000000.0; # must be float64.
sp = 1.0; # must be float64.
sims = 2; # take a single trajectory. Must be integer.

# simulate
@time data=SSAdt(sims, par_set, τ, N₀, total_time, sp);

# take the individual trajectories of each species.
n1_traj = data[2,1,:]; n2_traj = data[4,1,:];
n_traj = n1_traj + n2_traj;

writedlm("test_data/data_cc1.csv",n_traj[1:end])

function hist_prob(data::Vector{Float64})
        N = Int(floor(maximum(data)));
        mod_bins = LinRange(0.0,N,N);
        mid_pts = LinRange(0.0,N,N);
        bin_vals = normalize(fit(Histogram, data, mod_bins), mode=:probability).weights;
        return (mid_pts, bin_vals)
end

# plot first 1000 sample points
time = [sp*i for i in 1:1000];
plt = plot(time,n_traj[1:1000], label = L"$N$", legend = :none, grid = false);
xlabel!(L"\mathrm{time}/s"); ylabel!(L"\mathrm{Molecule}\quad \#s");
savefig(plt,"test_figs/fig_cc1.svg")

# plot histogram
barplt = bar(hist_prob(n_traj[1000:end]), legend = :none, grid = false)
xlabel!(L"\mathrm{Molecule}\quad \#s"); ylabel!(L"P(n)")
savefig(barplt,"test_figs/bar_cc1.svg")

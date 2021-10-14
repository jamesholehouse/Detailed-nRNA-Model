# in this file we have the code for the cell cycle model.

module SSAcellcycle

export SSAdt,mean,var

"""
propensity: function to output the propensities given the state of the system.

args:
- n: the current state vector.
- pars: the system parameters.
"""
function propensity(n::Vector{Int64}, pars::Vector{Float64}, replicated::Bool)
    # Reaction rates
    ρ = pars[1]; σoff = pars[2]; σon = pars[3]; κ = pars[4]; k = pars[5];
    f_r = zeros(7);

    if replicated == false
        f_r[1] = ρ*n[1];
        f_r[2] = σoff*n[1];
        f_r[3] = σon*(1-n[1]);
        f_r[4] = 0.0;
        f_r[5] = 0.0;
        f_r[6] = 0.0;
        f_r[7] = k;
    elseif replicated == true
        f_r[1] = κ*ρ*n[1];
        f_r[2] = σoff*n[1];
        f_r[3] = σon*(1-n[1]);
        f_r[4] = κ*ρ*n[3];
        f_r[5] = σoff*n[3];
        f_r[6] = σon*(1-n[3]);;
        f_r[7] = k;
    else
        error("'replicated' must be true or false!")
    end
    return f_r::Vector{Float64}
end

"""
SSAdt: function to perform the SSA (with the delayed degradation of nRNA) given the necessary parameters.

args:
- S_time: the number of individual simulations in the ensemble (set to 1 for single trajectory).
- pars: the parameters for the ensemble simulations.
- τ: the delayed degradation time.
- tol_time: the total simulation time to run for.
- sp: the storage time period, i.e., if sp = 1.0 then final state vector stored every 1.0s.

returns:
- the state vector at the specified times.
"""
function SSAdt(S_time::Int, pars::Array{Float64,1}, τ::Float64, N₀::Int64, tol_time::Float64, sp::Float64)
    N₀ >= 2 || error("The number of stages of the cell cycle needs to be at least 2.")
    sp <= tol_time || error("The storage time period must be less than the total simulation time!")

    # M = Number of reactions, N = Number of reactants + cell cycle stage (CCS) counter.
    M = 9::Int; N=5::Int;

    # Define stoichiometry matrix
    S_mat = zeros(M,N);
    # order is [G₁,N₁,G₂,N₂,CCS]
    S_mat[1,1:N] = [0,1,0,0,0]; # production of N₁
    S_mat[2,1:N] = [-1,0,0,0,0]; # turning on G₁
    S_mat[3,1:N] = [1,0,0,0,0]; # turning off G₁
    S_mat[4,1:N] = [0,0,0,1,0]; # production of N₂
    S_mat[5,1:N] = [0,0,-1,0,0]; # turning on G₂
    S_mat[6,1:N] = [0,0,1,0,0]; # turning off G₂
    S_mat[7,1:N] = [0,0,0,0,1]; # go to next cell-cycle stage
    S_mat[8,1:N] = [0,-1,0,0,0]; # this is the delayed reaction for N₁ degradation
    S_mat[9,1:N] = [0,0,0,-1,0]; # this is the delayed reaction for N₂ degradation

    times = convert(Array{Float64,1},LinRange(tol_time,0.0,floor(Int,tol_time/sp)+1));

    # Define reactants trjatory vector
    n = zeros(N,S_time,length(times));

    for sim in 1:S_time
        n_temp = [0;0;0;0;1]; # gene starts in G* (off) state with zero nascent in CC stage 1 (of N).
        T = 0;
        delay_times_1 = []; # delay times for N₁
        delay_times_2 = []; # delay times for N₂
        sim_times = copy(times);

        # define counter m for updating storage.
        m = 1;
        # initialise replicated to zero
        replicated = false;
        while T < tol_time
            # Step 1: Calculate propensity
            # println(n_temp)
            f_r = propensity(n_temp, pars, replicated); # propensity of each reaction.
            lambda = sum(f_r); # total propensity of any reaction.

            # Step 2: Calculate tau and mu using random number genrators
            r1 = rand(2,1);
            tau = (1/lambda)*log(1/r1[1]);
            next_r = findfirst(x -> x>=r1[2]*lambda,cumsum(f_r)); # chooses next gillespie reac.
            # println(next_r)

            # Step 3: check if there are any delayed reactions scheduled in [T,T+tau).
            fire_delay = 0; # if equal to 0 no delay reac, if 1 delay for N₁, if 2 delay for N₂.
            if length(delay_times_1)>0 || length(delay_times_2)>0
                length(delay_times_1)>0 ? nt1 = minimum(delay_times_1) : nt1 = Inf;
                length(delay_times_2)>0 ? nt2 = minimum(delay_times_2) : nt2 = Inf;
                nt, delay_reac_spec = minimum([nt1,nt2]), argmin([nt1,nt2]);
                if nt<T+tau && nt>T
                    fire_delay = delay_reac_spec; # equal to either 1 or 2.
                end
            end

            # Step 4: fire the next reaction (either delay or gillespie)
            # If no delay scheduled.
            if fire_delay == 0

                # If N production reaction fired then store the delay time.
                if next_r == 1
                    prepend!(delay_times_1,T+τ); # prepend so that next delay is on the end
                end
                # If N production reaction fired then store the delay time.
                if next_r == 4
                    prepend!(delay_times_2,T+τ); # prepend so that next delay is on the end
                end

                while T+tau >= sim_times[end]
                    n[1:N,sim,m] = n_temp; # m used here.
                    pop!(sim_times);
                    m += 1;
                    if length(sim_times) == 0
                        break
                    end
                end

                # update the system time
                T += tau;
                # update the state vector
                prod = S_mat[next_r,1:N];
                n_temp = convert(Vector{Int64}, n_temp + prod)

                # need to check for updates on the CC stage
                # first check is for replication
                if n_temp[5] == floor(Int,N₀/2)
                    replicated = true
                end
                # second check is for cell cycle end
                if n_temp[5] == N₀+1
                    replicated = false
                    # need to decide at random which gene copy to follow.
                    follow_n = rand([1,2])
                    if follow_n == 1
                        n_temp = [n_temp[1],n_temp[2],0,0,1] # inherit gene state and nascent from 1.
                        delay_times_2 = []; # reset delay times 2
                    else # follow_n == 2
                        n_temp = [n_temp[3],n_temp[4],0,0,1] # inherit gene state and nascent from 2.
                        delay_times_1 = delay_times_2; # reassign delay times 2-->1.
                        delay_times_2 = []; # reset delay times 2
                    end
                end

            # else delay is next scheduled reaction for N₁
            elseif fire_delay == 1
                deg_T = pop!(delay_times_1); # take the delay time and remove it from delay_times
                tau = deg_T - T; # time diff to the delayed reaction.
                next_r = 8; # set the next_r to the delayed reaction.

                while T+tau >= sim_times[end]
                    n[1:N,sim,m] = n_temp;
                    pop!(sim_times);
                    m += 1;
                    if length(sim_times) == 0
                        break
                    end
                end

                # update the system time
                T += tau;
                # update the state vector
                prod = S_mat[next_r,1:N];
                n_temp = convert(Vector{Int64}, n_temp + prod)
                # else delay is next scheduled reaction
            elseif fire_delay == 2
                deg_T = pop!(delay_times_2); # take the delay time and remove it from delay_times
                tau = deg_T - T; # time diff to the delayed reaction.
                next_r = 9; # set the next_r to the delayed reaction.

                while T+tau >= sim_times[end]
                    n[1:N,sim,m] = n_temp;
                    pop!(sim_times);
                    m += 1;
                    if length(sim_times) == 0
                        break
                    end
                end

                # update the system time
                T += tau;
                # update the state vector
                prod = S_mat[next_r,1:N];
                n_temp = convert(Vector{Int64}, n_temp + prod)
            else
                error("Reaction fired is either of type delay or Gillespie!")
            end

        end
        if mod(sim,1000) == 0
            println(sim)
        end
    end
    return n::Array{Float64, 3}
end

# data is a 2-D array:
# first dim is constant time, second dim diff. times. (cols of constant time)
function mean(data)
    means = zeros(size(data)[2])
    for m in 1:length(means)
        sum_all = sum(data[:,m]);
        avg = sum_all / length(data[:,m]);
        means[m] = avg;
    end
    return means
end

# data is a 2-D array:
# first dim is constant time, second dim diff. times.
function var(data,means)
    var = zeros(length(data[1,:]))
    for v in 1:length(var)
        m = means[Int(v)]
        squares = zeros(length(data[:,v]))
        for s in 1:length(squares)
            squares[s] = data[s,v] ^ 2;
        end
        sum_all_sq = sum(squares);
        va = ((sum_all_sq) / length(data[:,v])) - m^2;
        var[v] = va;
    end
    return var
end

end

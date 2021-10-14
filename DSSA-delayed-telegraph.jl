module SSAdelayedtelegraph

export SSAdt,mean,var

# define the propensity function.
function propensity(n::Vector{Int64},pars::Vector{Float64})
    # Reaction rates
    ρ = pars[1]; σoff = pars[2]; σon = pars[3];
    f_r = zeros(length(pars));

    f_r[1] = ρ*n[1];
    f_r[2] = σoff*n[1];
    f_r[3] = σon*(1-n[1]);
    return f_r
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
function SSAdt(S_time::Int, pars::Array{Float64,1}, τ::Float64, tol_time::Float64, sp::Float64)

    sp <= tol_time || error("The storage time period must be less than the total simulation time!")

    # M = Number of reactions, N = Number of reactants
    M = 4::Int; N=2::Int;

    # Define stoichiometry matrix
    S_mat = zeros(M,N);
    S_mat[1,1:N] = [0,1];
    S_mat[2,1:N] = [-1,0];
    S_mat[3,1:N] = [1,0];
    S_mat[4,1:N] = [0,-1]; # this is the delayed reaction.

    times = convert(Array{Float64,1},LinRange(tol_time,0.0,floor(Int,tol_time/sp)+1));

    # Define reactants trjatory vector
    n = zeros(N,S_time,length(times));

    for sim in 1:S_time
        n_temp = [0;0]; # gene starts in G* (off) state with zero nascent.
        T = 0;
        delay_times = [];
        sim_times = copy(times);

        # define counter m for updating storage.
        m = 1;
        while T < tol_time
            # Step 1: Calculate propensity
            f_r = propensity(n_temp,pars); # propensity of each reaction.
            lambda = sum(f_r); # total propensity of any reaction.

            # Step 2: Calculate tau and mu using random number genrators
            r1 = rand(2,1);
            tau = (1/lambda)*log(1/r1[1]);
            next_r = findfirst(x -> x>=r1[2]*lambda,cumsum(f_r));

            # Step 3: check if there are any delayed reactions scheduled in [T,T+tau).
            fire_delay = 0; # boolean var. set to 1 if delay is to be scheduled instead.
            if length(delay_times)>0
                nt = delay_times[end];
                if nt<T+tau && nt>T
                    fire_delay = 1;
                end
            end

            # Step 4: fire the next reaction (either delay or gillespie)
            # If no delay scheduled.
            if fire_delay == 0

                # Step 4.5: if N production reaction fired then store the delay time.
                if next_r == 1
                    prepend!(delay_times,T+τ); # prepend so that next delay is on the end
                end

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
                for i in 1:N
                    n_temp[i] += prod[i]
                end

            # else delay is next scheduled reaction
            elseif fire_delay == 1
                deg_T = pop!(delay_times); # take the delay time and remove it from delay_times
                tau = deg_T - T; # time diff to the delayed reaction.
                next_r = 4; # set the next_r to the delayed reaction.

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
                for i in 1:N
                    n_temp[i] += prod[i]
                end
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

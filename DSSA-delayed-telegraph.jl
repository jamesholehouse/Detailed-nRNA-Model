module SSAdelayedtelegraph

export SSAdt,mean,var

# define the propensity function.
function propensity(n,pars)
    # Reaction rates
    ρ = pars[1]; σoff = pars[2]; σon = pars[3];
    f_r = zeros(length(pars));

    f_r[1] = ρ*n[1];
    f_r[2] = σoff*n[1];
    f_r[3] = σon*(1-n[1]);
    return f_r
end


function SSAdt(S_time::Int, pars::Array{Float64,1}, τ::Float64, tol_time::Float64, sp::Float64)
    # M = Number of reactions
    # N = Number of reactants
    M = 4::Int; N=2::Int;
    # Define stoichiometry matrix
    S_mat = zeros(M,N);
    S_mat[1,1:N] = [0,1];
    S_mat[2,1:N] = [-1,0];
    S_mat[3,1:N] = [1,0];
    S_mat[4,1:N] = [0,-1];# this is the delayed reaction.
    # Define reactants trjatory vector
    n = zeros(N,S_time,Int(floor(tol_time/sp)));

    for sim in 1:S_time
        n_temp = [0;0]; # gene starts in U** state with zero nascent or mature.
        T = 0;
        delay_times = [];

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

                T += tau;

                # Update the trajectory vector
                if T <= tol_time
                    for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                        n[1:N,sim,t+1] = n_temp;
                    end
                    # Fire reaction next_r
                    prod = S_mat[next_r,1:N];
                    for i in 1:N
                        n_temp[i] = n_temp[i] + prod[i];
                    end
                else
                    for t in Int(ceil((T-tau)/sp)) : Int(floor(tol_time/sp)-1)
                        n[1:N,sim,t+1] = n_temp;
                    end
                end
            # else delay is scheduled and fire the next delayed reaction.
            elseif fire_delay == 1
                prev_T = T;
                T = pop!(delay_times); # take the delay time and remove it from delay_times
                next_r = 4; # set the next_r to the delayed reaction.

                # Update the trajectory vector
                if T <= tol_time
                    for t in Int(ceil((prev_T)/sp)) : Int(floor(T/sp))
                        n[1:N,sim,t+1] = n_temp;
                    end
                    # Fire reaction next_r
                    prod = S_mat[next_r,1:N];
                    for i in 1:N
                        n_temp[i] = n_temp[i] + prod[i];
                    end
                else
                    for t in Int(ceil((prev_T)/sp)) : Int(floor(tol_time/sp)-1)
                        n[1:N,sim,t+1] = n_temp;
                    end
                end
            else
                print("Error")
            end

        end
        if mod(sim,1000) == 0
            println(sim)
        end
    end
    return n
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

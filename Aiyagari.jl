# Bibliotecas
using Distributions
using Plots
using DataFrames
using Random
using ForwardDiff
using LinearAlgebra
using Interpolations
using DataFrames
using Optim
Random.seed!(1234)
using NBInclude

# Importando código 
@nbinclude("utils.ipynb") 
@nbinclude("RBC.ipynb")  
@nbinclude("EGM.ipynb")   
@nbinclude("SteadyState.ipynb") 
@nbinclude("GenBKM.ipynb") 

# Estado Estácionario
# Encontrando o valor do estado estacionario do capital
p = Params() 
H_ss = 1.0
@time oo = optimize(K -> eq_K(K,p), 10, 100, Brent()) 
K_star = oo.minimizer;
println("O valor do capital o estado estacionario é = $(K_star)")
g_star, c_star, g_low_star, g_high_star, success_flag= solve_EGM(x->log(x), x->log(x), R(K_star, H_ss, p), W(K_star, H_ss, p), p); #solve for policy functions
t_star = make_trans_mat(g_star, p)   
d_star = get_stationary_dist(t_star);

# Visualizando a convergência em direção ao estado estacionário
K_grid = collect(range(oo.minimizer-0.5, stop=oo.minimizer+0.5, length=20))
K_implied_grid = similar(K_grid)
R_grid = similar(K_grid)

for (K_index, K_value) in enumerate(K_grid)
    R_s, W_s = R(K_value, H_ss, p), W(K_value, H_ss, p) 
    gg, c_star, g_low, g_high, success_flag= solve_EGM(x -> log(x), x -> 2*log(x), R_s, W_s, p) 
    tt = make_trans_mat(gg, p)  
    dd = get_stationary_dist(tt) 
    K_implied = aggregate_K(dd, p) 
    R_grid[K_index] = R_s 
    K_implied_grid[K_index] = K_implied 
    K_grid[K_index] = K_value 
end


# Gerando os graficos da oferta e demanda do capital 

plot(K_grid, R_grid, label = "Demanda", ylabel="Taxa de Juros", xlabel="Capital")
plot!(K_implied_grid, R_grid, label = "Oferta")
savefig("/home/alexandre/Documentos/Mestrado/Macroeconomia/AiyagariAggregateUncertainty")
# Choque mit 
# Atualização posterior e anterior 
# Função Atualização anterior
function backward_update(g_low_ss::Function, g_high_ss::Function, K_path_guess::Array{Float64,1}, H_path::Array{Float64,1}, p::Params)

    
    nT = length(H_path)
    g_low_path = Array{Function}(undef,nT) 
    g_high_path =  Array{Function}(undef,nT)
    g_low_path[nT] = g_low_ss
    g_high_path[nT] = g_high_ss
    a_path = zeros(p.nI, p.grid_size, nT)
    R_path = zeros(nT) 
    W_path = zeros(nT)


    for t=nT:-1:2 
        R_path[t], W_path[t] = R(K_path_guess[t], H_path[t], p), W(K_path_guess[t], H_path[t], p)
        R_path[t-1], W_path[t-1] = R(K_path_guess[t-1], H_path[t-1], p), W(K_path_guess[t-1], H_path[t-1], p)
        a_path[:,:,t-1], c_new, g_low_path[t-1], g_high_path[t-1] = euler_back(g_low_path[t], g_high_path[t], R_path[t-1], W_path[t-1], R_path[t], W_path[t], p)
    end

    return a_path, g_low_path, g_high_path
end

# Função Atualização Posterior
function forward_update(K_star::Float64, a_path::Array{Float64,3}, d_ss::Array{Float64,1}, p::Params)
    nT = length(H_path)
    K_path_forward = zeros(nT)
    K_path_forward[1] = K_star
    dd_path_forward = zeros(size(d_ss,1), nT)
    dd_path_forward[:,1] = d_ss
    for t=2:nT
        tt = make_trans_mat(a_path[:,:,t-1], p) 
        dd_path_forward[:,t] = tt*dd_path_forward[:,t-1]
        K_path_forward[t] = aggregate_K(dd_path_forward[:,t], p)
    end

    return dd_path_forward, K_path_forward
end


# Encontrando o caminho de transição
function solve_mit!(K_path, g_low_ss::Function, g_high_ss::Function, d_ss::Array{Float64,1},
                    K_ss::Float64, H_path::Array{Float64,1}, p::Params; convex_combination::Float64 = 0.2,
                    shrink_factor::Float64 = 0.5, expand_factor::Float64 = 1.05,
                    max_iter::Int64 = 1000, tol::Float64=1e-6, verbose::Bool=true, display_iter::Int64 = 20)
    diff = Inf 
    diff_old = Inf 
    convergence_flag = 0 
    damp = convex_combination 
    for i_mit=1:max_iter

       
        a_path, g_low_path, g_high_path = backward_update(g_low_ss, g_high_ss, K_path[i_mit], H_path, p);

        dd_path_forward, K_path_forward = forward_update(K_ss, a_path, d_star, p);

        diff = maximum(abs.(K_path_forward - K_path[i_mit]))

        if verbose==true
            if mod(i_mit,display_iter) == 0
                println("Interação $(i_mit). diferença = $(diff)")
            end
        end

        if diff < tol
            if verbose==true
                println("Convergence rejeitada após $(i_mit) Interações.")
            end
            convergence_flag = 1
            break
        else
            
            if diff > diff_old
                damp = max(min(damp * shrink_factor, 1.0-eps()), eps())
            else
                damp = max(min(damp * expand_factor, 1.0-eps()), eps())
            end
            if mod(i_mit, 10) == 0
                if verbose==true
                    println("damp = $(damp); diff = $(diff)")
                end
            end
            push!(K_path, damp.*K_path_forward .+ (1.0 - damp).*K_path[i_mit])
            diff_old = diff

        end
    end

    return K_path, convergence_flag
end


# Encontrando o caminho para {K_t} para um choque de 1 padrão 
max_t = 300 
H_ss = 1.0  
H_shock = 2.0 

H_path = ones(max_t)
H_path[1] = H_ss*H_shock 


for t_index=2:max_t
    H_path[t_index] = H_path[t_index-1]^p.rho
end


K_path = []
push!(K_path, repeat([K_star], max_t))

@time K_path, convergence_flag = solve_mit!(K_path, g_low_star, g_high_star, d_star, K_star, H_path, p, convex_combination=0.25)


R_path = zeros(length(H_path)) 
W_path = zeros(length(H_path)) 
for t=length(H_path):-1:1 
    R_path[t], W_path[t] = R(K_path[end][t], H_path[t], p), W(K_path[end][t], H_path[t], p)
end
# Visualizando a convergência do caminho de transição

p0 = plot(1:max_t, K_path[1], label= "Interação 0", title="Convergencia de K(t)")
plot!(p0, 2:max_t, K_path[2][2:end], label = "Interação 1")
show_every = 5 
for k in 2:length(K_path)
    if mod(k,show_every) == 0
        plot!(p0, 2:max_t, K_path[k][2:end], xlabel="t", label = "Interação $(k)", title="Convergencia de K(t)", legend=:best)
    end
end

p0
savefig("/home/alexandre/Documentos/Mestrado/Macroeconomia/AiyagariAggregateUncertainty/2.png")


p1 = plot(1:max_t, K_path[end], label= "K_t Aiyagari", title="IRF Aiyagari versus RBC")
plot!(p1, xx[RBCp.iK,2:end] .+ K_star, label = "K_t RBC", color = "black", xlabel="t")
savefig("/home/alexandre/Documentos/Mestrado/Macroeconomia/AiyagariAggregateUncertainty/3.png")


p1 = plot(1:max_t, H_path./H_path[end] .-1, label = "H(t)", xlabel= "t")
p2 = plot(1:max_t, R_path./R_path[end] .-1, label= "R(t)", xlabel= "t")
p3 = plot(1:max_t, W_path./W_path[end].-1 , label= "W(t)", xlabel= "t")
p4 = plot(1:max_t, K_path[end]./K_path[end][end] .-1, label= "K(t)", xlabel= "t" )

p5 = plot(p1, p2, p3, p4)
savefig("/home/alexandre/Documentos/Mestrado/Macroeconomia/AiyagariAggregateUncertainty/5.png")
# checando a linearidade do modelo
max_t = 300 
H_ss = 1.0  

array_sigma = collect(range(-0.75, stop=0.75, step=0.25))

array_sigma = array_sigma[array_sigma .!= 0.]

x_mit_scaled_sigma = zeros(max_t, length(array_sigma))

H_path_sigma = zeros(max_t, length(array_sigma))
H_path_sigma_dev = zeros(max_t, length(array_sigma))

for (index_sigma, sigma) in enumerate(array_sigma)

    H_path = ones(max_t)
    H_path[1] = H_ss + H_ss*sigma

    for t_index=2:max_t
        H_path[t_index] = H_path[t_index-1]^p.rho
    end

    K_path = []
    push!(K_path, repeat([K_star], max_t))

    @time K_path, convergence_flag = solve_mit!(K_path, g_low_star, g_high_star, d_star, K_star, H_path, p, convex_combination=0.2, verbose=false);

    if convergence_flag!=1
        error("No convergence for H(1) = $(H_path[1]).")
    end

    H_path_sigma[:, index_sigma] = H_path

    H_path_sigma_dev[:, index_sigma] = H_path./H_ss .- 1.0

    # Scaled IRF: how a percentage deviation in z_t from its steady-state results in a % deviation of k_t
    x_mit_scaled_sigma[:, index_sigma] = (K_path[end]./K_star .- 1.0)./H_path_sigma_dev[1, index_sigma]

end

p0 = plot()
p1 = plot()
p2 = plot()
p3 = plot()

for (index_sigma, sigma) in enumerate(array_sigma)
    if index_sigma == 1
        p0 = plot(100 .*H_path_sigma_dev[:, index_sigma], label="H(t) com H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%")
        p1 = plot(H_path_sigma[:, index_sigma], label="H(t) H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%")
        p2 = plot(x_mit_scaled_sigma[:, index_sigma], label="H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%")
        p3 = plot(sign(H_path_sigma[1, index_sigma] - H_ss)*x_mit_scaled_sigma[:, index_sigma], label="H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%")
    else
        plot!(p0, 100 .*H_path_sigma_dev[:, index_sigma], label="H(t) with H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%", title = "% Desvio Padrão do Capital Humano", xlabel="t")
        plot!(p1, H_path_sigma[:, index_sigma], label="H(t) with H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%", title = "Agragação do Capital Humano", xlabel="t")
        plot!(p2, x_mit_scaled_sigma[:, index_sigma], label="H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%", title = "Impact of shock's size on scaled IRF: x", xlabel="t")
        plot!(p3, sign(H_path_sigma[1, index_sigma] - H_ss)*x_mit_scaled_sigma[:, index_sigma], label="H(1) = $(round(100 .*H_path_sigma_dev[1, index_sigma], digits=2))%", title = "Impact of shock's sign on scaled IRF: sign(x - x_ss)*x", xlabel="t")
    end
end

plot(p2,p3, fg_legend = :transparent, legend=:best, layout=(2,1))


# Dinamica fora do estado estacionario
max_t = 2000
shocks_t = rand(Normal(0,0.005), max_t) # Series of aggregate shocks
# Let's generate a path for the aggregate shock
H_path = ones(max_t)
H_path[1] = H_ss

# Evolution of aggregate productivity in level:
for t_index=2:max_t
    H_path[t_index] = H_path[t_index-1]^p.rho + shocks_t[t_index]
end

# Calculation of GenBKM path:
XT_GenBKM = zeros(max_t);# Initialization
@time GenBKM_path!(XT_GenBKM, max_t, x_mit_scaled_sigma, H_path./H_ss .- 1.0, array_sigma)


p1 = plot(100 .*(H_path./H_ss .- 1.0), label="H_t", xlabel="t")
p2 = plot(100 .*XT_GenBKM, label = "K_t", xlabel="t")
plot(p1,p2, fg_legend = :transparent, legend=:best, layout=(2,1))

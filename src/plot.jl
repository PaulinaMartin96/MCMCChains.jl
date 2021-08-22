@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@userplot EnergyPlot
#@shorthands energyplot

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end
#struct _EnergyPlot; marginal_energy; energy_transition; p_type; n_chains; end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot,
                        :pooleddensity => _DensityPlot,
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner, :energyplot)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains, i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false,
    plot_type = :density
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> string.(names(c))
        val = c.value[:, :, i]
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> map(x -> "Chain $x", MCMCChains.chains(c))
        val = c.value[:, i, :]
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity || st == :pooleddensity
        discrete = indiscretesupport(c, barbounds)
        st = if colordim == :chain
            discrete[i] ? :histogram : :density
        else
            # NOTE: It might make sense to overlay histograms and density plots here.
            :density
        end
        seriestype := st
    end

    if st == :autocorplot
        lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(range(c)))) : maxlag)
        ac = autocor(c; sections = nothing, lags = lags)
        ac_mat = convert(Array, ac)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    #elseif st == :energyplot
    #    p_type = plot_type
    #    energy_section = get(c, :hamiltonian_energy)
    #    #@show energy_section
    #    #@show params.hamiltonian_energy
    #    n_chains = (append_chains ? 1 : size(c, 3))
    #    energy_data = (append_chains ? vec(energy_section.hamiltonian_energy.data) : energy_section.hamiltonian_energy.data)
    #    mean_energy = vec(mean(energy_data, dims = 1))
    #    marginal_energy = [energy_data[:,i] .- mean_energy[i] for i in 1:n_chains]
    #    energy_transition = [energy_data[2:end,i] .- energy_data[1:end-1,i] for i in 1:n_chains]
    #    _EnergyPlot(marginal_energy, energy_transition, p_type, n_chains)
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        range(c), val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    bins --> 25
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), cummean(p.val)
end

@recipe function f(p::_AutocorPlot)
    seriestype := :path
    xaxis --> "Lag"
    yaxis --> "Autocorrelation"
    p.lags, p.val
end

@recipe function f(p::_TracePlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Sample value"
    range(p.c), p.val
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{Symbol};
    colordim = :chain
)
    colordim != :chain &&
        error("Symbol names are interpreted as parameter names, only compatible with ",
              "`colordim = :chain`")

    ret = indexin(parameters, names(chains))
    any(y === nothing for y in ret) && error("Parameter not found")

    return chains, Int.(ret)
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{<:Integer} = Int[];
    sections = _default_sections(chains),
    width = 500,
    height = 250,
    colordim = :chain,
    append_chains = false
)
    _chains = isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes*width, N*height)
        legend --> false

        multiple_plots = N * ntypes > 1
        if multiple_plots
            layout := (N, ntypes)
        end

        i = 0
        for par in parameters
            for ptype in ptypes
                i += 1

                @series begin
                    if multiple_plots
                        subplot := i
                    end
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        ntypes > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, names(c)[parameters])
    end
end

struct Corner
    c
    parameters
end

@recipe function f(corner::Corner)
    label --> permutedims(corner.parameters)
    compact --> true
    size --> (600, 600)
    ar = collect(Array(corner.c.value[:, corner.parameters,i]) for i in chains(corner.c))
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end

#function compute_energy(
#    chains::Chains,
#    combined = false,
#    plot_type = :density
#)
#    st = get(plotattributes, :seriestype, :traceplot)
#
#    if st == :energyplot
#        p_type = plot_type
#        params = get(chains, :hamiltonian_energy)
#        n_chains = (combined ? 1 : size(chains, 3))
#        energy_data = (combined ? vec(params.hamiltonian_energy.data) : params.hamiltonian_energy.data)
#        mean_energy = vec(mean(energy_data, dims = 1))
#        marginal_energy = energy_data[:,i] .- mean_energy[i]
#        energy_transition = energy_data[2:end,i] .- energy_data[1:end-1,i]
#        _EnergyPlot(marginal_energy, energy_transition, p_type, n_chains)
#    else
#
#    end
#end

#@recipe function f(
#    chains::Chains;
#    plot_type = :density,
#    append_chains = false
#)
#
#    st = get(plotattributes, :seriestype, :traceplot)
#    if st == :energyplot
#        p_type = plot_type
#        energy_section = get(chains, :hamiltonian_energy)
#        #@show energy_section
#        #@show params.hamiltonian_energy
#        n_chains = (append_chains ? 1 : size(chains, 3))
#        energy_data = (append_chains ? vec(energy_section.hamiltonian_energy.data) : energy_section.hamiltonian_energy.data)
#        mean_energy = vec(mean(energy_data, dims = 1))
#        marginal_energy = [energy_data[:,i] .- mean_energy[i] for i in 1:n_chains]
#        energy_transition = [energy_data[2:end,i] .- energy_data[1:end-1,i] for i in 1:n_chains]
#        _EnergyPlot(marginal_energy, energy_transition, p_type, n_chains)
#    elseif st ∈ supportedplots
#        translationdict[st](c, val)
#    end
#end

function compute_energy(
    chains::Chains,
    combined = false,
    plot_type = :density
)
        p_type = plot_type
        params = get(chains, :hamiltonian_energy)
        isempty(params) && error("EnergyPlot receives a Chains object containing only the
            :internals section. Please use Chains(chain, [:internals]) to create it")
        n_chains = (combined ? 1 : size(chains, 3))
        energy_data = (combined ? vec(params.hamiltonian_energy.data) : params.hamiltonian_energy.data)
        mean_energy = vec(mean(energy_data, dims = 1))
        marginal_energy = [energy_data[:,i] .- mean_energy[i] for i in 1:n_chains]
        energy_transition = [energy_data[2:end,i] .- energy_data[1:end-1,i] for i in 1:n_chains]
        return marginal_energy, energy_transition, p_type, n_chains
    end

@recipe function f(
    p::EnergyPlot;
    combined = false,
    plot_type = :density
    )

    c = p.args[1]
    #p_type = plot_type
    #params = get(c, :hamiltonian_energy)
    #isempty(params) && error("EnergyPlot receives a Chains object containing only the
    #    :internals section. Please use Chains(chain, [:internals]) to create it")
    #n_chains = (combined ? 1 : size(c, 3))
    #energy_data = (combined ? vec(params.hamiltonian_energy.data) : params.hamiltonian_energy.data)
    #mean_energy = vec(mean(energy_data, dims = 1))
    #marginal_energy = [energy_data[:,i] .- mean_energy[i] for i in 1:n_chains]
    #energy_transition = [energy_data[2:end,i] .- energy_data[1:end-1,i] for i in 1:n_chains]
    marginal_energy, energy_transition, p_type, n_chains = compute_energy(c, combined, plot_type)
    k = 0
    for i in 1:n_chains
        k += 1
        title --> "Chain $(MCMCChains.chains(c)[i])"
        subplot := i
        @series begin
            seriestype := p_type
            label --> "Marginal energy"
            marginal_energy[i]
        end

        @series begin
            seriestype := p_type
            label --> "Energy transition"
            energy_transition[i]
        end
    end
end

#@recipe function f(p::_EnergyPlot)
#
#    k = 0
#    for i in 1:p.n_chains
#        k = 1
#        @series begin
#            subplot := i
#            seriestype := p.p_type
#            label --> "Marginal energy"
#            p.marginal_energy[i]
#        end
#
#        @series begin
#            subplot := i
#            seriestype := p.p_type
#            label --> "Energy transition"
#            p.energy_transition[i]
#        end
#    end
#end
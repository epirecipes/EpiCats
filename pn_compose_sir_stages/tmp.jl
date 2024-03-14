using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab
using LabelledArrays
using OrdinaryDiffEq
using Plots

nstages = 4
sub(i::Int) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))
sub(x::String,i::Int) = x*sub(i);

epi_lpn = LabelledPetriNet(
  [:Pop],
  :infection=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :progression=>(:Pop=>:Pop),
  :recovery=>(:Pop=>:Pop)
);

# S-I infection transitions

# programatically generate the wiring diagram for infection from stages
function si_stages_uwd(n)
    uwd = RelationDiagram(repeat([:Pop], n+1))
    junctions = Dict(begin
        if i <= n
            variable = Symbol(sub("I",i))
        else
            variable = :S
        end 
        junction = add_junction!(uwd, :Pop, variable=variable)
        set_junction!(uwd, port, junction, outer=true)
        variable => junction
    end for (i, port) in enumerate(ports(uwd, outer=true)))
    for i in 1:n
        box_wires = [:S, Symbol(sub("I",i)),Symbol(sub("I",1)),Symbol(sub("I",i))]
        box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in box_wires], name=:infection)
        for (i, port) in enumerate(ports(uwd, box))
            set_junction!(uwd, port, junctions[box_wires[i]])
        end
    end
    return uwd
end

si_uwd = si_stages_uwd(nstages)

to_graphviz(si_uwd, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:size=>"6", :ratio=>"fill"))

si_uwd_hand = @relation (S, I₁, I₂, I₃, I₄) where (S::Pop, I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop) begin
    infection(S,I₁,I₁,I₁)
    infection(S,I₂,I₁,I₂)
    infection(S,I₃,I₁,I₃)
    infection(S,I₄,I₁,I₄)
end

@assert length(isomorphisms(si_uwd, si_uwd_hand)) > 0

betas = Symbol.([sub("β",i) for i=1:nstages])
si_acst = oapply_typed(epi_lpn, si_uwd, betas)
si_lpn = dom(si_acst)
to_graphviz(si_lpn)

# I-R progression transitions

# programatically generate the wiring diagram for infection from stages
function ir_stages_uwd(n)
    uwd = RelationDiagram(repeat([:Pop], n+1))
    junctions = Dict(begin
        variable = i <= n ? Symbol(sub("I",i)) : :R
        junction = add_junction!(uwd, :Pop, variable=variable)
        set_junction!(uwd, port, junction, outer=true)
        variable => junction
    end for (i, port) in enumerate(ports(uwd, outer=true)))
    for i in 1:n
        box_wires = i < n ? [Symbol(sub("I",i)), Symbol(sub("I",i+1))] : [Symbol(sub("I",i)), :R] 
        box_name = i < n ? :progression : :recovery
        box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in box_wires], name=box_name)
        for (i, port) in enumerate(ports(uwd, box))
            set_junction!(uwd, port, junctions[box_wires[i]])
        end
    end
    return uwd
end

ir_uwd = ir_stages_uwd(nstages)

to_graphviz(ir_uwd, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:size=>"6", :ratio=>"fill"))

ir_uwd_hand = @relation (I₁, I₂, I₃, I₄, R) where (I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop, R::Pop) begin
    progression(I₁,I₂)
    progression(I₂,I₃)
    progression(I₃,I₄)
    recovery(I₄,R)
end

@assert length(isomorphisms(ir_uwd, ir_uwd_hand)) > 0

deltas = Symbol.([sub("δ",i) for i=1:nstages])
ir_acst = oapply_typed(epi_lpn, ir_uwd, deltas)
ir_lpn = dom(ir_acst)
to_graphviz(ir_lpn)

# the whole shebang
function sir_stages_uwd(n)
    uwd = RelationDiagram(repeat([:Pop], n+2))
    states = [[Symbol(sub("I",i)) for i in 1:4]; :S; :R]
    junctions = Dict(begin
        junction = add_junction!(uwd, :Pop, variable=state)
        set_junction!(uwd, port, junction, outer=true)
        state => junction
    end for (state, port) in zip(states, ports(uwd, outer=true)))
    # add si box
    box_wires = states[[n+1;1:n]]
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in box_wires], name=:si)
    for (i, port) in enumerate(ports(uwd, box))
        set_junction!(uwd, port, junctions[box_wires[i]])
    end
    # add ir box
    box_wires = states[[1:n;end]]
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in box_wires], name=:ir)
    for (i, port) in enumerate(ports(uwd, box))
        set_junction!(uwd, port, junctions[box_wires[i]])
    end
    return uwd
end

sir_uwd = sir_stages_uwd(nstages)

to_graphviz(sir_uwd, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:size=>"6", :ratio=>"fill"))

sir_uwd_hand = @relation (S, I₁, I₂, I₃, I₄, R) where (S::Pop, I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop, R::Pop) begin
    si(S, I₁, I₂, I₃, I₄)
    ir(I₁, I₂, I₃, I₄, R)
end

@assert length(isomorphisms(sir_uwd, sir_uwd_hand)) > 0

sir_smc = oapply(sir_uwd, Dict(
    :si => Open(si_lpn),
    :ir => Open(ir_lpn),
))
sir_lpn = apex(sir_smc)
to_graphviz(sir_lpn)
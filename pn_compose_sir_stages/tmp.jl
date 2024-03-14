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
        si_pair = [:S, Symbol(sub("I",i)),Symbol(sub("I",1)),Symbol(sub("I",i))]
        box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in si_pair], name=:infection)
        for (i, port) in enumerate(ports(uwd, box))
            set_junction!(uwd, port, junctions[si_pair[i]])
        end
    end
    return uwd
end

si_uwd = si_stages_uwd(nstages)

to_graphviz(si_uwd, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:size=>"6", :ratio=>"fill"))

betas = Symbol.([sub("β",i) for i=1:nstages])
si_acst = oapply_typed(epi_lpn, si_uwd, betas)
si_lpn = dom(si_acst)
to_graphviz(si_lpn)



n = 4
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
    si_pair = [:S, Symbol(sub("I",i)),Symbol(sub("I",1)),Symbol(sub("I",i))]
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in si_pair], name=:infection)
    for (i, port) in enumerate(ports(uwd, box))
        set_junction!(uwd, port, junctions[si_pair[i]])
    end
end
    
to_graphviz(uwd, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:size=>"6", :ratio=>"fill"))





si_uwd = @relation (S, I₁, I₂, I₃, I₄) where (S::Pop, I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop) begin
    infection(S,I₁,I₁,I₁)
    infection(S,I₂,I₁,I₂)
    infection(S,I₃,I₁,I₃)
    infection(S,I₄,I₁,I₄)
end
betas = Symbol.([sub("β",i) for i=1:nstages])
si_acst = oapply_typed(epi_lpn, si_uwd, betas)
si_lpn = dom(si_acst)
to_graphviz(si_lpn)



si_acst1 = oapply_typed(epi_lpn, uwd, betas)
si_lpn1 = dom(si_acst1)
to_graphviz(si_lpn1)


ir_uwd = @relation (I₁, I₂, I₃, I₄, R) where (I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop, R::Pop) begin
    progression(I₁,I₂)
    progression(I₂,I₃)
    progression(I₃,I₄)
    recovery(I₄,R)
end
deltas = Symbol.([sub("δ",i) for i=1:nstages])
ir_acst = oapply_typed(epi_lpn, ir_uwd, deltas)
ir_lpn = dom(ir_acst)
to_graphviz(ir_lpn)

sir_uwd = @relation (S, I₁, I₂, I₃, I₄, R) where (S::Pop, I₁::Pop, I₂::Pop, I₃::Pop, I₄::Pop, R::Pop) begin
    si(S, I₁, I₂, I₃, I₄)
    ir(I₁, I₂, I₃, I₄, R)
end
sir_smc = oapply(sir_uwd, Dict(
    :si => Open(si_lpn),
    :ir => Open(ir_lpn),
))
sir_lpn = apex(sir_smc)
to_graphviz(sir_lpn)
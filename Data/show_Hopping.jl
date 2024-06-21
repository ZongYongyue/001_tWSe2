include("/Users/zongyy/Library/Mobile Documents/com~apple~CloudDocs/Clone/Fork/master/QuantumClusterTheories.jl/src/tools.jl")
using Plots


hps = loadData("./using/TBAMIM/Hoppings.jls")

f=plot()

cls = [:red, :blue, :green, :black]
ji = [1,3,5,7]
ou = [2,4,6,8]
for i in 1:4
    plot!(f,[hp[ji[i]] for hp in hps]; color=cls[i], legend=false) 
end

for i in 1:4
    plot!(f,[hp[ou[i]] for hp in hps];ls=:dash, color=cls[i],legend=false) 
end
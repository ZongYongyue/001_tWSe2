include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
using Distributed
#开启多进程
spawn(11)
@everywhere begin
include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
include("/fsa/home/jxl_zongyy/script/PaperData/L12MIM.jl")
end 

Vim = range(15.5, 15.6, 11)

musim = range(11.681270903010033,11.581270903010033, 11)

vcas = pmap(i->moireVCA(4.2, Vim[i], 4.7), [1,3,5,7,9,11])

sdfg = pmap(i->SDF(vcas[i], 11.581270903010033), [1,3,5,7,9,11])

saveData(sdfg, "L12IM_SU4p7_SDFG.jls")
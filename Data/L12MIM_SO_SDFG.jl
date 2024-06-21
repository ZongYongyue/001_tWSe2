include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
using Distributed
#开启多进程
spawn(21)
@everywhere begin
include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
include("/fsa/home/jxl_zongyy/script/PaperData/L12MIM.jl")
end 


vcas = pmap(i->moireVCA(4.2, V[i], 4.7, sp[i]), 1:length(sp))

sdfg = pmap(i->SDFG(vcas[i], mus[i]), 1:length(mus))

saveData(sdfg, "L12MIM_SO_SU4p7_SDFG.jls")

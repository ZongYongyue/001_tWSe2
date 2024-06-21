include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
using Distributed
#开启多进程
spawn(7)
@everywhere begin
include("/fsa/home/jxl_zongyy/mycode/QuantumClusterTheories.jl/src/tools.jl")
include("/fsa/home/jxl_zongyy/script/PaperData/L12MIM.jl")
    function VUD(vca, μ)
        unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
        rz = ReciprocalZone(reciprocals(unitcell.vectors); length=200)
        ω_range = 4.755330863071197*range(-10, 10, length=500)
        GG = singleParticleGreenFunction(:f, vca, rz, ω_range; η=0.02*4.755330863071197, μ=μ)
        D = densityofstates(GG)
        return D
    end
end 

us = [3,4,5,6,7,8,9]
vcas = pmap(i->moireVCA(4.2, 9, us[i]), 1:7)
musu = us/2 .+ 1
vud = pmap(i->VUD(vcas[i], musu[i]), 1:7)

saveData(vud, "L12MIM_VUD.jls")
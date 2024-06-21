include("/Users/zongyy/Library/Mobile Documents/com~apple~CloudDocs/Clone/Fork/master/QuantumClusterTheories.jl/src/tools.jl")
using Plots

colorbar = :gnuplot2#:bukavu
colorbarg = cgrad([:yellow, :red, :LightCoral, :white, :DeepSkyBlue, :blue, :MidnightBlue],[0.1,0.25,0.45,0.55,0.75,0.9],rev = true)
sdfg = loadData("./using/PaperData/L12IM_SU4p7_SDFG.jls")
so_sdfg = loadData("./using/PaperData/L12IM_SO_SU4p7_SDFG.jls")

unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
dz = ReciprocalZone(reciprocals([[0.5,0],[0,0.5]]); length=1000)
k_path = ReciprocalPath(reciprocals(unitcell.vectors), hexagon"K-Γ-M₂-K, 120°", length=300)
ω_range = 4.755330863071197*range(-10, 10, length=500)

for i in 1:21
    s = plot(k_path, ω_range, sdfg[i][1]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function",clims=(0, 0.3))
    savefig(s, "./result/PaperData_IM/L12IM_Spec_$i.pdf")

    su = plot(k_path, ω_range, sdfg[i][2]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function spin-up",clims=(0, 0.3))
    savefig(su, "./result/PaperData_IM/L12IM_Spec_spinup_$i.pdf")

    sd = plot(k_path, ω_range, sdfg[i][3]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function spin-down",clims=(0, 0.3))
    savefig(sd, "./result/PaperData_IM/L12IM_Spec_spindown_$i.pdf")

    d = plot(ω_range, sdfg[i][4]; xlabel="ω" , title="Density of states")
    savefig(d, "./result/PaperData_IM/L12IM_Dens_$i.pdf")

    # 定义正六边形的顶点坐标
    x = π*[2/√3, 0, -2/√3, -2/√3, 0, 2/√3]
    y = π*[2/3, 4/3, 2/3, -2/3, -4/3, -2/3]

    # 将首尾两个点连接起来形成闭合的多边形
    push!(x, x[1])
    push!(y, y[1])

    cls = (0, 0.3)

    f = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),sdfg[i][5], ratio=1, color=colorbar,clims=cls,title="Density of states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(f, "./result/PaperData_IM/L12IM_Fermi_$i.pdf")

    fu = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),sdfg[i][6], ratio=1, color=colorbar,clims=cls,title="Density of spin up states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(fu, "./result/PaperData_IM/L12IM_Fermi_spinup_$i.pdf")

    fd = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),sdfg[i][7], ratio=1, color=colorbar,clims=cls,title="Density of spin down states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(fd, "./result/PaperData_IM/L12IM_Fermi_spindown_$i.pdf")
    
    cls2 = (-0.3, 0.3)

    gsu = plot(k_path, ω_range, sdfg[i][8]; xlabel="k", ylabel="ω", color=colorbarg, title="Spectral Function spin-up",clims=cls2)
    savefig(gsu, "./result/PaperData_IM/L12IM_Gre_Spec_spinup_$i.pdf")

    gsd = plot(k_path, ω_range, sdfg[i][9]; xlabel="k", ylabel="ω", color=colorbarg, title="Spectral Function spin-down",clims=cls2)
    savefig(gsd, "./result/PaperData_IM/L12IM_Gre_Spec_spindown_$i.pdf")

    gfu = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),sdfg[i][10], ratio=1, color=colorbarg,clims=cls2,title="Gre of spin up states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(gfu, "./result/PaperData_IM/L12IM_Gre_Fermi_spinup_$i.pdf")

    gfd = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),sdfg[i][11], ratio=1, color=colorbarg,clims=cls2,title="Gre of spin down states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(gfd, "./result/PaperData_IM/L12IM_Gre_Fermi_spindown_$i.pdf")
end
#=
for i in 1:21
    s = plot(k_path, ω_range, so_sdfg[i][1]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function",clims=(0, 0.3))
    savefig(s, "./result/PaperData_IM_SO/SO_L12IM_Spec_$i.pdf")

    su = plot(k_path, ω_range, so_sdfg[i][2]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function spin-up",clims=(0, 0.3))
    savefig(su, "./result/PaperData_IM_SO/SO_L12IM_Spec_spinup_$i.pdf")

    sd = plot(k_path, ω_range, so_sdfg[i][3]; xlabel="k", ylabel="ω", color=colorbar, title="Spectral Function spin-down",clims=(0, 0.3))
    savefig(sd, "./result/PaperData_IM_SO/SO_L12IM_Spec_spindown_$i.pdf")

    d = plot(ω_range, so_sdfg[i][4]; xlabel="ω" , title="Density of states")
    savefig(d, "./result/PaperData_IM_SO/SO_L12IM_Dens_$i.pdf")

    # 定义正六边形的顶点坐标
    x = π*[2/√3, 0, -2/√3, -2/√3, 0, 2/√3]
    y = π*[2/3, 4/3, 2/3, -2/3, -4/3, -2/3]

    # 将首尾两个点连接起来形成闭合的多边形
    push!(x, x[1])
    push!(y, y[1])

    cls = (0, 0.3)

    f = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),so_sdfg[i][5], ratio=1, color=colorbar,clims=cls,title="Density of states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(f, "./result/PaperData_IM_SO/SO_L12IM_Fermi_$i.pdf")

    fu = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),so_sdfg[i][6], ratio=1, color=colorbar,clims=cls,title="Density of spin up states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(fu, "./result/PaperData_IM_SO/SO_L12IM_Fermi_spinup_$i.pdf")

    fd = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),so_sdfg[i][7], ratio=1, color=colorbar,clims=cls,title="Density of spin down states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(fd, "./result/PaperData_IM_SO/SO_L12IM_Fermi_spindown_$i.pdf")
    
    cls2 = (-0.3, 0.3)

    gsu = plot(k_path, ω_range, so_sdfg[i][8]; xlabel="k", ylabel="ω", color=colorbarg, title="Spectral Function spin-up",clims=cls2)
    savefig(gsu, "./result/PaperData_IM_SO/SO_L12IM_Gre_Spec_spinup_$i.pdf")

    gsd = plot(k_path, ω_range, so_sdfg[i][9]; xlabel="k", ylabel="ω", color=colorbarg, title="Spectral Function spin-down",clims=cls2)
    savefig(gsd, "./result/PaperData_IM_SO/SO_L12IM_Gre_Spec_spindown_$i.pdf")

    gfu = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),so_sdfg[i][10], ratio=1, color=colorbarg,clims=cls2,title="Gre of spin up states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(gfu, "./result/PaperData_IM_SO/SO_L12IM_Gre_Fermi_spinup_$i.pdf")

    gfd = heatmap(range(-2pi,2pi,1000),range(-2pi,2pi,1000),so_sdfg[i][11], ratio=1, color=colorbarg,clims=cls2,title="Gre of spin down states on fermisurface")
    plot!(x, y, seriestype=:path, linecolor=:white, linestyle=:dot, linewidth=2.5, legend=false, aspect_ratio=1)
    savefig(gfd, "./result/PaperData_IM_SO/SO_L12IM_Gre_Fermi_spindown_$i.pdf")
end
=#
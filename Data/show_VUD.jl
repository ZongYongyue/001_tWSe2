include("/Users/zongyy/Library/Mobile Documents/com~apple~CloudDocs/Clone/Fork/master/QuantumClusterTheories.jl/src/tools.jl")
using Plots


# vud = loadData("./using/PaperData/L12MIM_VUD_2.jls")

# ω_range = range(-60, 60, length=600)
# us = [4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5]
# for i in 1:length(us)
#     d = plot(ω_range, vud[i]; xlabel="ω" , title="Density of states")
#     savefig(d, "./result/PaperData_VUD/L12I9_2_Dens_$i.pdf")
# end

# gap = [3.000713931043201
# 5.863614356340376
# 12.62163159587175
# 16.332626160164313
# 20.18211909184686]

# scatter([20,24,28,32,36], gap, )

gp =   [2.6206901104513496
4.111158585291278
5.61929429389652
7.160522963225172
8.84302573046355
10.627984162199668
12.437943990557455
14.23997203640727
16.17235698253297
18.03885807053302
20.020433726679133
22.039110370679353]

U= [19.021323452284786
 21.398988883820383
 23.776654315355984
 26.15431974689158
 28.531985178427178
 30.909650609962778
 33.28731604149838
 35.664981473033976
 38.04264690456957
 40.42031233610517
 42.797977767640766
 45.17564319917637]


 pl = plot(U, gp, xlims=(18, 46), ylims=(0,25), color=parse(RGB,"#1E90FF"))
 scatter!(pl, U, gp,xlims=(18, 46), ylims=(0,25), color=:red, frame=:box)

 savefig("./result/PaperFig/linearGap.pdf")
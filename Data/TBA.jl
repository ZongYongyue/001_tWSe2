using MoireSuperlattices
using Plots
using QuantumLattices
using TightBindingApproximation
using DelimitedFiles: readdlm
#=
parameters = (a₀=3.28, m=0.45, θ=4.2, Vᶻ=0, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4); parameters=parameters, map=bltmdmap)

emin, emax = -200.0, 100.0
recipls = bltmd.frontend.reciprocallattice.translations
lattice = MoireTriangular(6, reciprocals(recipls))
hilbert = Hilbert(1=>Fock{:f}(1, 2))
brillouinzone = BrillouinZone(recipls, 24)
tba = Algorithm(:tba, TBA(lattice, hilbert, terms(bltmd, lattice, brillouinzone; atol=10^-6)))
#@test all(map((x, y)->isapprox(x, y, atol=10^-6), collect(tba.parameters), [-2.7598267, -4.3678292, -1.3035002, 0.0, 0.2447067, -0.6094541, -0.2700574, -0.3888003, 0.0260020, -0.0308282, -0.2999706, 0.0, 10.2102190]))
plt = plot(grid = false, frame=:box)
plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100))), ylim=(emin, emax), color=parse(RGB,"#1E90FF"),grid = false, title="", lw=2)
plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100))), ylim=(emin, emax), color="black",grid = false, title="", ls=:dot, lw=2)
plot!(plt, tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100))), ylim=(emin, emax), ls=:dash, color="red", grid = false, size=(400, 300), title="", lw=2)
#savefig("WeSe₂-4p2V0.pdf")
=#

names = ["./using/PaperData/dft0.csv", "./using/PaperData/dft0p3.csv"]

f = plot()
parameters = (a₀=3.28, m=0.45, θ=5.08, Vᶻ=0.0, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4); parameters=parameters, map=bltmdmap)

emin, emax = -150.0, 50.0
recipls = bltmd.frontend.reciprocallattice.translations
lattice = MoireTriangular(6, reciprocals(recipls))
hilbert = Hilbert(1=>Fock{:f}(1, 2))
brillouinzone = BrillouinZone(recipls, 24)
tba = Algorithm(:tba, TBA(lattice, hilbert, terms(bltmd, lattice, brillouinzone; atol=10^-6)))
path = ReciprocalPath(recipls, hexagon"Γ-K₁-M₁", length=100)
plot!(f, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁", length=100))), ylim=(emin, emax), color=parse(RGB,"#1E90FF"),grid = false, title="", lw=2)
plot!(f, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄", length=100))), ylim=(emin, emax), color="black",grid = false, title="", ls=:dot, lw=2)
data1 = readdlm(names[1], ','; comments=true, comment_char='x')
data1[:, 1] .= data1[:, 1]./(150/distance(path))
plot!(f, data1[:, 1], data1[:,2:end]; ls=:dash, color=:red, frame=:box, grid = false)

g = plot()
parameters = (a₀=3.28, m=0.45, θ=5.08, Vᶻ=14.13, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4); parameters=parameters, map=bltmdmap)

plot!(g, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁", length=100))), ylim=(emin, emax), color=parse(RGB,"#1E90FF"),grid = false, title="", lw=2)
plot!(g, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄", length=100))), ylim=(emin, emax), color="black",grid = false, title="", ls=:dot, lw=2)
data2 = readdlm(names[2], ','; comments=true, comment_char='x')
data2[:, 1] .= data2[:, 1]./(150/distance(path))
data2[:, 2:end] .+= 11.5 
plot!(g, data2[:, 1], data2[:,2:end]; ls=:dash, color=:red,frame=:box, grid = false)

savefig(f, "./result/PaperFig/dft0.pdf")
savefig(g,"./result/PaperFig/dft0p3.pdf")
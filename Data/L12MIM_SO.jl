function moireGP(θ, Vᶻ, u, x, μ)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])#Lattice([0, 0]; vectors = [[1, 0], [-1/2, √3/2]])
    cluster = Lattice([0,0],[√3/2, 1/2], [-√3/2, 1/2], [0, 1], [√3/2, 3/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [-√3, 2], [-√3/2, 5/2], [0, 3]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)#Lattice([0.0, 0.0], [1, 0], [-1/2, √3/2], [1/2, √3/2], [3/2, √3/2], [-1, √3], [0, √3], [1, √3], [2, √3], [-1/2, 3√3/2], [1/2, 3√3/2], [3/2, 3√3/2]; vectors = [[3, √3], [-3, √3]])
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
    bs = Sector(hilbert, ParticleNumber(12))
    t = moiretuning(θ, Vᶻ)
    U = Hubbard(:U, 4.755330863071197*u)
    origiterms = (t..., U)
    r = Onsite(:r, Complex(x),coupling_by_spinrotation(π/2, 2π/3); amplitude=belongs([1,5,7,12]))
    g = Onsite(:g, Complex(x),coupling_by_spinrotation(π/2, 0.0); amplitude=belongs([2,3,8,10]))
    b = Onsite(:b, Complex(x),coupling_by_spinrotation(π/2, 4π/3); amplitude=belongs([4,6,9,11]))
    referterms = (t..., U, r, g, b)
    rz = ReciprocalZone(reciprocals(cluster.vectors); length=100)
    vca = VCA(:S, unitcell, cluster, hilbert, origiterms, referterms, bs; modelname="Moire-MIM", cachepath="/fsa/home/jxl_zongyy/data/VCA_Moire_MIM/")
    return GrandPotential(:f, vca, rz, μ)
end

function moireVCA(θ, Vᶻ, u, x)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])#Lattice([0, 0]; vectors = [[1, 0], [-1/2, √3/2]])
    cluster = Lattice([0,0],[√3/2, 1/2], [-√3/2, 1/2], [0, 1], [√3/2, 3/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [-√3, 2], [-√3/2, 5/2], [0, 3]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)#Lattice([0.0, 0.0], [1, 0], [-1/2, √3/2], [1/2, √3/2], [3/2, √3/2], [-1, √3], [0, √3], [1, √3], [2, √3], [-1/2, 3√3/2], [1/2, 3√3/2], [3/2, 3√3/2]; vectors = [[3, √3], [-3, √3]])
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(cluster))
    bs = Sector(hilbert, ParticleNumber(12))
    t = moiretuning(θ, Vᶻ)
    U = Hubbard(:U, 4.755330863071197*u)
    origiterms = (t..., U)
    r = Onsite(:r, Complex(x),coupling_by_spinrotation(π/2, 2π/3); amplitude=belongs([1,5,7,12]))
    g = Onsite(:g, Complex(x),coupling_by_spinrotation(π/2, 0.0); amplitude=belongs([2,3,8,10]))
    b = Onsite(:b, Complex(x),coupling_by_spinrotation(π/2, 4π/3); amplitude=belongs([4,6,9,11]))
    referterms = (t..., U, r, g, b)
    return VCA(:S, unitcell, cluster, hilbert, origiterms, referterms, bs; modelname="Moire-MIM", cachepath="/fsa/home/jxl_zongyy/data/VCA_Moire_MIM/")
end
moireGPP(x) = moireGP(x[1], x[2], x[3], x[4],x[5])
options = Optim.Options(x_tol=1e-4, f_tol=5e-7, iterations=20, show_trace=false)
V = range(0, 20, 21)

mus = [14.285714285714286
14.24561403508772
14.155388471177945
14.045112781954888
15.889724310776943
15.468671679197994
14.967418546365915
14.37593984962406
13.74436090225564
13.152882205513784
12.711779448621554
12.551378446115288
12.601503759398497
12.731829573934837
12.81203007518797
10.80220300751879
10.08203007518797
9.98203007518797
9.97203007518797
9.92203007518797
9.90203007518797]
sp = loadData("/fsa/home/jxl_zongyy/script/L12MIM_SO_OP_6/L12MIM_OP_6_U4p7_sp.jls")
cluster = Lattice([0,0],[√3/2, 1/2], [-√3/2, 1/2], [0, 1], [√3/2, 3/2], [-√3, 1], [-√3/2, 3/2], [0, 2], [√3/2, 5/2], [-√3, 2], [-√3/2, 5/2], [0, 3]; vectors = [[-√3, 3], [√3, 3]], name=:Triangle_L12)
rz = ReciprocalZone(reciprocals(cluster.vectors); length=200)
so = (
    Onsite(:r, Complex(1.0),coupling_by_spinrotation(π/2, 2π/3); amplitude=belongs([1,5,7,12])),
    Onsite(:g, Complex(1.0),coupling_by_spinrotation(π/2, 0.0); amplitude=belongs([2,3,8,10])),
    Onsite(:b, Complex(1.0),coupling_by_spinrotation(π/2, 4π/3); amplitude=belongs([4,6,9,11]))
)

function mindens(vca, arr)
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    rz = ReciprocalZone(reciprocals(unitcell.vectors); length=200)
    dens = Vector{Float64}(undef, length(arr))
    for i in eachindex(arr)
        dens[i] = densityofstates(singleParticleGreenFunction(:f, vca, rz, 0:1:0; η=0.08*4.755330863071197, μ=arr[i]))[1]
    end
    dmin_index = argmin(dens)
    return arr[dmin_index]
end

function SDFG(vca, μ)
    function specre(gfpathv::AbstractVector; select::AbstractVector=Vector(1:size(gfpathv[1][1],2)))
        A = zeros(Float64, length(gfpathv), length(gfpathv[1]))
            for i in eachindex(gfpathv)
                for j in eachindex(gfpathv[i])
                    A[i, j] = (det(gfpathv[i][j][select, select])).re
                end
            end
            return A
    end
    unitcell = Lattice([0, 0]; vectors = [[√3/2, 1/2], [0, 1]])
    rz = ReciprocalZone(reciprocals(unitcell.vectors); length=200)
    dz = ReciprocalZone(reciprocals([[0.5,0],[0,0.5]]); length=1000)
    k_path = ReciprocalPath(reciprocals(unitcell.vectors), hexagon"K-Γ-M₂-K, 120°", length=300)
    ω_range = 4.755330863071197*range(-10, 10, length=500)
    G = singleParticleGreenFunction(:f, vca, k_path, ω_range; η=0.08*4.755330863071197, μ=μ)
    S = spectrum(G)
    Sd = spectrum(G; select=[1])
    Su = spectrum(G; select=[2])
    GG = singleParticleGreenFunction(:f, vca, rz, ω_range; η=0.08*4.755330863071197, μ=μ)
    D = densityofstates(GG)
    GGG = singleParticleGreenFunction(:f, vca, dz, 0:1:0; η=0.08*4.755330863071197, μ=μ)
    F = reshape(spectrum(GGG), Int(√(length(dz))), Int(√(length(dz))))
    Fd = reshape(spectrum(GGG; select=[1]), Int(√(length(dz))), Int(√(length(dz))))
    Fu = reshape(spectrum(GGG; select=[2]), Int(√(length(dz))), Int(√(length(dz))))
    GSd = specre(G; select=[1])
    GSu = specre(G; select=[2])
    GFd = reshape(specre(GGG; select=[1]), Int(√(length(rz))), Int(√(length(rz))))
    GFu = reshape(specre(GGG; select=[2]), Int(√(length(rz))), Int(√(length(rz))))
    return S, Su, Sd, D, F, Fu, Fd, GSu, GSd, GFu, GFd
end



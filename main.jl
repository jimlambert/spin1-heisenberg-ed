# main function for heisenberg diagonalization

include("Spin1Utls.jl")

using LinearAlgebra

# == β RANGE == #
#n = 2
#D = 0.08
#h = 0.0
#
#βmin = -0.4
#βmax = -0.3
#βstp = 100
#βinc = (βmax-βmin)/βstp
#
#for i=1:βstp
#  
#  β = βmin + βinc*i
#  magblock = Spin1Utls.inith(true, n, D, h, β)
#
#  println(β, '\t', eigvals(magblock)[1]/n)
#
#end


# == CHAIN LATTICE == #
n = 4
B = 0.0 # biquadratic term
β = 16  # inverse temperature
D = -0.3
h = 0.0

# Construct lattice geometry
chainbonds = Spin1Utls.chain(true, n)

# Construct observables
sz = -Spin1Utls.initσz(n, 1)
sx = -Spin1Utls.initσx(n, 1)

dz = sz*sz

for i=2:n
  global sz += ((-1)^i)*Spin1Utls.initσz(n,i)
  global sx += ((-1)^i)*Spin1Utls.initσx(n,i)
  global dz += Spin1Utls.initσz(n,i)^2
end


# Get eigenvalues and eigenvectors
H = Spin1Utls.inith(chainbonds, n, D, h, B)
eigvas = eigvals(H)
eigvcs = eigvecs(H)
gs = eigvcs[:,1]

println("Stag. Mag. Z= ", '\t', adjoint(gs)*sz*gs)
println("Stag. Mag. Z Var. = ", '\t' , adjoint(gs)*sz*sz*gs)
println("Stag. Mag. X= ", '\t', adjoint(gs)*sx*gs)
println("Stag. Mag. X Var. = ", '\t' , adjoint(gs)*sx*sx*gs)

# convert everything to by diagonal in the eigenbasis of the Hamiltonian
H = adjoint(eigvcs)*H*eigvcs
sz = adjoint(eigvcs)*sz*eigvcs
dz = adjoint(eigvcs)*dz*eigvcs


# construct density matrix
ρ = exp(-β*H) / tr(exp(-β*H))


# compute disconnected variance
totvarsz = tr(sz^2*ρ)
totvardz = tr(dz^2*ρ)

# compute thermal variance
αmin = 0.0
αmax = 1.0
αstp = 1000
dα = (αmax - αmin) / αstp

thrmvarsz = 0.0
thrmvardz = 0.0
for i=1:αstp
  global thrmvarsz += tr(sz*ρ^(i*dα)*sz*ρ^(1-i*dα))  
  global thrmvardz += tr(dz*ρ^(i*dα)*dz*ρ^(1-i*dα))
end

thrmvarsz = thrmvarsz / αstp
thrmvardz = thrmvardz / αstp

println("SM Z Quantum Var = ", '\t', totvarsz - thrmvarsz)
println("DZ Quantum Var = ", '\t', totvardz - thrmvardz)

println("D = ", '\t', D)
println("β = ", '\t', β)
println("e_0 = ", '\t', eigvas[1]/n)
#println(β, '\t', eigs[2])
#println(β, '\t', eigs[3])
#println(β, '\t', eigs[4])
#println(β, '\t', eigs[5])


# == SQUARE LATTICE == #
#lx = 2
#ly = 2
#bcx = false
#bcy = false
#nsites = lx*ly
#
#sqbonds = Spin1Utls.isosqlatt(bcx,bcy,lx,ly)
#
#println(sqbonds)
#
#H = Spin1Utls.inith(sqbonds, nsites, D, h, β)
#println(β, '\t', eigvals(H)[1]/(lx*ly))


# == HONEYCOMB LATTICE == #
#β = 0.0
#D = 0.0
#h = 0.0
#lx = 2
#ly = 2
#bcx = true
#bcy = true
#nsites = 2*lx*ly
#
#honeycomb = Spin1Utls.honeycombLattice(bcx,bcy,lx,ly)
#
#println(honeycomb)
#
#H = Spin1Utls.inith(honeycomb, nsites, D, h, β)
#println(β, '\t', eigvals(H)[1]/nsites)


#n = 4
#β = -0.5
#D = 0.0
#h = 0.0
#
#magblock = Spin1Utls.inith(false, n, D, h, β)

#println(β, '\t', eigvals(magblock)[1]/(n))



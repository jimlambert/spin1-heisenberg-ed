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
#n = 4
#β = -1.0/3.0
#D = 0.0
#h = 0.0
#
#chainbonds = Spin1Utls.chain(true, n)
#
#
#H = Spin1Utls.inith(chainbonds, n, D, h, β)
#println(β, '\t', eigvals(H)[1]/n)


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


# == SQUARE LATTICE == #
β = 0.0
D = 0.0
h = 0.0
lx = 2
ly = 2
bcx = true
bcy = true
nsites = 2*lx*ly

honeycomb = Spin1Utls.honeycombLattice(bcx,bcy,lx,ly)

println(honeycomb)

H = Spin1Utls.inith(honeycomb, nsites, D, h, β)
println(β, '\t', eigvals(H)[1]/nsites)


#n = 4
#β = -0.5
#D = 0.0
#h = 0.0
#
#magblock = Spin1Utls.inith(false, n, D, h, β)

#println(β, '\t', eigvals(magblock)[1]/(n))



# main function for heisenberg diagonalization

include("Spin1Utls.jl")

using LinearAlgebra

n = 2
D = 0.0
h = 0.0

βmin = -1.0
βmax = 0.0
βstp = 100
βinc = (βmax-βmin)/βstp

for i=1:βstp
  
  β = βmin + βinc*i
  magblock = Spin1Utls.inith(true, n, D, h, β)

  println(β, '\t', eigvals(magblock)[1]/n)

end

#n = 4
#β = 0.0
#D = 0.0
#h = 0.0

#chainbonds = Spin1Utls.chain(true, n)


#H = Spin1Utls.inith(chainbonds, n, D, h, β)
#println(β, '\t', eigvals(H)[1]/n)

#lx = 2
#ly = 2
#bcx = false
#bcy = false
#
#sqbonds = Spin1Utls.isosqlatt(bcx,bcy,lx,ly)
#
#println(sqbonds)
#
#H = Spin1Utls.inith(sqbonds, lx*ly, D, h, β)
#println(β, '\t', eigvals(H)[1]/(lx*ly))


#n = 4
#β = -0.5
#D = 0.0
#h = 0.0
#
#magblock = Spin1Utls.inith(false, n, D, h, β)

#println(β, '\t', eigvals(magblock)[1]/(n))



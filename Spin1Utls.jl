module Spin1Utls

export initMagBlock


"""
# chain - 1 member function

## Methods

  * chain(bc::Bool, n::Int)

### Arguments

"""
function chain(bc::Bool, n::Int)
  
  bonds = []

  bc ? m=n : m=n-1

  for i=1:m
    push!(bonds,[i, (i%n) + 1])
  end

  return bonds

end


"""
# isosqlatt -  1 member function

"""
function isosqlatt(bcx::Bool, bcy::Bool, lx::Int, ly::Int)

  bonds = []

  bcx ? mx = lx : mx = lx-1
  bcy ? my = ly : my = ly-1

  for y=1:ly
    for x=1:mx
      
      xmin = (y-1)*lx
      xmax = y*lx
      push!(bonds,[x+xmin, (x+xmin)%xmax + 1 + (xmin*((x)÷lx))])

    end
  end

  for y=1:my
    for x=1:lx
      xmin = (y-1)*lx
      push!(bonds,[x+xmin, (x+xmin+lx)%(lx*ly+1) + (x+xmin+lx)÷(lx*ly+1)])
    end
  end

  return bonds
end


"""
# honeycomb lattice - 1 member function
"""
function honeycombLattice(bcx::Bool, bcy::Bool, lx::Int, ly::Int)
  
  chlen = 2*lx
  nchains = ly
  nsites = 2*lx*ly
 
  bcx ? mch=nchains : mch = nchains - 1
  bcy ? mcl=chlen   : mcl = chlen - 1

  bonds = []
  
  for chain=0:(nchains-1)
    
    offset = chain * chlen

    # x and y bonds
    for r=1:mcl
      rstrt = offset+1
      rend = offset + chlen

      rpos1::Int64 = r+offset
      rpos2::Int64 = 0 # this is just a placeholder value
    
      # periodic boundary conditions
      if (rpos1==rend) rpos2=rstrt else rpos2 = rpos1+1 end
     
      # consistently choose the first position to be on sublattice A
      if ((rpos1%2)==0)
        bond = [rpos1, rpos2]
        push!(bonds, bond)
      else
        bond = [rpos2, rpos1]
        push!(bonds, bond)
      end
    
    end

  end
  
  # z bonds
  for chain=0:(mch-1)
    offset = chain * chlen
    for r=2:2:chlen 
      rpos1::Int64 = r + offset
      rpos2::Int64 = (rpos1+chlen-1) % (nsites+1) + (rpos1+chlen-1)÷(nsites+1)
      bond = [rpos1, rpos2]
      push!(bonds, bond)
    end
  end

  return bonds
end


"""
# state2id - 1 method function

## Methods

  * state2id(Array{Int, 1})

### Arguments

  * state::Array{Int64, 1}

### Returns

  * state id

## Descriptions

Accepts an integer containing a particular state of an integer spin chain and
maps that state to an identifying integer. The identifying integer is chosen by
representing that integer in ternary. Each bit of the ternary number corresponds
to the spin state at this position plus one. Thus the spin state -1 maps to 0, 0
to 1, +1 to 2. The spin pair (-1, 1) maps to ternary number (0, 2) which would
be 6.

For convenience we begin the state identifiers at 1.
  
"""
function state2id(state::Array{Int64, 1})
  sum = 0
  len = length(state)
  for (site, spin) in enumerate(state)
    sum += (spin+1)*3^(site-1)
  end
  return sum+1
end


"""
# spin - 1 method function

## Methods 

  * spin(Int, Int, Int)

### Arguments

  * id::Int   -   spin state in the form of integer id
  * n::Int    -   number of spins in state
  * i::Int    -   spin index

### Return

  * spin value

## Description

Returns the spin contained at site, i, of the spin state identified by the
integer, id.  
"""
function spin(id::Int, n::Int, i::Int)
  return digits(id-1, base=3, pad=n)[i] - 1
end

"""
# spinflips - 1 method function

## Methods

  * spinflips(Int, Int, Int, Int)

### Arguments

  * id::Int   -   integer identifying spin state
  * n::Int    -   system size (for padding purposes)
  * i::Int    -   first site on bond
  * j::Int    -   second site on bond

### Returns

  * Array{Int, 1}   -   list of states coupled to this state by spinflip
                        operators

## Description

This function accepts a spin state through an id and two site indices indicated
which bond is to be flipped. It then find all connected spin states and returns
those spin ids in the form of an array. This function is especially useful for
determining matrix elements. 
"""
function spinflips(id::Int, n::Int, i::Int, j::Int)
  # recover state from id number
  state = broadcast(-, digits(id-1, base=3, pad=n), 1)
  coupstates = Int64[]
  if state[i]<1 && state[j] > -1
    state[i] += 1
    state[j] -= 1
    push!(coupstates, state2id(state))
  end
  state = broadcast(-, digits(id-1, base=3, pad=n), 1)
  if(state[i]>-1 && state[j] < 1)
    state[i] -= 1
    state[j] += 1
    push!(coupstates, state2id(state))
  end
  return coupstates
end


function doubspinflips(id::Int, n::Int, i::Int, j::Int)
  
  state = broadcast(-, digits(id-1, base=3, pad=n), 1)
  
  coupstates = Int64[]

  if state[i]==1 && state[j]==-1
    state[i] -= 2
    state[j] += 2
    push!(coupstates, state2id(state))
  end
  
  state = broadcast(-, digits(id-1, base=3, pad=n), 1)

  if state[i]==-1 && state[j]==1
    state[i] += 2
    state[j] -= 2
    push!(coupstates, state2id(state))
  end
   
  return coupstates
end


"""
# inith - 2 method function

## Methods

  * inith(Bool, Int)
  * inith(Bool, Int, Float64)

### Arugments

  * bc    -   boundary conditions on the lattice
  * n     -   number of sites on lattice
  * D     -   uniaxial anisotropy term

### Returns

  * Array{Float, 2}   -   Hamiltonian of spin1 chain

## Description 

This function create the Hamiltonian of a spin-1 AFM Heisenberg chain with OBC
or PBC.
"""
function inith(bc::Bool, n::Int)
  H = zeros(3^n, 3^n) # initialize empty Hamiltonian
  for a=1:3^n
    bc ? m=n : m=n-1 # select number of bonds based on boundary conditions
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      H[a, a] += spin(a, n, i)*spin(a, n, j)
      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        H[a, b] += 1
      end
    end
  end
  return H
end

function inith(bc::Bool, n::Int, D::Float64)
  H = zeros(3^n, 3^n) # initialize empty Hamiltonian
  for a=1:3^n
    bc ? m=n : m=n-1 # select number of bonds based on boundary conditions
    for i=1:m
      
      j = (i%n)+1 # determine nearest neighbour
      s1 = spin(a, n, i) 
      s2 = spin(a, n, j) 
      H[a, a] += s1*s2 + (D/2.0)*(s1^2 + s2^2)
      
      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        H[a, b] += 1
      end

    end
  end
  return H
end


function inith(bc::Bool, n::Int, D::Float64, h::Float64, β::Float64)
  H = zeros(3^n, 3^n) # initialize empty Hamiltonian
  for a=1:3^n
    bc ? m=n : m=n-1 # select number of bonds based on boundary conditions
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      s1 = spin(a, n, i) 
      s2 = spin(a, n, j) 
      H[a, a] += s1*s2 - β*(s1^2)*(s2^2) + (D/2.0)*(s1^2 + s2^2) + (h/2.0)*(s1+s2)
      
      if (s1==s2 && abs(s1) == 0) 
        H[a, a] -= 2*β
      end
      
      if (s1==-s2 && s1 !=0)
        H[a,a] -= β
      end

      if(s1==0 && abs(s2)==1) || (abs(s1)==1 && s2==0)
        H[a,a] -= β
      end

      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        sb1 = spin(b, n, i)
        sb2 = spin(b, n, j)
        H[a, b] += 1 - β*sb1*sb2-β*s1*s2
      end

      doubcoupstates = doubspinflips(a, n, i, j)
      for b in doubcoupstates
        H[a,b] += -β
      end

    end
  end
  return H
end


function inith(bondlst::Array{Any,1}, n::Int, D::Float64, h::Float64, β::Float64)
  H = zeros(3^n, 3^n) # initialize empty Hamiltonian
  for a=1:3^n
    for bond in bondlst
      
      i = bond[1]
      j = bond[2]

      s1 = spin(a, n, i) 
      s2 = spin(a, n, j) 
      H[a, a] += s1*s2 - β*(s1^2)*(s2^2) + (D/2.0)*(s1^2 + s2^2) + (h/2.0)*(s1+s2)
      
      if (s1==s2 && abs(s1) == 0) 
        H[a, a] -= 2*β
      end
      
      if (s1==-s2 && s1 !=0)
        H[a,a] -= β
      end

      if(s1==0 && abs(s2)==1) || (abs(s1)==1 && s2==0)
        H[a,a] -= β
      end

      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        sb1 = spin(b, n, i)
        sb2 = spin(b, n, j)
        H[a, b] += 1 - β*sb1*sb2-β*s1*s2
      end

      doubcoupstates = doubspinflips(a, n, i, j)
      for b in doubcoupstates
        H[a,b] += -β
      end

    end
  end
  return H
end

"""
# initMagBlock - 2 method function

## Methods

  * inith(Bool, Int, Int)
  * inith(Bool, Int, Int, Float64)

### Arugments

  * bc    -   boundary conditions on the lattice
  * n     -   number of sites on lattice
  * mz    -   z-direction magnetization to be constructed
  * D     -   uniaxial anisotropy term

### Returns

  * Array{Float, 2}   -   Hamiltonian of spin1 chain

## Description 

This function create the Hamiltonian of a spin-1 AFM Heisenberg chain with OBC
or PBC.
"""
function initMagBlock(bc::Bool, n::Int, mz::Int)
  # select all states with the given magnetization to construct their subblock
  mz_states = Int64[]
  for a=1:3^n
    tot::Int = 0
    for i=1:n
      tot += spin(a, n, i)
    end
    if tot == mz
      push!(mz_states, a)
    end
  end
  H = zeros(length(mz_states), length(mz_states))
  for a=1:length(mz_states)
    bc ? m=n : m=n-1  # select appropriate boundary conditions
    # compute matrix elements for this state
    for i=1:m
      j = (i%n)+1
      s1 = spin(mz_states[a], n, i)
      s2 = spin(mz_states[a], n, j)
      H[a, a] += s1*s2
      coupstates = spinflips(mz_states[a], n, i, j)
      for bstate in coupstates
        b = findfirst(mz_states, bstate)
        H[a, b] += 1
      end
    end
  end
  return H
end

function initMagBlock(bc::Bool, n::Int, mz::Int, D::Float64)
  # select all states with the given magnetization to construct their subblock
  mz_states = Int64[]
  for a=1:3^n
    tot::Int = 0
    for i=1:n
      tot += spin(a, n, i)
    end
    if tot == mz
      push!(mz_states, a)
    end
  end
  H = zeros(length(mz_states), length(mz_states))
  for (i, a) in enumerate(mz_states)
    bc ? m=n : m=n-1  # select appropriate boundary conditions
    # compute matrix elements for this state
    for i=1:m
      j = (i%n)+1
      s1 = spin(a, n, i)
      s2 = spin(a, n, j)
      H[i, i] += s1*s2 + (D/2.0)*(s1^2 + s2^2)
      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        j = findfirst(mz_states, b)
        H[i, j] += 1
      end
    end
  end
  return H
end


"""
initσz - initialize the Pauli z matrix

"""
function initσz(n::Int, i::Int)

  id = [[1,0,0] [0,1,0] [0,0,1]]
  σz = [[1,0,0] [0,0,0] [0,0,-1]]

  if i==1
    tot = σz
    for i=2:n
      tot = kron(tot,id)
    end
    return tot
  else
    tot = id
    for j=2:(i-1)
      tot = kron(tot,id)
    end
    tot = kron(tot,σz)
    for j=(i+1):n
      tot = kron(tot,id)
    end
    return tot
  end
end


function initσx(n::Int, i::Int)

  e = sqrt(2)/2

  id = [[1,0,0] [0,1,0] [0,0,1]]
  σx = [[0,e,0] [e,0,e] [0,e,0]]

  if i==1
    tot = σx
    for i=2:n
      tot = kron(tot,id)
    end
    return tot
  else
    tot = id
    for j=2:(i-1)
      tot = kron(tot,id)
    end
    tot = kron(tot,σx)
    for j=(i+1):n
      tot = kron(tot,id)
    end
    return tot
  end
end


"""
initmagop - Initialize magnetization operator
"""
function initmagop(k::Int, n::Int)
  tot=0
  for r=1:n

  end
end

end # module Spin1Utls

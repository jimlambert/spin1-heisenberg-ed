
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
  coupstates::Array{Int64, 1}
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

"""
# inith - 1 method function

#
"""
function inith(bc::Bool, n::Int)
  H = zeros(3^n, 3^n) # initialize empty Hamiltonian
  for a=1:3^n
    bc ? m=n : m=n-1 # select number of bonds based on boundary conditions
    for i=1:m
      j = (1%n)+1 # determine nearest neighbour
      H[a, a] += spin(a, n, i)*spin(a, n, j)
      coupstates = spinflips(a, n, i, j)
      for b in coupstates
        H[a, b] += 1
      end
    end
  end
  return H
end



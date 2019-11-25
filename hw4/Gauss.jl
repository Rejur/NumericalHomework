struct EndPt
    label::Char
end

const neither = EndPt('N')
const left    = EndPt('L')
const right   = EndPt('R')
const both    = EndPt('B')

maxiterations = Dict(Float32 => 30, Float64 => 30, BigFloat => 40)
custom_gauss_rule(-1.0, 1.0, 2.0, both)

function custom_gauss_rule(lo::T, hi::T, 
    a::Array{T,1}, b::Array{T,1}, muzero::T, endpt::EndPt,
    maxits::Integer=maxiterations[T]) where {T <: AbstractFloat}
#
# On entry:
#
# a, b hold the coefficients (as given, for instance, by
# legendre_coeff) in the three-term recurrence relation
# for the orthonormal polynomials p_0, p_1, p_2, ... , that is,
#
#    b[j] p (x) = (x-a[j]) p   (x) - b[j-1] p   (x).
#          j                j-1              j-2
#      
# muzero holds the zeroth moment of the weight function, that is
#
#              / hi
#             |
#    muzero = | w(x) dx.
#             |
#             / lo
#
# On return:
#
# x, w hold the points and weights.
#
n = length(a)
@assert length(b) == n
if endpt == left 
   if n == 1
       a[1] = lo
   else
       a[n] = solve(n, lo, a, b) * b[n-1]^2 + lo
   end
elseif endpt == right
   if n == 1
       a[1] = hi
   else
       a[n] = solve(n, hi, a, b) * b[n-1]^2 + hi
   end
elseif endpt == both
   if n == 1 
       error("Must have at least two points for both ends.")
   end 
   g = solve(n, lo, a, b)
   t1 = ( hi - lo ) / ( g - solve(n, hi, a, b) )
   b[n-1] = sqrt(t1)
   a[n] = lo + g * t1
end
w = zero(a)
steig!(a, b, w, maxits)
for i = 1:n
   w[i] = muzero * w[i]^2
end
idx = sortperm(a)
return a[idx], w[idx]
end
function steig!(d::Array{T,1}, e::Array{T,1}, 
    z::Array{T,1}, maxits::Integer) where {T <: AbstractFloat}
#
# Finds the eigenvalues and first components of the normalised
# eigenvectors of a symmetric tridiagonal matrix by the implicit
# QL method.
#
# d[i]   On entry, holds the ith diagonal entry of the matrix. 
#        On exit, holds the ith eigenvalue.
#
# e[i]   On entry, holds the [i+1,i] entry of the matrix for
#        i = 1, 2, ..., n-1.  (The value of e[n] is not used.)
#        On exit, e is overwritten.
#
# z[i]   On exit, holds the first component of the ith normalised
#        eigenvector associated with d[i].
#
# maxits The maximum number of QL iterations.
#
# Martin and Wilkinson, Numer. Math. 12: 377-383 (1968).
# Dubrulle, Numer. Math. 15: 450 (1970).
# Handbook for Automatic Computation, Vol ii, Linear Algebra, 
#        pp. 241-248, 1971.
#
# This is a modified version of the Eispack routine imtql2.
#
n = length(z)
z[1] = 1
z[2:n] = 0
e[n] = 0

if n == 1 # Nothing to do for a 1x1 matrix.
return
end
for l = 1:n
for j = 1:maxits
# Look for small off-diagonal elements.
m = n
for i = l:n-1
if abs(e[i]) <= eps(T) * ( abs(d[i]) + abs(d[i+1]) )
m = i
break   
end
end
p = d[l]
if m == l
continue
end
if j == maxits
msg = ("No convergence after %d iterations", j)
msg *= " (try increasing maxits)"
error(msg)
end
# Form shift
g = ( d[l+1] - p ) / ( 2 * e[l] )
r = hypot(g, one(T))
g = d[m] - p + e[l] / ( g + copysign(r, g) )
s = one(T)
c = one(T)
p = zero(T)
for i = m-1:-1:l
f = s * e[i]
b = c * e[i]
if abs(f) <  abs(g)
s = f / g
r = hypot(s, one(T))
e[i+1] = g * r
c = one(T) / r
s *= c
else
c = g / f
r = hypot(c, one(T))
e[i+1] = f * r
s = one(T) / r
c *= s
end 
g = d[i+1] - p
r = ( d[i] - g ) * s + 2 * c * b
p = s * r
d[i+1] = g + p
g = c * r - b
# Form first component of vector.
f = z[i+1]
z[i+1] = s * z[i] + c * f
z[i]   = c * z[i] - s * f
end # loop over i
d[l] -= p
e[l] = g
e[m] = zero(T)
end # loop over j
end # loop over l
end

function solve(n::Integer, shift::T, 
                            a::Array{T,1}, b::Array{T,1}) where {T <: AbstractFloat}
#
# Perform elimination to find the nth component s = delta[n]
# of the solution to the nxn linear system
#
#     ( J_n - shift I_n ) delta = e_n,
#
# where J_n is the symmetric tridiagonal matrix with diagonal
# entries a[i] and off-diagonal entries b[i], and e_n is the nth
# standard basis vector.
#
t = a[1] - shift
for i = 2:n-1
   t = a[i] - shift - b[i-1]^2 / t
end
return one(t) / t
end
function orthonormal_poly(x::Array{T,1}, 
    a::Array{T,1}, b::Array{T,1}, muzero::T) where {T<:AbstractFloat}
# p[i,j] = value at x[i] of orthonormal polynomial of degree j-1.
m = length(x)
n = length(a)
p = zeros(T, m, n+1)
c = one(T) / sqrt(muzero)
rb = one(T) / b[1]
for i = 1:m
    p[i,1] = c
    p[i,2] = rb * ( x[i] - a[1] ) * c
end 
for j = 2:n
    rb = one(T) / b[j]
    for i = 1:m
    p[i,j+1] = rb * ( (x[i]-a[j]) * p[i,j] - b[j-1] * p[i,j-1] )
    end
end
return p
end
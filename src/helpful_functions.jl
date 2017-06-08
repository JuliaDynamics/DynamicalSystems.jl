export is1D

"""
    A_mul_signB!(A, B)
Perform `A = A * diagm(sign(diag(B)))` in-place for A without allocations.
"""
A_mul_signB!() = nothing
@inbounds function A_mul_signB!(A, B)
  @simd for i = 1:size(A,2)
    B[i,i] < 0 && (A[:, i] .*= -1.0)
  end
  return A
end

"""
    signB_mul_A!(B, A)
Perform `A = diagm(sign(diag(B))) * A` in-place for A without allocations.
"""
signB_mul_A!() = nothing
@inbounds function signB_mul_A!(B, A)
  @simd for i = 1:size(A,1)
    B[i,i] < 0 && (A[i, :] .*= -1.0)
  end
  return A
end

is1D(::Any) = false
is1D(::AbstractArray) = length(size(u0)) == 1

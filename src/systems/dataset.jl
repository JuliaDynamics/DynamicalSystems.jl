# A dataset will be a vector of vectors
# you can always move quickly from a
# vector of static arrays to a matrix using reinterpret though:
<< v = [rand(SVector{3}) for i in 1:1000];

<< @time V = reinterpret(Float64, v, (3, 1000));
  0.000009 seconds (7 allocations: 288 bytes)

<< typeof(V)
>> Array{Float64,2}

# and then back,
<< reinterpret(SVector{3, Float64}, V, (1000,))
>> 1000-element Array{SVector{3,Float64},1}:
 [0.511705, 0.880423, 0.502126]
...

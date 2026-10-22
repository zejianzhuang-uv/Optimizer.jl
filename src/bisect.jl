


function bisectv(f::Function, init_x::Vector{Float64}, n::Int64)
    x0 = init_x[1]
    solv = Float64[]
    for x in init_x[2:end]
        xx = f(x0) * f(x)
        if xx < 0
            sol = bisect(f, x0, x)
            if isnan(sol) == false
                push!(solv, sol)
                if length(solv) == n
                    break
                end
                x0 = x
            else
                x0 = x
            end
        else
            continue
        end
    end
    return solv
end



















"""
- Written by Charles Harris charles.harris@sdl.usu.edu
- Source code: https://github.com/willb/scipy/blob/v1.2.1/scipy/optimize/Zeros/bisect.c
- xtol=2e-12, rtol=1e-10
"""
function bisect(f, xa::Real, xb::Real; iter::Int=1000, xtol=1e-8, rtol=1e-8, ftol=1e-8)
    fa = f(xa)
    abs(fa) <= ftol && return float(xa)
    fb = f(xb)
    abs(fb) <= ftol && return float(xb)
    signbit(fa) == signbit(fb) && error("f(xa) and f(xb) have the same sign.")

    dm = xb - xa
    xm = xa
    for _ in 1:iter
        dm *= 0.5
        xm = xa + dm
        fm = f(xm)
        signbit(fm) == signbit(fa) && (xa = xm)
        (abs(fm) <= ftol || abs(dm) < xtol + rtol * abs(xm)) && return float(xm)
    end
    @warn "bisect: reached maxiter" iter
    return float(xm)
end

# function bisect(f::Function, xa::Float64, xb::Float64; iter::Int64=1000, xtol=2e-12, rtol=1e-10, ftol=1e-6)::Float64
#     # ftol = 1e-6
#     i::Int64 = 1
#     dm = 0e0
#     xm = 0e0
#     fm = 0e0
#     fa = 0e0
#     fb = 0e0
#     fa = f(xa)
#     fb = f(xb)
#     if abs(fa) <= ftol
#         return xa
#     end
#     if abs(fb) <= ftol
#         return xb
#     end
#     if signbit(fa) == signbit(fb)
#         error("f(xa) and f(xb) has the same sign.")
#         return 0e0
#     end
#     dm = xb - xa
#     while i < iter
#         dm *= 0.5
#         xm = xa + dm
#         fm = f(xm)
#         if signbit(fm) == signbit(fa)
#             xa = xm
#         end
#         if (abs(fm) <= ftol) || (abs(dm) < xtol + rtol*abs(xm))
#             # return xm
#             xa = xm
#             break
#         end
#         i += 1
#     end
#     if i >= iter
#         println("iteration reaches its maximum.")
#     end
#     xa = (abs(xa) < 1e-6 ? 0e0 : xa) 
#     if abs(f(xa)) > 0.1
#         return NaN64
#     else
#         return xa
#     end
# end

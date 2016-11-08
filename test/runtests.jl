
module TWPBVP_Test
	using TWPBVP
	const ϵ = 1e-2

	function f(x, y, dy)
		dy[1] = y[2]
		dy[2] = 1/ϵ*(-exp(y[1]) + 0.5π*sin(0.5π*x)*exp(2y[1]))
	end

	function df(x, y, jac)
		jac[1,1] = 0.0
		jac[1,2] = 1.0
		jac[2,1] = 1/ϵ*(-exp(y[1])+π*sin(0.5π*x)*exp(2y[1]))
		jac[1,2] = 0.0
	end

	function g(i, y)
		y[1]
	end

	function dg(i, y, jac)
		jac[1] = 1.0
		jac[2] = 0.0
	end

	# Other constants & values
	const fxpoint = Nullable{Vector{Float64}}()
	const ltol = [1]
	const tol = [1e-2]
	const u = Array(Float64, 2, 1000)
	const xx = Array(Float64, 1000)
	const wrk = Array(Float64, 10000)
	const iwrk = Array(Int64, 6000)
	const num_left_bc = 1
	const left = 0.0
	const right = 1.0
	const linear = false
	const giveu = false
	const givmsh = false
	const nmesh = Ref(0)
	const nmax = Ref(0)
	const κ1 = Ref(1.0)
	const κ = Ref(1.0)
	const γ1 = Ref(1.0)
	const result = Ref(0)

	perform() =
		twpbvpc(num_left_bc, left, right, fxpoint, ltol, tol, linear, givmsh, giveu, nmesh, xx, u, nmax, wrk, iwrk, f, df, g, dg, κ1, γ1, κ, result)

end

using Base.Test
import TWPBVP_Test

TT = TWPBVP_Test

TT.perform()

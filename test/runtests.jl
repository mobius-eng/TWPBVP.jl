
module TWPBVP_Test1
	using TWPBVP
	const ϵ = 1e-2

	function f(x, y, dy)
		println("f: input x = $x, y = $y")
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

module TWPBVP_Test2

	using TWPBVP

	function f(t,y, dy)
		dy[1] = y[2]
		dy[2] = -y[1]
	end

	function df(t, y, jac)
		jac[1,1] = 0.0
		jac[1,2] = 1.0
		jac[2,1] = -1.0
		jac[2,2] = 0.0
	end

	g(i, y) = i == 1 ? y[1] - 1.0 : y[1] + 1.0

	function dg(i, t, y, grad)
		gard[:] = [1.0, 0.0]
	end

	function perform()

		fixpoints = Nullable{Vector{Float64}}()
		local ltol :: Vector{FInt} = [2]
		tol = [1e-6]
		u = Array(Float64, 2, 1000)
		xx = Array(Float64, 1000)

		wrk = Array(Float64, 10000)
		iwrk = Array(FInt, 6000)
		local num_left_bc :: FInt = 1
		left = 0.0
		right = float(π)
		linear = false
		giveu = true

		givmsh = true

		local nmesh0 :: FInt = 20

		nmesh = Ref{FInt}(nmesh0)

		nmax = Ref{FInt}(0)
		κ1 = Ref(1.0)
		κ = Ref(1.0)
		γ1 = Ref(1.0)
		result = Ref{FInt}(-100)

		xx[1:nmesh0] = linspace(0.0, π, nmesh0)

		for i=1:nmesh0
			u[:, i] = [-2/π*xx[i] + 1, -2.0/π]
		end
		println("xx0 = $(xx[1:nmesh0+2])")

		twpbvpc(num_left_bc, left, right, fixpoints, ltol, tol, linear, givmsh, giveu, nmesh, xx, u, nmax, wrk, iwrk, f, df, g, dg, κ1, γ1, κ, result)

		println("Result is $(result[])")
		println("nmesh = $(nmesh[])")

	end

end

using Base.Test
import TWPBVP_Test2

TT = TWPBVP_Test2

TT.perform()

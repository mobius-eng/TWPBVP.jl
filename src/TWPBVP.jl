__precompile__()
module TWPBVP

export twpbvpc

const libtwpbvpc = joinpath(Pkg.dir("TWPBVP"), "deps", "libtwpbvpc.so")


"""
Thin layer of compatibility
"""
function twpbvpc(nlbc :: Int64,
		aleft :: Float64, aright :: Float64,
		fixpnt :: Nullable{Vector{Float64}},
		ltol :: Vector{Int64}, tol :: Vector{Float64},
		linear :: Bool, givmsh :: Bool, giveu :: Bool, nmsh :: Ref{Int64},
		xx :: Vector{Float64}, u :: Array{Float64, 2}, nmax :: Ref{Int64},
		wrk :: Vector{Float64}, iwrk :: Vector{Int64},
		fsub :: Function, dfsub :: Function,
		gsub :: Function, dgsub :: Function,
		ckappa1 :: Ref{Float64}, gamma1 :: Ref{Float64},
		ckappa :: Ref{Float64},
		# rpar :: Vector{Float64},
		# ipar :: Vector{Int64},
		iflbvp :: Ref{Int64})

	# As Fortran 77 has no clausers, rpar and ipar are used to capture external
	# parameters. In Julia: use clausers
	function fsub_par(n :: Int64, x :: Float64, y :: Vector{Float64}, dy :: Vector{Float64}, rpar :: Vector{Float64}, ipar :: Vector{Float64})
		fsub(x, y, dy)
	end

	function dfsub_par(n :: Int64, x :: Float64, y :: Vector{Float64}, jac :: Array{Float64, 2}, rpar :: Vector{Float64}, ipar :: Vector{Int64})
		dfsub(x, y, jac)
	end

	function gsub_par(i :: Int64, n :: Int64, u :: Vector{Float64}, g :: Ref{Float64}, rpar :: Vector{Float64}, ipar :: Vector{Int64})
		g.x = gsub(i, u)
	end

	function dgsub_par(i :: Int64, n :: Int64, u :: Vector{Float64}, dg :: Vector{Float64}, rpar :: Vector{Float64}, ipar :: Vector{Int64})
		dgsub(i, u, dg)
	end

	cf_fsub = cfunction(fsub_par, Void, (Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}))
	cf_dfsub = cfunction(dfsub_par, Void, (Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}))
	cf_gsub = cfunction(gsub_par, Void, (Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}))
	cf_dgsub = cfunction(dgsub_par, Void, (Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}))

	# Fake external parameters
	# Can't have it 0-length as it would be Any[0] and not Float64[0]
	local rpar :: Vector{Float64} = [0.0]
	local ipar :: Vector{Int64} = [0]

	# No need to pass these parameters
	# u is a matrix for the solution only!
	ncomp, nucol = size(u)
	# Get the maximum of xx
	nxxdim = length(xx)
	# max for mesh points must be the same as the number of column points of u
	assert(nucol == nxxdim)

	# Sizes of work arrays
	lwrkfl = length(wrk)
	lwrkin = length(iwrk)

	# Number of fixed mesh points
	if isnull(fixpnt)
		nfxpnt = 0
		fixpnt_v = [0.0]
	else
		fixpnt_v = get(fixpnt)
		nfxpnt = length(fixpnt_v)
	end

	# Size of tolerance vector ≤ ncomp
	ntol = length(ltol)

	ccall((:twpbvpc_, libtwpbvpc), Void,
		(Ref{Int64}, Ref{Int64},					# ncomp, nlbc,
		Ref{Float64}, Ref{Float64},					# aleft, aright
		Ref{Int64}, Ref{Float64}, 					# nfxpnt, fixpnt
		Ref{Int64}, Ref{Int64}, Ref{Float64},		# ntol, ltol, tol
		Ref{Int64}, Ref{Int64}, Ref{Int64},			# linear, givmsh, giveu
		Ref{Int64}, Ref{Int64},						# nmsh, nxxdim
		Ref{Float64}, Ref{Int64},					# xx, nudim
		Ref{Float64}, Ref{Int64},					# u, nmax
		Ref{Int64}, Ref{Float64},					# lwrkfl, wrk
		Ref{Int64}, Ref{Int64},						# lwrkin, iwrk
		Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void},	# fsub, dfsub, gsub, dgsub
		Ref{Float64}, Ref{Float64},					# ckappa1, gamma1
		Ref{Float64}, Ref{Float64}, Ref{Int64},		# ckappa, rpar, ipar
		Ref{Int64}),								# iflbvp
		ncomp, nlbc, aleft, aright,
		nfxpnt, fixpnt_v, ntol, ltol, tol,
		linear, givmsh, giveu, nmsh,
		nxxdim, xx, nucol, u, nmax,					# nudim = nucol
		lwrkfl, wrk, lwrkin, iwrk,
		cf_fsub, cf_dfsub, cf_gsub, cf_dgsub,
		ckappa1,gamma1,ckappa,rpar,ipar,iflbvp)
end

end # module

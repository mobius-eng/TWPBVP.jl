__precompile__()
module TWPBVP

export twpbvpc

const libtwpbvpc = joinpath(Pkg.dir("TWPBVP"), "deps", "libtwpbvpc.so")

immutable TWPBVPCProblem
	fsub :: Function
	dfsub :: Function
	gsub :: Function
	dgsub :: Function
end

"""
Provides the conversion betweeen Julia's function (including closures)
and Fortran's subroutine. The closure is passed in `rpar` parameter
that Fortran believes to be an array of Float64, it never looks into it
anyway and only passes it around.
"""
function unsafe_fsub(rn :: Ref{Int64}, rx :: Ref{Float64}, py :: Ptr{Float64}, pdy :: Ptr{Float64}, rpar :: Ptr{Float64}, ipar :: Ptr{Int64}) :: Void
	x = rx[]
	n = rn[]
	y = unsafe_wrap(Array, py, n)
	dy = unsafe_wrap(Array, pdy, n)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	problem.fsub(x, y, dy)
	return nothing
end

function unsafe_dfsub(rn :: Ref{Int64}, rx :: Ref{Float64}, ry :: Ptr{Float64}, rjac :: Ptr{Float64}, rpar :: Ptr{Float64}, ipar :: Ptr{Int64}) :: Void
	x = rx[]
	n = rn[]
	y = unsafe_wrap(Array, py, n)
	jac = unsafe_wrap(Array, pjac, (n,n))
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	problem.dfsub(x, y, jac)
	return nothing
end

function unsafe_gsub(ri :: Ref{Int64}, rn :: Ref{Int64}, py :: Ptr{Float64}, pg :: Ref{Float64}, rpar :: Ptr{Float64}, ipar :: Ptr{Int64}) :: Void
	i = ri[]
	n = rn[]
	y = unsafe_wrap(Array, py, n)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	rg[] = problem.gsub(i, y)
	return nothing
end

function unsafe_dgsub(ri :: Ref{Int64}, rn :: Ref{Int64}, py :: Ptr{Float64}, pdg :: Ptr{Float64}, rpar :: Ptr{Float64}, ipar :: Ptr{Int64}) :: Void
	i = ri[]
	n = rn[]
	y = unsafe_wrap(Array, py, n)
	dg = unsafe_wrap(Array, pdg, n)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	rg[] = problem.dgsub(i, y, dg)
	return nothing
end

const fsub_ptr = cfunction(unsafe_fsub, Void, (Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}))

const dfsub_ptr = cfunction(unsafe_dfsub, Void, (Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}))

const gsub_ptr = cfunction(unsafe_gsub, Void, (Ref{Int64}, Ref{Int64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Int64}))

const dgsub_ptr = cfunction(unsafe_dgsub, Void, (Ref{Int64}, Ref{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}))

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

	# Keep problem functions in rpar
	# HACK!
	rpar = TWPBVPCProblem(fsub, dfsub, gsub, dgsub)

	# Fake external parameters
	# Can't have it 0-length as it would be Any[0] and not Float64[0]
	# local rpar :: Vector{Float64} = [0.0]
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
		Ref{Int64}, Ptr{Float64}, 					# nfxpnt, fixpnt
		Ref{Int64}, Ptr{Int64}, Ptr{Float64},		# ntol, ltol, tol
		Ref{Int64}, Ref{Int64}, Ref{Int64},			# linear, givmsh, giveu
		Ref{Int64}, Ref{Int64},						# nmsh, nxxdim
		Ptr{Float64}, Ref{Int64},					# xx, nudim
		Ptr{Float64}, Ref{Int64},					# u, nmax
		Ref{Int64}, Ptr{Float64},					# lwrkfl, wrk
		Ref{Int64}, Ptr{Int64},						# lwrkin, iwrk
		Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void},	# fsub, dfsub, gsub, dgsub
		Ref{Float64}, Ref{Float64},					# ckappa1, gamma1
		Ref{Float64}, Ptr{Void}, Ptr{Int64},		# ckappa, rpar, ipar
		Ref{Int64}),								# iflbvp
		ncomp, nlbc, aleft, aright,
		nfxpnt, fixpnt_v, ntol, ltol, tol,
		linear, givmsh, giveu, nmsh,
		nxxdim, xx, nucol, u, nmax,					# nudim = nucol
		lwrkfl, wrk, lwrkin, iwrk,
		fsub_ptr, dfsub_ptr, gsub_ptr, dgsub_ptr,
		ckappa1,gamma1,ckappa,pointer_from_objref(rpar),ipar,iflbvp)
end

end # module

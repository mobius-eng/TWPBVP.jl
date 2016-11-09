__precompile__()
module TWPBVP

export twpbvpc, FInt

FInt = Int64

type TWPBVPCProblem
	fsub :: Function
	dfsub :: Function
	gsub :: Function
	dgsub :: Function
end

function __init__()
	global const libtwpbvpc = joinpath(Pkg.dir("TWPBVP"), "deps", "libtwpbvpc.so")
	global const fsub_ptr = cfunction(unsafe_fsub, Void, (Ref{FInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{FInt}))

	global const dfsub_ptr = cfunction(unsafe_dfsub, Void, (Ref{FInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{FInt}))

	global const gsub_ptr = cfunction(unsafe_gsub, Void, (Ref{FInt}, Ref{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{FInt}))

	global const dgsub_ptr = cfunction(unsafe_dgsub, Void, (Ref{FInt}, Ref{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{FInt}))
end


"""
Provides the conversion betweeen Julia's function (including closures)
and Fortran's subroutine. The closure is passed in `rpar` parameter
that Fortran believes to be an array of Float64, it never looks into it
anyway and only passes it around.
"""
function unsafe_fsub(n :: FInt, x :: Float64, py :: Ptr{Float64}, pdy :: Ptr{Float64}, rpar :: Ptr{Void}, ipar :: Ptr{FInt}) :: Void
	y = unsafe_wrap(Array, py, n)
	dy = unsafe_wrap(Array, pdy, n)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	problem.fsub(x, y, dy)
	return nothing
end

function unsafe_dfsub(n :: FInt, x :: Float64, py :: Ptr{Float64}, pjac :: Ptr{Float64}, rpar :: Ptr{Void}, ipar :: Ptr{FInt}) :: Void
	y = unsafe_wrap(Array, py, n)
	jac = unsafe_wrap(Array, pjac, (n,n))
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	problem.dfsub(x, y, jac)
	return nothing
end

function unsafe_gsub(i :: FInt, n :: FInt, py :: Ptr{Float64}, pg :: Ptr{Float64}, rpar :: Ptr{Void}, ipar :: Ptr{FInt}) :: Void
	y = unsafe_wrap(Array, py, n)
	g = unsafe_wrap(Array, pg, 1)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	g[1] = problem.gsub(i, y)
	return nothing
end

function unsafe_dgsub(i :: FInt, n :: FInt, py :: Ptr{Float64}, pdg :: Ptr{Float64}, rpar :: Ptr{Void}, ipar :: Ptr{FInt}) :: Void
	y = unsafe_wrap(Array, py, n)
	dg = unsafe_wrap(Array, pdg, n)
	problem = unsafe_pointer_to_objref(rpar) :: TWPBVPCProblem
	problem.dgsub(i, y, dg)
	return nothing
end

"""
Thin layer of compatibility
"""
function twpbvpc(nlbc :: FInt,
		aleft :: Float64, aright :: Float64,
		fixpnt :: Nullable{Vector{Float64}},
		ltol :: Vector{FInt}, tol :: Vector{Float64},
		linear :: Bool, givmsh :: Bool, giveu :: Bool, nmsh :: Ref{FInt},
		xx :: Vector{Float64}, u :: Array{Float64, 2}, nmax :: Ref{FInt},
		wrk :: Vector{Float64}, iwrk :: Vector{FInt},
		fsub :: Function, dfsub :: Function,
		gsub :: Function, dgsub :: Function,
		ckappa1 :: Ref{Float64}, gamma1 :: Ref{Float64},
		ckappa :: Ref{Float64},
		# rpar :: Vector{Float64},
		# ipar :: Vector{FInt},
		iflbvp :: Ref{FInt})

	# Keep problem functions in rpar
	# HACK!
	rpar = TWPBVPCProblem(fsub, dfsub, gsub, dgsub)

	# Fake external parameters
	# Can't have it 0-length as it would be Any[0] and not Float64[0]
	# local rpar :: Vector{Float64} = [0.0]
	local ipar :: Vector{FInt} = [0]

	# No need to pass these parameters
	# u is a matrix for the solution only!
	ncomp, nucol = size(u)
	nudim = ncomp
	# Get the maximum of xx
	nxxdim = length(xx)
	# max for mesh points must be the same as the number of columns of u
	assert(nucol == nxxdim)

	# Sizes of work arrays
	lwrkfl = length(wrk)
	lwrkin = length(iwrk)

	# Number of fixed mesh points
	if isnull(fixpnt)
		nfxpnt = 0
		fixpnt_v = [0.0]
	else
		fixpnt_v = fixpnt[]
		nfxpnt = length(fixpnt_v)
	end

	# Size of tolerance vector ≤ ncomp
	ntol = length(ltol)

	println("tol[1] = $(tol[1])")
	println("xx[2] = $(xx[2])")
	println("aright = $aright")

	varight = [aright]

	ccall((:twpbvpc_, libtwpbvpc), Void,
		(Ref{FInt}, Ref{FInt},					# ncomp, nlbc,
		Ref{Float64}, Ptr{Float64},					# aleft, aright
		Ref{FInt}, Ptr{Float64}, 					# nfxpnt, fixpnt
		Ref{FInt}, Ptr{FInt}, Ptr{Float64},		# ntol, ltol, tol
		Ref{FInt}, Ref{FInt}, Ref{FInt},			# linear, givmsh, giveu
		Ref{FInt}, Ref{FInt},						# nmsh, nxxdim
		Ptr{Float64}, Ref{FInt},					# xx, nudim
		Ptr{Float64}, Ref{FInt},					# u, nmax
		Ref{FInt}, Ptr{Float64},					# lwrkfl, wrk
		Ref{FInt}, Ptr{FInt},						# lwrkin, iwrk
		Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void},	# fsub, dfsub, gsub, dgsub
		Ref{Float64}, Ref{Float64},					# ckappa1, gamma1
		Ref{Float64}, Ptr{Void}, Ptr{FInt},		# ckappa, rpar, ipar
		Ref{FInt}),								# iflbvp
		ncomp, nlbc,
		aleft, varight,
		nfxpnt, fixpnt_v,
		ntol, ltol, tol,
		linear, givmsh, giveu, nmsh,
		nxxdim, xx,
		nudim, u, nmax,					# nudim = ncomp
		lwrkfl, wrk,
		lwrkin, iwrk,
		fsub_ptr, dfsub_ptr, gsub_ptr, dgsub_ptr,
		ckappa1, gamma1, ckappa,
		pointer_from_objref(rpar), ipar,
		iflbvp)
end

end # module

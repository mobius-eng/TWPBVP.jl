cd(joinpath(Pkg.dir("TWPBVP"), "deps"))
pic = @windows ? "" : "-fPIC"
run(`gfortran -m$WORD_SIZE -ffixed-form $pic -shared -O2 -o libtwpbvpc.so twpbvpc.f`)
# -fdefault-integer-8
# -fdefault-real-8
# -fdefault-integer-8

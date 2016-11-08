cd(joinpath(Pkg.dir("TWPBVP"), "deps"))
pic = @windows ? "" : "-fPIC"
run(`gfortran -m$WORD_SIZE -fdefault-real-8 -fdefault-integer-8 -ffixed-form $pic -shared -O3 -o libtwpbvpc.so twpbvpc.f`)

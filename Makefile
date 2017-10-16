
all: pgcalc

pgcalc: deps
	julia juliac.jl -ve pgcalc.jl

deps:
	julia deps.jl

clean:
	rm -r builddir

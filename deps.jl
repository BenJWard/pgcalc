#!/usr/bin/env julia

macro pkgcheck(pkgname)
    quote
        if !($pkgname ∈ keys(Pkg.installed()))
            Pkg.add($pkgname)
        end
    end
end

Pkg.update()

@pkgcheck "ArgParse"
@pkgcheck "BuildExecutable"
@pkgcheck "BioSequences"
@pkgcheck "GeneticVariation"
@pkgcheck "NaturalSelection"

Pkg.update()

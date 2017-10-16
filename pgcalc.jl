module pgcalc

using ArgParse
using BioSequences

using GeneticVariation
using NaturalSelection

macro header_generator(name, opts, fixed)
    optionals = [Expr(:call, :ifelse, :((flags & $(0x01 << (i - 1))) > 0x00), opts.args[i], "") for i in eachindex(opts.args)]
    println_expr = Expr(:call, :println, :out, optionals..., fixed)
    return Expr(:function, :($name(flags::UInt8, out::IOStream)) , println_expr)
end

function parsecmdline(args)
    s = ArgParse.ArgParseSettings()
    s.add_help = true

    ArgParse.@add_arg_table s begin
        "--outdir", "-o"
            help = "Specify an output directory for table files."
            arg_type = String
            default = "./"
        "input"
            help = "A FASTA input file containing sequence alignment."
            arg_type = String
            required = true
        "--sinput", "-I"
            help = "The FASTA input file containing similated sequence alignments."
            arg_type = String
            default = ""
        "sample"
            help = "Compute statistics that apply to the whole sample."
            action = :command
        "pairwise"
            help = "Compute pairwise statistics for your sequence sample."
            action = :command
        "dndscodon"
            help = "Compute pairwise dNdS according to Nei & Gojoborei's method, on a per codon basis."
            action = :command
    end

    ArgParse.@add_arg_table s["sample"] begin
        "--seg", "-s"
            help = "Compute the number of segregating sites."
            action = :store_true
        "--nucdiv", "-n"
            help = "Compute the nucleotide diversity of the sample."
            action = :store_true
        "--tajima", "-d"
            help = "Compute Tajima's D statistic for the sample."
            action = :store_true
    end

    ArgParse.@add_arg_table s["pairwise"] begin
        "--dnds", "-d"
            help = "Compute the pairwise dNdS according to Nei & Gojoborei's method."
            action = :store_true
        "--dndscor", "-D"
            help = "Compute a different dNdS."
            action = :store_true
        "--mut", "-m"
            help = "Compute the number of mutations between each pair of sequences."
            action = :store_true
    end

    return ArgParse.parse_args(args, s, as_symbols = true)
end

function extract_op_flags(args::Dict{Symbol,Any}, flags::Vararg{Symbol})
    bitflags = 0x00
    @inbounds for i in eachindex(flags)
        flag = args[flags[i]]::Bool
        bitflags |= reinterpret(UInt8, flag) << (i - 1)
    end
    return bitflags
end

function compute_pairwise_stats(seqs)
    mcounts = count_pairwise(Mutated, seqs...)
    dnds = NaturalSelection.pairwise_dNdS_NG86(seqs)

    return mcounts, dnds
end

function fill_from_records!(names, seqs, records)
    @inbounds for i in eachindex(seqs)
        names[i] = FASTA.identifier(records[i])
        seqs[i] = FASTA.sequence(DNASequence, records[i])
    end
end

function process_sims(input::IO, sout::IO, pout::IO, records::Vector{FASTA.Record}, seq_names::Vector{String}, seq_store::Vector{DNASequence})
    simrep = 1
    while !eof(input)
        @inbounds for i in 1:endof(records)
            read!(input, records[i])
        end
        fill_from_records!(snames, seqs, recs)

        seg, pi, td = compute_sample_stats(seq_store)
        mcounts, dnds = compute_pairwise_stats(seq_store)

        write_rep_to_file(output, seq_names, dnds, mcounts, pi, seg, td, simrep)
    end
end

function base_path(dir::String, file::String)
    return abspath(expanduser(joinpath(dir,splitext(basename(file))[1])))
end

function extract_sequences(records::Vector{FASTA.Record})
    seq_names = [FASTA.identifier(rec) for rec in records]
    seq_store = [FASTA.sequence(DNASequence, rec) for rec in records]
    return seq_names, seq_store
end

function head_generator(io::IOStream, flags::UInt8, head_options::Vararg{String})
    @inbounds for i in 1:(endof(head_options) - 1)
        if flags & (0x01 << (i - 1)) > 0x00
            print(io, head_options[i], ", ")
        end
    end
    if flags & (0x01 << (endof(head_options) - 1)) > 0x00
        print(io, head_options[endof(head_options)])
    end
end

function compute_sample_stats(seqs::Vector{BioSequences.DNASequence},
                              out::IOStream, comp::UInt8, rep::Int)
    if comp == 0x00
        error("You must specify a sample statistic to calculate!")
    end
    if comp & 0x01 > 0x00
        seg = count(GeneticVariation.Segregating, seqs)
        print(out, seg[1], ", ")
    end
    if comp & 0x02 > 0x00
        print(out, GeneticVariation.NL79(seqs), ", ")
    end
    if comp & 0x04 > 0x00
        print(out, NaturalSelection.tajimad(seqs), ", ")
    end
    println(out, rep)
end

function compute_pairwise_stats(names::Vector{String}, seqs::Vector{BioSequences.DNASequence},
                                out::IOStream, comp::UInt8, dndscor::Bool, rep::Int)
    if comp == 0x00
        error("You must specify a sample statistic to calculate!")
    end
    if comp & 0x01 > 0x00
        mut = BioSequences.count_pairwise(GeneticVariation.Mutated, seqs...)
    end
    if comp & 0x02 > 0x00
        dnds = NaturalSelection.pairwise_dNdS_NG86(seqs, 1.0, BioSequences.ncbi_trans_table[1], dndscor)
    end

    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        print(out, names[i], ", ", names[j], ", ")
        if comp & 0x01 > 0x00
            print(out, mut[i, j][1], ", ")
        end
        if comp & 0x02 > 0x00
            dN, dS = dnds[i, j]
            print(out, dN, ", ", dS, ", ")
        end
        println(out, rep)
    end
end

function compute_codon_stats(names::Vector{String}, seqs::Vector{BioSequences.DNASequence},
                                out::IOStream, rep::Int)

    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        xcdns, ycdns = NaturalSelection.aligned_codons(seqs[i], seqs[j])

        pos = 1

        for (x, y) in zip(xcdns, ycdns)
            sx, nx = NaturalSelection.S_N_NG86(x, 1.0, BioSequences.ncbi_trans_table[1])
            sy, ny = NaturalSelection.S_N_NG86(y, 1.0, BioSequences.ncbi_trans_table[1])
            S = (sx + sy) / 2.0
            N = (nx + ny) / 2.0
            DS, DN = NaturalSelection.DS_DN_NG86(x, y, BioSequences.ncbi_trans_table[1])
            pS = DS / S
            pN = DN / N

            println(out, names[i], ", ", names[j], ", ", pos, ", ",
                    S, ", ", N, ", ",
                    DS, ", ", DN, ", ",
                    pS, ", ", pN, ", ", rep)

            pos += 1

        end
    end

end

function pgcalc_main(args::Vector{String})
    args = parsecmdline(args)

    # Process the inputs and generate the filestreams.
    input_name::String = args[:input]
    simulation_name::String = args[:sinput]
    cmdname::String = args[:_COMMAND_]
    output_dir::String = args[:outdir]
    process_simulations = simulation_name != ""

    input_records = open(realpath(input_name), "r") do input_stream
        collect(FASTA.Reader(input_stream))
    end

    output_base = base_path(output_dir, input_name)
    output_path = string(output_base, "_", cmdname, ".csv")
    output_stream = open(output_path, "w")

    seq_names, seq_store = extract_sequences(input_records)

    if cmdname == "sample"

        sample_args::Dict{Symbol, Any} = args[:sample]
        flags = extract_op_flags(sample_args, :seg, :nucdiv, :tajima)
        head_generator(output_stream, flags, "n_segregating", "nuc_div", "tajima_d")
        println(output_stream, ", simrep")
        compute_sample_stats(seq_store, output_stream, flags, 0)

        if process_simulations
            sim_stream = FASTA.Reader(open(realpath(simulation_name), "r"))
            simrep = 1
            while !eof(sim_stream)
                @inbounds for i in 1:endof(input_records)
                    read!(sim_stream, input_records[i])
                end
                fill_from_records!(seq_names, seq_store, input_records)
                compute_sample_stats(seq_store, output_stream, flags, simrep)
                simrep += 1
            end
        end
    end

    if cmdname == "pairwise"

        pairwise_args::Dict{Symbol, Any} = args[:pairwise]
        flags = extract_op_flags(pairwise_args, :mut, :dnds)
        dnds_addone = pairwise_args[:dndscor]
        print(output_stream, "sequence_1, sequence_2, ")
        head_generator(output_stream, flags, "n_mutations", "dN, dS")
        println(output_stream, ", simrep")
        compute_pairwise_stats(seq_names, seq_store, output_stream, flags, dnds_addone, 0)

        if process_simulations
            sim_stream = FASTA.Reader(open(realpath(simulation_name), "r"))
            simrep = 1
            while !eof(sim_stream)
                @inbounds for i in 1:endof(input_records)
                    read!(sim_stream, input_records[i])
                end
                fill_from_records!(seq_names, seq_store, input_records)
                compute_pairwise_stats(seq_names, seq_store, output_stream, flags, dnds_addone, simrep)
                simrep += 1
            end
        end
    end

    if cmdname == "dndscodon"

        println(output_stream, "sequence_1, sequence_2, codon_i, S, N, DS, DN, pS, pN, simrep")
        compute_codon_stats(seq_names, seq_store, output_stream, 0)

        if process_simulations
            sim_stream = FASTA.Reader(open(realpath(simulation_name), "r"))
            simrep = 1
            while !eof(sim_stream)
                @inbounds for i in 1:endof(input_records)
                    read!(sim_stream, input_records[i])
                end
                fill_from_records!(seq_names, seq_store, input_records)
                compute_codon_stats(seq_names, seq_store, output_stream, simrep)
                simrep += 1
            end
        end
    end
end

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    pgcalc_main(ARGS)
    return 0
end

end # Module

using PackageCompiler
create_sysimage(:Channel,
  sysimage_path=joinpath(@__DIR__,"..","Channel.so"),
  precompile_execution_file=joinpath(@__DIR__,"..","src", "Channel.jl"))
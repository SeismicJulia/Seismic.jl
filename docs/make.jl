using Documenter, Seismic

makedocs(
    modules = [Seismic],
    doctest = false)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "https://github.com/fercarozzi/Seismic.jl.git",
    julia  = "0.6",
    osname = "osx")

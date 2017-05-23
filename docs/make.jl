using Documenter, Seismic

makedocs(
    modules = [Seismic],
    doctest = true)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/msacchi/Seismic.jl.git",
    julia  = "0.5",
    osname = "osx")

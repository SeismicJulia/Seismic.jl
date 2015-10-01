using Seismic,Lexicon

save("docs/user-docs/docs.md",Seismic)
run(`mkdocs build --clean`)


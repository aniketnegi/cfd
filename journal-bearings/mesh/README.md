# Surface meshing

Generate cylindrical or conical quadrilateral housing surfaces, either unported
or fitted around a circular drilling. See [generation commands and geometry
options](2d/README.md#cli-and-examples), or run:

```sh
julia --project=. mesh/2d/generate.jl --help
julia --project=. mesh/tests/runtests.jl
```

Run commands from `journal-bearings/`. Use `-v` for detailed diagnostics or
`-q` for errors only. Select tests by filename:

```sh
julia --project=. mesh/tests/runtests.jl mesh_cli.jl cylindrical.jl
```

- `2d/`: construction, quality and fitted-inlet placement.
- `vtk.jl`: mesh export.
- `examples/`: basic API examples and specialist optimization studies.
- `tests/`: automated checks and explicit snapshot comparisons.

Start with `examples/mesh_unported.jl`, `mesh_adjusted_unported.jl` or
`mesh_inlet.jl`. Count/geometry comparisons, fixed-budget optimizer diagnostics
and the long-run queue are [specialist studies](2d/README.md#optional-numerical-diagnostics),
not prerequisites for generation.

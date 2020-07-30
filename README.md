# DynLoco optimization of dynamic locomotion
# in Julia

Uses JuMP optimization package with rimless wheel, demonstrating walking on a ramp or slope, and short walks.

1. Install Julia, recommended also Jupyter notebook and/or Juno+Atom editing environment
2. Install `DynLoco` as a package.
  * If using command line (REPL), enter the package manager with "]" and
enter `dev https://bitbucket.org/hbcl/dynloco`.
  * Or, if inside a notebook, do the following:
  ```julia
    using Pkg
    Pkg.develop(PackageSpec(url="https://bitbucket.org/hbcl/dynloco"))
```    
3. The file `shortwalks.jl` can be opened in editor, and executed in chunks. Or paste the code into a
Jupyter notebook, and run as cells.

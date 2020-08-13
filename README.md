# DynLoco optimization of dynamic locomotion
# in Julia

Uses JuMP optimization package with rimless wheel, demonstrating walking on a ramp or slope, and short walks.

1. Install Julia, recommended also Jupyter notebook and/or Juno+Atom editing environment
2. Install `DynLoco` as a package.
  * If using command line (REPL), enter the package manager with "]" and enter 
    ```
    develop https://bitbucket.org/hbcl/dynloco
    activate .
    ```
    Then use delete/backspace key to exit package manager.

  * Or, if inside a notebook, do the following:
    ```julia
    using Pkg
    Pkg.develop(PackageSpec(url="https://bitbucket.org/hbcl/dynloco"))
    Pkg.activate(".")
    ```   
    `Pkg` manages dependencies, where `develop` retrieves code from the repository, and `activate` makes sure the appropriate dependencies are loaded and available.

3. The file `shortwalks.jl` can be opened in editor, and executed in chunks. Or paste the code into a
Jupyter notebook, and run as cells. Or in REPL,
    ```julia
    include("shortwalks.jl")
    ```

4. Recommended Julia workflow with [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/). You should develop code within Modules, and use `Revise` to automatically track modifications. This avoids issues with redefining functions within a single REPL session. To install, enter package manager with "]" and enter `add Revise` (then backspace/delete to exist package manager). It is helpful to start `Revise` from every Julia session. Do this by adding `using Revise` to your [`startup.jl`](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default-1) file.


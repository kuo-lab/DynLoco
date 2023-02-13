# DynLoco optimization of dynamic locomotion in Julia

Uses JuMP optimization package with rimless wheel, demonstrating walking on a ramp or slope, and short walks.

* Install Julia, recommended also Jupyter notebook and/or [VS Code](https://code.visualstudio.com/) editing environment
* Install `DynLoco` as a package. (Note: Need to insert an Enter after each command. Bitbucket is not rendering the newlines; you can view the Mardown source to see where the newlines are.)
    * If using command line (REPL), enter the package manager with "]" and enter 

        ```
        develop https://bitbucket.org/hbcl/dynloco
        ```


        ```
        activate .
        ```

    Then use delete/backspace key to exit package manager.
    
    * Or, if inside a notebook, do the following:

        ```
        using Pkg
        ```

        ```
        Pkg.develop(PackageSpec(url="https://bitbucket.org/hbcl/dynloco"))
        ```

        ```
        Pkg.activate(".")
        ```

    `Pkg` manages dependencies, where `develop` retrieves code from the repository, and `activate` makes sure the appropriate dependencies are loaded and available.

* The file `shortwalks.jl` can be opened in editor, and executed in chunks. Or paste the code into a
Jupyter notebook, and run as cells. Or in REPL,

    ```
    include("shortwalks.jl")
    ```

* Recommended Julia workflow with [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/). You should develop code within Modules, and use `Revise` to automatically track modifications. This avoids issues with redefining functions within a single REPL session. To install, enter package manager with "]" and enter `add Revise` (then backspace/delete to exist package manager). It is helpful to start `Revise` from every Julia session. Do this by adding `using Revise` to your [`startup.jl`](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default-1) file.

* Jupyter notebooks provide a fast, accessible interface for demos. Jupyter stands for Julia-Python-R, and requires Python. In package manager, `add IJulia`. Then in Julia REPL,

    ```
    using IJulia
    ```

    ```
    notebook()
    ```

    If you already have a Python installation, you may prefer to install manually, e.g. using Conda. See [documentation](https://github.com/JuliaLang/IJulia.jl).
    
    Jupyter notebooks run in a browser, but a desktop app experience is also available, with [nteract](https://nteract.io/desktop).

## References
Kuo, Arthur D. 2002. “[Energetics of Actively Powered Locomotion Using the Simplest Walking Model](https://www.ncbi.nlm.nih.gov/pubmed/11871597).”
*Journal of Biomechanical Engineering* 124 (1): 113–20.

This repository: [![DOI](https://zenodo.org/badge/594234010.svg)](https://zenodo.org/badge/latestdoi/594234010)

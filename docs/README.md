# Documentation for DynLoco

This documentation was built with `Documenter.jl`. The main editable files are contained in `src/`, from which the `docs/*html` files are generated. The remaining files in `docs/` are the result of this process, and not generally considered editable.

## How to generate documentation
From the `docs/` directory, and `using Documenter` within Julia, enter `makedocs(sitename="DynLoco Documentation")`. This generates html pages starting with `build/index.html`. These are intended for viewing as GitHub pages, so the links may not work. Also, do not commit anything from the `build` directory to the main branch; the commits can go to a separate gh-pages branch.
# EpiCats
Applied category theory applied to epidemiological models 

## Introduction

This is a repository for applying category theory to epidemiological models. The goal is to provide a framework for epidemiologists to build models that are more modular, composable, and reusable. This makes use of the [AlgebraicJulia ecosystem of packages](https://algebraicjulia.org), and complements the examples in the AlgebraicJulia documentation with more detail and explanations targeted predominantly at modellers who are unfamiliar with applied category theory (ACT).

## Examples

### Petri nets
- [SIR model composed from infection and recovery submodels](https://github.com/epirecipes/EpiCats/blob/main/pn_compose_sir/pn_compose_sir.ipynb)
- [The linear chain trick](https://github.com/epirecipes/EpiCats/blob/main/pn_compose_sir_stages/pn_compose_sir_stages.ipynb)
- [SIR model stratified by risk](https://github.com/epirecipes/EpiCats/blob/main/pn_stratify_multigroup/pn_stratify_multigroup.ipynb)

## Glossary of applied category theory applied to epidemiological models

ACT has a lot of terminology that may be unfamiliar; in the table below, terminology in ACT is mapped to concepts in epidemiological models.

| *Term*                                    | *Definition*                                                                                                                                    | *Example in epidemiology*                                                                                                                                                |
| --------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Category                                | Includes “objects” representing entities or states and “arrows” (morphisms) representing transitions or relationships between entities.       | Objects: different disease states of an individual (Susceptible, Infected, Recovered).                                                                                 |
| Morphism                                | A mapping, function or transformation from one object to another.                                                                             | The transformation of a person from a susceptible state to an infected state.                                                                                          |
| Domain                                  | The object from which the arrow or morphism originates.                                                                                       | For an arrow representing transmission, the domain would be the susceptible state.                                                                                     |
| Codomain                                | The object to which the arrow points.                                                                                                         | For an arrow representing transmission, the codomain would be the infected state.                                                                                      |
| Functor                                 | A mapping between categories that preserves the structure, by mapping objects to objects and arrows to arrows.                                | The effort of a vaccination program, where the functor maps the pre-vaccination states and transitions to the post-vaccination states and transitions.                 |
| Monad                                   | A structure that encapsulates side effects.                                                                                                   | Stochastic effects such as the randomness of disease transmission.                                                                                                     |
| Product                                 | The product of two objects in a category is a new object that captures the "essence" of both objects simultaneously.                          | The combination of age and vaccination classes to create ‘young/unvaccinated’ classes etc.                                                                             |
| Coproduct                               | Represents the "sum" or "choice" of two objects                                                                                               | A grouping of all vaccinated individuals across age classes.                                                                                                           |
| Pullback  (generalization of product)  | Represents the most general object that maps to two given objects in a way that is compatible with a specified mapping between those objects. | Can represent interactions between different populations or compartments. Unlike a product, specific interactions can be selected, rather than all.                    |
| Pushout (generalization of coproduct) | Represents the most general object that two given objects map to, in a way that is compatible with a specified mapping between those objects. | Can be used to represent the fusion or aggregation of different disease states or populations. Unlike a coproduct, specific states can be chosen to be glued together. |
| Span                                    | A diagram consisting of two morphisms with a common domain                                                                                    | Leaving a susceptible state through infection or vaccination.                                                                                                          |
| Cospan                                  | A diagram consisting of two morphisms with a common codomain                                                                                  | Entering an immune state through either recovery or vaccination.                                                                                                       |
| Hom-set                                 | The set of all morphisms from one object to another in a category.                                                                            | All possible transitions between a pair of states; in an SIR model, there would be a hom-set between S and I composed of the morphism representing transmission.       |

## Contributing

We welcome contributions to this repository. Each example should be placed in its own directory, as a Quarto (`.qmd`) file. To assist in conversion to Pluto notebooks, please avoid multiple definitions within a single cell (wrapping a cell with `begin` and `end` if this is necessary).

The Quarto file should be rendered to (at least) a Jupyter notebook and placed in the same directory, to allow the input and output to be displayed in GitHub. The Jupyter notebook should be named to match the directory name. For example, the directory `pn_compose_sir` contains the Quarto file `pn_compose_sir.qmd` and the Jupyter notebook `pn_compose_sir.ipynb`. The `README.md` in the root directory should contain a link to the rendered Jupyter notebook.

To render Quarto files to Jupyter notebooks, install Quarto and run the following command in the example directory (where `{FILENAME_OF_QMD}` is the name of the Quarto file):

```bash
quarto render {FILENAME_OF_QMD} --to ipynb
```

In addition, the Jupyter notebooks can be rendered to Pluto notebooks using [this tool](https://observablehq.com/@olivier_plas/pluto-jl-jupyter-conversion) to create a plain Julia script. This script may need to be edited (e.g. to avoid multiple definitions within a single cell) in order to run properly.

## License

This repository is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgements

This repository borrows heavily from the AlgebraicJulia documentation, and exchanges with the AlgebraicJulia community. We are grateful to the AlgebraicJulia community for their support and encouragement.
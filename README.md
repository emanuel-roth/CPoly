<p align="center">
  <img src="cpolylogo.svg#gh-light-mode-only" width="300" alt="logo" />
  <img src="cpolylogodark.svg#gh-dark-mode-only" width="300" alt="logo" />
</p>
 

*QuiverCombinatoricsTools* is a *SageMath* package that adds combinatorial functions to *QuiverTools* to calculate symplectic leaves of quiver varieties. It adds on to the *QuiverTools* package written by Pieter Belmans, Hans Franzen and Gianni Petrella, as seen [here](https://sage.quiver.tools/) and [here](https://github.com/QuiverTools/QuiverTools). The website is hosted [here](https://github.com/QuiverCombinatoricsTools/quivercombinatoricstools.github.io).

# Authors

* Tudor-Ioan Caba (University of Edinburgh)
* Mia Lam (University of Edinburgh)
* [Emanuel Roth](https://emanuel-roth.github.io/) (University of Edinburgh)

We were supervised by [Gwyn Bellamy](https://sites.google.com/view/gwynbellamy/home) (University of Glasgow), as part of an [AGQ](https://www.agq-cdt.org/) computing project. The project is licensed under the *MIT license*. We are open to improvements and suggestions.

# Instructions

To install it, make sure you have both QuiverTools and QuiverCombinatoricsTools
```
sage --pip install git+https://github.com/QuiverTools/QuiverTools.git
sage --pip install git+https://github.com/QuiverCombinatoricsTools/QuiverCombinatoricsTools.git
```
and then you can simply run
```
from quiver import *
from quivercombinatorics import *
```
to get started.

Alternatively, you can run it from your browser in a notebook here using [binder](https://mybinder.org/v2/gh/QuiverCombinatoricsTools/binder/HEAD).

# Documentation

You can read the documentation of *QuiverCombinatoricsTools* as:

* [a webpage](https://quivercombinatoricstools.github.io/)
* [a pdf](https://raw.githubusercontent.com/QuiverCombinatoricsTools/quivercombinatoricstools.github.io/main/docs/_static/quivercombinatoricstools.pdf)

This builds on *QuiverTools*, and you can read the documentation of *QuiverTools* as:

* [a webpage](https://sage.quiver.tools)
* [a pdf](https://sage.quiver.tools/documentation.pdf)

# How to cite QuiverCombinatoricsTools

If you have used this code in any way, please consider citing it in the following way.
```
@software{quivercombinatoricstools,
  author = {Caba, Tudor-Ioan and Lam, Mia and Roth, Emanuel},
  title = {QuiverCombinatoricsTools},
  url = {https://quivercombinatoricstools.github.io},
}
```
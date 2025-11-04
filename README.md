# storm-lists

[![PyPI badge](https://badge.fury.io/py/storm-lists.svg)](https://badge.fury.io/py/storm-lists)
[![Zenodo badge](https://zenodo.org/badge/1088748199.svg)](https://doi.org/10.5281/zenodo.17521945)

This package compiles Python implementations of various geomagnetic storm lists.

## Citations

When using this software, please cite [the Zenodo record](https://doi.org/10.5281/zenodo.17521945) as well as the paper corresponding to the list you are using:

- [Noora Partamies et al. (2013), Statistical properties of substorms during different storm and solar cycle phases, _Annales Geophysicae_, 31(2), 349–358. https://doi.org/10.5194/angeo-31-349-2013](https://doi.org/10.5194/angeo-31-349-2013)
- [Maria-Theresia Walach and Adrian Grocott (2019). Superdarn observations during geomagnetic storms, geomagnetically active times, and enhanced solar wind driving, _Journal of Geophysical Research: Space Physics_, 124(7), 5828–5847. https://doi.org/10.1029/2019JA026816](https://doi.org/10.1029/2019JA026816)

## Funding

John C Coxon was supported during this work by Science and Technology Facilities Council (STFC) Ernest Rutherford Fellowship ST/V004883/1.

## Installation

To install and test the model, simply type:

```
pip install storm-lists
```

## Usage

For examples of how to use the models, see the Jupyter notebooks in the `notebooks` directory.

## Model-specific notes

### Walach and Grocott (2019)

The Walach and Grocott (2019) implementation in this module is directly based on the IDL code used to generate the list for that paper.

### Partamies et al. (2013)

The Partamies et al. (2013) implementation in this module is based on the description in the paper, and it makes assumptions outlined below.

1. The rate of change at a timestamp is calculated between Dst at the timestamp and Dst at the previous timestamp.
2. The recovery phase start is set to the last timestamp at which the rate of change is still under the threshold.
3. If there is only one hour where the rate of change is below the threshold, the main phase and recovery phase will be set to the same timestamp.
4. When a local minimum in Dst is not associated with a rate of change below the threshold, it is ignored and the algorithm runs again neglecting that minimum Dst. (This implies that Dst can theoretically still continue to decrease during the recovery phase before starting to recover.)  
5. When looking back to the point at which Dst exceeded the start threshold, if that point has already been labelled as belonging to the storm, the start of the storm under consideration is set to one hour after the previous storm end and the storm is labelled as a type 2 storm.
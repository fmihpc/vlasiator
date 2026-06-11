
![](https://github.com/fmihpc/vlasiator/blob/master/doc/artwork/logo_black.png?raw=true)

Vlasiator - ten letters you can count on

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10600112.svg)](https://doi.org/10.5281/zenodo.10600112)

## Code description
Space weather is a term used to describe the variable environmental effects within near-Earth space, caused by the Sun emitting solar wind, a stream of charged particles carrying the solar electromagnetic field. Space weather can be caused by solar high-energy particles or by dynamic variations of the solar wind that can cause extended periods of major disturbances on ground and space, affecting technological systems (e.g., telecommunication and weather spacecraft at geostationary orbit, and ground-based power grids).

In Vlasiator, ions are represented as velocity distribution functions, while electrons are magnetohydrodynamic fluid, enabling a self-consistent global plasma simulation that can describe multi-temperature plasmas to resolve non-MHD processes that currently cannot be self-consistently described by the existing global space weather simulations. The novelty is that by modelling ions as velocity distribution functions the outcome will be numerically noiseless. 

Due to the multi-dimensional approach at ion scales, Vlasiator's computational challenges are immense. We use advanced high performance computing techniques to allow massively parallel computations on tens of thousands of cores.

## git submodules
We are transferring to use `git submodules` for the dependent libraries. Some of the header libraries have already been moved to this framework. Thus, we recommend to use the `--recurse-submodules` option when pulling or checking out branches.

For first-time cloning, the following is required in order to initialize submodules correctly:
```
git clone --recurse-submodules https://github.com/fmihpc/vlasiator
git checkout <branch>
git submodule update --init --recursive
```

## Documentation
See the [wiki](https://github.com/fmihpc/vlasiator/wiki) for build instructions and general advice.

## Paper references
- [Palmroth, M., Ganse, U., Pfau-Kempf, Y., Battarbee, M., Turc, L., Brito, T., Grandin, M., Hoilijoki, S., Sandroos, A. & von Alfthan, S. (2018). Vlasov methods in space physics and astrophysics. Living Reviews in Computational Astrophysics, 4, 1. doi: 10.1007/s41115-018-0003-2](https://link.springer.com/article/10.1007/s41115-018-0003-2)
- [von Alfthan, S., Pokhotelov, D., Kempf, Y., Hoilijoki, S., Honkonen, I., Sandroos, A. & Palmroth, M. (2014). Vlasiator: First global hybrid-Vlasov simulations of Earth's foreshock and magnetosheath . Journal of Atmospheric and Solar-Terrestrial Physics , 120, 24 - 35. doi: http://dx.doi.org/10.1016/j.jastp.2014.08.012](http://www.sciencedirect.com/science/article/pii/S1364682614001916)

[<img src="https://github.com/zulip/zulip/blob/main/static/images/logo/zulip-icon-128x128.png" width="32"/>](https://zulip.com) Sponsored by Zulip, an open-source modern team chat app designed to keep both live and asynchronous conversations organized.

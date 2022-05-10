**MATsint**
======

## Description

**MATsint** is a Matlab package for the interpolating sparse geophysical data. This package is a temporary developing package in the initial stage.

## Reference
Chen, Y., Chen, X., Wang, Y. and Zu, S., 2019. The interpolation of sparse geophysical data. Surveys in Geophysics, 40(1), pp.73-105.

BibTeX:

	@article{sint,
	  title={The interpolation of sparse geophysical data},
	  author={Yangkang Chen and Xiaohong Chen and Yufeng Wang and Shaohuan Zu},
	  journal={Surveys in Geophysics},
	  volume={40},
	  number={1},
	  pages={73-105},
	  year={2019}
	}

-----------
## Copyright
    Initial version: Yangkang Chen (chenyk2016@gmail.com), 2021-2022
	Later version: MATsint developing team, 2022-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    svn co https://github.com/chenyk1990/reproducible_research/trunk/nonMada/sint ./
	addpath(genpath('./sint'));

-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of MATsint. 

-----------
## Dependence Packages
* seistr (https://github.com/chenyk1990/seistr)

-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------
## Modules
    yc_sint2d.m -> 2D interpolation for sparse data

-----------
## Counting lines
Counting lines of the Package:

    port install cloc

Using the following comman in src/MATsint/ or in the main directory to get a report

    cloc ./

-----------


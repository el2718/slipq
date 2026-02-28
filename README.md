# Slip-Squashing Factors

Slip-squashing factors $Q_\mathrm{sf}$ and $Q_\mathrm{sb}$ ([Titov_2009_ApJ_693_1029](https://iopscience.iop.org/article/10.1088/0004-637X/693/1/1029)) are defined by two field line mappings and two boundary flow mappings between two instants; their large values define the surfaces that border of the reconnected or to-be-reconnected magnetic flux tubes for a given period of time during the magnetic evolution. 

For the case of static boundaries ($\vec{u}(t) = \vec{0}$ at the bottom boundary), we can compute the slip-squashing factors using the coordinate mapping provided by FastQSL. Following the initial coordinate mapping within the first magnetic field, the resulting mapped coordinates can be served as a seed grid for applying FastQSL to the second magnetic field. 

This program can be downloaded with the command
```
git clone https://github.com/el2718/slipq
```

And FastQSL2 should be downloaded first
```
git clone https://github.com/el2718/FastQSL2
```

Please address comments and suggestions to [Dr. Chen, Jun (陈俊)](mailto:chenjun@pmo.ac.cn)

If your markdown reader can not can not render the formulae in README\.md, please read README.html directly.

-----------------------------
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This program is licensed under a [CC BY-NC-SA 4.0 License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

-----------------------------
## Cite as

* Jun Chen*, Thomas Wiegelmann, Li Feng*, Bernhard Kliem, Chaowei Jiang, and Rui Liu. FastQSL 2: A Comprehensive Toolkit for Magnetic Connectivity Analysis.  2026, SCIENCE CHINA Physics, Mechanics & Astronomy, submitted

-----------------------------
## Parameters
* **bx0**, **by0**, **bz0**: the three components of magnetic field at $t_0$

* **bx1**, **by1**, **bz1**: the three components of magnetic field at $t_1$

* **xa**, **ya**, **za**, **spherical**, **xreg**, **yreg**, **factor**, **lon_delta**, **lat_delta**, **preview**, **fname** are the same as those in FastQSL2

-----------------------------
## Products
Both par2solarwind\.pro and par2solarwind\.py are functions, their return is a slip-squashing factor. If $t_0 < t_1$, $Q_\mathrm{sf}$ is returned; If $t_1 < t_0$, $Q_\mathrm{sb}$ is returned.

**qsl0** is the product of FastQSL 2 at for the field **bx0**, **by0**, **bz0** at $t_0$. If use slipq\.pro, it is return by the keyword **qsl0**. If use slipq\.py, $t_0 < t_1$, the results are returned by the tuple `(qsf, qsl0)`, see the example of demo_Btitov2009\.py

-----------------------------
## Demos

### If use fastqsl\.pro
```idl
IDL> .r demo_Btitov2009.pro
```
### If use fastqsl\.py
```python
python3 demo_Btitov2009.py
```

-----------------------------
## History
* Feb 26, 2026 Jun Chen, version 1.0

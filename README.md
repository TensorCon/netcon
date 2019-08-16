# netcon
A mirror of source provided in [***Faster identification of optimal contraction sequences for tensor networks***](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.033315)

Work is provided in accordance with [Physical Review Transfer of Copyright Agreement](https://journals.aps.org/authors/transfer-of-copyright-agreement). More can be found [here](https://journals.aps.org/pre/authors/supplemental-materials-journals).

Thank you to Robert Pfeifer and the [American Physical Society](https://www.aps.org/), who have granted re-use permissions.

You can cite this work with:
```
@article{pfeifer2014faster,
  title={Faster identification of optimal contraction sequences for tensor networks},
  author={Pfeifer, Robert NC and Haegeman, Jutho and Verstraete, Frank},
  journal={Physical Review E},
  volume={90},
  number={3},
  pages={033315},
  year={2014},
  publisher={APS}
}
```

You may find the preprint of this work on [arXiv](https://arxiv.org/abs/1304.6112).


### Examples
Can be found near the top of `netcon.m`

### Dependencies

* `gcc`
* `octave` or `Matlab`

### Setting up with `octave`

Be sure to create a `.mex` file from `netcon_nondisj_cpp.cpp` to ensure the C++ source is utilized:

```
mkoctfile --mex -O3 netcon_nondisj_cpp.cpp
```

### Running
To run the algorithm outside of the octave interactive mode:
```
$ octave --no-gui-libs --no-gui --path <path/to/netcon.m> --eval "netcon({[-1 1 2 3],[2 4 5 6],[1 5 7 -3],[3 8 4 9],[6 9 7 10],[-2 8 11 12],[10 11 12 -4]},1,1,1,1);"
 
Looking for solutions with maximum cost of 2
Looking for solutions with maximum cost of 64
Looking for solutions with maximum cost of 128
Looking for solutions with maximum cost of 256
Looking for solutions with maximum cost of 512
Looking for solutions with maximum cost of 1024
 
Best sequence:  11 12 2 1 5 3 4 6 7 9 8 10
Cost:           896
```

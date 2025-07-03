# Femur Stress Simulation (FEniCSx + ParaView)

This project simulates the stress distribution in a human femur using a scalar finite element model implemented in FEniCSx. The output includes von Mises stress approximations and deformation visualized in ParaView.

> ⚠️ Due to Docker constraints, scalar PDEs were used instead of full vector elasticity.

## Installation

Make sure you have [FEniCSx](https://docs.fenicsproject.org/dolfinx/main/python/) installed (via Docker or conda). Place your mesh at:

/home/fenics/shared/femur.msh


## Usage

```bash
python vonmises_fea.py
```
## Features

Scalar FEA using FEniCSx

Custom region-based force loading (radial + vertical control)

NumPy-based von Mises stress computation

ParaView visualization (warp and stress fields)

Clean workaround for vector displacement using scalar stacking

## Result

Open results in ParaView and apply:

WarpByVector using vec

Color by col (von Mises stress)

## Limitations

Scalar approximation does not capture realistic deformation physics

Downward force results in upward visual displacement due to Poisson PDE behavior

No traction boundary conditions (only internal force distribution)

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)

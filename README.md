# Femur Stress Simulation (FEniCSx + ParaView)

This project simulates the stress distribution in a human femur using a scalar finite element model implemented in FEniCSx. The output includes von Mises stress approximations and deformation visualized in ParaView.

> ⚠️ Due to Docker constraints, scalar PDEs were used instead of full vector elasticity.

## Installation

Make sure you have [FEniCSx](https://docs.fenicsproject.org/dolfinx/main/python/) installed (via Docker or conda). Place your mesh at:

/home/fenics/shared/femur.msh


## Usage

```bash
python vonmises_fea.py

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)

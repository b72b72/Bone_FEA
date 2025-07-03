print("Running FEniCSx scalar FEA model...", flush=True)

from petsc4py import PETSc
PETSc.Sys.popErrorHandler()

from dolfinx import mesh, fem
from dolfinx.io import gmshio
import ufl
import numpy as np
from mpi4py import MPI


# Load mesh
domain, _, _ = gmshio.read_from_msh("/home/fenics/shared/femur.msh", MPI.COMM_WORLD)
dx = ufl.dx(domain)


# Create scalar function space (Lagrange P1)
V = fem.functionspace(domain, ("Lagrange", 1))


# Function to solve scalar PDE
def solve_scalar_component():
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Tabulate coordinates of DOFs
    x_dof = V.tabulate_dof_coordinates()
    x_dof = x_dof.reshape((-1, 3))
    x_vals, y_vals, z_vals = x_dof[:, 0], x_dof[:, 1], x_dof[:, 2]

    # Z thresholds
    z_max = np.max(z_vals)
    z_min = np.min(z_vals)
    height = z_max - z_min
    top_thresh = z_min + 0.90 * height
    top_cutoff = z_max - 0.01

    # Center of top region
    top_mask = (z_vals > top_thresh)
    center_x = np.mean(x_vals[top_mask])
    center_y = np.mean(y_vals[top_mask])

    # Radial distance from center of femoral head
    r = np.sqrt((x_vals - center_x)**2 + (y_vals - center_y)**2)

    # Final region: top Z + close to center (radius < r_max)
    r_max = 2.0  # adjust this as needed based on femoral head size
    region_mask = (z_vals > top_thresh) & (z_vals < top_cutoff) & (r < r_max)

    # Cosine force taper inside region
    force_vals = np.zeros_like(z_vals)
    force_vals[region_mask] = -5000 * 0.5 * (
        1 + np.cos(np.pi * (z_vals[region_mask] - top_thresh) / (top_cutoff - top_thresh))
    )

    f = fem.Function(V)
    f.x.array[:] = force_vals

    # Variational problem
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
    L = f * v * dx

    # Fixed bottom
    facets_bottom = mesh.locate_entities_boundary(
        domain, domain.topology.dim - 1,
        lambda x: x[2] < z_min + 1e-3
    )
    dofs_bottom = fem.locate_dofs_topological(V, domain.topology.dim - 1, facets_bottom)
    bc = fem.dirichletbc(PETSc.ScalarType(0.0), dofs_bottom, V)

    uh = fem.Function(V)
    from dolfinx.fem.petsc import LinearProblem
    problem = LinearProblem(a, L, bcs=[bc], u=uh)
    problem.solve()

    return uh.x.array, f



# Run solves for x, y, z directions
_, f_dummy = solve_scalar_component()  # no x
_, f_dummy = solve_scalar_component()  # no y
u_z, f = solve_scalar_component()      # actual z-displacement + keep f
u_x = np.zeros_like(u_z)
u_y = np.zeros_like(u_z)



# Combine into vector displacement (node-wise)
u_combined = np.stack([u_x, u_y, -u_z], axis=-1)

if MPI.COMM_WORLD.rank == 0:
    mags = np.linalg.norm(u_combined, axis=1)
    print("Max displacement:", np.max(mags))
    print("Min displacement:", np.min(mags))

# Material properties
E = 1e3
nu = 0.3
mu = E / (2 * (1 + nu))
lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))

coordinates = domain.geometry.x

domain.topology.create_connectivity(domain.topology.dim, 0)
cells = domain.topology.connectivity(domain.topology.dim, 0).array.reshape(-1, 4)  # each cell = 4 nodes
von_mises_stress = np.zeros(cells.shape[0])  # ✅ one value per cell


for i in range(len(cells)):
    node_ids = cells[i]                     # 4 node indices (tetrahedron)
    X = coordinates[node_ids]               # shape (4, 3)
    U = u_combined[node_ids]                # shape (4, 3)

    X0 = X[0]
    DX = X[1:] - X0                         # (3, 3)
    DU = U[1:] - U[0]                       # (3, 3)

    try:
        J = np.linalg.inv(DX.T) @ DU.T      # grad_u ≈ inv(DX^T) @ DU^T
        grad_u = J.T

        strain = 0.5 * (grad_u + grad_u.T)
        stress = lmbda * np.trace(strain) * np.eye(3) + 2 * mu * strain
        dev = stress - (1/3) * np.trace(stress) * np.eye(3)
        vm = np.sqrt(1.5 * np.sum(dev**2))
        von_mises_stress[i] = vm
    except np.linalg.LinAlgError:
        von_mises_stress[i] = 0.0

if MPI.COMM_WORLD.rank == 0:
    print("✅ von Mises stress (NumPy) computed.")
    print("Max σ_vm:", np.max(von_mises_stress))
    print("Min σ_vm:", np.min(von_mises_stress))



from dolfinx.io import XDMFFile
from dolfinx.fem import Function, FunctionSpace

from dolfinx.fem import FunctionSpace, Function, locate_dofs_topological

# Define DG0 function space (1 value per cell)
S = fem.functionspace(domain, ("DG", 0))
stress_fem = Function(S)

# ⚠️ Ensure connectivity is built
domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim)

stress_fem.x.array[:] = von_mises_stress



with XDMFFile(MPI.COMM_WORLD, "von_mises.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(stress_fem)

    
f.name = "applied_force"  # Optional: name that appears in ParaView

with XDMFFile(MPI.COMM_WORLD, "applied_force.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(f)

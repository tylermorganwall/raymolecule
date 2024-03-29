% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_bond_scene.R
\name{generate_bond_scene}
\alias{generate_bond_scene}
\title{Build Scene (bonds only)}
\usage{
generate_bond_scene(
  model,
  x = 0,
  y = 0,
  z = 0,
  scale = 1,
  center = TRUE,
  force_single_bonds = FALSE,
  pathtrace = TRUE,
  material = rayrender::glossy,
  material_vertex = material_list(diffuse = "grey33", ambient = "grey33", type = "phong",
    ambient_intensity = 0.3)
)
}
\arguments{
\item{model}{Model extracted from a PDB or SDF file.}

\item{x}{Default `0`. X offset, applied after centering.}

\item{y}{Default `0`. Y offset, applied after centering.}

\item{z}{Default `0`. Z offset, applied after centering.}

\item{scale}{Default `1`. Amount to scale the interatom spacing.}

\item{center}{Default `TRUE`. Centers the bounding box of the model.}

\item{force_single_bonds}{Default `FALSE`. Whether to force all bonds to show as a single connection.}

\item{pathtrace}{Default `TRUE`. If `FALSE`, the `rayvertex` package will be used to render the scene.}

\item{material}{Default `rayrender::glossy`. Rayrender material to use when `pathtrace = TRUE`. Must be either `glossy`, `diffuse`, or `dielectric`.}

\item{material_vertex}{Default `material_list(diffuse="grey33",ambient="grey33",type="phong", ambient_intensity=0.3)`.
Material to use for the bonds when `pathtrace = FALSE`.}
}
\value{
Rayrender/rayvertex scene containing only the connections between atoms in a molecule/protein.
}
\description{
Reads an SDF file and extracts the 3D molecule model
}
\examples{
#Generate a scene with benzene molecule with just the atoms
\donttest{
get_example_molecule("benzene") \%>\%
  read_sdf() \%>\%
  generate_bond_scene() \%>\%
  render_model(lights = "both", samples=256,sample_method="sobol_blue")

#Force single bonds to just show the shape of the molecule
get_example_molecule("benzene") \%>\%
  read_sdf() \%>\%
  generate_bond_scene(force_single_bonds = TRUE) \%>\%
  render_model(lights = "both", samples=256,sample_method="sobol_blue")

#Generate a scene with PFOA, reducing the inter-atom spacing
get_example_molecule("pfoa") \%>\%
  read_sdf() \%>\%
  generate_bond_scene(scale=0.3,force_single_bonds = TRUE) \%>\%
  render_model(lights = "both", samples=256,sample_method="sobol_blue")
}
}

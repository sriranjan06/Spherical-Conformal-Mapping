# Spherical Conformal Mapping

## Overview

This project implements **Spherical Conformal Mapping** on Stanford's brain model. The goal is to map the brain's surface onto a spherical domain while preserving the conformal structure, making it useful for applications in neuroanatomy, medical imaging, and other fields requiring complex surface analysis.

## Project Structure

The project is structured into several key components, each playing a vital role in the process of spherical conformal mapping:

### Core Components

- **Geometric Entities:**
  - `Edge.cpp/h`, `Face.cpp/h`, `HalfEdge.cpp/h`, `Vertex.cpp/h`, `Point.cpp/h`: These files define the fundamental geometric components of the brain model's mesh, including vertices, edges, faces, and half-edges.
  - `Solid.cpp/h`, `SolidDelegate.cpp/h`: These files manage the overall structure of the mesh, handling higher-level operations and interactions between the geometric entities.

- **Traits and Traits Management:**
  - `CTrait.h`, `EditTrait.h`, `Trait.cpp/h`: These files define and manage traits associated with geometric entities, such as attributes or metadata that may be necessary for the conformal mapping process.

- **File Parsing and Exception Handling:**
  - `OBJFileReader.cpp/h`, `Parser.cpp/h`, `StringTokenizer.cpp/h`, `string_token_iterator.h`: These files handle the parsing of input files (e.g., OBJ files), which contain the 3D model data. The `StringTokenizer` and related iterator assist in breaking down and interpreting file content.
  - `FException.cpp/h`, `TopologyException.cpp/h`: These files manage exceptions, handling errors related to file parsing, geometric processing, or topological inconsistencies.

- **Data Structures and Algorithms:**
  - `avltree.h`, `DList.h`, `SList.h`, `iterators.h`: These files implement core data structures (e.g., AVL trees, linked lists) and iterators for traversing and manipulating the mesh data.

### Data and Output Files

- **Data Input:**
  - `brain.m`: A MATLAB script potentially used for pre-processing or analyzing the brain model data.
  - `brain.obj`: The primary 3D brain model in OBJ format, used as the input for spherical conformal mapping.
  - `brain.off`: An alternative mesh representation of the brain model in OFF format.

- **Output:**
  - `sphericalConfOp_harmonic.obj`: The output file containing the brain model after spherical conformal mapping. The suffix `_harmonic` indicates the use of harmonic functions or optimizations in the mapping process.

### Build and Execution Instructions

The project is built using CMake, a cross-platform build system, that ensures easy compilation across different environments.

#### Build Instructions

1. Navigate to the root directory of the project.
2. Run the following commands:

   ```bash
   $ rm -r build
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make
   ```

   **Note:** The program has been tested on Ubuntu 16.04 with g++ 5.4.0.

#### Usage

After building the project, you can generate the spherical conformal mapping output by running:

```bash
$ ./conformalMap ../Data/brain.obj ../Output/sphericalConfOp
```

This command will process the input brain model (`brain.obj`) and produce the mapped output (`sphericalConfOp_harmonic.obj`).

## Conclusion

This project successfully implements spherical conformal mapping on a 3D brain model, providing a valuable tool for researchers and professionals in fields such as neuroanatomy and medical imaging. The comprehensive handling of mesh data, combined with robust file parsing and error management, ensures a reliable and efficient mapping process.

## Acknowledgments

Special thanks to Stanford for providing the brain model and to the developers of the tools and libraries that made this project possible.

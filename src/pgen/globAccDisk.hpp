#ifndef GLOB_ACC_DISK
#define GLOB_ACC_DISK

// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"

#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"

extern void (*userFieldBvals) (MeshBlock*);

#endif // GLOB_ACC_DISK

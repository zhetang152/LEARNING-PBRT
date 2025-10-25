#pragma once

// Float Type Definitions
#ifdef PBRT_FLOAT_AS_DOUBLE
using Float = double;
#else
using Float = float;
#endif

template<typename T>
class Vector3;
template<typename T>
class Point3;
template<typename T>
class Normal3;
#pragma once

#include "vecmath.h"

class ray{
public:
    Vector3f direction;
    Point3f origin;
    Float time = 0;
    
    ray() = default;
    ray(Point3f o, Vector3f d, Float t = 0):origin(o), direction(d), time(t) {}

    bool HasNaN() const { return (origin.HasNaN() || direction.HasNaN()); }

    Point3f operator()(Float t) const{
        return origin + direction * t;
    }
};

class rayDifferential : public ray {
public:
    bool hasDifferentials = false;
    Point3f rxOrigin, ryOrigin;
    Vector3f rxDirection, ryDirection;

    rayDifferential() = default;
    rayDifferential(Point3f o, Vector3f d, Float t =0): ray(o, d, t) {}
    explicit rayDifferential(const ray& r): ray(r) {}

    bool HasNaN() const {
        return ray::HasNaN() ||
               (hasDifferentials && (rxOrigin.HasNaN() || ryOrigin.HasNaN() ||
                                     rxDirection.HasNaN() || ryDirection.HasNaN()));
    }

    void ScaleDifferentials(Float s) {
        rxOrigin = origin + (rxOrigin - origin) * s;
        ryOrigin = origin + (ryOrigin - origin) * s;
        rxDirection = direction + (rxDirection - direction) * s;
        ryDirection = direction + (ryDirection - direction) * s;
    }
};
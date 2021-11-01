#ifndef _HITABLE_H_
#define _HITABLE_H_
#include "ray.h"

class material;

struct hit_record{
    float t_ = 0; // 方程解t
    vec3  p_ = vec3(0,0,0); // point
    vec3  normal_= vec3(0,0,0); //法线
    material *mate_ptr;//Ch8 add:object's material attribution
};

// 抽象类：适应多物体，多光线情况
// 给定t的区间 tmin , tmax
class hitable{
    public:
    virtual bool hit(const ray&r, float t_min, float t_max,hit_record& rec) const = 0;
};


#endif
#ifndef _HITABLE_H_
#define _HITABLE_H_
#include "ray.h"

class material;

struct hit_record{
    float t_ = 0; // 方程解t
    vec3  p_ = vec3(0,0,0); // hit point
    vec3  normal_= vec3(0,0,0); //法线
    material *mate_ptr;//Ch8 add:object's material attribution
    bool frontFace = true;//靠近光源一侧为front
    // 这里的normal都是从球心指向球面,一束光贯穿球体
    inline void setFaceNormal(const ray& r,const vec3& outward_normal){
        frontFace = dot(r.direction(), outward_normal) < 0;
        normal_ = frontFace ? outward_normal : (-outward_normal);
    }
};

// 抽象类:适应多物体,多光线情况
// 给定t的区间 tmin,tmax
class hitable{
    public:
    virtual bool hit(const ray&r, float t_min, float t_max,hit_record& rec) const = 0;
};


#endif

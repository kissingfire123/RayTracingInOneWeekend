#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "vec3.h"
#include "ray.h"
#include "hitable.h"

class material{
    public:
    virtual bool scatter(const ray& r_in,const hit_record& reco, vec3& attenuation, ray & scattered) const = 0;
};

class lambertian : public material{
public:
    lambertian(const vec3 & a): albedo_(a){}
    virtual bool  scatter(const ray& r_in,const hit_record& reco,vec3& attenuation, ray & scattered) const {
        vec3 target = reco.p_  + reco.normal_ + random_in_unit_sphere();
        scattered = ray(reco.p_ , target - reco.p_);
        attenuation = albedo_;
        return true;
    }
private:
    vec3 albedo_;//反射率,入射光量/散射光量
};


class metal : public material{
public:
    metal(const vec3& a): albedo_(a){}
    virtual bool  scatter(const ray& r_in,const hit_record& reco,vec3& attenuation, ray & scattered) const {
        vec3 reflected = metal::_reflect(unit_vector(r_in.direction()), reco.normal_);
        scattered = ray(reco.p_ , reflected);
        attenuation = albedo_;
        return dot(scattered.direction(), reco.normal_) > 0;
    }
    static vec3 _reflect(const vec3 &v ,const vec3 &n){
        return v - 2 * dot(v, n)* n;
    }
private:

    vec3 albedo_;
};

#endif
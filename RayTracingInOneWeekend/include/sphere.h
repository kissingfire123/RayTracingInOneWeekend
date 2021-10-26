#ifndef _SPHERE_H_
#define _SPHERE_H_

#include "hitable.h"

class sphere:public hitable{
public:
    sphere(){}
    sphere(vec3 cen,float r):center_(cen),radius_(r){}
    virtual bool hit(const ray& r, float tmin,float tmax,hit_record& rec)const override; 
private:
    vec3 center_;
    float radius_;
};

bool sphere::hit(const ray& r, float tmin,float tmax,hit_record& reco) const{
    // 计算delta
    vec3 oc = r.origin() - center_;
    float a = dot(r.direction() , r.direction());
    float b = 2.0 * dot(oc,r.direction());
    float c = dot(oc,oc) - radius_ * radius_;
    float discriminant = b * b - 4 * a * c;
    auto has_resolve_t = [&](float t)->bool{
        if(t < tmax && t > tmin){
            reco.t_ = t;
            reco.p_ = r.at_Parameter(reco.t_);
            reco.normal_ = (reco.p_ - center_) / radius_;
            return true;
        }
        return false;
    };
    if(discriminant > 0){
        bool hitable = has_resolve_t((-b - sqrt(discriminant))/(2.0 * a));
        if(hitable) return true;
        hitable =  has_resolve_t((-b + sqrt(discriminant))/(2.0 * a));
        if(hitable) return true;
    }
    return false;
}


#endif
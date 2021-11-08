#ifndef _SPHERE_H_
#define _SPHERE_H_

#include "hitable.h"

class lambertian;

class sphere:public hitable{
public:
    sphere(){}
    explicit sphere(vec3 cen, float r, const std::shared_ptr<material>& material = nullptr) :
        center_(cen),
        radius_(r) ,
        material_(material){}
    virtual bool hit(const ray& r, float tmin,float tmax,hit_record& rec)const override; 
private:
    vec3 center_;
    float radius_;
    shared_ptr<material> material_;
};

bool sphere::hit(const ray& r, float tmin,float tmax,hit_record& reco) const{
    // 计算delta
    vec3 oc = r.origin() - center_;
    float a = dot(r.direction() , r.direction());
    float half_b = dot(oc,r.direction());//简化运算,统一除以2
    float c = dot(oc,oc) - radius_ * radius_;
    float discriminant = half_b * half_b - a * c;
    auto has_resolve_t = [&](float t)->bool{
        if(t < tmax && t > tmin){
            reco.t_ = t;
            reco.p_ = r.at_Parameter(reco.t_);
            vec3 outward_normal = (reco.p_ - center_) / radius_;
            reco.setFaceNormal(r,outward_normal);
            reco.mate_ptr = material_.get();
            return true;
        }
        return false;
    };
    if (discriminant < 0) {
        return false;//方程无解, 光线与球没有交点
    }
    auto sqrtd = sqrt(discriminant);
    
    bool hitable = has_resolve_t((-half_b - sqrtd)/ a);
    if(hitable) return true;
    //hitable =  has_resolve_t((-half_b + sqrtd)/ a); //可能造成无限递归
    //if(hitable) return true;
    
    return false;
}


#endif
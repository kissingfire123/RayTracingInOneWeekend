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
    virtual bool  scatter(const ray& r_in,const hit_record& reco,vec3& attenuation, ray & scattered) const override {
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
    virtual bool  scatter(const ray& r_in,const hit_record& reco,vec3& attenuation, ray & scattered) const override {
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


class dielectric : public material {
public:
    dielectric(float ri) : ref_idx_(ri) {}
    virtual bool  scatter(const ray& r_in, const hit_record& reco, vec3& attenuation, ray & scattered) const override{
        attenuation = color(1.0, 1.0, 1.0);//透明物质,例如玻璃，不吸收光线
        float ni_over_nt = reco.frontFace ? (1.0 / ref_idx_) : ref_idx_;
        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = fmin(dot(-unit_direction, reco.normal_), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        bool hasRefract = ni_over_nt * sin_theta <= 1.0;
        
        float dt = dot(unit_direction, reco.normal_);
        float discriminant = 1.0 - ni_over_nt * ni_over_nt *(1 - dt * dt);
        float outProb = 0;
        //if (discriminant > 0) {
            outProb = _schlick(cos_theta, ni_over_nt) ;
        //}
        //else {
           // outProb = 1.0;
        //}
        vec3 outDirection;
        if ((!hasRefract) || outProb > random_double()) {
            outDirection = reflect(unit_direction,reco.normal_);
        }
        else {
            outDirection = refract(unit_direction,reco.normal_,ni_over_nt);
        }

        scattered = ray(reco.p_,outDirection);
        return true;
    }
private:
    float ref_idx_ ;


    //  入射光线v, 法线n, 相对折射率ni_over_nt
    bool _refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) const {
        vec3 uv = unit_vector(v);
        float dt = dot(uv, n);
        float discriminant = 1.0 - ni_over_nt * ni_over_nt *(1 - dt * dt);
        if (discriminant > 0) {
            refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
            return true;
        } 
        else {
            return false;// 此时没有折射光
        }
    };

    // 反射光和折射光的比重关系讨论
    // 菲涅尔方程可以很好地描述,但计算太复杂,使用Schlick近似
    // 可参考文章:https://zhuanlan.zhihu.com/p/31534769
    float _schlick(float cosine, float ref_idx) const{
        float r0 = (1 - ref_idx) / (1 + ref_idx);//r0为材质提供光线垂直反射时的反射光占比
        r0 = r0 * r0;
        return r0 + (1 - r0)*std::pow((1-cosine),5);
    }
    
};
#endif
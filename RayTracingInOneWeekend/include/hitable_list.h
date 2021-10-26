#ifndef _HITABLE_LIST_H_
#define _HITABLE_LIST_H_

#include "hitable.h"

class hitable_list: public hitable{
public:
    hitable_list(){} 
    hitable_list(std::shared_ptr<hitable> *l, int n){
        list_ = l;
        list_size_ = n;
    }
    virtual bool hit(const ray&r, float t_min, float t_max,hit_record& reco) const override;
private:
    std::shared_ptr<hitable> * list_;
    int list_size_;
};

bool hitable_list::hit(const ray&r, float t_min, float t_max,hit_record& reco) const{
    hit_record tmp_reco;
    bool hit_anything = false;
    double closest_so_far = t_max;//类似深度缓冲
    for(int i = 0; i < list_size_; ++i){
        if(list_[i] -> hit(r,t_min,closest_so_far,tmp_reco)){
            hit_anything = true;
            closest_so_far = tmp_reco.t_;
            reco = tmp_reco;
        }
    }
    return hit_anything;
}

#endif
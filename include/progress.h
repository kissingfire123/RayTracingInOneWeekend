#ifndef _RTW_PROGRESS_H_
#define _RTW_PROGRESS_H_
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cassert>

// 因为生成图片比较耗时,实现一个简单的进度记录类
class RtwProgress{
public:
    RtwProgress(std::string imageFilePath, int imageLine):imgTotalLine_(imageLine){
        begin_ = std::chrono::high_resolution_clock::now();
        std::cout << "Start to produce image: <<" << imageFilePath << ">> ,please wait...\n";
    }

    void Refresh(int currLine){
        assert(imgTotalLine_ > 0);
        //std::cout << static_cast<double>(100.0 * currLine/imgTotalLine_) << " %"<< std::endl;
        print_progress(currLine,imgTotalLine_);
    }

    ~RtwProgress(){
        end_  = std::chrono::high_resolution_clock::now();
        auto cost = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - begin_).count();
        std::cout << "\nIt's done! Spend time : " << cost / 1000.0 << " seconds\n";
    }


private :
    std::chrono::high_resolution_clock::time_point begin_;
    std::chrono::high_resolution_clock::time_point end_;
    int imgTotalLine_ = 0;

    //this function 'print_progress' is copied from Github:https://gist.github.com/juliusikkala/946f505656ed3c35f6c2741f29f26080
    void print_progress(int p, int total, int width = 80)
    {
        std::string total_str = std::to_string(total);
        std::string p_str = std::to_string(p);
        int bar_width = width - total_str.size() * 2 - 4;

        std::cout << '\r';
        if (bar_width > 1)
        {
            int fill_width = bar_width * p / total;
            std::cout << "[";
            for (int i = 0; i < bar_width; ++i)
            {
                char c = ' ';
                if (i < fill_width || p == total) c = '=';
                else if (i == fill_width) c = '>';

                std::cout << c;
            }
            std::cout << "] ";
        }
        std::cout << std::setfill(' ') << std::setw(total_str.size())
            << p_str << "/" << total_str << std::flush;
        if (p == total) std::cout << std::endl;
    }
};



#endif 
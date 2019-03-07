#ifndef _MISC
#define _MISC
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
// useful functions that c++ lacks but everyone else provides

namespace misc{
    inline bool file_exists (const std::string& name) {
        std::ifstream f(name);
        return f.good();
    }

    inline std::string basename(const std::string& pathname){
        return {std::find_if(pathname.rbegin(), pathname.rend(),
                             [](char c) { return c == '/'; }).base(),
                pathname.end()};
    }

    inline std::string trim(const std::string& str)
    {
        size_t first = str.find_first_not_of(' ');
        if (std::string::npos == first)
        {
            return str;
        }
        size_t last = str.find_last_not_of(' ');
        return str.substr(first, (last - first + 1));
    }

    

    // inline double exp_fast(double x){
    //     // WARNING fails if |x| > 1024
    //     //https://codingforspeed.com/using-faster-exponential-approximation/
    //     x = 1 + x/1024;
    //     x *= x; x *= x; x *= x; x *= x;
    //     x *= x; x *= x; x *= x; x *= x;
    //     x *= x; x *= x;
    //     return x;
    // }

    // inline double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}



}


#endif

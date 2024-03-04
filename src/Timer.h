//
// Created by Edward on 2022/4/22.
//

#ifndef HGVSGO_TIMER_H
#define HGVSGO_TIMER_H

#include <iostream>
#include <chrono>


class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float> duration;

public:
    Timer() {
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        float ms = duration.count();
        std::cout << "Duration: " << ms << "s" << "\n";
    }
};


#endif //HGVSGO_TIMER_H

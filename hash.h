#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

typedef std::size_t HASH_INDEX_T;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    MyStringHash(bool debug = true)
    {
        if(false == debug){
            generateRValues();
        }
    }
    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {
    
    const int MAX_GROUPS = 5;
    const int GROUP_SIZE = 6;

    unsigned long long w[MAX_GROUPS] = {0, 0, 0, 0, 0};

    int len = k.length();
    int groupIndex = MAX_GROUPS - 1;

    for(int i = len; i > 0 && groupIndex >= 0; i -= GROUP_SIZE)
    {
        unsigned long long value = 0;
        unsigned long long base = 1;

        int start = std::max(0, i - GROUP_SIZE);

        for(int j = i - 1; j >= start; --j)
        {
            value += letterDigitToNumber(k[j]) * base;
            base *= 36;
        }

        w[groupIndex] = value;
        groupIndex--;
    }

    unsigned long long hash = 0;
    for(int i = 0; i < MAX_GROUPS; ++i)
    {
        hash += rValues[i] * w[i];
    }

    return hash;


    }

    // A likely helper function is to convert a-z,0-9 to an integral value 0-35
    HASH_INDEX_T letterDigitToNumber(char letter) const
    {
        if(letter >= 'A' && letter <= 'Z') {
            letter = letter - 'A' + 'a';
        }
        if(letter >= 'a' && letter <= 'z') {
            return letter - 'a';
        }
        else if(letter >= '0' && letter <= '9') {
            return letter - '0' + 26;
        }
        return 0;

    }

    // Code to generate the random R values
    void generateRValues()
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator (seed);  // mt19937 is a standard random number generator

        // Simply call generator() [it has an operator()] to get another random number
        for(int i{ 0 }; i < 5; ++i)
        {
            rValues[i] = generator();
        }
    }
};

#endif

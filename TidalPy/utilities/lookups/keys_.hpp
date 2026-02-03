#pragma once

#include <cstdint>

// Define KeyType
typedef uint64_t RefKeyType;

// Correct OFFSET (int32_t) to prevent overflow
const int32_t OFFSET = 32768;

// Helper to pack keys
inline RefKeyType convert_4key(int16_t a, int16_t b, int16_t c, int16_t d)
{
    // Helper to pack the key
    return (static_cast<RefKeyType>(a + OFFSET) << 48) | 
           (static_cast<RefKeyType>(b + OFFSET) << 32) | 
           (static_cast<RefKeyType>(c + OFFSET) << 16) | 
           static_cast<RefKeyType>(d + OFFSET);
};

inline RefKeyType convert_3key(int16_t a, int16_t b, int16_t c)
{
    // Helper to pack the key
    return (static_cast<RefKeyType>(a + OFFSET) << 32) | 
           (static_cast<RefKeyType>(b + OFFSET) << 16) | 
           static_cast<RefKeyType>(c + OFFSET);
};

inline RefKeyType convert_2key(int16_t a, int16_t b)
{
    // Helper to pack the key
    return (static_cast<RefKeyType>(a + OFFSET) << 16) | 
           static_cast<RefKeyType>(b + OFFSET);
};

inline RefKeyType convert_1key(int16_t a)
{
    // Helper to pack the key
    return static_cast<RefKeyType>(a + OFFSET);
};

class c_Key4
{
public:
    int16_t a;
    int16_t b;
    int16_t c;
    int16_t d;
    RefKeyType reference;

    c_Key4() :
        a(0),
        b(0),
        c(0),
        d(0),
        reference(0)
    {
    }

    c_Key4(int16_t a_, int16_t b_, int16_t c_, int16_t d_) : 
        a(a_),
        b(b_),
        c(c_),
        d(d_),
        reference(convert_4key(a_, b_, c_, d_))
    {
    }

};

class c_Key3
{
public:
    int16_t a;
    int16_t b;
    int16_t c;
    RefKeyType reference;

    c_Key3() :
        a(0),
        b(0),
        c(0),
        reference(0)
    {
    }

    c_Key3(int16_t a_, int16_t b_, int16_t c_) : 
        a(a_),
        b(b_),
        c(c_),
        reference(convert_3key(a_, b_, c_))
    {
    }

};

class c_Key2
{
public:
    int16_t a;
    int16_t b;
    RefKeyType reference;

    c_Key2() :
        a(0),
        b(0),
        reference(0)
    {
    }

    c_Key2(int16_t a_, int16_t b_) : 
        a(a_),
        b(b_),
        reference(convert_2key(a_, b_))
    {
    }

};

class c_Key1
{
public:
    int16_t a;
    RefKeyType reference;

    c_Key1() :
        a(0),
        reference(0)
    {
    }

    c_Key1(int16_t a_) : 
        a(a_),
        reference(convert_1key(a_))
    {
    }

};

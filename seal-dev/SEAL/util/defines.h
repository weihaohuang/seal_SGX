#pragma once

// For security reasons one should never throw when decoding fails due to overflow, but in
// some cases this might help in diagnosing problems.
#undef THROW_ON_DECODER_OVERFLOW

// Microsoft Visual Studio 2012 or newer
#if (_MSC_VER >= 1700)

// X64
#ifdef _M_X64

// Use compiler intrinsics for better performance
#define ENABLE_INTRIN

#include <intrin.h>

#pragma intrinsic(_addcarry_u64)
#define ADD_CARRY_UINT64(operand1, operand2, carry, result) _addcarry_u64(          \
    carry,                                                                          \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_subborrow_u64)
#define SUB_BORROW_UINT64(operand1, operand2, borrow, result) _subborrow_u64(       \
    borrow,                                                                         \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_BitScanReverse64)
#define MSB_INDEX_UINT64(result, value) _BitScanReverse64(result, value)

#pragma intrinsic(_umul128)
#define MULTIPLY_UINT64(operand1, operand2, result128) {                            \
    result128[0] = _umul128(                                                        \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(result128 + 1));                      \
}
#define MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                            \
    _umul128(                                                                       \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(hw64));                               \
}

#else //_M_X64

#undef ENABLE_INTRIN

#endif //_M_X64

#endif //_MSC_VER


// GNU GCC/G++
#if (__GNUC__ >= 5) && defined(__cplusplus)

// Read in config.h to disable unavailable intrinsics
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Are intrinsics enabled?
#ifdef ENABLE_INTRIN
#include <x86intrin.h>

#ifdef ENABLE___BUILTIN_CLZLL
// Builtin
#define MSB_INDEX_UINT64(result, value) {                                           \
    *result = 63 - __builtin_clzll(value);                                          \
}
#endif //ENABLE___BUILTIN_CLZLL

#ifdef ENABLE___INT128
// Builtin
#define MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                            \
    *hw64 = static_cast<uint64_t>((static_cast<unsigned __int128>(operand1)         \
            * static_cast<unsigned __int128>(operand2)) >> 64);                     \
}
// Builtin
#define MULTIPLY_UINT64(operand1, operand2, result128) {                            \
    unsigned __int128 product = static_cast<unsigned __int128>(operand1) * operand2;\
    result128[0] = static_cast<uint64_t>(product);                                  \
    result128[1] = product >> 64;                                                   \
}
#endif //ENABLE___INT128

#ifdef ENABLE__ADDCARRY_U64
#define ADD_CARRY_UINT64(operand1, operand2, carry, result) _addcarry_u64(          \
    carry,                                                                          \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))
#endif //ENABLE__ADDCARRY_U64

#ifdef ENABLE__SUBBORROW_U64
// Warning: Note the inverted order of operand1 and operand2
#define SUB_BORROW_UINT64(operand1, operand2, borrow, result) _subborrow_u64(       \
    borrow,                                                                         \
    static_cast<unsigned long long>(operand2),                                      \
    static_cast<unsigned long long>(operand1),                                      \
    reinterpret_cast<unsigned long long*>(result))
#endif //ENABLE__SUBBORROW_U64

#endif //ENABLE_INTRIN

#endif //(defined(__GNUC__) || defined(__GNUG__))


// Use generic functions as (slower) fallback
#ifndef ADD_CARRY_UINT64
#define ADD_CARRY_UINT64(operand1, operand2, carry, result) add_uint64_generic(operand1, operand2, carry, result)
//#pragma message("ADD_CARRY_UINT64 not defined. Using add_uint64_generic (see util/defines.h)")
#endif

#ifndef SUB_BORROW_UINT64
#define SUB_BORROW_UINT64(operand1, operand2, borrow, result) sub_uint64_generic(operand1, operand2, borrow, result)
//#pragma message("SUB_BORROW_UINT64 not defined. Using sub_uint64_generic (see util/defines.h).")
#endif

#ifndef MULTIPLY_UINT64
#define MULTIPLY_UINT64(operand1, operand2, result128) {                           \
    multiply_uint64_generic(operand1, operand2, result128);                        \
}
//#pragma message("MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef MULTIPLY_UINT64_HW64
#define MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                           \
    multiply_uint64_hw64_generic(operand1, operand2, hw64);                        \
}
//#pragma message("MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef MSB_INDEX_UINT64
#define MSB_INDEX_UINT64(result, value) get_msb_index_generic(result, value)
//#pragma message("MSB_INDEX_UINT64 not defined. Using get_msb_index_generic (see util/defines.h).")
#endif

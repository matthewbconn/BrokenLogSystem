//
// Created by matth on 9/8/2020.
//

/*
 * This file is the intermediary between your test bench, which likes real values
 * (at least the ones we can represent with C++ supported data types) and the
 * Log Engine, which uses a fixed point representation.  This file will not
 * synthesize, but should also be #include in Vivado under the testbench files
 *
 * */

#include <cmath>
#include "lognum.h"

// Forward Declarations Here

/*
 * Input: real number
 * Output: signbit part of a lognum
 *         True if input, a real number, is non negative
 * */
template<typename T>
bool getSign(T);

/*
 * Returns '+' for positives, '-' for negatives
 * You shouldn't need this now that you can use conversion functions
 * */
template<typename T>
char getSignSymbol(T);


/*
 * Input: real number
 * Output: fixed point part of the lognum
 * */
template<typename T>
// see warning about a change you made while you were tired...
double getlog2val(T);


/*
 * Input: a real number
 * Output: it's lognum representation
 * */
template<typename T>
lognum toLogNum(T);

// For some reason these forward declarations weren't working

///////////////////////////////////////////////////////////////////

// Function Definitions here
template <typename T>
bool getSign(T item) {
    return (item>0);
}

template <typename T>
char getSignSymbol(T item) {
    if (item>0){return ('+');}
    return('-');
}

// take in a real number, spit out a log number
template <typename T>
// changed return type from int to double while you were tired
double getlogval(T item) {

    if (item > 0) {
        return log(item)/log(BASE);
        //return log2(item);
    } else if (item < 0){
        return log(-1*item)/log(BASE);
        return log2(-1*item);
    } else { // item == 0
        return INT16_MIN;
    }
}

template<typename T>
lognum toLogNum(T realnum) {
    return lognum(getSign(realnum),getlogval(realnum));
}
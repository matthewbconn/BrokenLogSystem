//
// Created by matth on 9/7/2020.
//

#include "lognum.h"

// with no better information, generate for x(real) = 1
// so sgn(x) = 1, log2(x) = 0
lognum::lognum(): signBit(true), logval(0) {}

// copy ctor...needed for MAC
lognum::lognum(const lognum &original): signBit(original.getSignBit()),logval(original.logval) {}

// ctor to instantiate for a specific value (someone else's job to convert)
lognum::lognum(bool sign, fixedtype fixedNum):
        signBit(sign), logval(fixedNum){}

lognum lognum::addReals(lognum v1, lognum v2) {
    lognum result;

    result.signBit = v1.logval > v2.logval ? v1.signBit : v2.signBit;

    fixedtype d = v1.logval > v2.logval ?
            v1.logval - v2.logval : v2.logval - v1.logval;

    fixedtype larger = v1.logval > v2.logval ? v1.logval : v2.logval;

    result.logval  = v1.signBit == v2.signBit ?
            (larger + deltaPlus(d)): (larger + deltaMinus(d));

    return result;
}

lognum lognum::multiplyReals(lognum v1, lognum v2) {
    return lognum(! (v1.signBit xor v2.signBit), v1.logval + v2.logval);
}

void lognum::MAC(lognum A, lognum B) {
    lognum mult = multiplyReals(A,B);
    (*this)=(addReals(*this,mult));
}

// Q: why are these functions defined?
// A: for use in copy ctor
bool lognum::getSignBit() const  {
    return signBit;
}

double lognum::getLogval() {
    return logval.to_double();
}

/*
 * the bitshift will be performed after
 * d is implicitly cast to an integer
 * using the FLOOR round scheme.
 * Note that d >= 0 from how we call this func
 * */
fixedtype deltaPlus(fixedtype d) {
    fixedtype bitshift(1);
    bitshift = (bitshift >> d);
    return bitshift;
}

fixedtype deltaMinus(fixedtype d) {
/*  old implementation: shift then multiply
    fixedtype bitshift(1.5);
    bitshift = -1*(bitshift >> d);
*/

// new implementation: start at -1.5, then shift
    fixedtype bitshift(-1.5);
    bitshift = (bitshift >> d);
    // need to 'OR' mask to set MSB == '1' to preserve arithmetic shift
    return (bitshift | MIN_VAL);
}

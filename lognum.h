//
// Created by matth on 9/7/2020.
//

#ifndef VIVADOLOGSYSTEM_LOGNUM_H
#define VIVADOLOGSYSTEM_LOGNUM_H

#define WBITS 8
#define INTBITS 4
#define FRACBITS (WBITS-INTBITS)
#define NBITS 1
#define BASE 2

#include "ap_fixed.h"
#include "cmath"    // ill advised to use this
//#include "hls_math.h" // for some reason this causes a dumb issue

typedef ap_fixed<WBITS,INTBITS,AP_RND_ZERO,AP_SAT,NBITS> fixedtype;

// Get a mask that is "100000.000000"
static fixedtype MIN_VAL = fixedtype(-1*pow(BASE,INTBITS-1));

class lognum {
public:
    // Default ctor, copy ctor, component-wise ctor
    lognum();
    lognum(const lognum &original);
    lognum(bool sign, fixedtype);

    // basic operations
    static lognum addReals(lognum v1, lognum v2);
    static lognum multiplyReals (lognum v1, lognum v2);
    void MAC(lognum A, lognum B);

    // access to components
    bool getSignBit() const;
    double getLogval();

private:
    // defining fields
    bool signBit;
    fixedtype logval;

};

// Delta function operating on a fixed type
fixedtype deltaPlus(fixedtype);
fixedtype deltaMinus(fixedtype);

#endif //VIVADOLOGSYSTEM_LOGNUM_H

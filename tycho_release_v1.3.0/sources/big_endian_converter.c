/*
 * big_endian_converter.c
 *
 * Author: Wolfgang Kapferer
 * Thanks to 
 * http://www.linuxquestions.org/questions/programming-9/c-function-to-reverse-the-byte-order-in-a-double-825743/
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "variables_global.h"
#include "prototypes.h"

/*!
 converts double values to BigEndian 
 VTK binaries are in BigEndian format
 as well as AMIRA MESH Files
 */
long long int ntohll(const long long int data) {

    enum {
        TYP_INIT, TYP_SMLE, TYP_BIGE
    };

    union {
        long long int ull;
        uint8_t c[8];
    } x;

    // Test if on BigEndian machine
    static int typ = TYP_INIT;

    if (typ == TYP_INIT) {
        x.ull = 0x01;
        typ = (x.c[7] == 0x01) ? TYP_BIGE : TYP_SMLE;
    }

    // If system is BigEndian; return data as is.
    if (typ == TYP_BIGE) {
        return data;
    }

    // convert data to Big Endian
    x.ull = data;

    int8_t c = 0;
    c = x.c[0];
    x.c[0] = x.c[7];
    x.c[7] = c;
    c = x.c[1];
    x.c[1] = x.c[6];
    x.c[6] = c;
    c = x.c[2];
    x.c[2] = x.c[5];
    x.c[5] = c;
    c = x.c[3];
    x.c[3] = x.c[4];
    x.c[4] = c;

    return x.ull;
}

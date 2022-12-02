#ifndef TRAILS_H
#define TRAILS_H
#include <random>
#include <iostream>
#include <stdint.h>


#if !defined SIZE
#error SIZE is not defined.
#endif

using namespace std;

static void conversion(uint8_t alpha[30][16],uint32_t trail_round,int size);
static bool sanityCheck(uint8_t alpha[30][16],uint32_t trail_round,int size);
#if SIZE == 4
void SK_4_2021_1179_1(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2021_1179_2(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2021_1179_3(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2019_49_9(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2019_49_12(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2019_49_13(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2018_390_Table4(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_4_2018_390_Table6(uint8_t alpha[30][16], uint32_t& trail_rounds);
#elif SIZE == 8
void SK_8_12_2019_25(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_13_2019_25(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_21_2019_25(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_2019_49_21(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_2018_390_Table10(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_2018_390_Table15(uint8_t alpha[30][16], uint32_t& trail_rounds);
void SK_8_2018_390_Table16(uint8_t alpha[30][16], uint32_t& trail_rounds);
#endif
#endif
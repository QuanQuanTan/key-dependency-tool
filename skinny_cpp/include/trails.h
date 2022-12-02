#ifndef TRAILS_H
#define TRAILS_H

#include <iostream>

using namespace std;

#if !defined SIZE
#error SIZE is not defined.
#endif

static void conversion(uint8_t alpha[20][4][4],uint8_t trail_rounds,uint8_t key_diff[4][4][4]);
static bool sanityCheck(uint8_t alpha[20][4][4],uint8_t trail_rounds,uint8_t key_diff[4][4][4]);

#if SIZE == 4
void TK3_2016_1108(uint8_t alpha[20][4][4], uint8_t key_diff_1[4][4][4] ,uint32_t &upper_trail_rounds, int &TK_NUM, uint8_t gamma[20][4][4], uint8_t key_diff_2[4][4][4] ,uint32_t &lower_trail_rounds);
void TK2_2016_1108(uint8_t alpha[20][4][4], uint8_t key_diff_1[4][4][4] ,uint32_t &upper_trail_rounds, int &TK_NUM, uint8_t gamma[20][4][4], uint8_t key_diff_2[4][4][4] ,uint32_t &lower_trail_rounds);
// Differential characteristics for SKINNY-64 //
void Corrected_TK1_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void SK_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void TK1_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void TK1_2020_1402_CP(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void TK2_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void TK3_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_1_17_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_1_17_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_1_18_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_1_19_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_2_17_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_2_17_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_2_18_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_2_19_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_22_L_BMD1_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_22_U_BMD1_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_22_L_BMD2_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_22_U_BMD2_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_22_U_BMD3_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_23_U_BMD1_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK3_1_23_U_BMD2_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_4TK2_18_L_2021_656(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM, uint8_t u, uint8_t v);
void get_4TK3_22_L_2021_656(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM, uint8_t u, uint8_t v);


// Differential characteristics of SKINNY-128
#elif SIZE == 8
void get_8TK2_19_L_2021_656(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM, uint8_t u, uint8_t v, uint8_t w);
void get_8TK3_22_L_2021_656(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM, uint8_t u, uint8_t v);
void get_8TK2_1_18_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_18_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_19_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_19_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_20_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_20_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_21_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_1_21_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void SK_tosc_2017_i4_99_129(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void TK1_8_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void TK2_8_2020_1402(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_18_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_18_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_19_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_19_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_20_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_20_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK2_2_21_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_22_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_22_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_23_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_23_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_24_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_24_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_25_L_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_8TK3_1_25_U_2020_1317(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_2021_856_Table20_8TK3_U(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_2021_856_Table20_8TK3_L(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_2021_856_Table22_8TK2_U(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM);
void get_2021_856_Table22_8TK2_L(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4] ,uint32_t& trail_rounds, int &TK_NUM, uint8_t u, uint8_t v, uint8_t w1, uint8_t w2, uint8_t w3);
#endif
#endif


#ifndef SKINNY_H
#define SKINNY_H

#include <stdint.h>
#include <iostream>
#include <cstring>

using namespace std;

#ifndef SIZE
#error SIZE is not defined.
#endif



void perm(uint8_t key_state[4][4]);
void invPerm(uint8_t key_state[4][4]);
void AddConstant(uint8_t state[4][4], int r);
void AddRoundTweakey(uint8_t state[4][4],uint8_t TK[4][4][4]);
void reset_key(uint8_t key_tmp[4][4][4], uint8_t key[4][4][4]);
void reset_plaintext(uint8_t p_tmp[4][4], uint8_t p[4][4]);
void SR(uint8_t state[4][4]);
void invSR(uint8_t state[4][4]);
void MC(uint8_t state[4][4]);
void invMC(uint8_t state[4][4]);
uint8_t getConstants(uint8_t val);
uint8_t getPermSchedule(uint8_t val);
uint8_t getInvPermSchedule(uint8_t val);
// for 4 bit skinny
#if SIZE==4
void computeDDT(int DDT[16][16]);
uint8_t LFSR(uint8_t n, int v);
uint8_t invLFSR(uint8_t n, int v);
void Substitution(uint8_t state[4][4]);
void invSubstitution(uint8_t state[4][4]);
void key_schedule_round_function(uint8_t key_states[4][4][4]);
void inv_key_schedule_round_function(uint8_t key_states[4][4][4]);
void RoundFunction(uint8_t state[4][4], uint8_t TK[4][4][4], int r);
void InvRoundFunction(uint8_t state[4][4], uint8_t TK[4][4][4], int r);
void Skinny_Enc(uint8_t state[4][4], uint8_t TK[4][4][4], int nr);
void Skinny_Dec(uint8_t state[4][4], uint8_t TK[4][4][4], int nr);
uint8_t getSbox(uint8_t val);
uint8_t getInvSbox(uint8_t val);


// for the 8 bit skinny
#elif SIZE==8
void computeDDT8(int DDT[256][256]);
uint8_t getSbox8(uint8_t val);
uint8_t getInvSbox8(uint8_t val);
uint8_t LFSR8(uint8_t n, int v);
uint8_t invLFSR8(uint8_t n, int v);
void key_schedule_round_function8(uint8_t key_states[4][4][4]);
void inv_key_schedule_round_function8(uint8_t key_states[4][4][4]);
void Substitution8(uint8_t state[4][4]);
void invSubstitution8(uint8_t state[4][4]);
void RoundFunction8(uint8_t state[4][4], uint8_t TK[4][4][4], int r);
void InvRoundFunction8(uint8_t state[4][4], uint8_t TK[4][4][4], int r);
void Skinny_Enc8(uint8_t state[4][4], uint8_t TK[4][4][4], int nr);
void Skinny_Dec8(uint8_t state[4][4], uint8_t TK[4][4][4], int nr);
#endif
#endif
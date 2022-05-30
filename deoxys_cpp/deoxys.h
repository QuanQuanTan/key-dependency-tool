#ifndef DEOXYS_H
#define DEOXYS_H

uint8_t getSbox(uint8_t val);
uint8_t getInvSbox(uint8_t val);
uint8_t getConstants(uint8_t val);
uint8_t getPermSchedule(uint8_t val);
void transpose(uint8_t a[4][4]);
uint8_t mul2(uint8_t val);
uint8_t mul3(uint8_t val);
uint8_t LFSR(uint8_t n, int v);
uint8_t invLFSR(uint8_t n, int v);
void perm(uint8_t key_state[4][4]);
void invPerm(uint8_t key_state[4][4]);
void keyScheduleRoundFunction(uint8_t key_states[3][4][4]);
void Substitution(uint8_t state[4][4]);
void invSubstitution(uint8_t state[4][4]);
void SR(uint8_t state[4][4]);
void invSR(uint8_t state[4][4]);
void MC(uint8_t state[4][4]);
void invMC(uint8_t state[4][4]);
void AddRoundTweakey(uint8_t state[4][4], uint8_t TK[3][4][4]);
void AddConstant(uint8_t state[4][4], int r);
void RoundFunction(uint8_t state[4][4], uint8_t TK[3][4][4], int r);
void InvRoundFunction(uint8_t state[4][4], uint8_t TK[3][4][4], int r);
void deoxys_encrypt(uint8_t state[4][4], uint8_t TK[3][4][4], int nr);
void deoxys_decrypt(uint8_t state[4][4], uint8_t TK[3][4][4], int nr);
void computeDDT(uint32_t DDT[256][256]);
void computeXDDT(uint8_t x, uint8_t y, uint32_t xddt[256]);
void computeYDDT(uint8_t x, uint8_t y, uint32_t yddt[256]);
#endif
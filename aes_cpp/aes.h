#ifndef AES_H
#define AES_H

uint8_t getSbox(uint8_t val);
uint8_t getInvSbox(uint8_t val);
uint8_t getConstants(uint8_t val);
void transpose(uint8_t a[4][4]);
uint8_t mul2(uint8_t val);
uint8_t mul3(uint8_t val);
void keyScheduleRoundFunction(uint8_t key_states[3][4][4]);
void Substitution(uint8_t state[4][4]);
void invSubstitution(uint8_t state[4][4]);
void SR(uint8_t state[4][4]);
void invSR(uint8_t state[4][4]);
void MC(uint8_t state[4][4]);
void invMC(uint8_t state[4][4]);
void AddRoundKey(uint8_t state[4][4], uint8_t key[4][4]);
void RoundFunction(uint8_t state[4][4], uint8_t key[4][4]);
void InvRoundFunction(uint8_t state[4][4], uint8_t key[4][4]);
void deoxys_encrypt(uint8_t state[4][4], uint8_t key[4][8], int nr, int keySize);
void deoxys_decrypt(uint8_t state[4][4], uint8_t key[4][8], int nr, int keySize);
void computeDDT(uint32_t DDT[256][256]);
void computeXDDT(uint8_t x, uint8_t y, uint32_t xddt[256]);
void computeYDDT(uint8_t x, uint8_t y, uint32_t yddt[256]);
#endif
#ifndef GIFT_H
#define GIFT_H
uint8_t getSbox(uint8_t val);
uint8_t getInvSbox(uint8_t val);
uint8_t getConstants(uint8_t val);
uint8_t getPerm(uint8_t val);
uint8_t getPerm8(uint8_t val);
uint8_t getinvPerm(uint8_t val);
uint8_t getinvPerm8(uint8_t val);
void Substitution(uint8_t state[16]);
void invSubstitution(uint8_t state[16]);
void Permutation(uint8_t state[16]);
void invPermutation(uint8_t state[16]);
void AddRoundKey(uint8_t state[16], uint16_t key[2]);
void AddConstant(uint8_t state[16], uint8_t constant);
void keyScheduleRoundFunction(uint16_t masterKey[8]);
void RoundFunction(uint8_t state[16], uint16_t key[2], int nr);
void invRoundFunction(uint8_t state[16], uint16_t key[2], int nr);
void gift_encrypt(uint8_t state[16], uint16_t masterKey[8]);
void gift_decrypt(uint8_t state[16], uint16_t masterKey[8]);
void Substitution8(uint8_t state[16]);
void invSubstitution8(uint8_t state[16]);
void Permutation8(uint8_t state[16]);
void invPermutation8(uint8_t state[16]);
void AddRoundKey8(uint8_t state[16], uint16_t key[4]);
void AddConstant8(uint8_t state[16], uint8_t constant);
void RoundFunction8(uint8_t state[16], uint16_t key[4], int nr);
void invRoundFunction8(uint8_t state[16], uint16_t key[4], int nr);
void gift_encrypt8(uint8_t state[16], uint16_t masterKey[8]);
void gift_decrypt8(uint8_t state[16], uint16_t masterKey[8]);
void computeDDT(uint32_t DDT[16][16]);
void computeXDDT(uint8_t x, uint8_t y, uint32_t xddt[16]);
void computeYDDT(uint8_t x, uint8_t y, uint32_t yddt[16]);
#endif
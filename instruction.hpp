#ifndef RISC_V_INSTRUCTION_HPP
#define RISC_V_INSTRUCTION_HPP

#include <iostream>
#include <cstring>

typedef unsigned int u32;
typedef unsigned char u8;

const int REG_SIZE = 32;
const int MEM_SIZE = 1048576;
const int RS_SIZE = 64;
const int ROB_SIZE = 64;
const int LSB_SIZE = 64;

const char *OPNAME[] = {
    "NONE",   //                                                                     0
    "HALT",   //                                                                     1 
    "LUI",    // U    Load Upper Immediate                                           2
    "AUIPC",  // U    Add Upper Immediate to PC
    "JAL",    // J    Jump & Link
    "JALR",   // I    Jump & Link Register                                           5
    "BEQ",    // B    Branch Equal                                                   6
    "BNE",    // B    Branch Not Equal
    "BLT",    // B    Branch Less Than
    "BGE",    // B    Branch Greater than or Equal
    "BLTU",   // B    Branch Less than Unsigned
    "BGEU",   // B    Branch Greater than or Equal Unsigned                          11
    "ADDI",   // I    ADD Immediate                                                  12
    "SLTI",   // I    Set Less than Immediate
    "SLTIU",  // I    Set Less than Immediate Unsigned
    "XORI",   // I    XOR Immediate
    "ORI",    // I    OR Immediate
    "ANDI",   // I    AND Immediate
    "SLLI",   // I    Shift Left Immediate
    "SRLI",   // I    Shift Right Immediate
    "SRAI",   // I    Shift Right Arith Immediate                                    20
    "ADD",    // R    ADD                                                            21
    "SUB",    // R    SUBtract
    "SLL",    // R    Shift Left
    "SLT",    // R    Set Less than
    "SLTU",   // R    Set Less than Unsigned
    "XOR",    // R    XOR
    "SRL",    // R    Shift Right
    "SRA",    // R    Shift Right Arithmetic
    "OR",     // R    OR
    "AND",    // R    AND                                                            30
    "LB",     // I    Load Byte                                                      31
    "LH",     // I    Load Halfword
    "LW",     // I    Load Word
    "LBU",    // I    Load Byte Unsigned
    "LHU",    // I    Load Halfword Unsigned                                         35
    "SB",     // S    Store Byte                                                     36
    "SH",     // S    Store Halfword
    "SW"      // S    Store Word                                                     38
};

enum insMnemonic {
    NONE,   //                                                                     0
    HALT,   //                                                                     1 

    LUI,    // U    Load Upper Immediate                                           2
    AUIPC,  // U    Add Upper Immediate to PC
    JAL,    // J    Jump & Link
    JALR,   // I    Jump & Link Register                                           5

    BEQ,    // B    Branch Equal                                                   6
    BNE,    // B    Branch Not Equal
    BLT,    // B    Branch Less Than
    BGE,    // B    Branch Greater than or Equal
    BLTU,   // B    Branch Less than Unsigned
    BGEU,   // B    Branch Greater than or Equal Unsigned                          11

    ADDI,   // I    ADD Immediate                                                  12
    SLTI,   // I    Set Less than Immediate
    SLTIU,  // I    Set Less than Immediate Unsigned
    XORI,   // I    XOR Immediate
    ORI,    // I    OR Immediate
    ANDI,   // I    AND Immediate
    SLLI,   // I    Shift Left Immediate
    SRLI,   // I    Shift Right Immediate
    SRAI,   // I    Shift Right Arith Immediate                                    20

    ADD,    // R    ADD                                                            21
    SUB,    // R    SUBtract
    SLL,    // R    Shift Left
    SLT,    // R    Set Less than
    SLTU,   // R    Set Less than Unsigned
    XOR,    // R    XOR
    SRL,    // R    Shift Right
    SRA,    // R    Shift Right Arithmetic
    OR,     // R    OR
    AND,    // R    AND                                                            30

    LB,     // I    Load Byte                                                      31
    LH,     // I    Load Halfword
    LW,     // I    Load Word
    LBU,    // I    Load Byte Unsigned
    LHU,    // I    Load Halfword Unsigned                                         35

    SB,     // S    Store Byte                                                     36
    SH,     // S    Store Halfword
    SW      // S    Store Word                                                     38
};

inline u32 getRange(u32 x, const int &l, const int &r) {
    x <<= 31 - r, x >>= 31 - r + l;
    return x;
}

inline u32 signExtend(const u32 &x, const int &l, const int &r) {
    return (((x >> (r - l)) & 1) ? ((x << l) | (~((1u << r) - 1))) : (x << l));
}

inline u32 zeroExtend(const u32 &x, const int &l) {
    return x << l;
}

struct InsNode {
    u32 imm, rs1, rs2, rd;
    insMnemonic type;
};

struct Memory {
    u8 x[MEM_SIZE];
    Memory() { memset(x, 0, sizeof x); }
    void scan() {
        std::string str;
        int pos = 0;
        while (std::cin >> str) {
            if (str[0] == '@') {
                char *p;
                pos = strtoul(str.substr(1, 8).c_str(), &p, 16);
            } else {
                char *p;
                x[pos++] = strtoul(str.c_str(), &p, 16);
            }
        }
    }
    u32 getWord(u32 p) const {
        return x[p] + (u32(x[p + 1]) << 8) + (u32(x[p + 2]) << 16) + (u32(x[p + 3]) << 24);
    }
    u32 getHalfwordUnsigned(u32 p) const {
        return x[p] + (u32(x[p + 1]) << 8);
    }
    u32 getByteUnsigned(u32 p) const {
        return x[p];
    }
    u32 getHalfword(u32 p) const {
        return signExtend(x[p] + (u32(x[p + 1]) << 8), 0, 15);
    }
    u32 getByte(u32 p) const {
        return signExtend(x[p], 0, 7);
    }
    void setWord(u32 p, u32 val) {
        x[p] = val & 255;
        x[p + 1] = (val >> 8) & 255;
        x[p + 2] = (val >> 16) & 255;
        x[p + 3] = (val >> 24) & 255;
    }
    void setHalfword(u32 p, u32 val) {
        x[p] = val & 255;
        x[p + 1] = (val >> 8) & 255;
    }
    void setByte(u32 p, u32 val) {
        x[p] = val & 255;
    }
};

inline void instructionDecode(const Memory &mem, u32 &pc, InsNode &now) {
    u32 ins = mem.getWord(pc);
    u32 opcode = getRange(ins, 0, 6);
    u32 func3 = getRange(ins, 12, 14);
    u32 func7 = getRange(ins, 25, 31);

    u32 &rs1 = now.rs1;
    u32 &rs2 = now.rs2;
    u32 &rd = now.rd;
    u32 &imm = now.imm;
    insMnemonic &insType = now.type;
    rs1 = getRange(ins, 15, 19);
    rs2 = getRange(ins, 20, 24);
    rd = getRange(ins, 7, 11);
    imm = 0;

    if (ins == 0xff00513)
        return insType = HALT, void();

    insType = NONE;

    if (opcode == 0b0110011) { // R
        if (func3 == 0b000) {
            if (func7 == 0b0000000)
                insType = ADD;
            if (func7 == 0b0100000)
                insType = SUB;
        }
        if (func3 == 0b001)
            insType = SLL;
        if (func3 == 0b010)
            insType = SLT;
        if (func3 == 0b011)
            insType = SLTU;
        if (func3 == 0b100)
            insType = XOR;
        if (func3 == 0b101) {
            if (func7 == 0b0000000)
                insType = SRL;
            if (func7 == 0b0100000)
                insType = SRA;
        }
        if (func3 == 0b110)
            insType = OR;
        if (func3 == 0b111)
            insType = AND;
    }
    
    if (opcode == 0b0010011) { // I
        if (func3 == 0b000) {
            insType = ADDI;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b010) {
            insType = SLTI;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b011) {
            insType = SLTIU;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b100) {
            insType = XORI;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b110) {
            insType = ORI;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b111) {
            insType = ANDI;
            imm = signExtend(getRange(ins, 20, 31), 0, 11);
        }
        if (func3 == 0b001) {
            insType = SLLI;
            imm = getRange(ins, 20, 24);
        }
        if (func3 == 0b101) {
            if (func7 == 0b0000000) {
                insType = SRLI;
                imm = getRange(ins, 20, 24);
            }
            if (func7 == 0b0100000) {
                insType = SRAI;
                imm = getRange(ins, 20, 24);
            }
        }
    }

    if (opcode == 0b0100011) { // S
        imm = zeroExtend(getRange(ins, 7, 11), 0) | signExtend(getRange(ins, 25, 31), 5, 11);
        if (func3 == 0b000)
            insType = SB;
        if (func3 == 0b001)
            insType = SH;
        if (func3 == 0b010)
            insType = SW;
    }

    if (opcode == 0b0000011) { // I
        imm = signExtend(getRange(ins, 20, 31), 0, 11);
        if (func3 == 0b000)
            insType = LB;
        if (func3 == 0b001)
            insType = LH;
        if (func3 == 0b010)
            insType = LW;
        if (func3 == 0b100)
            insType = LBU;
        if (func3 == 0b101)
            insType = LHU;
    }

    if (opcode == 0b1100011) { // B
        imm = zeroExtend(getRange(ins, 8, 11), 1) |
                zeroExtend(getRange(ins, 25, 30), 5) |
                zeroExtend(getRange(ins, 7, 7), 11) |
                signExtend(getRange(ins, 31, 31), 12, 12);
        if (func3 == 0b000)
            insType = BEQ;
        if (func3 == 0b001)
            insType = BNE;
        if (func3 == 0b100)
            insType = BLT;
        if (func3 == 0b101)
            insType = BGE;
        if (func3 == 0b110)
            insType = BLTU;
        if (func3 == 0b111)
            insType = BGEU;
    }

    if (opcode == 0b1100111) {
        insType = JALR;
        imm = signExtend(getRange(ins, 20, 31), 0, 11);
    }

    if (opcode == 0b1101111) {
        insType = JAL;
        imm = zeroExtend(getRange(ins, 21, 30), 1) |
                zeroExtend(getRange(ins, 20, 20), 11) |
                zeroExtend(getRange(ins, 12, 19), 12) |
                signExtend(getRange(ins, 31, 31), 20, 20);
    }

    if (opcode == 0b0010111) {
        insType = AUIPC;
        imm = signExtend(getRange(ins, 12, 31), 12, 31);
    }

    if (opcode == 0b0110111) {
        insType = LUI;
        imm = signExtend(getRange(ins, 12, 31), 12, 31);
    }
}

#endif
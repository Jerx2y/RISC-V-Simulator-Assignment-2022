#ifndef RISC_V_SIMULATOR_HPP
#define RISC_V_SIMULATOR_HPP

#include <iostream>
#include <cstring>
#include <cassert>
#include <bitset>
#include <queue>

#include "instruction.hpp"

using std::cout;
using std::endl;

struct Register {
    u32 vpre[REG_SIZE], vnow[REG_SIZE];
    int rpre[REG_SIZE], rnow[REG_SIZE];
    Register() {
        memset(vpre, 0, sizeof vpre);
        memset(vnow, 0, sizeof vnow);
        memset(rpre, -1, sizeof rpre);
        memset(rnow, -1, sizeof rnow);
    }
    void tick() {
        vnow[0] = 0, rnow[0] = -1;
        memcpy(vpre, vnow, sizeof vpre);
        memcpy(rpre, rnow, sizeof rpre);
    }
    u32 get(int x) { return vpre[x]; }
    int getr(int x) { return rpre[x]; }
    void set(int x, u32 v) { vnow[x] = v; rnow[x] = -1; }
    void setr(int x, int name) { rnow[x] = name; }
    void clean() {
        memset(rpre, -1, sizeof rpre);
        memset(rnow, -1, sizeof rnow);
    }
};

template<class T, unsigned int SIZE>
class queue {
public:
    int hd, tl, sz;
    T q[SIZE];
    queue() { clear(); }
    void clear() { hd = tl = sz = 0; }
    bool empty() { return sz == 0; }
    int size() { return sz; }
    bool full() { return sz == SIZE; }
    T &front() { return q[(hd + 1) % SIZE]; }
    T &back() { return q[tl]; }
    int push_back(const T &x) {
        tl = (tl + 1) % SIZE;
        q[tl] = x, ++sz;
        return tl;
    }
    void pop_front() {
        if (empty()) return ;
        hd = (hd + 1) % SIZE, --sz;
    }
    void pop_back() {
        if (empty()) return ;
        tl = (tl + SIZE - 1) % SIZE, --sz;
    }
    int getidx() { return (tl + 1) % SIZE; }
};

class ReservationStation {
public:
    struct node {
        bool busy;
        insMnemonic op;
        u32 Vj, Vk;
        int Qj, Qk, name;
        node() { busy = 0, Qj = Qk = -1, Vj = Vk = 0; }
    } pre[RS_SIZE], now[RS_SIZE];
    void tick() {
        for (int i = 0; i < RS_SIZE; ++i)
            pre[i] = now[i];
    }
    void push(const node &tmp) {
        for (int i = 0; i < RS_SIZE; ++i) {
            if (now[i].busy) continue;
            now[i] = tmp;
            return ;
        }
        assert(0); // TODO
    }
    void clean() {
        for (int i = 0; i < RS_SIZE; ++i)
            pre[i].busy = now[i].busy = false;
    }
    bool full() {
        for (int i = 0; i < RS_SIZE; ++i)
            if (!now[i].busy) return false;
        return true;    
    }
};

class LoadStoreBuffer {
public:
    struct node {
        insMnemonic op;
        bool busy, ready;
        u32 Vj, Vk;
        int Qj, Qk;
        node() { busy = ready = false, Qj = Qk = -1; }
    };
    queue<node, LSB_SIZE> pre, now;
    bool full() { return now.full(); }
    void tick() { pre = now; }
    int push(node x) { x.busy = 1; return now.push_back(x); }
    void pop() { now.pop_front(), now.q[now.hd].busy = 0; }
    void setready(int idx) { now.q[idx].ready = 1; }
    void clean() {
        while (!now.empty() && !now.back().ready)
            now.q[now.tl].busy = 0, now.pop_back();
        for (int i = 0; i < LSB_SIZE; ++i) {
            if (now.q[i].busy && 31 <= now.q[i].op && now.q[i].op <= 35)
                now.q[i].busy = 0, now.q[i].ready = 0;
            if (now.q[i].busy && 36 <= now.q[i].op && now.q[i].op <= 38 && !now.q[i].ready)
                now.q[i].busy = 0;
        }
        while (!now.empty() && !now.front().busy)
            pop();
        pre = now;
    }
};

class ReorderBuffer {
public:
    struct node {
        bool ready = 0;
        int dest, npc = 0;
        u32 value;
    };
    queue<node, ROB_SIZE> pre, now;
    void tick() { pre = now; }
    int push(const node &x) { return now.push_back(x); }
    bool full() { return now.full(); }
    int getidx() { return now.getidx(); }
    void clean() { pre.clear(), now.clear(); }
};

class ArithmeticLogicUnit {
public:
    struct node {
        insMnemonic op;
        u32 x, y;
        int name;
        bool stall;
    } pre, now;
    void tick() { pre = now; }
    void set(u32 x, u32 y, insMnemonic op, int name) {
        now.x = x, now.y = y, now.op = op, now.name = name, now.stall = 0;
    }
    void stall() { now.stall = 1; }
    void clean() { pre.stall = now.stall = 1; }
};

struct MemoryController {
    insMnemonic op;
    int clk;
    u32 pos, val;
    bool stall = 1;
    void tick() { if (!stall) --clk; }
    void setstall() { stall = 1;}
    void set(insMnemonic op_, u32 pos_, u32 val_) {
        op = op_, pos = pos_, val = val_, clk = 2, stall = 0;
    }
    void clean() {
        if (31 <= op && op <= 35)
            setstall();
    }
};

struct CommonDataBus {
    int name = -1; u32 val;
};

struct BranchPredictor {
    char cnt[65536];
    BranchPredictor() { memset(cnt, 0, sizeof cnt); }
    void calc(u32 pos, int val) {
        if (val) { if (cnt[pos & 65535] < 3) ++cnt[pos]; }
        else { if (cnt[pos & 65535] > 0) --cnt[pos]; }
    }
    bool predict(u32 pos) {
        return cnt[pos & 65535] & 2;
    }
};

class Simulator {
private:
    Register reg;
    Memory mem;
    ReservationStation RS;
    LoadStoreBuffer LSB;
    ReorderBuffer RoB;
    ArithmeticLogicUnit ALU;
    MemoryController MC;
    CommonDataBus alubus, lsbus;
    BranchPredictor BP;
    u32 clk = 0, pc = 0;
    int rpc = -1;
public:

/*
 * 故事的小黄花
 * 从出生那年就飘着
 * 童年的荡秋千
 * 随 记忆一直晃到现在
 * 2557176 5677776765
 * 吹着前奏望着天空
 * 我想起花瓣试着掉落
 * 为你翘课的那一天
 * 风吹的那一天
 * 教室的那一间
 * 我怎么看不见
 * 消失的下雨天
 * 我好想再淋一遍
 */
    
    void scanmem() { mem.scan(); }

    void run() {
        while (++clk) {
            update();
            ex(); // ALU and ALUbus
            run_lsb(); // LSbus
            run_rs();
            run_rob(); // commit
            issue(); 
            /*
             * commit must be in front of issue, or :
             * 1. issue changed the now of a register-r to name;
             *    commit change the now of a register-r to -1;
             *    then the register-r == -1, meaning that the issue is useless.
             *    in a word, the name of new issue must cover -1
             * 2. branch predictor will clean all the thing issue have written.
             */
        }
    }

    void update() {
        reg.tick();
        RS.tick();
        LSB.tick();
        RoB.tick();
        ALU.tick();
        MC.tick();
        alubus.name = lsbus.name = -127;
    }

    void ex() {
        const ArithmeticLogicUnit::node &alu = ALU.pre;
        CommonDataBus &bus = alubus;
        if (alu.stall) return bus.name = -127, void();
        bus.name = alu.name;

        const u32 &x = alu.x;
        const u32 &y = alu.y;
        const u32 &op = alu.op;

        if (op == BEQ) bus.val = (x == y);
        if (op == BNE) bus.val = (x != y);
        if (op == BLT) bus.val = ((int)x < (int)y);
        if (op == BGE) bus.val = ((int)x >= (int)y);
        if (op == BLTU) bus.val = (x < y);
        if (op == BGEU) bus.val = (x >= y);
        if (op == JALR) bus.val = (x + y) & (~1);
        if (op == ADDI) bus.val = x + y;
        if (op == SLTI) bus.val = ((int)x < (int)y);
        if (op == SLTIU) bus.val = (x < y);
        if (op == XORI) bus.val = (x ^ y);
        if (op == ORI) bus.val = (x | y);
        if (op == ANDI) bus.val = (x & y);
        if (op == SLLI) bus.val = x << y;
        if (op == SRLI) bus.val = x >> y;
        if (op == SRAI)  {
            u32 tmp = ((x >> 31) & 1) ? (~((1u << (32 - y)) - 1)) : 0;
            bus.val = (x >> y) | tmp;
        }
        if (op == ADD) bus.val = x + y;
        if (op == SUB) bus.val = x - y;
        if (op == SLL) bus.val = x << (y & 31);
        if (op == SLT) bus.val = ((int)x < (int)y) ;
        if (op == SLTU) bus.val = (x < y);
        if (op == XOR) bus.val = x ^ y;
        if (op == SRL) bus.val = x >> (y & 31);
        if (op == SRA) {
            u32 tmp = ((x >> 31) & 1) ? (~((1u << (32 - y)) - 1)) : 0;
            bus.val = (x >> (y & 31)) | tmp;
        }
        if (op == OR) bus.val = x | y;
        if (op == AND) bus.val = x & y;

        listen(bus);
    }

    void listen(const CommonDataBus &bus) {
        // pc

        if (bus.name == -127)
            return ;

        if (bus.name == -1) {
            rpc = -1, pc = bus.val;
            return ;
        }

        // RoB
        if (0 <= bus.name && bus.name < ROB_SIZE) {
            RoB.now.q[bus.name].ready = 1;
            RoB.now.q[bus.name].value = bus.val;
        }

        // RS
        for (int i = 0; i < RS_SIZE; ++i) if (RS.pre[i].busy) {
            if (RS.pre[i].Qj == bus.name)
                RS.now[i].Qj = -1, RS.now[i].Vj = bus.val;
            if (RS.pre[i].Qk == bus.name)
                RS.now[i].Qk = -1, RS.now[i].Vk = bus.val;
        }

        // LSB
        for (int i = 0; i < LSB_SIZE; ++i) if (LSB.pre.q[i].busy) {
            if (LSB.pre.q[i].Qj == bus.name)
                LSB.now.q[i].Qj = -1, LSB.now.q[i].Vj = bus.val;
            if (LSB.pre.q[i].Qk == bus.name)
                LSB.now.q[i].Qk = -1, LSB.now.q[i].Vk = bus.val;
        }
    }

    void run_rob() { // commit
        if (RoB.pre.empty()) return ;
        const ReorderBuffer::node &tmp = RoB.pre.front();
        // cout << std::dec << "@" << RoB.pre.hd + 1 << " " << tmp.dest  << " " << tmp.ready << endl;
        if (!tmp.ready) return ;
        if (tmp.dest == -1) { // halt
            cout << std::dec << (reg.get(10) & 255u) << endl;
            // std::cerr << "clock: " << clk << endl;
            exit(0);
        } else if (tmp.dest < 0) { // branch
            // cout << " !! "<< RoB.pre.hd + 1 << " " <<tmp.npc << " " << tmp.value << endl;
            if ((tmp.value && tmp.npc >= 0) || (!tmp.value && tmp.npc < 0)) { // need to jump: predicted wrong
                pc = tmp.npc >= 0 ? tmp.npc : -tmp.npc;
                BP.calc(-tmp.dest, tmp.value);
                LSB.clean();
                reg.clean();
                RS.clean();
                RoB.clean();
                ALU.clean();
                MC.clean();
                rpc = -1;
                // cout << " @@@ " << lsbus.name << " " << alubus.name << endl;
                lsbus.name = -127;
                alubus.name = -127;
            }
            RoB.now.pop_front();
        } else if (tmp.dest == 32) { // JALR
            // have updated pc in listen();
            RoB.now.pop_front();
        } else if (tmp.dest >= 256) { // store mem
            LSB.setready(tmp.dest - 256);
            if (LSB.pre.q[tmp.dest - 256].Qj == -1 && LSB.pre.q[tmp.dest - 256].Qk == -1) RoB.now.pop_front();
        } else { // register
            reg.vnow[tmp.dest] = tmp.value;
            // if (tmp.dest == 15)
                // cout << std::dec<< tmp.value << " @@ " << RoB.pre.hd + 1 << endl;
            if (reg.getr(tmp.dest) == (RoB.pre.hd + 1) % ROB_SIZE) // only equal can update
                reg.rnow[tmp.dest] = -1;
            RoB.now.pop_front();
        }
    }

    void run_rs() {
        for (int i = 0; i < RS_SIZE; ++i)
            if (RS.pre[i].busy && RS.pre[i].Qj == -1 && RS.pre[i].Qk == -1) {
                ALU.set(RS.pre[i].Vj, RS.pre[i].Vk, RS.pre[i].op, RS.pre[i].name);
                RS.now[i].busy = false;
                return ;
            }
        ALU.stall();
    }

    void ex_loadstore() {
        if (MC.op == SB)
            return mem.setByte(MC.pos, MC.val);
        if (MC.op == SH)
            return mem.setHalfword(MC.pos, MC.val);
        if (MC.op == SW)
            return mem.setWord(MC.pos, MC.val);
        CommonDataBus &bus = lsbus;
        bus.name = MC.val;
        if (MC.op == LB)
            bus.val = mem.getByte(MC.pos);
        if (MC.op == LH)
            bus.val = mem.getHalfword(MC.pos);
        if (MC.op == LW)
            bus.val = mem.getWord(MC.pos);
        if (MC.op == LBU)
            bus.val = mem.getByteUnsigned(MC.pos);
        if (MC.op == LHU)
            bus.val = mem.getHalfwordUnsigned(MC.pos);
        listen(bus);
    }

    void run_lsb() {
        lsbus.name = -127;
        if (!MC.stall && MC.clk == 1)
            ex_loadstore();
        if (!MC.stall && MC.clk > 0) return ; 
        while (!LSB.pre.empty() && !LSB.pre.front().busy)
            LSB.pre.pop_front(), LSB.now.pop_front();
        if (!LSB.pre.empty() && LSB.pre.front().ready && LSB.pre.front().Qj == -1 && LSB.pre.front().Qk == -1) {
            MC.set(LSB.pre.front().op, LSB.pre.front().Vj, LSB.pre.front().Vk);
            LSB.pop();
        } else MC.setstall();
    }

    void issue() {

        // cout << std::hex << pc << " "<< std::dec << rpc << " "<< RoB.full() << " " << RS.full() << " " << LSB.full() << endl;

        if (rpc != -1 || RoB.full() || RS.full() || LSB.full()) return ;

        InsNode ins;
        instructionDecode(mem, pc, ins);
        int name = RoB.getidx();
        u32 delta_pc = 4;

        if (ins.type == NONE) // 0
            return ;

        if (ins.type == HALT) { // 1
            ReorderBuffer::node now;
            now.dest = -1;
            rpc = -2;
            now.ready = 1;
            RoB.push(now);
        }

        if (ins.type == LUI) { // 2
            ReorderBuffer::node now;
            now.dest = ins.rd;
            now.ready = 1;
            now.value = ins.imm;

            reg.setr(ins.rd, name);

            RoB.push(now);
        }

        if (ins.type == AUIPC) { // 3
            ReservationStation::node tmp;
            tmp.busy = 1, tmp.op = ADD, tmp.name = name;
            tmp.Vj = pc;
            tmp.Vk = ins.imm;
            RS.push(tmp);
            ReorderBuffer::node now;
            now.dest = ins.rd;
            reg.setr(ins.rd, name);
            RoB.push(now);
        }

        if (ins.type == JAL) { // 4
            ReorderBuffer::node now;
            now.dest = ins.rd;
            now.ready = 1;
            now.value = pc + 4;
            reg.setr(ins.rd, name);
            RoB.push(now);
            delta_pc = ins.imm;
        }

        if (ins.type == JALR) { // 5
            ReservationStation::node tmp;
            tmp.busy = 1, tmp.op = JALR, tmp.name = -1;
            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) tmp.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) tmp.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vj = lsbus.val;
            else if (rpos == alubus.name) tmp.Vj = alubus.val;
            else tmp.Qj = rpos;
            tmp.Vk = ins.imm;
            // if (pc == 0x1144) {
            //     cout << " @ " << " " << tmp.Vj << " " << tmp.Qj << " " << tmp.Vk << endl;
            //     cout << RoB.pre.hd + 1 << " @ " << RoB.pre.tl << endl;
            //     for (int i = 0; i < 64; ++i)
            //         if (RS.pre[i].busy && RS.pre[i].name == tmp.Qj) cout << " # # # " << endl;
            // }
            RS.push(tmp);
            ReorderBuffer::node now;
            now.dest = ins.rd;
            now.ready = 1;
            now.value = pc + 4;
            reg.setr(ins.rd, name);
            RoB.push(now);
            rpc = name, delta_pc = 0;
        }

        if (6 <= ins.type && ins.type <= 11) { // Branch
            ReservationStation::node tmp;
            tmp.busy = 1, tmp.op = ins.type, tmp.name = name;
            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) tmp.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) tmp.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vj = lsbus.val;
            else if (rpos == alubus.name) tmp.Vj = alubus.val;
            else tmp.Qj = rpos;
            rpos = reg.getr(ins.rs2);
            if (rpos == -1) tmp.Vk = reg.get(ins.rs2);
            else if (RoB.pre.q[rpos].ready) tmp.Vk = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vk = lsbus.val;
            else if (rpos == alubus.name) tmp.Vk = alubus.val;
            else tmp.Qk = rpos;
            RS.push(tmp);

            // cout << name << " ## " << endl;

            ReorderBuffer::node now;
            now.dest = -pc;
            if (BP.predict(pc)) // predict to jump
                delta_pc = ins.imm, now.npc = -(pc + 4);
            else now.npc = pc + ins.imm;
            RoB.push(now);
        }

        if (12 <= ins.type && ins.type <= 20) { // Arith : reg and imm
            ReservationStation::node tmp;
            tmp.busy = 1, tmp.op = ins.type, tmp.name = name;

            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) tmp.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) tmp.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vj = lsbus.val;
            else if (rpos == alubus.name) tmp.Vj = alubus.val;
            else tmp.Qj = rpos;

            tmp.Vk = ins.imm;

            RS.push(tmp);
            ReorderBuffer::node now;
            now.dest = ins.rd;
            reg.setr(ins.rd, name);
            RoB.push(now);
        }

        if (21 <= ins.type && ins.type <= 30) { // Arith : reg and reg
            ReservationStation::node tmp;
            tmp.busy = 1, tmp.op = ins.type, tmp.name = name;

            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) tmp.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) tmp.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vj = lsbus.val;
            else if (rpos == alubus.name) tmp.Vj = alubus.val;
            else tmp.Qj = rpos;

            rpos = reg.getr(ins.rs2);
            if (rpos == -1) tmp.Vk = reg.get(ins.rs2);
            else if (RoB.pre.q[rpos].ready) tmp.Vk = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) tmp.Vk = lsbus.val;
            else if (rpos == alubus.name) tmp.Vk = alubus.val;
            else tmp.Qk = rpos;


            RS.push(tmp);
            ReorderBuffer::node now;
            now.dest = ins.rd;
            reg.setr(ins.rd, name);
            RoB.push(now);
        }

        if (31 <= ins.type && ins.type <= 35) { // Mem : load
            ReservationStation::node rs;
            rs.busy = 1, rs.op = ADD, rs.name = name + 256; // TODO: ATTENTION: one name, the same as lsb.Qj

            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) rs.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) rs.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) rs.Vj = lsbus.val;
            else if (rpos == alubus.name) rs.Vj = alubus.val;
            else rs.Qj = rpos;

            rs.Vk = ins.imm;

            RS.push(rs);
            LoadStoreBuffer::node lsb;
            lsb.op = ins.type;
            lsb.Qj = name + 256;
            lsb.ready = 1;
            lsb.Vk = name;
            LSB.push(lsb);
            ReorderBuffer::node rob;

            rob.dest = ins.rd;
            reg.setr(ins.rd, name);
            RoB.push(rob);
        }

        if (36 <= ins.type && ins.type <= 38) { // Mem : store
            ReservationStation::node rs;
            rs.busy = 1, rs.op = ADD, rs.name = name + 256;

            int rpos = reg.getr(ins.rs1);
            if (rpos == -1) rs.Vj = reg.get(ins.rs1);
            else if (RoB.pre.q[rpos].ready) rs.Vj = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) rs.Vj = lsbus.val;
            else if (rpos == alubus.name) rs.Vj = alubus.val;
            else rs.Qj = rpos;
            rs.Vk = ins.imm;
            RS.push(rs);

            LoadStoreBuffer::node lsb;
            lsb.op = ins.type;
            lsb.Qj = name + 256;

            rpos = reg.getr(ins.rs2);
            if (rpos == -1) lsb.Vk = reg.get(ins.rs2);
            else if (RoB.pre.q[rpos].ready) lsb.Vk = RoB.pre.q[rpos].value;
            else if (rpos == lsbus.name) lsb.Vk = lsbus.val;
            else if (rpos == alubus.name) lsb.Vk = alubus.val;
            else lsb.Qk = rpos;

            int lsb_name = LSB.push(lsb);
            ReorderBuffer::node rob;
            rob.dest = lsb_name + 256;
            rob.ready = 1;
            RoB.push(rob);
        }

        pc += delta_pc;
    }

};

#endif
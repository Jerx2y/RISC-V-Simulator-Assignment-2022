#include "simulator.hpp"

int main() {

    static Simulator cpu;

    cpu.scanmem();

    cpu.run();

    return 0;
}
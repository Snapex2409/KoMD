#include <iostream>

#include "Registry.h"
#include "IO/Logging.h"

static void init();

int main(int argc, char** argv) {
    init();

    Log::finalize();
    return 0;
}

static void init() {
    Registry::instance = std::make_unique<Registry>();
    Log::init();
}

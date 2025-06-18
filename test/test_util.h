#pragma once

#define ASSERT(cond) do { if (!(cond)) { std::exit(__LINE__); } } while(0)

#include <gtest/gtest.h>

#include "tst_frequency.h"
#include "tst_block_frequency.h"
#include "tst_runs.h"
#include "tst_longest_run_of_ones.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

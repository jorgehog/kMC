#pragma once

#include <unittest++/UnitTest++.h>


#define TESTCORE(which, ...) \
    std::cout << "Running test " << UnitTestSuite::GetSuiteName() << "::" << #which << std::endl; \
    \
    testBed::timer.tic(); \
    testBed::test##which(__VA_ARGS__); \
    std::cout << "Done (" << testBed::timer.toc() << " s)" << std::endl; \
    \
    testBed::solver->reset()


#ifdef focusTest
#define isFocusTest(which)  (strcmp(#which, focusTest) == 0)
#else
#define isFocusTest(which) true
#endif

#define TESTWRAPPER(which, ...) TEST(which){ if (isFocusTest(which)) { TESTCORE(which, ##__VA_ARGS__); }}



#define SUITE_IMPL(runner, which) runner.RunTestsIf(UnitTest::Test::GetTestList(), which, UnitTest::True(), 0)

#ifdef focusSuite
#define RUNSUITE(runner, which) \
    (strcmp(which, focusSuite) == 0) \
    ? SUITE_IMPL(runner, which) \
    : 0
#else
#define RUNSUITE SUITE_IMPL
#endif

#pragma once

#include <unittest++/UnitTest++.h>
#include <sstream>
#include <iomanip>


#define TESTCORE(which, ...) \
    std::stringstream s;\
    s << "Running test " << UnitTestSuite::GetSuiteName() << "::" << #which; \
    \
    std::cout << std::left << std::setfill(' ') << std::setw(65); \
    std::cout << s.str(); \
    \
    BADAssBool(Site::boundariesIsInitialized(), "Boundaries are not initialized."); \
    BADAss(Site::_refCount(), !=, 0, "Sites need to be initialized."); \
    \
    testBed::timer.tic(); \
    testBed::test##which(__VA_ARGS__); \
    std::cout << "Done (" << std::setprecision(1) << std::fixed << testBed::timer.toc() << " s)" << std::endl; \
    \
    testBed::solver->reset(); \
    testBed::initSimpleSystemParameters()


#ifdef focusTest
#define isFocusTest(which)  (strcmp(#which, focusTest) == 0)
#else
#define isFocusTest(which) true
#endif

#define TESTWRAPPER(which, ...) \
TEST(which) \
{ \
    if (isFocusTest(which)) \
    { \
        TESTCORE(which, ##__VA_ARGS__); \
     } \
} \
\


#define SUITE_IMPL(runner, which) runner.RunTestsIf(UnitTest::Test::GetTestList(), which, UnitTest::True(), 0)

#ifdef focusSuite
#define RUNSUITE(runner, which) \
    (strcmp(which, focusSuite) == 0) \
    ? SUITE_IMPL(runner, which) \
    : 0
#else
#define RUNSUITE SUITE_IMPL
#endif

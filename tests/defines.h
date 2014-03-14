#pragma once

#define TESTCORE_IMPL(which, ...) \
    std::cout << "Running test " << #which << std::endl; \
    \
    testBed::timer.tic(); \
    testBed::test##which(__VA_ARGS__); \
    std::cout << "Done (" << testBed::timer.toc() << " s)" << std::endl; \
    \
    if (UnitTest::CurrentTest::Results()->GetFailureCount() != 0) \
    { \
        KMCDebugger_DumpFullTrace(); \
    } \
    \
    testBed::solver->reset()

#define isFocusSuite()      (strcmp(UnitTestSuite::GetSuiteName(), focusSuite) == 0)
#define isFocusTest(which)  (strcmp(#which, focusTest) == 0)

#define skipTest static_cast<void>(0)

#define selectByReq(req, which, ...) \
    if (req) { TESTCORE_IMPL(which, ##__VA_ARGS__); } else { skipTest; }

#define requireSuite(which, ...) \
    selectByReq(isFocusSuite(), which, ##__VA_ARGS__)

#define requireTest(which, ...) \
    selectByReq(isFocusTest(which), which, ##__VA_ARGS__)

#define requireSuiteAndTest(which, ...) \
    selectByReq(isFocusSuite() && isFocusTest(which), which, ##__VA_ARGS__)


#ifdef focusSuite
    #ifdef focusTest
        //both focus are defined
        #define TESTCORE(which, ...) requireSuiteAndTest(which, ##__VA_ARGS__)
    #else
        //only suite defined.
        #define TESTCORE(which, ...) requireSuite(which, ##__VA_ARGS__)
    #endif
#else
    #ifdef focusTest
        //only test defined.
        #define TESTCORE(which, ...) requireTest(which, ##__VA_ARGS__)
    #else
        //no focus defined.
        #define TESTCORE(which, ...) TESTCORE_IMPL(which, ##__VA_ARGS__)
    #endif
#endif



//Defined in one line to made unittest++ file line match.
#define TESTWRAPPER(which, ...) TEST(which) {TESTCORE(which, ##__VA_ARGS__);}

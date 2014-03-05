#include "kmcdebugger_class.h"

#ifndef KMC_NO_DEBUG

#include "../../site.h"
#include "../../reactions/reaction.h"

#include <fstream>
#include <sys/time.h>

#include "intrinsicmacros.h"

bool KMCDebugger::enabled = true;

std::vector<std::string> KMCDebugger::reactionTraceBefore;
std::vector<std::string> KMCDebugger::reactionTraceAfter;
std::vector<std::string> KMCDebugger::implicationTrace;
std::vector<double>      KMCDebugger::timerData;

std::set<Site*>          KMCDebugger::affectedUnion;


std::string KMCDebugger::implications;
std::string KMCDebugger::reactionString;

Reaction* KMCDebugger::currentReaction;
Reaction* KMCDebugger::lastCurrentReaction;

uint KMCDebugger::traceCount;
uint KMCDebugger::implicationCount;

std::string KMCDebugger::traceFileName = "";
std::string KMCDebugger::traceFilePath = "";

wall_clock KMCDebugger::timer;


void KMCDebugger::setFilename(const string & filename)
{
    KMCDebugger::traceFileName = filename;
}

void KMCDebugger::setFilepath(const string & filepath)
{
    KMCDebugger::traceFilePath = filepath;
}

void KMCDebugger::setEnabledTo(bool state)
{
    enabled = state;

    if (enabled)
    {
        initialize();
    }

    else
    {
        reset();
    }
}

void KMCDebugger::pushTraces()
{
    if (!enabled) return;

    timerData.push_back(timer.toc());

    stringstream s;

    s << setupAffectedUnion();

    s << addFlagsToImplications();

    string full_str = s.str();

    implications += full_str;

    implicationTrace.push_back(implications);

    reactionTraceBefore.push_back(reactionString);

    if (currentReaction != NULL)
    {
        reactionTraceAfter.push_back(_KMCDebugger_SITE_STR(currentReaction->reactionSite()));
    }

    else
    {
        reactionTraceAfter.push_back("");
    }

    currentReaction = NULL;
    implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;
    implicationCount = 0;
    reactionString = "No Reaction Selected";
    traceCount++;
    affectedUnion.clear();
    timer.tic();

}

void KMCDebugger::pushImplication(Site * site, const char * _pre, const char * _new)
{

    if (!enabled) return;

    using namespace std;

    stringstream s;

    s << "[Implication"; \
    if (currentReaction == NULL) s << " (standalone)";
    s << " " << implicationCount << "]\n";
    s << _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new);

    s << setupAffectedUnion();

    s << "[End of implication";
    if (currentReaction == NULL) s << " (standalone)";
    s << " " << implicationCount << "]\n\n";
    implications += s.str();

    implicationCount++;

}

void KMCDebugger::markPartialStep(const char * msg)
{

    if (!enabled) return;

    using namespace std;

    stringstream s;

    s << "##### " << msg << "  " << "prev. imp.: " << implicationCount << " #####\n\n";

    s << addFlagsToImplications();

    string full_string = s.str();

    implications += full_string;

    implicationCount = 0;

}

void KMCDebugger::setActiveReaction(Reaction *reaction)
{

    if (!enabled) return;

    currentReaction = reaction;
    lastCurrentReaction = reaction;
    reactionString = _KMCDebugger_REACTION_STR();
}

void KMCDebugger::initialize()
{
    if (!enabled) return;

    reset();

    currentReaction = NULL;
    lastCurrentReaction = NULL;

    implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;
    reactionString = _KMCDebugger_INITIAL_REACTION_STR;

    traceCount = 0;
    implicationCount = 0;

    (void)timer.toc();

    timer.tic();
}

std::string KMCDebugger::setupAffectedUnion()
{
    using namespace std;

    stringstream s;

    vector<Site*> intersect;

    bool currentSiteInPrev;

    for (Site * currentSite : Site::affectedSites())
    {

        currentSiteInPrev = false;
        for (Site * previousSite : affectedUnion)
        {
            if (currentSite == previousSite)
            {
                currentSiteInPrev = true;
                break;
            }
        }

        if (!currentSiteInPrev)
        {
            intersect.push_back(currentSite);
        }
    }

    int X, Y, Z;
    for (Site * site : intersect)
    {
        s << "   -" << site->str();

        if (currentReaction != NULL)
        {
            currentReaction->reactionSite()->distanceTo(site, X, Y, Z);
            s << " [" << X << ", " << Y  << ", " << Z << "]";
        }

        s << "\n";

        affectedUnion.insert(site);
    }

    if (intersect.size() == 0)
    {
        return "";
    }

    s << "Total: " << intersect.size() << endl;
    intersect.clear();

    KMCDebugger_Assert(affectedUnion.size(), ==, Site::affectedSites().size());

    return "New affected site(s):\n" + s.str();

}

std::string KMCDebugger::addFlagsToImplications()
{

#ifdef KMC_VERBOSE_DEBUG
    stringstream s, sitess, reacss;

    for (Site * site : affectedUnion)
    {
        bool addSite = false;

        sitess << "   " << site->str() << "\n";
        for (Reaction * r : site->siteReactions())
        {

            reacss << "    .    " << r->str() << " Flags: ";

            bool addReac = false;

            for (int flag : r->updateFlags())
            {
                reacss << flag << " ";
            }

            if (!r->isNotBlocked())
            {
                reacss << " (blocked)";
            }

            addReac = true;
            addSite = true;


            if (addReac)
            {
                sitess << reacss.str() << "\n";
            }

            reacss.str(string());

        }

        if (addSite)
        {
            s << sitess.str();
            sitess.str(string());
        }
    }

    if (s.str().empty())
    {
        return "";
    }


    return "New non-default flags given to affected site reactions:\n\n" + s.str();

#else
    return "";
#endif

}

void KMCDebugger::dumpFullTrace(int line, const char * filename, const string additionalInfo)
{
    if (!enabled) return;


    using namespace std;

    ofstream file;
    stringstream path;


    if (currentReaction != NULL || (traceCount == 0 && !implications.empty()))
    {
        pushTraces();
    }

    if (!traceFilePath.empty())
    {
        path << traceFilePath << "/";
    }

    path << "trace_";

    if (!traceFileName.empty())
    {
        path << traceFileName;
    }

    else
    {
        path << time(NULL);
    }

    path << ".txt";

    file.open(path.str().c_str());
    file << fullTrace(line, filename, additionalInfo);
    file.close();

    cout << "Trace saved to " << path.str() << endl;

    exit(1);

}

void KMCDebugger::dumpPartialTrace(const int &i)
{

    if (!enabled) return;

    using namespace std;

    assert(i >= 0);

    cout << partialTrace(i);
}


string KMCDebugger::fullTrace(int line, const string filename, const string additionalInfo)
{
    using namespace std;

    stringstream s;


    s << "=====================\n TRACING REACTIONS \n=====================\n" << endl;


    for (uint i = 0; i < traceCount; ++i)
    {
        s << partialTrace(i);
    }


    s << "Full trace initiated at line " << line;
    s << "\nin " << filename << "\n";
    s << "=====================\n  TRACE FINISHED  \n=====================\n" << endl;

    if (!additionalInfo.empty())
    {
        s << "\nFinalizing debug info:\n============================\n\n";

        s << additionalInfo;
    }


    s << endl;

    return s.str();

}

string KMCDebugger::partialTrace(const uint &i)
{
    using namespace std;

    stringstream s;

    s << "---[" << i << " / " << traceCount-1 << " ] " << timerData.at(i)*1000 << " ms" << endl;

    s << reactionTraceBefore.at(i) << endl;

    s << implicationTrace.at(i);

#ifdef KMC_VERBOSE_DEBUG
    if (!reactionTraceAfter.at(i).empty())
    {
        s << "\nEnd of reaction look from initial reaction site view: " << endl;
        s << reactionTraceAfter.at(i);
    }
#endif

    s << endl;

    return s.str();

}

void KMCDebugger::reset()
{

    reactionTraceBefore.clear();
    reactionTraceAfter.clear();
    implicationTrace.clear();
    timerData.clear();

    affectedUnion.clear();

}

#endif

#include "defs.hpp"
#include "ioroutines.hpp"
#include "graphmodel.hpp"
#include "hidden_cascade.hpp"

#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices_nan.txt"
#define TRANSACTFILE "table_transacitions.txt"
#define ADIR "/opt/lintula/worktmp/hansen/TAU_epidemic_modelling_for_insiders/announcements"

#define REFDATE 20101230

int main()
{
    // Initialize a metagraph.      
    IoR zed;
    zed.SetAnnouncementDirectory(ADIR);
    zed.SetCompanyDictionaryFile(CDFILE);
    zed.SetNodeDictionaryFile(NDFILE);
    zed.ReadCompanyDictionary();
    zed.ReadNodeDictionary();
    zed.ReadAnnouncements();
    MetaGraph foo(NWDIR);
    MonoGraph* bar = foo.GetGraph(REFDATE);
    double p = 0.2;
    double tp = 0.1;
    double fp = 0.08;
    HiddenCascade cas(bar,p,fp,tp);

    return 0;
}
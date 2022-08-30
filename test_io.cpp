#include "utils/ioroutines.hpp"
#include "data_gen/stat_tester.hpp"

#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices_nan.txt"
#define TRANSACTFILE "table_transacitions.txt"

using namespace util;
using namespace data;

int main()
{
    
    StatTester foo;
    std::cerr << "init.\n";
    
    foo.CreateInsideDayWindows();
    std::cerr << "inside windows calculated\n";
    foo.CreateProfitWindows();
    std::cerr << " profit windows done\n";
    
    //foo.TestHyperG();
    //std::cerr << " hypergeometric tests done\n";
    // foo.PrintHGTest();
    
    // foo.TestGraphIntegrity();
    //std::cerr << " graph tests done \n";
    
    foo.GenerateSmallDataMatrix();
    std::cerr << " Datamatrix generated \n";
    foo.GenerateCSV();
    std::cerr << " Datamatrix sent to files\n";
    // std::cerr << " Printhin isin-map to standard output\n";
    // foo.PrintIsinMap();
    return 0;
    
    /* IoR foo;
    foo.ReadReasonsTable();
    return 0;
    */  
}
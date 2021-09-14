#include "ioroutines.hpp"
#include "stat_tester.hpp"

#define TABLEDIR "/worktmp/hansen/TAU_epidemic_modelling_for_insiders/raw_tables/" 
#define ANFILE "table_announcements.txt"
#define INSFILE "table_insiderships.txt"
#define PRICEFILE "table_prices.txt"
#define TRANSACTFILE "table_transacitions.txt"

int main()
{
    StatTester foo;
    std::cerr << "init.\n";
    foo.CreateInsideDayWindows();
    std::cerr << "inside windows calculated\n";
    foo.CreateProfitWindows();
    std::cerr << " profit windows done\n";
    foo.TestHyperG();
    std::cerr << " hypergeometric tests done\n";
    foo.PrintHGTest();
    foo.GenerateSmallDataMatrix();
    std::cerr << " Datamatrix generated \n";
    foo.GenerateCSV();
    std::cerr << " Datamatrix sent to files\n";
}
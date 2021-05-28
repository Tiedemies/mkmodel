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
    foo.CreateProfitWindows();
}
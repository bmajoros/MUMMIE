#!/bin/sh
rm Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var1.5 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var1.25 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.75 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.5 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.375 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.25 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.20 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.15 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.1 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.075 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.05 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.01 >> Mummie.results.raw
eval-shuffle.pl fc/shuffling/parallel-var0.005 >> Mummie.results.raw

#eval-shuffle.pl fc/shuffling/parallel-PD0.0 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.1 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.2 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.3 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.4 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.5 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.6 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.7 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.8 >> Mummie.results.raw
#eval-shuffle.pl fc/shuffling/parallel-PD0.9 >> Mummie.results.raw




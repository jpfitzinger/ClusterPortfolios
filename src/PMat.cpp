
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix PMat(const NumericMatrix& x) {
size_t n = x.rows();
size_t p = x.cols();
NumericMatrix result(n, n);

for (size_t i = 0; i < n; i++) {
for (size_t j = 0; j < i; j++) {
bool rows_equal = 1;

for (size_t k = 0; k < p; k++) {
if (x(i, k) != x(j, k)) {
rows_equal = 0;
break;
}
}

result(i, j) = rows_equal;
result(j, i) = rows_equal;
}
result(i, i) = 1;
}

return result;
}

// Independent integer check of Corollary 5.5: F^(2d-1) < 12/5 and chi_* < 7/2 at the explicit pairs,
// and chi_* < 7/4 at (7, 3/20). With t = p/q: F = (q^2 + p^2) / (q^2 - p^2) =: a / b.
// Run with Node.js 16 or later: node verify_explicit_balance.mjs
const pairs = [[12,3,25],[12,11,100],[12,23,200],[12,1,8],[12,13,100],[11,13,100],[11,3,25],
  [10,27,200],[10,13,100],[9,7,50],[9,29,200],[9,3,20],[9,31,200],[9,4,25],[8,3,20],[7,31,200],[7,3,20]];
let ok = true;
for (const [d, p, q] of pairs) {
  const a = BigInt(q*q + p*p), b = BigInt(q*q - p*p), n = BigInt(2*d - 1);
  const an = a ** n, bn = b ** n;
  // F^(2d-1) < 12/5  <=>  5 a^n < 12 b^n
  const yOk = 5n * an < 12n * bn;
  // chi_* = 1 + (a^(n+1)/b^(n+1) - 1) / (3 - a^n/b^n) = 1 + (a^(n+1) - b^(n+1)) / (b (3 b^n - a^n))
  const num = a * an - b * bn, den = b * (3n * bn - an);
  const chiOk = den > 0n && 2n * (den + num) < 7n * den;            // chi_* < 7/2
  const sevenOk = !(d === 7 && p === 3) || 4n * (den + num) < 7n * den; // chi_* < 7/4 at (7, 3/20)
  ok &&= yOk && chiOk && sevenOk;
  console.log(`d=${d} t=${p}/${q}: F^(2d-1)<12/5 ${yOk}, chi_*<7/2 ${chiOk}${d === 7 && p === 3 ? `, chi_*<7/4 ${sevenOk}` : ''}`);
}
console.log(ok ? 'PASS' : 'FAIL');
process.exit(ok ? 0 : 1);

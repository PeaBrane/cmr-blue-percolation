/* Independent literal brute force for the Path-1 Bernoulli boxes (p = 1/2, no diamonds).
 * usage: jbrute1 KIND LAB X Y
 *   KIND = side  : box [0,3]x[-1,1] around the horizontal grid segment (0,0)-(3,0); 10 boundary vertices
 *          sideo : same box, origin (1,0) interior, appended as boundary position 10 (own block)
 *          int   : box [-1,1]^2 around the intersection (0,0); 8 boundary vertices
 * Weighted counts v = (Gc, Gm, P, T): sum_e w_e #{z : z_e = 0, e (+)-pivotal for X<->Y}; weights from the family
 * design recomputed here: side: segment edges 120 (Gc/Gm), P-edges 120 (in 1 box), stubs 60 (in 2 boxes), others 0;
 * int: the 4 edges at the centre 120 (Gc), the 8 stubs 120.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int nv, ne, vx[32], vy[32], eu[32], ew[32], w[32][4];
static int nb, bvert[16], blk[16], X, Y;
static int m3(int a) { return ((a % 3) + 3) % 3; }
static int cand(int x, int y) { return m3(x) != 0 && m3(y) != 0; }
static int findv(int x, int y) {
    for (int i = 0; i < nv; i++) if (vx[i] == x && vy[i] == y) return i;
    return -1;
}
static int par[64];
static int fr(int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; }
static int conn(uint32_t z) {
    for (int i = 0; i < nv + 16; i++) par[i] = i;
    for (int k = 0; k < nb; k++) { int a = fr(bvert[k]), b = fr(nv + blk[k]); if (a != b) par[a] = b; }
    for (int k = 0; k < ne; k++) if (z >> k & 1u) { int a = fr(eu[k]), b = fr(ew[k]); if (a != b) par[a] = b; }
    return fr(nv + X) == fr(nv + Y);
}
int main(int argc, char **argv) {
    if (argc != 5) return 1;
    const char *kind = argv[1];
    int side = strncmp(kind, "side", 4) == 0, orig = strcmp(kind, "sideo") == 0;
    int x0 = side ? 0 : -1, x1 = side ? 3 : 1, y0 = -1, y1 = 1;
    for (int y = y0; y <= y1; y++) for (int x = x0; x <= x1; x++) { vx[nv] = x; vy[nv] = y; nv++; }
    for (int i = 0; i < nv; i++) {
        int j = findv(vx[i] + 1, vy[i]);
        if (j >= 0) { eu[ne] = i; ew[ne] = j; ne++; }
        j = findv(vx[i], vy[i] + 1);
        if (j >= 0) { eu[ne] = i; ew[ne] = j; ne++; }
    }
    for (int k = 0; k < ne; k++) {
        int a = eu[k], b = ew[k];
        int ca = cand(vx[a], vy[a]), cb = cand(vx[b], vy[b]);
        int atI = (m3(vx[a]) == 0 && m3(vy[a]) == 0) || (m3(vx[b]) == 0 && m3(vy[b]) == 0);
        memset(w[k], 0, sizeof w[k]);
        if (ca && cb) { w[k][2] = side ? 120 : 0; continue; }
        if (ca || cb) { w[k][3] = side ? 60 : 120; continue; }
        if (side) {  /* grid edge: counted iff on the segment y=0, 0<=x<=3 */
            if (vy[a] == 0 && vy[b] == 0) w[k][atI ? 0 : 1] = 120;
        } else {     /* int box: counted iff incident to the centre */
            if ((vx[a] == 0 && vy[a] == 0) || (vx[b] == 0 && vy[b] == 0)) w[k][0] = 120;
        }
    }
    /* boundary cyclic order */
    for (int x = x0; x <= x1; x++) bvert[nb++] = findv(x, y0);
    for (int y = y0 + 1; y <= y1; y++) bvert[nb++] = findv(x1, y);
    for (int x = x1 - 1; x >= x0; x--) bvert[nb++] = findv(x, y1);
    for (int y = y1 - 1; y > y0; y--) bvert[nb++] = findv(x0, y);
    if (orig) bvert[nb++] = findv(1, 0);
    if ((int)strlen(argv[2]) != nb) { fprintf(stderr, "lab length %zu != %d\n", strlen(argv[2]), nb); return 2; }
    for (int k = 0; k < nb; k++) { char c = argv[2][k]; blk[k] = c <= '9' ? c - '0' : c - 'a' + 10; }
    X = atoi(argv[3]); Y = atoi(argv[4]);
    uint32_t NZ = 1u << ne;
    uint8_t *A = malloc(NZ);
    for (uint32_t z = 0; z < NZ; z++) A[z] = (uint8_t)conn(z);
    uint64_t v[4] = {0, 0, 0, 0};
    for (int k = 0; k < ne; k++) {
        uint32_t b = 1u << k;
        uint64_t cnt = 0;
        for (uint32_t z = 0; z < NZ; z++) if (!(z & b) && A[z | b] && !A[z]) cnt++;
        for (int c = 0; c < 4; c++) v[c] += (uint64_t)w[k][c] * cnt;
    }
    printf("%llu %llu %llu %llu\n", (unsigned long long)v[0], (unsigned long long)v[1], (unsigned long long)v[2],
           (unsigned long long)v[3]);
    return 0;
}

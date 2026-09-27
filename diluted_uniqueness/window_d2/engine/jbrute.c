/* Independent literal brute force for the exact local ratios of the sparse cell (and of a radius-1 box, not used).
 * Written from the definitions (isolated diamond; (+)-pivotal; s-pivotal), sharing no code with engine/lrm.c.
 * usage: jbrute cell LAB X Y     (4x4 cell, candidate {1,2}^2, boundary = 12 cell-perimeter vertices)
 *        jbrute box  LAB X Y     (B_1(0) plus 12 outer neighbours; boundary = outer neighbours; not used here)
 * LAB: 12 chars (0-9,a-f) = outside-connection block of each boundary vertex, in the cyclic boundary order.
 * Output (cell): D N0_P N0_T N0_G N1_P N1_T N1_G   [counts over 2^24 configs; pivot counts x2 per 23-config]
 * Output (box):  Nplus_e1 Nplus_e2 Ns              [same normalisation; Ns summed over the 4 inside plaquettes]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXV 32
static int nv, ne, eu[32], ew[32], ecls[32];
static int bvert[12];
static int blk[12];
static int X, Y;
static int npl, plc[4][4]; /* plaquette corners in cyclic order */
static uint32_t plm[4];    /* plaquette edge masks */
static int vx[MAXV], vy[MAXV];

static int vid(int x, int y) {
    for (int i = 0; i < nv; i++) if (vx[i] == x && vy[i] == y) return i;
    vx[nv] = x; vy[nv] = y; return nv++;
}
static int findv(int x, int y) {
    for (int i = 0; i < nv; i++) if (vx[i] == x && vy[i] == y) return i;
    fprintf(stderr, "no vertex %d %d\n", x, y); exit(3);
}
static void add_edge(int x1, int y1, int x2, int y2, int c) {
    eu[ne] = findv(x1, y1); ew[ne] = findv(x2, y2); ecls[ne] = c; ne++;
}
static int edge_idx(int a, int b) {
    for (int k = 0; k < ne; k++) if ((eu[k] == a && ew[k] == b) || (eu[k] == b && ew[k] == a)) return k;
    fprintf(stderr, "no edge\n"); exit(3);
}
static int par[MAXV + 16];
static int fr(int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; }
/* connectivity of block X and block Y given open edge set z */
static int conn(uint32_t z) {
    for (int i = 0; i < nv + 16; i++) par[i] = i;
    for (int k = 0; k < 12; k++) { int a = fr(bvert[k]), b = fr(nv + blk[k]); if (a != b) par[a] = b; }
    for (int k = 0; k < ne; k++) if (z >> k & 1u) { int a = fr(eu[k]), b = fr(ew[k]); if (a != b) par[a] = b; }
    return fr(nv + X) == fr(nv + Y);
}
static int deg(uint32_t z, int v) {
    int d = 0;
    for (int k = 0; k < ne; k++) if ((z >> k & 1u) && (eu[k] == v || ew[k] == v)) d++;
    return d;
}
/* Definition 2.1 of the diluted-model manuscript: four edges open and exactly one diagonal pair of degree-2 vertices */
static int is_diamond(uint32_t z, int q) {
    if ((z & plm[q]) != plm[q]) return 0;
    int p1 = deg(z, plc[q][0]) == 2 && deg(z, plc[q][2]) == 2;
    int p2 = deg(z, plc[q][1]) == 2 && deg(z, plc[q][3]) == 2;
    return p1 != p2;
}
static void set_plaq(int q, int x, int y) {
    plc[q][0] = findv(x, y); plc[q][1] = findv(x + 1, y); plc[q][2] = findv(x + 1, y + 1); plc[q][3] = findv(x, y + 1);
    plm[q] = 0;
    for (int i = 0; i < 4; i++) plm[q] |= 1u << edge_idx(plc[q][i], plc[q][(i + 1) % 4]);
}

int main(int argc, char **argv) {
    if (argc != 5 || strlen(argv[2]) != 12) { fprintf(stderr, "usage\n"); return 1; }
    int mode_cell = strcmp(argv[1], "cell") == 0;
    for (int k = 0; k < 12; k++) { char c = argv[2][k]; blk[k] = c <= '9' ? c - '0' : c - 'a' + 10; }
    X = atoi(argv[3]); Y = atoi(argv[4]);
    if (mode_cell) {
        for (int y = 0; y < 4; y++) for (int x = 0; x < 4; x++) vid(x, y);
        for (int y = 0; y < 4; y++) for (int x = 0; x < 3; x++) {
            int in1 = (x >= 1 && x <= 2 && y >= 1 && y <= 2), in2 = (x + 1 >= 1 && x + 1 <= 2 && y >= 1 && y <= 2);
            add_edge(x, y, x + 1, y, in1 + in2 == 2 ? 0 : (in1 + in2 == 1 ? 1 : 2));
        }
        for (int x = 0; x < 4; x++) for (int y = 0; y < 3; y++) {
            int in1 = (x >= 1 && x <= 2 && y >= 1 && y <= 2), in2 = (x >= 1 && x <= 2 && y + 1 >= 1 && y + 1 <= 2);
            add_edge(x, y, x, y + 1, in1 + in2 == 2 ? 0 : (in1 + in2 == 1 ? 1 : 2));
        }
        int bx[12] = {0, 1, 2, 3, 3, 3, 3, 2, 1, 0, 0, 0}, by[12] = {0, 0, 0, 0, 1, 2, 3, 3, 3, 3, 2, 1};
        for (int k = 0; k < 12; k++) bvert[k] = findv(bx[k], by[k]);
        npl = 1; set_plaq(0, 1, 1);
    } else {
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) vid(a, b);
        int tx[12] = {2, 2, 2, 1, 0, -1, -2, -2, -2, -1, 0, 1}, ty[12] = {-1, 0, 1, 2, 2, 2, 1, 0, -1, -2, -2, -2};
        for (int k = 0; k < 12; k++) bvert[k] = vid(tx[k], ty[k]);
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
            if (a < 1) add_edge(a, b, a + 1, b, 0);
            if (b < 1) add_edge(a, b, a, b + 1, 0);
        }
        for (int k = 0; k < 12; k++) {  /* edge from outer neighbour to its unique box neighbour */
            int x = tx[k], y = ty[k];
            int ax = x > 1 ? 1 : (x < -1 ? -1 : x), ay = y > 1 ? 1 : (y < -1 ? -1 : y);
            add_edge(ax, ay, x, y, 1);
        }
        npl = 4; set_plaq(0, -1, -1); set_plaq(1, 0, -1); set_plaq(2, -1, 0); set_plaq(3, 0, 0);
    }
    if (ne != 24) { fprintf(stderr, "ne=%d\n", ne); return 2; }
    uint32_t NZ = 1u << ne;
    uint8_t *A = malloc(NZ);
    for (uint32_t z = 0; z < NZ; z++) A[z] = (uint8_t)conn(z);
    if (mode_cell) {
        /* alpha=1: delete the four P-edges if the candidate is an isolated diamond */
        uint8_t *A1 = malloc(NZ);
        uint64_t Dc = 0;
        for (uint32_t z = 0; z < NZ; z++) {
            uint32_t e1 = is_diamond(z, 0) ? (z & ~plm[0]) : z;
            A1[z] = A[e1];
            if (A[z] && !A1[z]) Dc++;
        }
        uint64_t N0[3] = {0, 0, 0}, N1[3] = {0, 0, 0};
        for (int k = 0; k < ne; k++) {
            uint32_t b = 1u << k;
            for (uint32_t z = 0; z < NZ; z++) {
                if (z & b) continue;
                if (A[z | b] && !A[z]) N0[ecls[k]] += 2;
                if (A1[z | b] && !A1[z]) N1[ecls[k]] += 2;
            }
        }
        printf("%llu %llu %llu %llu %llu %llu %llu\n", (unsigned long long)Dc, (unsigned long long)N0[0],
               (unsigned long long)N0[1], (unsigned long long)N0[2], (unsigned long long)N1[0], (unsigned long long)N1[1],
               (unsigned long long)N1[2]);
    } else {
        int k1 = edge_idx(findv(0, 0), findv(1, 0)), k2 = edge_idx(findv(0, 0), findv(0, 1));
        uint64_t Np1 = 0, Np2 = 0, Ns = 0;
        for (uint32_t z = 0; z < NZ; z++) {
            if (!(z >> k1 & 1u) && A[z | 1u << k1] && !A[z]) Np1 += 2;
            if (!(z >> k2 & 1u) && A[z | 1u << k2] && !A[z]) Np2 += 2;
            if (A[z])
                for (int q = 0; q < npl; q++)
                    if (is_diamond(z, q) && !A[z & ~plm[q]]) Ns++;
        }
        printf("%llu %llu %llu\n", (unsigned long long)Np1, (unsigned long long)Np2, (unsigned long long)Ns);
    }
    return 0;
}

/*
 * lrm.c -- exact local-ratio engine for the d = 2 window certificate.
 * Written from the definitions of the sparse model and of the records, independently of the program that first
 * produced the certificate.
 *
 * A box is read from a spec file (geom.py): vertices, NE <= 24 inside edges with class c in {0=Gc,1=Gm,2=P,3=T}
 * and integer weight w, the NB boundary vertices in cyclic order, an optional interior origin, and an optional
 * candidate plaquette (corners in cyclic order; all its corners must be interior vertices of the box).
 *
 * Inside configuration z in {0,1}^NE. Candidate activation alpha in {0,1}. Inside diminished configuration
 *   Ein(z,alpha) = z minus the 4 candidate edges if alpha = 1 and the candidate is an isolated diamond of z
 *   (Def. 2.1: its 4 edges open and EXACTLY ONE diagonal pair consists of two vertices of degree 2);
 *   Ein(z,alpha) = z otherwise.
 * Inside partition part(E) = partition of the positions (boundary vertices, then the origin if present) induced by
 * connectivity in the open edge set E.
 *
 * Boundary condition beta = (pi, X, Y): pi a non-crossing partition of the NB boundary positions (cyclic order);
 * if the origin is present it is an extra singleton block (it has no outside edges). X, Y two distinct blocks
 * (unordered, X < Y as RGS labels); with an origin, Y = the origin block and X ranges over the other blocks.
 * Joint connectivity: X ~ Y in pi v part(E).
 *
 * Counts (all exact integers):
 *   N[alpha][c] = sum_{f in class c} w_f * #{ z in {0,1}^{E minus f} :  X~Y in pi v part(Ein(z+f,alpha))
 *                                                                      and not X~Y in pi v part(Ein(z,alpha)) }
 *   D           = #{ z in {0,1}^E : X~Y in pi v part(Ein(z,0)) and not X~Y in pi v part(Ein(z,1)) }
 * (N: (+)-pivotality of f; D: s-pivotality of the candidate.)
 *
 * Method: phase 1 enumerates all z once and aggregates, per transition key (i_before, i_after) of inside-partition
 * indices, the integer weights of all (z, f, alpha) with that transition (and the D transitions). Phase 2 loops over
 * the boundary partitions pi; for each inside partition i it computes the joint classes of the pi-blocks, and every
 * key adds its weights to the pairs (X,Y) that are joined after but not before. Keys whose after-partition is the
 * before-partition with exactly two blocks merged are processed by the product rule S1 x S2 (joins are
 * associative); all other keys by the literal definition.
 *
 * usage: lrm SPEC MODE OUT NTHREADS
 *   MODE = all        : all non-crossing partitions of the NB boundary positions (restricted growth strings)
 *          list:FILE  : partitions read from FILE (NB bytes per partition, RGS labels), each checked non-crossing
 * output OUT: records (uint32 pid, uint8 X, uint8 Y, uint16 0, uint64 cnt[9]) for EVERY beta, including zero ones;
 *   cnt = D, N[0][Gc,Gm,P,T], N[1][Gc,Gm,P,T].  OUT.parts: the processed partitions (NP bytes each).
 */
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXV 32
#define MAXE 24
#define MAXP 16
#define NCNT 9

static int NV, NE, NB, NP, ORIG = -1, HASC = 0;
static int ev[MAXE][2], ecls[MAXE], ew[MAXE];
static int posv[MAXP]; /* vertex of each position */
static int cand[4];
static uint32_t inc[MAXV]; /* incident edge mask per vertex */
static uint32_t Pmask = 0;

static void die(const char *m) {
    fprintf(stderr, "error: %s\n", m);
    exit(2);
}

static void read_spec(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) die("spec");
    char line[512];
    int nv = 0, ne = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char key[16];
        if (sscanf(line, "%15s", key) != 1) continue;
        if (!strcmp(key, "NV")) sscanf(line + 2, "%d", &NV);
        else if (!strcmp(key, "V")) nv++;
        else if (!strcmp(key, "NE")) sscanf(line + 2, "%d", &NE);
        else if (!strcmp(key, "E")) {
            if (ne >= MAXE) die("too many edges");
            sscanf(line + 1, "%d %d %d %d", &ev[ne][0], &ev[ne][1], &ecls[ne], &ew[ne]);
            ne++;
        } else if (!strcmp(key, "NB")) sscanf(line + 2, "%d", &NB);
        else if (!strcmp(key, "B")) {
            char *p = line + 1;
            for (int k = 0; k < NB; k++) {
                int n;
                if (sscanf(p, "%d%n", &posv[k], &n) != 1) die("B line");
                p += n;
            }
        } else if (!strcmp(key, "ORIGIN")) sscanf(line + 6, "%d", &ORIG);
        else if (!strcmp(key, "CAND")) {
            int a = -1, b = -1, c = -1, d = -1;
            sscanf(line + 4, "%d %d %d %d", &a, &b, &c, &d);
            if (a >= 0) { HASC = 1; cand[0] = a; cand[1] = b; cand[2] = c; cand[3] = d; }
        }
    }
    fclose(f);
    if (nv != NV || ne != NE || NV > MAXV || NB + 1 > MAXP) die("spec sizes");
    NP = NB;
    if (ORIG >= 0) posv[NP++] = ORIG;
    for (int e = 0; e < NE; e++) {
        inc[ev[e][0]] |= 1u << e;
        inc[ev[e][1]] |= 1u << e;
    }
    if (HASC) {
        for (int k = 0; k < 4; k++) {
            int a = cand[k], b = cand[(k + 1) % 4], found = 0;
            for (int e = 0; e < NE; e++)
                if ((ev[e][0] == a && ev[e][1] == b) || (ev[e][0] == b && ev[e][1] == a)) { Pmask |= 1u << e; found = 1; }
            if (!found) die("candidate edge missing");
        }
        /* all candidate corners must be interior (not positions), so their degrees are inside quantities */
        for (int k = 0; k < 4; k++)
            for (int q = 0; q < NP; q++)
                if (posv[q] == cand[k]) die("candidate corner on the boundary");
        /* and every edge at a candidate corner must be an inside edge: degree 4 inside */
        for (int k = 0; k < 4; k++)
            if (__builtin_popcount(inc[cand[k]]) != 4) die("candidate corner not interior");
    }
}

/* Def. 2.1, literally */
static int is_diamond(uint32_t z) {
    if ((z & Pmask) != Pmask) return 0;
    int d0 = __builtin_popcount(z & inc[cand[0]]), d1 = __builtin_popcount(z & inc[cand[1]]);
    int d2 = __builtin_popcount(z & inc[cand[2]]), d3 = __builtin_popcount(z & inc[cand[3]]);
    int pairA = (d0 == 2 && d2 == 2), pairB = (d1 == 2 && d3 == 2);
    return pairA != pairB;
}
static uint32_t ein(uint32_t z, int alpha) { return (alpha && HASC && is_diamond(z)) ? (z & ~Pmask) : z; }

static int fr(int *p, int x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}
/* canonical code of the partition of the NP positions: RGS labels, 4 bits each */
static uint64_t part_code(uint32_t E) {
    int p[MAXV];
    for (int i = 0; i < NV; i++) p[i] = i;
    for (int e = 0; e < NE; e++)
        if (E >> e & 1u) {
            int a = fr(p, ev[e][0]), b = fr(p, ev[e][1]);
            if (a != b) p[a] = b;
        }
    int lab[MAXV];
    for (int i = 0; i < NV; i++) lab[i] = -1;
    int nl = 0;
    uint64_t code = 0;
    for (int k = 0; k < NP; k++) {
        int r = fr(p, posv[k]);
        if (lab[r] < 0) lab[r] = nl++;
        code |= (uint64_t)lab[r] << (4 * k);
    }
    return code;
}

/* ---- hash: partition code -> index ---- */
static uint64_t *pkey;
static int32_t *pval;
static uint64_t *pcode;
static uint32_t npart = 0, PBITS = 20;
static uint32_t hmix(uint64_t k, uint32_t bits) { return (uint32_t)((k * 0x9E3779B97F4A7C15ull) >> (64 - bits)); }
static uint32_t pidx(uint64_t c) {
    uint32_t s = hmix(c, PBITS), m = (1u << PBITS) - 1;
    while (1) {
        if (pval[s] < 0) {
            if (npart >= (1u << (PBITS - 1))) die("partition table full");
            pkey[s] = c; pval[s] = (int32_t)npart; pcode[npart] = c;
            return npart++;
        }
        if (pkey[s] == c) return (uint32_t)pval[s];
        s = (s + 1) & m;
    }
}

/* ---- hash: transition key (ib, ia) -> weight counters ---- */
typedef struct {
    uint64_t key;
    uint64_t c[NCNT];
} Tr;
static Tr *tr;
static uint32_t TBITS = 22, ntr = 0;
static void tr_add(uint32_t ib, uint32_t ia, int slot, uint64_t w) {
    uint64_t key = ((uint64_t)ib << 32) | ia;
    uint32_t s = hmix(key, TBITS), m = (1u << TBITS) - 1;
    while (1) {
        if (tr[s].key == ~0ull) {
            if (ntr >= (1u << (TBITS - 1))) die("transition table full");
            tr[s].key = key;
            memset(tr[s].c, 0, sizeof tr[s].c);
            tr[s].c[slot] = w;
            ntr++;
            return;
        }
        if (tr[s].key == key) { tr[s].c[slot] += w; return; }
        s = (s + 1) & m;
    }
}

/* compact key list for phase 2 */
typedef struct {
    uint32_t ib, ia;
    int8_t k1, k2; /* single merge: positions of the two merged blocks; k1 < 0 for general keys */
    uint32_t c[NCNT];
} Key;
static int keycmp(const void *a, const void *b) {
    const Key *x = (const Key *)a, *y = (const Key *)b;
    if (x->ib != y->ib) return x->ib < y->ib ? -1 : 1;
    if (x->ia != y->ia) return x->ia < y->ia ? -1 : 1;
    return 0;
}
static Key *keys;
static uint32_t nkeys;

static void decode(uint64_t code, int *lab) {
    for (int k = 0; k < NP; k++) lab[k] = (int)((code >> (4 * k)) & 15);
}

static void classify_keys(void) {
    keys = malloc(sizeof(Key) * (size_t)ntr);
    nkeys = 0;
    uint32_t nsingle = 0;
    for (uint64_t s = 0; s < (1ull << TBITS); s++) {
        if (tr[s].key == ~0ull) continue;
        Key *K = &keys[nkeys++];
        K->ib = (uint32_t)(tr[s].key >> 32);
        K->ia = (uint32_t)(tr[s].key & 0xFFFFFFFFu);
        for (int t = 0; t < NCNT; t++) {
            if (tr[s].c[t] > 0xFFFFFFFFull) die("key counter overflow");
            K->c[t] = (uint32_t)tr[s].c[t];
        }
        K->k1 = K->k2 = -1;
        int lb[MAXP], la[MAXP];
        decode(pcode[K->ib], lb);
        decode(pcode[K->ia], la);
        int nbb = 0, nba = 0;
        for (int k = 0; k < NP; k++) {
            if (lb[k] + 1 > nbb) nbb = lb[k] + 1;
            if (la[k] + 1 > nba) nba = la[k] + 1;
        }
        /* coarsening check: every before-block inside one after-block */
        int map[MAXP], firstpos[MAXP], ok = 1;
        for (int b = 0; b < nbb; b++) { map[b] = -1; firstpos[b] = -1; }
        for (int k = 0; k < NP; k++) {
            if (firstpos[lb[k]] < 0) firstpos[lb[k]] = k;
            if (map[lb[k]] < 0) map[lb[k]] = la[k];
            else if (map[lb[k]] != la[k]) ok = 0;
        }
        if (ok && nba == nbb - 1) {
            int b1 = -1, b2 = -1;
            for (int x = 0; x < nbb && b1 < 0; x++)
                for (int y = x + 1; y < nbb; y++)
                    if (map[x] == map[y]) { b1 = x; b2 = y; break; }
            K->k1 = (int8_t)firstpos[b1];
            K->k2 = (int8_t)firstpos[b2];
            nsingle++;
        }
    }
    qsort(keys, nkeys, sizeof(Key), keycmp);
    free(tr);
    fprintf(stderr, "keys: %u (single merges %u, general %u)\n", nkeys, nsingle, nkeys - nsingle);
}

/* ---- non-crossing partitions ---- */
static uint8_t (*plist)[MAXP];
static uint64_t nplist = 0, cap_plist = 0;
static int is_nc(const uint8_t *lab, int n) {
    /* stack test: scanning left to right, a returning block must be the most recently opened open block */
    int last[MAXP], stack[MAXP], sp = 0;
    for (int b = 0; b < MAXP; b++) last[b] = -1;
    for (int k = 0; k < n; k++) last[lab[k]] = k;
    int seen[MAXP];
    memset(seen, 0, sizeof seen);
    for (int k = 0; k < n; k++) {
        int l = lab[k];
        if (!seen[l]) {
            seen[l] = 1;
            if (last[l] > k) stack[sp++] = l;
        } else {
            if (sp == 0 || stack[sp - 1] != l) return 0;
            if (last[l] == k) sp--;
        }
    }
    return 1;
}
static int is_nc_brute(const uint8_t *lab, int n) {
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            for (int c = b + 1; c < n; c++)
                for (int d = c + 1; d < n; d++)
                    if (lab[a] == lab[c] && lab[b] == lab[d] && lab[a] != lab[b]) return 0;
    return 1;
}
static uint8_t rgs[MAXP];
static uint64_t nrgs = 0;
static void push_part(const uint8_t *lab) {
    if (nplist == cap_plist) {
        cap_plist = cap_plist ? 2 * cap_plist : 1024;
        plist = realloc(plist, cap_plist * MAXP);
    }
    memcpy(plist[nplist++], lab, MAXP);
}
static void enum_rgs(int k, int mx) {
    if (k == NB) {
        nrgs++;
        int a = is_nc(rgs, NB);
        if (NB <= 12 && a != is_nc_brute(rgs, NB)) die("nc test mismatch");
        if (a) push_part(rgs);
        return;
    }
    for (int l = 0; l <= mx + 1; l++) {
        rgs[k] = (uint8_t)l;
        enum_rgs(k + 1, l > mx ? l : mx);
    }
}

/* ---- phase 2 ---- */
static FILE *outf;
static pthread_mutex_t omx = PTHREAD_MUTEX_INITIALIZER;
typedef struct { int tid, nth; uint64_t nrec, nnz; } Job;

static void *worker(void *arg) {
    Job *jb = (Job *)arg;
    uint16_t *jm = malloc(sizeof(uint16_t) * (size_t)npart * MAXP);
    typedef struct { uint32_t pid; uint8_t x, y; uint16_t pad; uint64_t c[NCNT]; } Rec;
    Rec *buf = malloc(sizeof(Rec) * 256);
    jb->nrec = jb->nnz = 0;
    for (uint64_t pid = (uint64_t)jb->tid; pid < nplist; pid += (uint64_t)jb->nth) {
        const uint8_t *pl = plist[pid];
        int pb[MAXP], nb = 0;
        for (int k = 0; k < NB; k++) { pb[k] = pl[k]; if (pl[k] + 1 > nb) nb = pl[k] + 1; }
        int nbt = nb;
        if (ORIG >= 0) pb[NB] = nbt++; /* origin: own block */
        /* joint classes of the pi-blocks for every inside partition */
        for (uint32_t i = 0; i < npart; i++) {
            int lab[MAXP];
            decode(pcode[i], lab);
            int p[MAXP], firstb[MAXP];
            for (int b = 0; b < nbt; b++) p[b] = b;
            for (int l = 0; l < NP; l++) firstb[l] = -1;
            for (int k = 0; k < NP; k++) {
                int l = lab[k];
                if (firstb[l] < 0) firstb[l] = pb[k];
                else {
                    int a = fr(p, pb[k]), b = fr(p, firstb[l]);
                    if (a != b) p[a] = b;
                }
            }
            uint16_t rm[MAXP];
            memset(rm, 0, sizeof rm);
            for (int b = 0; b < nbt; b++) rm[fr(p, b)] |= (uint16_t)(1u << b);
            for (int b = 0; b < nbt; b++) jm[(size_t)i * MAXP + b] = rm[fr(p, b)];
        }
        uint64_t M[MAXP][MAXP][NCNT];
        memset(M, 0, sizeof M);
        for (uint32_t q = 0; q < nkeys; q++) {
            const Key *K = &keys[q];
            const uint16_t *mb = &jm[(size_t)K->ib * MAXP];
            if (K->k1 >= 0) {
                uint16_t S1 = mb[pb[K->k1]], S2 = mb[pb[K->k2]];
                if (S1 == S2) continue;
                for (uint16_t a = S1; a; a &= (uint16_t)(a - 1)) {
                    int x = __builtin_ctz(a);
                    for (uint16_t b = S2; b; b &= (uint16_t)(b - 1)) {
                        int y = __builtin_ctz(b);
                        int lo = x < y ? x : y, hi = x < y ? y : x;
                        for (int t = 0; t < NCNT; t++) M[lo][hi][t] += K->c[t];
                    }
                }
            } else {
                const uint16_t *ma = &jm[(size_t)K->ia * MAXP];
                for (int x = 0; x < nbt; x++) {
                    uint16_t nw = (uint16_t)(ma[x] & ~mb[x]);
                    nw &= (uint16_t)~((2u << x) - 1u); /* y > x */
                    for (; nw; nw &= (uint16_t)(nw - 1)) {
                        int y = __builtin_ctz(nw);
                        for (int t = 0; t < NCNT; t++) M[x][y][t] += K->c[t];
                    }
                }
            }
        }
        int nr = 0;
        for (int x = 0; x < nbt; x++)
            for (int y = x + 1; y < nbt; y++) {
                if (ORIG >= 0 && y != nbt - 1) continue;
                Rec *r = &buf[nr++];
                r->pid = (uint32_t)pid; r->x = (uint8_t)x; r->y = (uint8_t)y; r->pad = 0;
                int nz = 0;
                for (int t = 0; t < NCNT; t++) { r->c[t] = M[x][y][t]; nz |= (M[x][y][t] != 0); }
                jb->nnz += (uint64_t)nz;
            }
        pthread_mutex_lock(&omx);
        fwrite(buf, sizeof(Rec), (size_t)nr, outf);
        pthread_mutex_unlock(&omx);
        jb->nrec += (uint64_t)nr;
    }
    free(jm);
    free(buf);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: lrm SPEC MODE OUT NTHREADS [TBITS]\n"); return 1; }
    read_spec(argv[1]);
    if (argc > 5) TBITS = (uint32_t)atoi(argv[5]);
    fprintf(stderr, "spec %s: NV=%d NE=%d NB=%d NP=%d origin=%d candidate=%d\n", argv[1], NV, NE, NB, NP, ORIG, HASC);
    pkey = malloc(8ull << PBITS);
    pval = malloc(4ull << PBITS);
    pcode = malloc(8ull << PBITS);
    memset(pval, 0xFF, 4ull << PBITS);
    tr = malloc(sizeof(Tr) << TBITS);
    for (uint64_t s = 0; s < (1ull << TBITS); s++) tr[s].key = ~0ull;

    /* ---- phase 1 ---- */
    uint64_t NZ = 1ull << NE, ndiam = 0, ncyc = 0;
    int na = HASC ? 2 : 1;
    for (uint64_t z = 0; z < NZ; z++) {
        uint32_t zz = (uint32_t)z;
        for (int alpha = 0; alpha < na; alpha++) {
            uint32_t ib = pidx(part_code(ein(zz, alpha)));
            for (int f = 0; f < NE; f++) {
                if (zz >> f & 1u) continue;
                uint32_t ia = pidx(part_code(ein(zz | (1u << f), alpha)));
                if (ia != ib) tr_add(ib, ia, 1 + 4 * alpha + ecls[f], (uint64_t)ew[f]);
            }
        }
        if (HASC) {
            if ((zz & Pmask) == Pmask) {
                int d0 = __builtin_popcount(zz & inc[cand[0]]), d1 = __builtin_popcount(zz & inc[cand[1]]);
                int d2 = __builtin_popcount(zz & inc[cand[2]]), d3 = __builtin_popcount(zz & inc[cand[3]]);
                if (d0 == 2 && d1 == 2 && d2 == 2 && d3 == 2) ncyc++;
            }
            if (is_diamond(zz)) {
                ndiam++;
                uint32_t ideleted = pidx(part_code(zz & ~Pmask)), ipresent = pidx(part_code(zz));
                if (ipresent != ideleted) tr_add(ideleted, ipresent, 0, 1);
            }
        }
    }
    fprintf(stderr, "phase 1: %llu configurations, inside partitions K=%u, transition keys %u, isolated diamonds %llu, "
            "isolated 4-cycles %llu\n", (unsigned long long)NZ, npart, ntr, (unsigned long long)ndiam,
            (unsigned long long)ncyc);
    classify_keys();

    /* ---- partitions ---- */
    if (!strcmp(argv[2], "all")) {
        enum_rgs(0, -1);
        fprintf(stderr, "RGS of length %d: %llu, non-crossing: %llu\n", NB, (unsigned long long)nrgs,
                (unsigned long long)nplist);
    } else if (!strncmp(argv[2], "list:", 5)) {
        FILE *f = fopen(argv[2] + 5, "rb");
        if (!f) die("list file");
        uint8_t lab[MAXP];
        memset(lab, 0, sizeof lab);
        while (fread(lab, 1, (size_t)NB, f) == (size_t)NB) {
            /* must be a restricted growth string and non-crossing */
            int mx = -1;
            for (int k = 0; k < NB; k++) {
                if (lab[k] > mx + 1) die("list entry not an RGS");
                if (lab[k] > mx) mx = lab[k];
            }
            if (!is_nc(lab, NB) || !is_nc_brute(lab, NB)) die("list entry not non-crossing");
            push_part(lab);
        }
        fclose(f);
        fprintf(stderr, "read %llu partitions from %s\n", (unsigned long long)nplist, argv[2] + 5);
    } else die("mode");
    char fn[1024];
    snprintf(fn, sizeof fn, "%s.parts", argv[3]);
    FILE *fp = fopen(fn, "wb");
    for (uint64_t i = 0; i < nplist; i++) {
        uint8_t lab[MAXP];
        memcpy(lab, plist[i], MAXP);
        if (ORIG >= 0) {
            int mx = -1;
            for (int k = 0; k < NB; k++) if (lab[k] > mx) mx = lab[k];
            lab[NB] = (uint8_t)(mx + 1);
        }
        fwrite(lab, 1, (size_t)NP, fp);
    }
    fclose(fp);

    /* ---- phase 2 ---- */
    outf = fopen(argv[3], "wb");
    int nth = atoi(argv[4]);
    pthread_t th[64];
    Job jobs[64];
    for (int t = 0; t < nth; t++) {
        jobs[t].tid = t;
        jobs[t].nth = nth;
        pthread_create(&th[t], NULL, worker, &jobs[t]);
    }
    uint64_t nrec = 0, nnz = 0;
    for (int t = 0; t < nth; t++) {
        pthread_join(th[t], NULL);
        nrec += jobs[t].nrec;
        nnz += jobs[t].nnz;
    }
    fclose(outf);
    printf("%s: partitions %llu, records %llu (nonzero %llu), K=%u, keys=%u\n", argv[1], (unsigned long long)nplist,
           (unsigned long long)nrec, (unsigned long long)nnz, npart, nkeys);
    return 0;
}

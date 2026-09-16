/*
 * frequent_pattern_dfd.cpp
 * =========================================================================
 * Implementation of Section 4 ("Frequent pattern discovery") of:
 *
 *   Tang, Yiu, Mouratidis, Zhang, Wang.
 *   "On discovering motifs and frequent patterns in spatial trajectories
 *    with discrete Frechet distance." GeoInformatica 26, 29-66 (2022).
 *
 * Implemented algorithms
 * -----------------------
 *   1. BruteFP   -- Algorithm 4  (brute-force baseline, Sec. 4)
 *   2. BFP       -- Algorithm 5  (bounding-based pruning using the Boolean
 *                    feasibility matrix beta, Sec. 4.1: cell-based,
 *                    cross-based and band-based pruning rules, Sec 4.1.2,
 *                    accelerated with the O(1)-amortized sliding-window
 *                    counting trick of Sec. 4.1.3)
 *   3. BFP*      -- space-efficient variant (Sec. 4.2): does not store the
 *                    O(n^2) matrices; computes beta/Brow/Bcol/Bband and the
 *                    exact DFD on the fly, trading a bit of speed for much
 *                    lower memory (useful for very long trajectories).
 *
 * Problem solved (Problem 2 in the paper)
 * ----------------------------------------
 * Given a trajectory S, a subtrajectory length xi, a distance threshold
 * dthres and a frequency threshold tthres, find every "frequent pattern"
 * (S_{i,i+xi-1}, Gamma_i) such that:
 *   (1) every member S_{j,j+xi-1} in Gamma_i satisfies
 *          dF(i,i+xi-1, j,j+xi-1) <= dthres          (discrete Frechet distance)
 *   (2) the two subtrajectories' timestamp/index intervals do not overlap
 *   (3) |Gamma_i| >= tthres
 *
 * Input format
 * ------------
 * The input file is a single (already concatenated) trajectory given as
 * one point per line:
 *
 *      x  y  trajectory_id
 *
 * As in the paper's multi-trajectory extension (Section 5, Fig. 13), when
 * several original trajectories are concatenated into one long sequence we
 * must forbid any candidate subtrajectory of length xi from straddling a
 * boundary between two different trajectory_id values -- such windows are
 * "invalid" and are never used as the start of a subtrajectory (a
 * representative S_{i,ie} nor a candidate S_{j,je}).
 *
 * Ground distance dG(i,j)
 * ------------------------
 * Implemented as plain Euclidean distance between (x,y) pairs. If your
 * x,y are (longitude, latitude) and you want the great-circle distance
 * used in the paper, replace groundDistance() with a haversine formula.
 *
 * Build
 * -----
 *   g++ -O2 -std=c++17 -o frequent_pattern_dfd frequent_pattern_dfd.cpp
 *
 * Usage
 * -----
 *   ./frequent_pattern_dfd <input_file> <xi> <dthres> <tthres> [algo] [-v]
 *
 *      xi      : subtrajectory length (number of points), integer > 1
 *      dthres  : distance threshold (same unit as x,y), double
 *      tthres  : frequency threshold (min. number of matches), integer >= 1
 *      algo    : "brute" | "bfp" | "bfpstar"   (default: bfp)
 *      -v      : verbose statistics (candidates pruned by each rule, etc.)
 *
 * Output
 * ------
 * For every discovered frequent pattern, prints the representative
 * subtrajectory's [start,end] point-index range together with the list of
 * matching subtrajectories' [start,end] ranges.
 * =========================================================================
 */

#include <bits/stdc++.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

// -------------------------------------------------------------------------
// Basic data structures
// -------------------------------------------------------------------------

struct Point {
    double x, y;
    long long tid;
};

static inline double groundDistance(const Point &a, const Point &b) {
    // dG(i,j): Euclidean distance. Swap for haversine if x,y are lon/lat.
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

struct FrequentPattern {
    int repStart;                  // representative subtrajectory start index i
    int repEnd;                    // representative subtrajectory end index  ie
    int num_members; // (start,end) of each matching subtrajectory
};


// -------------------------------------------------------------------------
// Input loading
// -------------------------------------------------------------------------

static vector<Point> loadTrajectory(const string &path) {
    ifstream in(path);
    if (!in) {
        cerr << "ERROR: cannot open input file: " << path << "\n";
        exit(1);
    }
    vector<Point> pts;
    double x, y;
    long long tid;
    while (in >> x >> y >> tid) {
        pts.push_back({x, y, tid});
    }
    if (pts.empty()) {
        cerr << "ERROR: no points read from " << path << "\n";
        exit(1);
    }
    return pts;
}

// valid[i] == true  iff  points[i .. i+xi-1] all belong to the same
// trajectory_id, i.e. the window does not cross a trajectory boundary
// (see Fig. 13 / Section 5 of the paper: "invalid subtrajectories").
static vector<char> computeValidStarts(const vector<Point> &pts, int xi) {
    int n = (int)pts.size();
    vector<char> valid(n, 0);
    if (n < xi) return valid;

    // runEnd[i] = last index of the maximal run of identical trajectory_id
    // that contains position i.
    vector<int> runEnd(n);
    runEnd[n - 1] = n - 1;
    for (int i = n - 2; i >= 0; --i) {
        runEnd[i] = (pts[i].tid == pts[i + 1].tid) ? runEnd[i + 1] : i;
    }
    for (int i = 0; i + xi - 1 < n; ++i) {
        valid[i] = (runEnd[i] >= i + xi - 1) ? 1 : 0;
    }
    return valid;
}

static inline bool intervalsOverlap(int i, int ie, int j, int je) {
    return !(je < i || ie < j);
}

// -------------------------------------------------------------------------
// Exact discrete Frechet distance between two *fixed-length* subtrajectories
// S[i..ie] and S[j..je]  (standard O(xi^2) dynamic programming, Section 2.1)
// -------------------------------------------------------------------------

static double exactDFD(const vector<Point> &pts, int i, int ie, int j, int je) {
    int li = ie - i + 1;
    int lj = je - j + 1;
    // rolling 2-row DP to save memory; dp[b] = ca(a,b)
    vector<double> prev(lj), cur(lj);
    for (int a = 0; a < li; ++a) {
        for (int b = 0; b < lj; ++b) {
            double d = groundDistance(pts[i + a], pts[j + b]);
            double val;
            if (a == 0 && b == 0) {
                val = d;
            } else if (a == 0) {
                val = std::max(d, cur[b - 1]);
            } else if (b == 0) {
                val = std::max(d, prev[b]);
            } else {
                val = std::max(d, std::min({prev[b], cur[b - 1], prev[b - 1]}));
            }
            cur[b] = val;
        }
        std::swap(prev, cur);
    }
    return prev[lj - 1];
}

// =========================================================================
// 1) Algorithm 4 : BruteFP  (baseline brute-force solution, Section 4)
// =========================================================================
//
// For every valid representative start i, try every valid, non-overlapping
// candidate start j and compute the exact DFD directly. O(n^2 * xi^2).

static vector<FrequentPattern> bruteFP(const vector<Point> &pts, int xi,
                                        double dthres, int tthres,
                                        bool verbose = false) {
    int n = (int)pts.size();
    vector<char> valid = computeValidStarts(pts, xi);
    vector<FrequentPattern> results;
    long long dfdCalls = 0;

    for (int i = 0; i + xi - 1 < n; ++i) {
        if (!valid[i]) continue;
        int ie = i + xi - 1;
        FrequentPattern pat;
        pat.repStart = i;
        pat.repEnd = ie;

        for (int j = 0; j + xi - 1 < n; ++j) {
            if (!valid[j]) continue;
            int je = j + xi - 1;
            if (intervalsOverlap(i, ie, j, je)) continue;

            double d = exactDFD(pts, i, ie, j, je);
            ++dfdCalls;
            if (d <= dthres) {
                pat.num_members++;
            }
        }
        if ((int)pat.num_members >= tthres) {
            
            results.push_back(std::move(pat));
        }
    }

    if (verbose) {
        cerr << "[BruteFP] exact DFD computations: " << dfdCalls << "\n";
    }
    return results;
}

// =========================================================================
// 2) Algorithm 5 : BFP (bounding-based pruning, Section 4.1)
// =========================================================================
//
// Steps (mirroring Sec. 4.1.1 - 4.1.4):
//   (a) Build the Boolean feasibility matrix beta(i,j) = [dG(i,j) <= dthres]
//       (Section 4.1.1, Observation 5).
//   (b) Build Brow(i,j), Bcol(i,j)      (cross-based rule, Section 4.1.2)
//       Brow(i,j) = OR_{i' in [i,i+xi-1]} beta(i', j+1)
//       Bcol(i,j) = OR_{j' in [j,j+xi-1]} beta(i+1, j')
//   (c) Build Brow_band(i,j), Bcol_band(i,j) (band-based rule, Sec. 4.1.2)
//       Brow_band(i,j) = AND_{j' in [j,j+xi-1]} Brow(i,j')
//       Bcol_band(i,j) = AND_{i' in [i,i+xi-1]} Bcol(i',j)
//   All of the above are built with the O(1)-amortized sliding-window /
//   prefix-sum technique of Section 4.1.3, so the whole construction is
//   O(n^2) total (instead of O(n^2 * xi) or O(n^2 * xi^2) naively).
//
// (d) Main loop: for each (i,j) pair
//       - cell-based pruning  : if !beta(i,j)                     -> skip
//       - cross-based pruning : if !(Brow(i,j) && Bcol(i,j))      -> skip
//       - band-based pruning  : if !(Brow_band && Bcol_band)      -> skip
//       - otherwise compute the exact DFD and test against dthres.
//
// Memory: O(n^2) bits per matrix, 5 matrices -> ~5*n^2 bits total.
// For very large n prefer BFP* (Section 4.2 / bfpStar() below).

struct BoolMatrix {
    int n;
    vector<vector<bool>> m;
    BoolMatrix() : n(0) {}
    void init(int n_) { n = n_; m.assign(n, vector<bool>(n, false)); }
    inline bool get(int i, int j) const {
        if (i < 0 || i >= n || j < 0 || j >= n) return false;
        return m[i][j];
    }
    inline void set(int i, int j, bool v) { m[i][j] = v; }
};

static vector<FrequentPattern> bfp(const vector<Point> &pts, int xi,
                                    double dthres, int tthres,
                                    bool verbose = false) {
    int n = (int)pts.size();
    vector<char> valid = computeValidStarts(pts, xi);

    // ---- (a) beta(i,j) --------------------------------------------------
    BoolMatrix beta; beta.init(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            beta.set(i, j, groundDistance(pts[i], pts[j]) <= dthres);

    // ---- (b1) Brow(i,j) via sliding row-window (Section 4.1.3 idea) ----
    // Brow(i,j) = OR_{i' in [i,i+xi-1]} beta(i', j+1)
    // RowTrueCount[c] = number of True among beta[i .. i+xi-1][c]
    BoolMatrix Brow; Brow.init(n);
    {
        vector<int> RowTrueCount(n, 0);
        int lastI = n - xi; // last valid i for a length-xi row band
        if (lastI >= 0) {
            for (int c = 0; c < n; ++c) {
                int cnt = 0;
                for (int r = 0; r < xi; ++r) cnt += beta.get(r, c) ? 1 : 0;
                RowTrueCount[c] = cnt;
            }
            for (int i = 0; i <= lastI; ++i) {
                for (int j = 0; j < n; ++j) {
                    bool v = (j + 1 < n) && (RowTrueCount[j + 1] > 0);
                    Brow.set(i, j, v);
                }
                if (i < lastI) {
                    // slide the row-window from [i,i+xi-1] to [i+1,i+xi]
                    for (int c = 0; c < n; ++c) {
                        RowTrueCount[c] += (beta.get(i + xi, c) ? 1 : 0)
                                          - (beta.get(i, c) ? 1 : 0);
                    }
                }
            }
        }
    }

    // ---- (b2) Bcol(i,j) via per-row sliding column-window ---------------
    // Bcol(i,j) = OR_{j' in [j,j+xi-1]} beta(i+1, j')
    BoolMatrix Bcol; Bcol.init(n);
    for (int i = 0; i + 1 < n; ++i) {
        int row = i + 1;
        // prefix sum of beta(row, *)
        vector<int> prefix(n + 1, 0);
        for (int c = 0; c < n; ++c)
            prefix[c + 1] = prefix[c] + (beta.get(row, c) ? 1 : 0);
        for (int j = 0; j + xi <= n; ++j) {
            bool v = (prefix[j + xi] - prefix[j]) > 0;
            Bcol.set(i, j, v);
        }
    }

    // ---- (c1) Brow_band(i,j) via per-row sliding window over Brow -------
    // Brow_band(i,j) = AND_{j' in [j, j+xi-2]} Brow(i,j')     (xi-1 terms)
    //
    // NOTE on the off-by-one vs. the paper's literal formula: Brow(i,j')
    // itself inspects column j'+1. If j' were allowed to range all the way
    // up to j+xi-1 (=je), the check would inspect column je+1, i.e. one
    // position *past* the last point of the fixed-length candidate
    // subtrajectory S_{j,je} -- a column that does not belong to the
    // candidate at all. That variant of the rule is sound for the
    // *variable*-length motif problem of Section 3 (where an actual
    // end je can be arbitrarily larger than j+xi), but for Section 4's
    // *fixed*-length candidates (je = j+xi-1 exactly) it must be clamped
    // to j' in [j, je-1] so every inspected column stays inside [j+1, je].
    // This keeps the rule a valid (no-false-negative) necessary condition.
    BoolMatrix BrowBand; BrowBand.init(n);
    for (int i = 0; i < n; ++i) {
        vector<int> prefix(n + 1, 0);
        for (int c = 0; c < n; ++c)
            prefix[c + 1] = prefix[c] + (Brow.get(i, c) ? 1 : 0);
        for (int j = 0; j + xi <= n; ++j) {
            int lo = j, hi = j + xi - 1; // exclusive upper bound = je (=j+xi-1)
            bool v = (prefix[hi] - prefix[lo]) == (hi - lo); // xi-1 terms
            BrowBand.set(i, j, v);
        }
    }

    // ---- (c2) Bcol_band(i,j) via sliding window over rows on Bcol -------
    // Bcol_band(i,j) = AND_{i' in [i, i+xi-2]} Bcol(i',j)     (xi-1 terms)
    // (symmetric clamp for the same reason as Brow_band above)
    BoolMatrix BcolBand; BcolBand.init(n);
    {
        vector<int> ColTrueCount(n, 0); // count over rows [i, i+xi-2]
        int lastI = n - xi;
        if (lastI >= 0 && xi >= 2) {
            for (int c = 0; c < n; ++c) {
                int cnt = 0;
                for (int r = 0; r < xi - 1; ++r) cnt += Bcol.get(r, c) ? 1 : 0;
                ColTrueCount[c] = cnt;
            }
            for (int i = 0; i <= lastI; ++i) {
                for (int j = 0; j < n; ++j) {
                    BcolBand.set(i, j, ColTrueCount[j] == (xi - 1));
                }
                if (i < lastI) {
                    // slide row-window [i, i+xi-2] -> [i+1, i+xi-1]
                    for (int c = 0; c < n; ++c) {
                        ColTrueCount[c] += (Bcol.get(i + xi - 1, c) ? 1 : 0)
                                          - (Bcol.get(i, c) ? 1 : 0);
                    }
                }
            }
        }
    }

    // ---- (d) main candidate loop -----------------------------------------
    long long total = 0, prunedCell = 0, prunedCross = 0, prunedBand = 0, dfdCalls = 0;
    vector<FrequentPattern> results;

    for (int i = 0; i + xi - 1 < n; ++i) {
        if (!valid[i]) continue;
        int ie = i + xi - 1;
        FrequentPattern pat;
        pat.repStart = i;
        pat.repEnd = ie;

        for (int j = 0; j + xi - 1 < n; ++j) {
            if (!valid[j]) continue;
            int je = j + xi - 1;
            if (intervalsOverlap(i, ie, j, je)) continue;
            ++total;

            // cell-based pruning (Section 4.1.2)
            if (!beta.get(i, j)) { ++prunedCell; continue; }

            // cross-based pruning
            if (!(Brow.get(i, j) && Bcol.get(i, j))) { ++prunedCross; continue; }

            // band-based pruning
            if (!(BrowBand.get(i, j) && BcolBand.get(i, j))) { ++prunedBand; continue; }

            // survived all pruning rules -> exact DFD test
            double d = exactDFD(pts, i, ie, j, je);
            ++dfdCalls;
            if (d <= dthres) pat.num_members++;
        }

        if ((int)pat.num_members >= tthres) {
            results.push_back(std::move(pat));
        }
    }

    if (verbose) {
        cerr << "[BFP] total candidates considered : " << total << "\n";
        cerr << "[BFP] pruned by cell-based rule    : " << prunedCell << "\n";
        cerr << "[BFP] pruned by cross-based rule   : " << prunedCross << "\n";
        cerr << "[BFP] pruned by band-based rule    : " << prunedBand << "\n";
        cerr << "[BFP] exact DFD computations       : " << dfdCalls << "\n";
        if (total > 0) {
            double prunedPct = 100.0 * (prunedCell + prunedCross + prunedBand) / total;
            cerr << "[BFP] overall pruning ratio        : " << prunedPct << "%\n";
        }
    }
    return results;
}

// =========================================================================
// 3) BFP*  -- space-efficient variant (Section 4.2)
// =========================================================================
//
// Does not materialize the O(n^2) beta / Brow / Bcol / Bband matrices.
// Instead, everything is (re)computed on demand with only O(xi) or O(1)
// extra memory per test, exactly as described in Sec. 4.2: "it computes
// dG, beta, Bcross, Bband and dF on-the-fly". This trades some CPU time
// for drastically lower memory (useful for very long trajectories where
// an O(n^2) matrix would not fit in RAM).

static inline bool betaOnFly(const vector<Point> &pts, int n, double dthres,
                              int i, int j) {
    if (i < 0 || i >= n || j < 0 || j >= n) return false;
    return groundDistance(pts[i], pts[j]) <= dthres;
}

// Brow(i,j) = OR_{i' in [i,i+xi-1]} beta(i', j+1)     -- O(xi)
static inline bool browOnFly(const vector<Point> &pts, int n, double dthres,
                              int xi, int i, int j) {
    if (j + 1 >= n) return false;
    for (int ip = i; ip < i + xi; ++ip) {
        if (betaOnFly(pts, n, dthres, ip, j + 1)) return true;
    }
    return false;
}

// Bcol(i,j) = OR_{j' in [j,j+xi-1]} beta(i+1, j')     -- O(xi)
static inline bool colOnFly(const vector<Point> &pts, int n, double dthres,
                             int xi, int i, int j) {
    if (i + 1 >= n) return false;
    for (int jp = j; jp < j + xi; ++jp) {
        if (betaOnFly(pts, n, dthres, i + 1, jp)) return true;
    }
    return false;
}

// Brow_band(i,j) = AND_{j' in [j, j+xi-2]} Brow(i,j')  -- O(xi^2) worst case
// (clamped to xi-1 terms; see the long comment above BrowBand in bfp() for
//  why the naive xi-term range from the paper's motif-section formula would
//  probe one column past the fixed-length candidate's last valid column)
static inline bool rowBandOnFly(const vector<Point> &pts, int n, double dthres,
                                 int xi, int i, int j) {
    for (int jp = j; jp < j + xi - 1; ++jp) {
        if (!browOnFly(pts, n, dthres, xi, i, jp)) return false;
    }
    return true;
}

// Bcol_band(i,j) = AND_{i' in [i, i+xi-2]} Bcol(i',j)  -- O(xi^2) worst case
static inline bool colBandOnFly(const vector<Point> &pts, int n, double dthres,
                                 int xi, int i, int j) {
    for (int ip = i; ip < i + xi - 1; ++ip) {
        if (!colOnFly(pts, n, dthres, xi, ip, j)) return false;
    }
    return true;
}

static vector<FrequentPattern> bfpStar(const vector<Point> &pts, int xi,
                                        double dthres, int tthres,
                                        bool verbose = false) {
    int n = (int)pts.size();
    vector<char> valid = computeValidStarts(pts, xi);
    vector<FrequentPattern> results;

    long long total = 0, prunedCell = 0, prunedCross = 0, prunedBand = 0, dfdCalls = 0;

    for (int i = 0; i + xi - 1 < n; ++i) {
        if (!valid[i]) continue;
        int ie = i + xi - 1;
        FrequentPattern pat;
        pat.repStart = i;
        pat.repEnd = ie;

        for (int j = 0; j + xi - 1 < n; ++j) {
            if (!valid[j]) continue;
            int je = j + xi - 1;
            if (intervalsOverlap(i, ie, j, je)) continue;
            ++total;

            if (!betaOnFly(pts, n, dthres, i, j)) { ++prunedCell; continue; }

            bool br = browOnFly(pts, n, dthres, xi, i, j);
            bool bc = br ? colOnFly(pts, n, dthres, xi, i, j) : false; // short-circuit
            if (!(br && bc)) { ++prunedCross; continue; }

            bool bandOk = rowBandOnFly(pts, n, dthres, xi, i, j) &&
                          colBandOnFly(pts, n, dthres, xi, i, j);
            if (!bandOk) { ++prunedBand; continue; }

            double d = exactDFD(pts, i, ie, j, je);
            ++dfdCalls;
            if (d <= dthres) pat.num_members++;
        }

        if ((int)pat.num_members >= tthres) {
            results.push_back(std::move(pat));
        }
    }

    if (verbose) {
        cerr << "[BFP*] total candidates considered : " << total << "\n";
        cerr << "[BFP*] pruned by cell-based rule    : " << prunedCell << "\n";
        cerr << "[BFP*] pruned by cross-based rule   : " << prunedCross << "\n";
        cerr << "[BFP*] pruned by band-based rule    : " << prunedBand << "\n";
        cerr << "[BFP*] exact DFD computations       : " << dfdCalls << "\n";
    }
    return results;
}

// -------------------------------------------------------------------------
// Output helpers
// -------------------------------------------------------------------------

static void printResults(const vector<FrequentPattern> &results) {
    if (results.empty()) {
        cout << "No frequent patterns found.\n";
        return;
    }
    cout << "Found " << results.size() << " frequent pattern(s):\n";
    for (size_t k = 0; k < results.size(); ++k) {
        const auto &p = results[k];
        cout << "Pattern #" << (k + 1)
             << "  representative = [" << p.repStart << ", " << p.repEnd << "]"
             << "  |Gamma| = " << p.num_members << "\n";
        cout << "    matches: ";
        /*
        for (size_t m = 0; m < p.members.size(); ++m) {
            cout << "[" << p.members[m].first << "," << p.members[m].second << "]";
            if (m + 1 < p.members.size()) cout << ", ";
        }        
        */
        cout << "\n";
    }
}

// -------------------------------------------------------------------------
// main
// -------------------------------------------------------------------------

static void usage(const char *prog) {
    cerr << "Usage: " << prog
         << " <input_file> <xi> <dthres> <tthres> [brute|bfp|bfpstar] [-v]\n"
         << "  input_file : lines of \"x y trajectory_id\"\n"
         << "  xi         : subtrajectory length (integer > 1)\n"
         << "  dthres     : distance threshold (double)\n"
         << "  tthres     : frequency threshold (integer >= 1)\n"
         << "  algo       : which algorithm to run (default: bfp)\n"
         << "  -v         : print pruning statistics to stderr\n";
}

int main(int argc, char **argv) {
    if (argc < 5) {
        usage(argv[0]);
        return 1;
    }
    string inputFile = argv[1];
    int xi = std::stoi(argv[2]);
    double dthres = std::stod(argv[3]);
    int tthres = std::stoi(argv[4]);
    string algo = (argc >= 6) ? string(argv[5]) : "bfp";
    bool verbose = false;
    for (int a = 6; a < argc; ++a) {
        if (string(argv[a]) == "-v") verbose = true;
    }
    // allow "-v" to be passed as the 5th positional arg too
    if (algo == "-v") { verbose = true; algo = "bfp"; }

    if (xi <= 1) {
        cerr << "ERROR: xi must be > 1\n";
        return 1;
    }

    vector<Point> pts = loadTrajectory(inputFile);
    cerr << "Loaded " << pts.size() << " points from " << inputFile << "\n";

    auto t0 = high_resolution_clock::now();
    vector<FrequentPattern> results;

    if (algo == "brute") {
        results = bruteFP(pts, xi, dthres, tthres, verbose);
    } else if (algo == "bfp") {
        results = bfp(pts, xi, dthres, tthres, verbose);
    } else if (algo == "bfpstar") {
        results = bfpStar(pts, xi, dthres, tthres, verbose);
    } else {
        cerr << "ERROR: unknown algorithm '" << algo << "'\n";
        usage(argv[0]);
        return 1;
    }

    auto t1 = high_resolution_clock::now();
    double secs = duration_cast<duration<double>>(t1 - t0).count();

    printResults(results);
    cerr << "[" << algo << "] response time: " << secs << " sec\n";
    return 0;
}

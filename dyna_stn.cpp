/**
 * Dyna-STN: Dynamic Discrete Topology-Oriented Wide-Area Routing
 * Li, Wu, Wang — IEEE/ACM Transactions on Networking, Vol.32, No.5, 2024
 *
 * Complete implementation with all improvements:
 *   1.  Algorithm 1 (SID/VID Binding Procedure) — full 10-step state machine,
 *       Boolean(LocSID,C), Inquire(SID,LocSID,t), all message types
 *   2.  RouterState: RIB, FIB, RouterConfig, status RUNNING/SUSPENDED/MIGRATING
 *   3.  Service migration: tcp/tdp/tl/Tm (data/control-plane split, Section III-A-3)
 *   4.  Q(t) = fn_nodes + fl_links  (eq.5, counted separately)
 *   5.  OSPF simulation: Hello, LSA flooding, LSDB sync, convergence, LSA overhead
 *   6.  Ground stations + traffic: throughput, PLR, E2E delay
 *   7.  Calibrated ISL quality (SNR+duration utility, eq.28 meaningful values)
 *   8.  Correct Keplerian satellite propagation (inclined polar orbit)
 *
 * Compile:  g++ -std=c++17 -O2 -o dyna_stn dyna_stn_full.cpp
 */

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <cmath>
#include <limits>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <functional>

// ================================================================
//  §0  CONSTANTS  (Table IV and simulation parameters)
// ================================================================
namespace C {
    constexpr double PI      = 3.141592653589793;
    constexpr double Re      = 6371.0;              // km
    constexpr double h_orb   = 780.0;               // km  Iridium altitude
    constexpr double R_orb   = Re + h_orb;           // km
    constexpr double incl    = 86.4 * PI / 180.0;   // rad  inclination
    constexpr int    M       = 6;                    // orbital planes
    constexpr int    N       = 11;                   // satellites per plane
    constexpr int    T_orb   = 6900;                // s   orbital period
    constexpr double polar_B = 70.0 * PI / 180.0;   // rad  polar blackout boundary

    // ── ISL quality (eqs.14-22, calibrated for Iridium distances) ──
    // U_SNR(L) = 1 – exp( –Pt_eff · (Lv/L)² )
    // At L = Lv ≈ 4054 km  →  U_SNR ≈ 0.950  (Pt_eff = 3.0)
    // At L = Lh ≈ 2035 km  →  U_SNR ≈ 0.9999
    constexpr double Pt_eff  = 3.0;   // calibrated SNR coefficient
    constexpr double ws      = 0.6;   // SNR utility weight       (eq.28)
    constexpr double wl      = 0.4;   // duration utility weight  (eq.28)

    // ── Coverage ──────────────────────────────────────────────────
    constexpr double theta_min = 10.0 * PI / 180.0;  // rad  min elevation (eqs.8-11)

    // ── Service-migration timing (calibrated: Tm_66nodes ≈ 0.13 s) ──
    // Step 2-3: control-plane migration
    constexpr double tcp_base      = 0.050;   // s  base cost (tunnel + image copy)
    constexpr double tcp_per_route = 0.0005;  // s  per RIB entry
    // Step 4: data-plane clone (FIB installation 100-500 µs per entry)
    constexpr double tdp_per_fib   = 0.00030; // s  per FIB entry
    // Step 5: asynchronous link migration
    constexpr double tl_per_link   = 0.010;   // s  per ISL

    // ── OSPF parameters ───────────────────────────────────────────
    constexpr int    hello_ivl  = 10;        // s
    constexpr int    dead_ivl   = 40;        // s
    constexpr double c_light    = 3.0e5;     // km/s (propagation speed)
    // Dijkstra calibration: 66 nodes → ~0.35 ms total for all RIBs
    constexpr double dijk_coeff = 1.5e-6;    // s per (V*(E+V)*logV) op

    // ── Traffic ───────────────────────────────────────────────────
    constexpr double bw_link   = 100.0;   // Mbps per ISL
    constexpr double queue_ms  = 0.5;     // ms queuing delay per hop
}

// ================================================================
//  §1  GEOMETRY
// ================================================================

struct Vec3 {
    double x{0}, y{0}, z{0};

    double norm() const { return std::sqrt(x*x + y*y + z*z); }

    double dist(const Vec3& o) const {
        double dx=x-o.x, dy=y-o.y, dz=z-o.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Geodetic latitude (radians) from ECEF
    double lat_rad() const {
        double r = norm();
        return (r > 0) ? std::asin(z / r) : 0.0;
    }

    std::string str() const {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(0)
           << "(" << x << "," << y << "," << z << ")";
        return ss.str();
    }
};

// Geographic (lat°, lon°) → ECEF km
Vec3 llToECEF(double lat_deg, double lon_deg, double r = C::Re) {
    double la = lat_deg * C::PI / 180.0;
    double lo = lon_deg * C::PI / 180.0;
    return { r * std::cos(la) * std::cos(lo),
             r * std::cos(la) * std::sin(lo),
             r * std::sin(la) };
}

// ── ISL lengths (eqs.1-2) ──────────────────────────────────────

// Intra-plane ISL length (eq.1)
double Lv_km() {
    return std::sqrt(2.0) * C::R_orb
           * std::sqrt(1.0 - std::cos(2.0 * C::PI / C::N));
}

// Inter-plane ISL length at latitude lat_rad, F=0 (eq.2)
double Lh_km(double lat_rad) {
    return std::sqrt(2.0) * C::R_orb
           * std::sqrt((1.0 - std::cos(C::PI / C::N)) * std::cos(lat_rad));
}

// Max communication time T (eqs.8-11)
double maxCommTime_s() {
    double theta = C::theta_min;
    double ratio = C::Re / C::R_orb;
    double gamma = std::acos(ratio * std::cos(theta)) - theta;    // eq.9
    double L     = 2.0 * C::R_orb * gamma;                        // eq.10
    double vs    = std::sqrt(3.986e5 / C::R_orb);                 // km/s satellite speed
    return L / vs;                                                 // eq.11
}

// ================================================================
//  §2  ISL QUALITY MODEL  (eqs.14-28, calibrated)
// ================================================================

namespace ISL {

static double Lv_ref = 0.0;   // set once in main

// ── U_SNR (eqs.14-22, calibrated) ────────────────────────────
// Replaces the raw Gaussian-beam formula that yields U≈0 at Iridium
// distances due to missing unit normalisation in the paper.
// We keep the functional form SNR ∝ L^-2, calibrated so that:
//   U_SNR(Lv) ≈ 0.95  (intra-plane ISL at reference distance)
double uSNR(double Lij) {
    if (Lij <= 0.0 || Lv_ref <= 0.0) return 0.0;
    double snr = C::Pt_eff * (Lv_ref / Lij) * (Lv_ref / Lij);  // eq.17 spirit
    return 1.0 - std::exp(-snr);                                  // eq.22
}

// ── U_dur (eqs.23-27) ─────────────────────────────────────────
// Intra-plane ISLs are always active (paper Section IV-B).
// Inter-plane ISLs are disrupted inside polar region (|lat| ≥ 70°).
double uDur(bool intraPlane, const Vec3& posA, const Vec3& posB) {
    if (intraPlane) return 1.0;              // intra-plane: always connected
    double latA = posA.lat_rad();
    double latB = posB.lat_rad();
    if (std::abs(latA) >= C::polar_B || std::abs(latB) >= C::polar_B)
        return 0.0;                          // polar blackout → ISL down
    // Linear decay as satellites approach the polar boundary
    double maxLat = std::max(std::abs(latA), std::abs(latB));
    return (C::polar_B - maxLat) / C::polar_B;   // eq.27 spirit
}

// ── U_total (eq.28) ───────────────────────────────────────────
double uTotal(double Lij, bool intraPlane,
              const Vec3& posA, const Vec3& posB) {
    double us = uSNR(Lij);
    double ul = uDur(intraPlane, posA, posB);
    if (ul <= 0.0) return 0.0;
    return std::pow(us, C::ws) * std::pow(ul, C::wl);   // eq.28
}

} // namespace ISL

// ================================================================
//  §3  DATA STRUCTURES
// ================================================================

using SID = std::string;
using VID = std::string;

// ── Router status ─────────────────────────────────────────────
enum class RouterStatus { RUNNING, SUSPENDED, MIGRATING };

inline std::string rsStr(RouterStatus s) {
    switch(s) {
        case RouterStatus::RUNNING:   return "RUNNING";
        case RouterStatus::SUSPENDED: return "SUSPENDED";
        case RouterStatus::MIGRATING: return "MIGRATING";
    }
    return "?";
}

// ── RouterConfig ──────────────────────────────────────────────
struct RouterConfig {
    int    helloIvl = C::hello_ivl;   // s
    int    deadIvl  = C::dead_ivl;    // s
    double linkBW   = C::bw_link;     // Mbps
    int    area     = 0;              // OSPF area
    int    lsaSeq   = 0;              // LSA sequence number base
};

// ── RIB / FIB entries ─────────────────────────────────────────
struct RIBEntry { VID nextHop; double metric; };
struct FIBEntry { VID outIface; double cost;  };

// ── RouterState (per virtual router / virtual node) ────────────
struct RouterState {
    RouterStatus              status     = RouterStatus::RUNNING;
    std::map<VID, RIBEntry>  rib;        // Routing Information Base
    std::map<VID, FIBEntry>  fib;        // Forwarding Information Base
    RouterConfig              cfg;
    double                    downStart  = -1.0;  // s (-1 = not suspended)
    double                    totalDown  =  0.0;  // accumulated downtime, s

    void suspend(double now) {
        if (status != RouterStatus::SUSPENDED) {
            status    = RouterStatus::SUSPENDED;
            downStart = now;
        }
    }

    void restart(double now) {
        if (downStart >= 0.0) {
            totalDown += now - downStart;
            downStart  = -1.0;
        }
        status = RouterStatus::RUNNING;
    }

    // Populate FIB directly from RIB
    void updateFIBFromRIB() {
        fib.clear();
        for (auto& [dst, e] : rib)
            fib[dst] = { e.nextHop, e.metric };
    }
};

// ── Satellite (Keplerian propagation) ─────────────────────────
struct Satellite {
    SID    sid;
    int    plane, idx;
    double raan;       // right ascension of ascending node, rad
    double initAnom;   // initial mean anomaly, rad
    Vec3   pos;        // current ECEF position, km

    // Update position for absolute time t (seconds)
    void update(double t) {
        double omega = 2.0 * C::PI / C::T_orb;
        double anom  = initAnom + omega * t;
        double inc   = C::incl;
        double R     = C::R_orb;
        // Standard Euler rotation for circular orbit with RAAN and inclination
        pos.x = R * (std::cos(raan)*std::cos(anom)
                    - std::sin(raan)*std::sin(anom)*std::cos(inc));
        pos.y = R * (std::sin(raan)*std::cos(anom)
                    + std::cos(raan)*std::sin(anom)*std::cos(inc));
        pos.z = R *  std::sin(anom)*std::sin(inc);
    }

    bool inPolar() const { return std::abs(pos.lat_rad()) >= C::polar_B; }
};

// ── Virtual link (labelled intra/inter plane) ──────────────────
struct VLink {
    VID  u, v;
    bool intra;   // true = intra-plane, false = inter-plane
};

// ── Virtual node (fixed service cube) ─────────────────────────
struct VirtualNode {
    VID         vid;
    int         plane, idx;
    Vec3        ctr;            // fixed centre of service cube
    bool        active   = false;
    SID         bound    = "";  // currently bound satellite
    SID         prevBound= "";  // binding at previous slot
    RouterState router;
};

// ── Ground station ─────────────────────────────────────────────
struct GS {
    std::string id;
    Vec3        pos;
    VID         nearVID;
    double      nearDist = 1e18;
};

// ================================================================
//  §4  DIJKSTRA + LSDB  (Section III-B: Routing Calculation)
// ================================================================

struct LSAEntry {
    VID                               origin;
    int                               seq = 0;
    std::vector<std::pair<VID,double>> links;  // (neighbour, cost)

    bool sameLinks(const LSAEntry& o) const { return links == o.links; }
};

using LSDB = std::map<VID, LSAEntry>;
using RIBMap = std::map<VID, RIBEntry>;

// Dijkstra shortest-path tree from `src` over the LSDB (eq.2 of Routing Calc)
RIBMap dijkstra(const LSDB& lsdb, const VID& src) {
    constexpr double INF = std::numeric_limits<double>::infinity();

    std::map<VID,double>  dist;
    std::map<VID,VID>     prev;
    for (auto& [v,_] : lsdb) { dist[v] = INF; prev[v] = ""; }
    dist[src] = 0.0;

    using PQE = std::pair<double,VID>;
    std::priority_queue<PQE, std::vector<PQE>, std::greater<PQE>> pq;
    pq.push({0.0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        auto it = lsdb.find(u);
        if (it == lsdb.end()) continue;
        for (auto& [v, w] : it->second.links) {
            if (w <= 0.0 || w >= INF) continue;
            double nd = dist[u] + w;
            if (nd < dist[v]) {
                dist[v] = nd;  prev[v] = u;
                pq.push({nd, v});
            }
        }
    }

    RIBMap rib;
    for (auto& [dst, _] : lsdb) {
        if (dst == src) continue;
        if (dist[dst] >= INF) { rib[dst] = {"UNREACHABLE", INF}; continue; }
        // Trace back to find first hop from src
        VID cur = dst;
        while (!prev[cur].empty() && prev[cur] != src) cur = prev[cur];
        rib[dst] = {cur, dist[dst]};
    }
    return rib;
}

// Reconstruct full path from per-node RIBs
std::vector<VID> getPath(
    const std::map<VID, RIBMap>& ribs,
    const VID& src, const VID& dst)
{
    std::vector<VID> path;
    if (ribs.find(src) == ribs.end()) return path;
    path.push_back(src);
    VID cur = src;
    std::set<VID> visited;
    while (cur != dst && !visited.count(cur)) {
        visited.insert(cur);
        auto rit = ribs.find(cur);
        if (rit == ribs.end()) break;
        auto dit = rit->second.find(dst);
        if (dit == rit->second.end() ||
            dit->second.nextHop == "UNREACHABLE") break;
        cur = dit->second.nextHop;
        path.push_back(cur);
    }
    return path;
}

// ================================================================
//  §5  OSPF SIMULATOR  (Section III-B)
// ================================================================

struct OSPFStats {
    double routeCalc_ms  = 0;   // Dijkstra time (all nodes)
    double converge_ms   = 0;   // flooding + route-calc
    int    lsaPkts       = 0;   // LSA packets generated this slot
    int    helloPkts     = 0;   // Hello packets
    int    changedLinks  = 0;   // number of link-state changes
};

class OSPFSim {
public:
    LSDB                      lsdb;         // global link-state database
    std::map<VID,int>         seqNum;       // per-node LSA sequence numbers
    std::map<VID, RIBMap>     ribs;         // per-node RIBs
    int                       totalLSA = 0; // cumulative LSA overhead

    // One OSPF epoch per time slot
    // Implements: topology-establishment, route-calculation, link-failure-response
    OSPFStats runEpoch(const std::vector<VirtualNode>& vns,
                       const std::vector<VLink>&       links,
                       const std::vector<Satellite>&   sats,
                       double /*slotTime*/)
    {
        OSPFStats st;
        LSDB newLSDB;
        int V = 0;

        // ── 1. Build fresh LSA for each active virtual node ──────
        for (const auto& vn : vns) {
            if (!vn.active) continue;
            V++;
            LSAEntry entry;
            entry.origin = vn.vid;
            entry.seq    = seqNum[vn.vid]++;

            const Satellite* sA = findSat(sats, vn.bound);

            for (const auto& lk : links) {
                // Which neighbour VID is this link's other endpoint?
                VID nb;
                bool intra = lk.intra;
                if      (lk.u == vn.vid) nb = lk.v;
                else if (lk.v == vn.vid) nb = lk.u;
                else continue;

                const VirtualNode* vnB = findVN(vns, nb);
                if (!vnB || !vnB->active) continue;

                // Skip counter-rotating seam inter-plane links
                if (!intra && isSeam(vn.plane, vnB->plane)) continue;

                const Satellite* sB = findSat(sats, vnB->bound);
                if (!sA || !sB) continue;

                double L = sA->pos.dist(sB->pos);
                double U = ISL::uTotal(L, intra, sA->pos, sB->pos);
                if (U <= 0.0) continue;   // ISL down → exclude from LSDB

                entry.links.push_back({nb, 1.0 / U});
            }
            newLSDB[vn.vid] = entry;
        }

        // ── 2. Detect topology changes vs previous slot ───────────
        for (auto& [vid, newE] : newLSDB) {
            auto it = lsdb.find(vid);
            if (it == lsdb.end() || !it->second.sameLinks(newE))
                st.changedLinks++;
        }

        // ── 3. Hello packets (Section III-B, Topology Establishment)
        // Each active node sends Hello to ≤4 neighbours per slot
        st.helloPkts = V * 4;

        // ── 4. LSA flooding (Section III-B)
        // Each changed LSA is flooded to all V nodes (OSPF flooding mechanism).
        // Dyna-STN: once the virtual overlay is established,
        // topology changes = 0 → near-zero ongoing LSA overhead.
        st.lsaPkts  = st.changedLinks * V;
        totalLSA   += st.lsaPkts;

        // ── 5. Convergence time estimate ──────────────────────────
        // diameter of Iridium grid+ ≈ M + N/2
        int  E          = 0;
        for (auto& [vid, e] : newLSDB) E += (int)e.links.size();
        E /= 2;   // undirected
        int  diam       = C::M + C::N / 2;
        double hopMs    = 2034.0 / C::c_light * 1000.0;   // ~6.8 ms
        double floodMs  = diam * hopMs;
        // Dijkstra for all V nodes: O(V·(V+E)·logV)
        double dijkOps  = (double)V * (V + E) * std::log2(V + 1.0);
        double calcMs   = dijkOps * C::dijk_coeff * 1000.0;
        st.routeCalc_ms = calcMs;
        st.converge_ms  = floodMs + calcMs;

        // ── 6. Update global LSDB and compute per-node RIBs ───────
        lsdb = newLSDB;
        ribs.clear();
        for (const auto& vn : vns)
            if (vn.active)
                ribs[vn.vid] = dijkstra(lsdb, vn.vid);

        return st;
    }

    std::pair<VID,double> route(const VID& src, const VID& dst) const {
        auto it = ribs.find(src);
        if (it == ribs.end()) return {"NO_RIB", -1};
        auto jt = it->second.find(dst);
        if (jt == it->second.end()) return {"NOT_FOUND", -1};
        return {jt->second.nextHop, jt->second.metric};
    }

private:
    const VirtualNode* findVN(const std::vector<VirtualNode>& vns,
                               const VID& vid) const {
        for (auto& vn : vns) if (vn.vid == vid) return &vn;
        return nullptr;
    }
    const Satellite* findSat(const std::vector<Satellite>& sats,
                              const SID& sid) const {
        for (auto& s : sats) if (s.sid == sid) return &s;
        return nullptr;
    }
    // Counter-rotating seam between plane 0 and plane M-1
    bool isSeam(int a, int b) const {
        return (a == 0 && b == C::M-1) || (a == C::M-1 && b == 0);
    }
};

// ================================================================
//  §6  SERVICE MIGRATION  (Section III-A-3, Fig.7)
// ================================================================

struct MigRes {
    SID    from, to;
    VID    vid;
    double tcp;      // control-plane migration time, s
    double tdp;      // data-plane clone time, s
    double tl;       // link migration time, s
    double Tm;       // total migration time (tcp+tdp+tl), s
    double downtime; // VID downtime (stall-and-copy phase only), s
};

// Data/control-plane split service migration
// nLinks = number of ISLs the satellite has (up to 4 for Iridium)
MigRes performMigration(const SID& from, const SID& to, const VID& vid,
                         const RouterState& oldRS, int nLinks,
                         double /*now*/)
{
    int nRoutes = (int)oldRS.rib.size();
    int nFIB    = std::max((int)oldRS.fib.size(), nRoutes);

    MigRes m;
    m.from = from;  m.to = to;  m.vid = vid;

    // Step 1 (tunnel setup) is included in tcp_base
    // Steps 2-3: pre-copy + stall-and-copy (router image + memory)
    m.tcp = C::tcp_base + nRoutes * C::tcp_per_route;

    // Step 4: data-plane clone via data-plane hypervisor
    m.tdp = nFIB * C::tdp_per_fib;

    // Step 5: asynchronous link migration
    m.tl  = nLinks * C::tl_per_link;

    m.Tm      = m.tcp + m.tdp + m.tl;
    m.downtime = m.tcp / 3.0;   // only stall-and-copy (Step 3-2) is "down"
    return m;
}

// ================================================================
//  §7  SID/VID MAPPING SERVER  (Table III)
// ================================================================

class MapServer {
public:
    std::map<int, std::map<SID,VID>> table;   // slot → {SID→VID}
    std::map<VID, bool>              actv;    // VID activation status

    void store(int slot, const SID& sid, const VID& vid) {
        table[slot][sid] = vid;
        actv[vid] = true;
    }
    void deactivate(const VID& vid) { actv[vid] = false; }
    bool isActive(const VID& vid) const {
        auto it = actv.find(vid);
        return it != actv.end() && it->second;
    }

    // eq.13: VID = Inquire(SID, LocSID, t)
    // Returns VID of the service cube containing LocSID (or "" if none)
    VID inquire(const Vec3& loc,
                const std::vector<VirtualNode>& vns,
                double cubeR) const
    {
        double best = cubeR;
        VID    bestVID;
        for (auto& vn : vns) {
            double d = vn.ctr.dist(loc);
            if (d < best) { best = d; bestVID = vn.vid; }
        }
        return bestVID;
    }

    void print(int maxSlots = 2) const {
        std::cout << "\n══ SID/VID Mapping Table (Table III, last "
                  << maxSlots << " slots) ══\n";
        std::cout << std::left << std::setw(6)  << "Slot"
                              << std::setw(16) << "SID"
                              << "VID\n"
                  << std::string(42, '-') << "\n";
        int cnt = 0;
        for (auto rit = table.rbegin();
             rit != table.rend() && cnt < maxSlots; ++rit, ++cnt)
        {
            for (auto& [sid, vid] : rit->second)
                std::cout << std::setw(6) << rit->first
                          << std::setw(16) << sid << vid << "\n";
        }
    }
};

// ================================================================
//  §8  ALGORITHM 1 — SID/VID BINDING PROCEDURE (full state machine)
//      Implements: eqs.12-13, Fig.5 (dynamic binding mechanism),
//                  Fig.6 (signalling process), Algorithm 1 pseudocode
// ================================================================

// Message types (from Figs.5-6 of the paper)
enum class MsgType {
    Inquire_FromSID,            // satellite → mapping server
    Response_FromMapSrv,        // mapping server → new satellite
    Shift_FromMapSrv,           // mapping server → old satellite
    Request_FromSID,            // old satellite → new satellite
    Confirm_FromSID,            // new satellite → old satellite
    Migrate_FromSID,            // old satellite → new satellite (carries RIB+config)
    Complete_FromSID            // new satellite → old satellite
};

inline std::string msgStr(MsgType t) {
    switch(t) {
        case MsgType::Inquire_FromSID:       return "Inquire";
        case MsgType::Response_FromMapSrv:   return "Response";
        case MsgType::Shift_FromMapSrv:      return "Shift";
        case MsgType::Request_FromSID:       return "Request";
        case MsgType::Confirm_FromSID:       return "Confirm";
        case MsgType::Migrate_FromSID:       return "Migrate";
        case MsgType::Complete_FromSID:      return "Complete";
    }
    return "?";
}

// A single signalling message
struct Msg {
    MsgType     type;
    std::string from, to;   // SID / "MapSrv"
    VID         vid;
    RouterState rs;          // carried by Migrate message
};

// Per-satellite agent (tracks Algorithm 1 state across slots)
struct Agent {
    SID         sid;
    VID         curVID;      // currently bound VID (empty = not bound)
    RouterState router;
};

// ── DDTM Plane (manages all SID/VID bindings and migrations) ────
class DDTMPlane {
public:
    MapServer              mapSrv;
    std::map<SID, Agent>   agents;             // per-satellite state
    std::map<VID, RouterState> vidRS;          // authoritative router state per VID
    std::map<VID, SID>     vidCurSID;          // VID → bound SID
    std::vector<MigRes>    migrations;         // migrations this slot
    std::vector<std::string> msgLog;           // Algorithm 1 trace

    int fn_sum = 0;   // eq.5 node-binding count
    int fl_sum = 0;   // eq.5 link-binding count

    // ── Run Algorithm 1 for every satellite, process migrations ──
    int runSlot(int slot, double now,
                std::vector<VirtualNode>& vns,
                const std::vector<Satellite>& sats,
                const std::vector<VLink>& vlinks,
                double cubeR)
    {
        fn_sum = fl_sum = 0;
        migrations.clear();
        msgLog.clear();

        // ── Phase 1: determine best satellite for each VN ────────
        // Implements eq.12  Boolean(LocSID, C)
        // For each satellite, check if it has entered a *new* service cube
        struct BestBind { VID vid; SID sid; double dist; };
        std::map<VID, BestBind> vnToBest;   // VID → closest in-cube satellite

        for (const auto& sat : sats) {
            // eq.12: Boolean(LocSID, C) — is satellite inside any service cube?
            VID inCubeVID = mapSrv.inquire(sat.pos, vns, cubeR);
            if (inCubeVID.empty()) continue;

            double d = 0.0;
            for (auto& vn : vns)
                if (vn.vid == inCubeVID) { d = vn.ctr.dist(sat.pos); break; }

            auto it = vnToBest.find(inCubeVID);
            if (it == vnToBest.end() || d < it->second.dist)
                vnToBest[inCubeVID] = {inCubeVID, sat.sid, d};
        }

        // ── Phase 2: Message queue (index-based, allows push_back) ─
        std::vector<Msg> Q;

        // Algorithm 1, line 2: if Boolean → send Inquire_FromSID_n
        for (auto& [vid, bb] : vnToBest) {
            auto& ag = agents[bb.sid];
            ag.sid = bb.sid;

            if (ag.curVID != vid) {
                // Satellite entering a new service cube
                Q.push_back({MsgType::Inquire_FromSID,
                             bb.sid, "MapSrv", vid, {}});
                msgLog.push_back("Inquire: " + bb.sid + " → MapSrv (VID=" + vid + ")");
            }
        }

        // ── Phase 3: Process message queue ───────────────────────
        for (size_t i = 0; i < Q.size(); i++) {
            // NOTE: Q may grow inside this loop; use index, not iterator
            const Msg msg = Q[i];   // copy to avoid dangling ref after push_back

            switch (msg.type) {

            // ── Step 2: Mapping server → Response + optional Shift ──
            case MsgType::Inquire_FromSID: {
                const SID& newSID = msg.from;
                const VID& vid    = msg.vid;

                // Algorithm 1 lines 4-6: Response_FromMappingServer
                Q.push_back({MsgType::Response_FromMapSrv,
                             "MapSrv", newSID, vid, {}});
                msgLog.push_back("Response: MapSrv → " + newSID);

                // Algorithm 1 lines 7-8: Shift to old satellite if VID occupied
                auto it = vidCurSID.find(vid);
                if (it != vidCurSID.end() && !it->second.empty()
                    && it->second != newSID)
                {
                    Q.push_back({MsgType::Shift_FromMapSrv,
                                 "MapSrv", it->second, vid, {}});
                    msgLog.push_back("Shift: MapSrv → " + it->second
                                     + " (VID=" + vid + ")");
                }
                break;
            }

            // ── Step 3: New satellite binds VID (lines 4-6) ────────
            case MsgType::Response_FromMapSrv: {
                const SID& newSID = msg.to;
                const VID& vid    = msg.vid;
                auto& ag          = agents[newSID];
                ag.curVID         = vid;
                mapSrv.store(slot, newSID, vid);
                break;
            }

            // ── Steps 4-5: Old satellite gets Shift → sends Request ─
            // Algorithm 1 lines 7-8
            case MsgType::Shift_FromMapSrv: {
                const SID& oldSID = msg.to;
                const VID& vid    = msg.vid;
                // Find the new satellite now bound to this VID
                SID newSID;
                for (auto& [sid, ag] : agents)
                    if (ag.curVID == vid && sid != oldSID)
                        { newSID = sid; break; }
                if (newSID.empty()) break;

                // old sends Request to new
                RouterState rsToSend = vidRS.count(vid) ? vidRS[vid]
                                                         : agents[oldSID].router;
                Q.push_back({MsgType::Request_FromSID,
                             oldSID, newSID, vid, rsToSend});
                msgLog.push_back("Request: " + oldSID + " → " + newSID);
                break;
            }

            // ── Step 6: New satellite suspends + sends Confirm ──────
            // Algorithm 1 lines 9-11
            case MsgType::Request_FromSID: {
                const SID& oldSID = msg.from;
                const SID& newSID = msg.to;
                const VID& vid    = msg.vid;
                // New satellite suspends its router
                agents[newSID].router.suspend(now);
                // Send Confirm to old satellite
                Q.push_back({MsgType::Confirm_FromSID,
                             newSID, oldSID, vid, {}});
                msgLog.push_back("Confirm: " + newSID + " → " + oldSID
                                 + " (suspended)");
                break;
            }

            // ── Step 7: Old satellite migrates RIB+config ───────────
            // Algorithm 1 lines 12-13
            case MsgType::Confirm_FromSID: {
                const SID& newSID = msg.from;   // sent the Confirm
                const SID& oldSID = msg.to;     // receives Confirm
                const VID& vid    = msg.vid;

                RouterState& rs = vidRS.count(vid) ? vidRS[vid]
                                                   : agents[oldSID].router;
                int nLinks = 4;  // Iridium: up to 4 ISLs
                MigRes mr = performMigration(oldSID, newSID, vid, rs, nLinks, now);
                migrations.push_back(mr);

                // Old satellite sends Migrate (carries RIB + config)
                Q.push_back({MsgType::Migrate_FromSID,
                             oldSID, newSID, vid, rs});
                msgLog.push_back("Migrate: " + oldSID + " → " + newSID
                                 + " Tm=" + std::to_string(mr.Tm).substr(0,5) + "s");
                break;
            }

            // ── Steps 8-9: New satellite updates RIB, restarts, sends Complete
            // Algorithm 1 lines 14-17
            case MsgType::Migrate_FromSID: {
                const SID& oldSID = msg.from;
                const SID& newSID = msg.to;
                const VID& vid    = msg.vid;

                // Update RIB and config files
                auto& ag    = agents[newSID];
                ag.router   = msg.rs;           // copy RIB + config from old
                ag.router.restart(now);          // Restart
                vidRS[vid]  = ag.router;

                // Send Complete to old satellite
                Q.push_back({MsgType::Complete_FromSID,
                             newSID, oldSID, vid, {}});
                msgLog.push_back("Complete: " + newSID + " → " + oldSID);
                break;
            }

            // ── Step 10: Old satellite receives Complete ─────────────
            // Algorithm 1 pseudocode step 10: old satellite acknowledges.
            // We must NOT clear curVID unconditionally — the old satellite
            // may have already been re-bound to a NEW service cube within
            // this same slot (Inquire→Response processed earlier in the
            // message queue). Only clear if curVID still points to the VID
            // being migrated away from.
            case MsgType::Complete_FromSID: {
                const SID& oldSID = msg.to;
                const VID& vid    = msg.vid;
                if (agents.count(oldSID) && agents[oldSID].curVID == vid)
                    agents[oldSID].curVID = "";
                break;
            }

            } // end switch
        } // end message queue processing

        // ── Phase 4: Commit bindings to virtual nodes ─────────────
        for (auto& vn : vns) {
            vn.prevBound = vn.bound;
            vn.active    = false;
            vn.bound     = "";
        }
        for (auto& [sid, ag] : agents) {
            if (ag.curVID.empty()) continue;
            for (auto& vn : vns) {
                if (vn.vid == ag.curVID) {
                    vn.bound  = sid;
                    vn.active = true;
                    vn.router = ag.router;
                    break;
                }
            }
        }

        // ── Phase 5: Update vidCurSID map ─────────────────────────
        vidCurSID.clear();
        for (auto& vn : vns)
            if (vn.active) vidCurSID[vn.vid] = vn.bound;

        // ── Phase 6: Compute Q(t) = fn + fl  (eq.5) ──────────────

        // fn(i,j): one per active node binding
        for (auto& vn : vns)
            if (vn.active) fn_sum++;

        // fl(i,j): one per active *physical* ISL bound to a virtual link
        std::set<std::pair<VID,VID>> seenLinks;
        for (auto& lk : vlinks) {
            auto key = std::make_pair(std::min(lk.u, lk.v),
                                      std::max(lk.u, lk.v));
            if (seenLinks.count(key)) continue;

            const VirtualNode* vnA = findVN(vns, lk.u);
            const VirtualNode* vnB = findVN(vns, lk.v);
            if (!vnA || !vnB || !vnA->active || !vnB->active) continue;

            const Satellite* sA = findSat(sats, vnA->bound);
            const Satellite* sB = findSat(sats, vnB->bound);
            if (!sA || !sB) continue;

            // Inter-plane: skip polar blackout and counter-rotating seam
            if (!lk.intra) {
                if (sA->inPolar() || sB->inPolar()) continue;
                if (isSeam(vnA->plane, vnB->plane)) continue;
            }

            seenLinks.insert(key);
            fl_sum++;
        }

        return fn_sum + fl_sum;  // Q(t)
    }

    double avgTm() const {
        if (migrations.empty()) return 0.0;
        double s = 0;
        for (auto& m : migrations) s += m.Tm;
        return s / migrations.size();
    }

private:
    const VirtualNode* findVN(const std::vector<VirtualNode>& vns,
                               const VID& vid) const {
        for (auto& vn : vns) if (vn.vid == vid) return &vn;
        return nullptr;
    }
    const Satellite* findSat(const std::vector<Satellite>& sats,
                              const SID& sid) const {
        for (auto& s : sats) if (s.sid == sid) return &s;
        return nullptr;
    }
    bool isSeam(int a, int b) const {
        return (a == 0 && b == C::M-1) || (a == C::M-1 && b == 0);
    }
};

// ================================================================
//  §9  TRAFFIC SIMULATION
// ================================================================

struct TrafficResult {
    std::string src, dst;
    double throughput_mbps = 0;
    double pktLoss_pct     = 100;
    double e2eDelay_ms     = 9999;
    int    hops            = 0;
};

TrafficResult simTraffic(const GS& gsS, const GS& gsD,
                          const OSPFSim& ospf,
                          const std::vector<VirtualNode>& vns,
                          const std::vector<Satellite>& sats,
                          const std::vector<VLink>& /*vlinks*/)
{
    TrafficResult res;
    res.src = gsS.id;  res.dst = gsD.id;

    if (gsS.nearVID.empty() || gsD.nearVID.empty()) return res;

    // Find path through virtual overlay
    std::vector<VID> path = getPath(ospf.ribs, gsS.nearVID, gsD.nearVID);
    if (path.size() < 2) return res;

    res.hops = (int)path.size() - 1;

    // ── Path quality: product of per-hop ISL utilities ─────────
    double pathQ  = 1.0;
    double pathKm = 0.0;

    for (int i = 0; i < res.hops; i++) {
        const VirtualNode* a = nullptr, *b = nullptr;
        for (auto& vn : vns) {
            if (vn.vid == path[i])     a = &vn;
            if (vn.vid == path[i+1])   b = &vn;
        }
        if (!a || !b) { pathQ *= 0.5; continue; }
        pathKm += a->ctr.dist(b->ctr);

        const Satellite* sA = nullptr, *sB = nullptr;
        for (auto& s : sats) {
            if (s.sid == a->bound) sA = &s;
            if (s.sid == b->bound) sB = &s;
        }
        // Determine intra vs inter plane for this hop
        bool intra = (a->plane == b->plane);
        double L   = (sA && sB) ? sA->pos.dist(sB->pos)
                                 : a->ctr.dist(b->ctr);
        Vec3   pA  = sA ? sA->pos : a->ctr;
        Vec3   pB  = sB ? sB->pos : b->ctr;
        double U   = ISL::uTotal(L, intra, pA, pB);
        if (U <= 0.0) U = 1e-4;
        pathQ *= U;
    }

    // ── Distances: add uplink + downlink ──────────────────────
    auto vnCtr = [&](const VID& vid) -> Vec3 {
        for (auto& vn : vns) if (vn.vid == vid) return vn.ctr;
        return {};
    };
    double uplinkKm   = gsS.pos.dist(vnCtr(gsS.nearVID));
    double downlinkKm = gsD.pos.dist(vnCtr(gsD.nearVID));
    double totalKm    = pathKm + uplinkKm + downlinkKm;

    // ── E2E delay (ms) ─────────────────────────────────────────
    double propMs = totalKm / C::c_light * 1000.0;
    res.e2eDelay_ms = propMs + res.hops * C::queue_ms;

    // ── Packet loss rate (%) ───────────────────────────────────
    res.pktLoss_pct = (1.0 - pathQ) * 100.0;

    // ── Throughput (Mbps) ──────────────────────────────────────
    double effBW = C::bw_link / std::max(1, res.hops);  // bottleneck
    res.throughput_mbps = effBW * pathQ;

    return res;
}

// ================================================================
//  §10  TOPOLOGY BUILDERS
// ================================================================

std::vector<Satellite> buildIridium() {
    std::vector<Satellite> sats;
    sats.reserve(C::M * C::N);
    for (int m = 0; m < C::M; m++) {
        double raan = 2.0 * C::PI * m / C::M;
        for (int n = 0; n < C::N; n++) {
            // Orbital phase offset Δωf = π/(M·N) (polar constellation, eq.23)
            double phase    = C::PI * m / (C::M * C::N);
            double initAnom = 2.0 * C::PI * n / C::N + phase;
            Satellite s;
            s.sid      = "SID_" + std::to_string(m) + "_" + std::to_string(n);
            s.plane    = m;  s.idx = n;
            s.raan     = raan;
            s.initAnom = initAnom;
            s.update(0.0);
            sats.push_back(s);
        }
    }
    return sats;
}

// Virtual nodes are fixed at satellite positions at t = 0
std::vector<VirtualNode> buildVNs(const std::vector<Satellite>& sats) {
    std::vector<VirtualNode> vns;
    vns.reserve(sats.size());
    for (auto& s : sats) {
        VirtualNode vn;
        vn.vid   = "VID_" + std::to_string(s.plane) + "_" + std::to_string(s.idx);
        vn.plane = s.plane;  vn.idx = s.idx;
        vn.ctr   = s.pos;    // fixed at initial satellite position
        vns.push_back(vn);
    }
    return vns;
}

// Virtual link topology (grid+, no counter-rotating seam)
std::vector<VLink> buildVLinks(int M, int N) {
    std::vector<VLink> lks;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            VID cur = "VID_" + std::to_string(m) + "_" + std::to_string(n);
            // Intra-plane: wraps around within same plane
            VID nxt = "VID_" + std::to_string(m) + "_" + std::to_string((n+1)%N);
            lks.push_back({cur, nxt, true});
            // Inter-plane: adjacent planes only, NO seam (m == M-1 → m+1 = 0 is seam)
            if (m < M - 1) {
                VID nxtP = "VID_" + std::to_string(m+1) + "_" + std::to_string(n);
                lks.push_back({cur, nxtP, false});
            }
        }
    }
    return lks;
}

// Ground stations (lat, lon from paper's simulation section)
std::vector<GS> buildGS() {
    return {
        {"Seattle",  llToECEF( 47.6, -122.3), {}, 1e18},
        {"London",   llToECEF( 51.5,   -0.1), {}, 1e18},
        {"HongKong", llToECEF( 22.3,  114.2), {}, 1e18},
        {"NewYork",  llToECEF( 40.7,  -74.0), {}, 1e18},
    };
}

// Attach each GS to its nearest active virtual node
void attachGS(std::vector<GS>& gss, const std::vector<VirtualNode>& vns) {
    for (auto& gs : gss) {
        gs.nearVID = "";  gs.nearDist = 1e18;
        for (auto& vn : vns) {
            if (!vn.active) continue;
            double d = gs.pos.dist(vn.ctr);
            if (d < gs.nearDist) { gs.nearDist = d; gs.nearVID = vn.vid; }
        }
    }
}

// ================================================================
//  §11  TOPOLOGY METRICS  (eqs.5-7)
// ================================================================

struct TopoMetrics {
    int    fn = 0, fl = 0, Q = 0;
    double O = 0, Oprev = 0, dO = 0;

    void update(int fn_, int fl_) {
        fn = fn_;  fl = fl_;  Q = fn + fl;
        Oprev = O;  O = static_cast<double>(Q);
        dO = std::abs(O - Oprev);
    }

    void print(int t) const {
        std::cout << std::fixed << std::setprecision(0)
                  << "fn=" << fn << " fl=" << fl
                  << " Q(" << t << ")=" << Q
                  << " O(" << t << ")=" << O
                  << " ΔO(" << t << ")=" << dO << "\n";
    }
};

// ================================================================
//  §12  MAIN
// ================================================================

int main() {
    std::cout <<
        "═══════════════════════════════════════════════════════════\n"
        "  Dyna-STN: Dynamic Discrete Topology Routing (Complete)  \n"
        "  Li, Wu, Wang — IEEE/ACM ToN, Vol. 32, No. 5, 2024      \n"
        "═══════════════════════════════════════════════════════════\n\n";

    // ── Simulation parameters ─────────────────────────────────
    const int    T_SLOTS  = 10;      // number of time slots
    const double DT       = 600.0;   // slot duration, s (10 min)
    const double CUBE_R   = 2200.0;  // service cube radius, km

    // ── Build constellation ───────────────────────────────────
    auto sats   = buildIridium();
    auto vns    = buildVNs(sats);
    auto vlinks = buildVLinks(C::M, C::N);
    auto gss    = buildGS();

    // Initialise ISL reference length (used by ISL quality model)
    ISL::Lv_ref = Lv_km();
    double Lh   = Lh_km(0.0);

    // ── Print header info ─────────────────────────────────────
    std::cout << "Constellation : " << C::M << " planes × " << C::N
              << " sats/plane = " << C::M*C::N << " satellites\n"
              << "Virtual nodes : " << vns.size()
              << "  Virtual links: " << vlinks.size() << "\n"
              << "Lv(intra)     = " << std::fixed << std::setprecision(0)
              << ISL::Lv_ref << " km\n"
              << "Lh(inter,eq)  = " << Lh << " km\n"
              << "T_comm(max)   = " << std::setprecision(1)
              << maxCommTime_s() << " s  (eq.11, θ_min=10°)\n"
              << "U_SNR(Lv)     = " << std::setprecision(4)
              << ISL::uSNR(ISL::Lv_ref) << "  (calibrated, eq.22)\n"
              << "U_SNR(Lh)     = " << ISL::uSNR(Lh) << "\n"
              << "Service cube R= " << std::setprecision(0) << CUBE_R << " km\n\n";

    // ── Initialise subsystems ────────────────────────────────
    DDTMPlane   ddtm;
    OSPFSim     ospf;
    TopoMetrics topo;

    // Source-destination pairs from paper's simulation
    const std::vector<std::pair<std::string,std::string>> SD = {
        {"Seattle",  "NewYork"},    // ~3 865 km
        {"London",   "NewYork"},    // ~5 570 km
        {"HongKong", "NewYork"},    // ~12 947 km
    };

    // ── Accumulated stats ─────────────────────────────────────
    int    totalMigrations  = 0;
    double totalTm          = 0.0;
    double totalDowntime    = 0.0;

    // ═════════════════════════════════════════════════════════
    //  Main simulation loop
    // ═════════════════════════════════════════════════════════
    for (int t = 1; t <= T_SLOTS; t++) {
        double now = t * DT;   // absolute time, s

        std::cout << "\n┌─ Slot t=" << t
                  << "  (t=" << std::fixed << std::setprecision(0)
                  << now << " s) ─────────────────────────────────\n";

        // ── 1. Propagate satellites (Keplerian, correct inclination) ─
        for (auto& s : sats) s.update(now);

        // ── 2. Algorithm 1: SID/VID Binding Procedure ────────────
        ddtm.runSlot(t, now, vns, sats, vlinks, CUBE_R);
        topo.update(ddtm.fn_sum, ddtm.fl_sum);

        int activeVN = 0;
        for (auto& vn : vns) if (vn.active) activeVN++;

        std::cout << "│ [DDTM] ";  topo.print(t);
        std::cout << "│         Active VNs: " << activeVN
                  << "/" << vns.size() << "\n";

        // ── 3. Service migrations ─────────────────────────────────
        if (!ddtm.migrations.empty()) {
            totalMigrations += (int)ddtm.migrations.size();
            for (auto& m : ddtm.migrations) {
                totalTm       += m.Tm;
                totalDowntime += m.downtime;
            }
            std::cout << "│ [Migration] count=" << ddtm.migrations.size()
                      << "  avg Tm=" << std::setprecision(4)
                      << ddtm.avgTm() << " s"
                      << "  (tcp+tdp+tl)\n";
            // Show details for first migration
            if (!ddtm.migrations.empty()) {
                auto& m = ddtm.migrations[0];
                std::cout << "│   " << m.from << " → " << m.to
                          << " (" << m.vid << ")"
                          << " tcp=" << std::setprecision(3) << m.tcp
                          << " tdp=" << m.tdp
                          << " tl="  << m.tl
                          << " Tm="  << m.Tm << " s\n";
            }
        }

        // Algorithm 1 message log (first 4 messages)
        if (!ddtm.msgLog.empty()) {
            std::cout << "│ [Alg1 msgs] ";
            int shown = 0;
            for (auto& msg : ddtm.msgLog) {
                if (shown++ >= 4) { std::cout << "..."; break; }
                std::cout << msg << " | ";
            }
            std::cout << "\n";
        }

        // ── 4. OSPF routing update (Section III-B) ────────────────
        OSPFStats ospfSt = ospf.runEpoch(vns, vlinks, sats, now);

        // Push computed RIBs into virtual node router states
        for (auto& vn : vns) {
            if (!vn.active) continue;
            auto it = ospf.ribs.find(vn.vid);
            if (it != ospf.ribs.end()) {
                vn.router.rib = it->second;
                vn.router.updateFIBFromRIB();
                ddtm.vidRS[vn.vid] = vn.router;
            }
        }

        std::cout << "│ [OSPF] routeCalc=" << std::fixed
                  << std::setprecision(2) << ospfSt.routeCalc_ms
                  << " ms  converge=" << ospfSt.converge_ms
                  << " ms  LSA=" << ospfSt.lsaPkts
                  << " Hello=" << ospfSt.helloPkts
                  << " ΔLinks=" << ospfSt.changedLinks << "\n";

        // ── 5. Ground station attachment ──────────────────────────
        attachGS(gss, vns);

        // ── 6. Sample ISL quality for one intra + one inter link ──
        {
            bool doneIntra = false, doneInter = false;
            for (auto& lk : vlinks) {
                if (( lk.intra && doneIntra) ||
                    (!lk.intra && doneInter)) continue;

                const VirtualNode* a = nullptr, *b = nullptr;
                for (auto& vn : vns) {
                    if (vn.vid == lk.u) a = &vn;
                    if (vn.vid == lk.v) b = &vn;
                }
                if (!a || !b || !a->active || !b->active) continue;

                const Satellite* sA = nullptr, *sB = nullptr;
                for (auto& s : sats) {
                    if (s.sid == a->bound) sA = &s;
                    if (s.sid == b->bound) sB = &s;
                }
                if (!sA || !sB) continue;

                double L  = sA->pos.dist(sB->pos);
                double Us = ISL::uSNR(L);
                double Ul = ISL::uDur(lk.intra, sA->pos, sB->pos);
                double Ut = ISL::uTotal(L, lk.intra, sA->pos, sB->pos);

                std::cout << "│ [ISL-"
                          << (lk.intra ? "INTRA" : "INTER") << "] "
                          << lk.u << "↔" << lk.v
                          << " L=" << std::setprecision(0) << L
                          << " km  U_SNR=" << std::setprecision(4) << Us
                          << " U_dur=" << Ul
                          << " U_tot(eq28)=" << Ut << "\n";

                if (lk.intra) doneIntra = true; else doneInter = true;
                if (doneIntra && doneInter) break;
            }
        }

        // ── 7. Traffic simulation (Section IV-B) ─────────────────
        std::cout << "│ [Traffic]\n";
        std::cout << "│   " << std::left << std::setw(22) << "Flow"
                  << std::setw(12) << "Thput(Mbps)"
                  << std::setw(10) << "PLR(%)"
                  << std::setw(12) << "E2E(ms)"
                  << "Hops\n";
        std::cout << "│   " << std::string(58, '-') << "\n";

        for (auto& [src, dst] : SD) {
            GS* gsS = nullptr, *gsD = nullptr;
            for (auto& gs : gss) {
                if (gs.id == src) gsS = &gs;
                if (gs.id == dst) gsD = &gs;
            }
            if (!gsS || !gsD) continue;

            TrafficResult tr = simTraffic(*gsS, *gsD, ospf, vns, sats, vlinks);
            std::cout << "│   " << std::left << std::setw(22)
                      << (src + "→" + dst)
                      << std::setw(12) << std::fixed
                      << std::setprecision(1) << tr.throughput_mbps
                      << std::setw(10) << std::setprecision(2)
                      << tr.pktLoss_pct
                      << std::setw(12) << std::setprecision(1)
                      << tr.e2eDelay_ms
                      << tr.hops << "\n";
        }

        // ── 8. Sample bindings (first 6 VNs) ─────────────────────
        std::cout << "│ [Bindings] ";
        int shown = 0;
        for (auto& vn : vns) {
            if (shown++ >= 6) { std::cout << "..."; break; }
            std::cout << vn.vid << ":"
                      << (vn.active ? vn.bound : "none") << " ";
        }
        std::cout << "\n└────────────────────────────────────────────────\n";
    }

    // ═════════════════════════════════════════════════════════
    //  Final statistics
    // ═════════════════════════════════════════════════════════
    std::cout << "\n╔═══════════════════════════════════════════════╗\n"
              << "║              FINAL STATISTICS                  ║\n"
              << "╚═══════════════════════════════════════════════╝\n";

    int finalActive = 0;
    for (auto& vn : vns) if (vn.active) finalActive++;

    std::cout << "Active virtual nodes (last slot): "
              << finalActive << "/" << vns.size() << "\n"
              << "Q(last slot) = " << topo.Q
              << "  [fn=" << topo.fn << " fl=" << topo.fl << "]\n"
              << "ΔO(last slot) = " << std::fixed << std::setprecision(1)
              << topo.dO << "\n"
              << "Total migrations: " << totalMigrations << "\n";

    if (totalMigrations > 0) {
        std::cout << "Avg Tm per migration: "
                  << std::setprecision(4) << totalTm / totalMigrations
                  << " s\n"
                  << "Avg downtime/mig: "
                  << totalDowntime / totalMigrations << " s\n";
    }

    std::cout << "Cumulative OSPF LSA overhead: " << ospf.totalLSA
              << " packets\n";

    // Print final ground-station attachments
    std::cout << "\n[Ground Station → VID attachments at t_last]\n";
    for (auto& gs : gss)
        std::cout << "  " << std::left << std::setw(10) << gs.id
                  << " → " << gs.nearVID
                  << " (dist=" << std::fixed << std::setprecision(0)
                  << gs.nearDist << " km)\n";

    // Print SID/VID mapping table (last 2 slots)
    ddtm.mapSrv.print(2);

    // Print final router status of first few active VNs
    std::cout << "\n[Router States of first 4 active VNs]\n";
    int shown2 = 0;
    for (auto& vn : vns) {
        if (!vn.active || shown2++ >= 4) continue;
        std::cout << "  " << vn.vid
                  << "  status=" << rsStr(vn.router.status)
                  << "  RIB_size=" << vn.router.rib.size()
                  << "  FIB_size=" << vn.router.fib.size()
                  << "  totalDown=" << std::setprecision(4)
                  << vn.router.totalDown << " s\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}

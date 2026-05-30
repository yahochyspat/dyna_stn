
//g++ dyna_stn.cpp -std=c++17 -o dyna_stn.exe; .\dyna_stn.exe
//g++ back.cpp -std=c++17 -o back.exe; .\back.exe
//python plot_dyna_stn.py

#include <iostream>
#include <fstream>
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
#include <numeric>


namespace C {
    constexpr double PI      = 3.141592653589793;
    constexpr double Re      = 6371.0;               
    constexpr double h_orb   = 780.0;                
    constexpr double R_orb   = Re + h_orb;              
    constexpr double incl    = 86.4 * PI / 180.0;     
    constexpr int    M       = 6;                     
    constexpr int    N       = 11;                  
    constexpr int    T_orb   = 6900;                
    constexpr double polar_B = 70.0 * PI / 180.0;  
    
    constexpr double Pt_eff  = 3.0;   
    constexpr double ws      = 0.6;  
    constexpr double wl      = 0.4;  
    
    constexpr double theta_min = 10.0 * PI / 180.0; 

    constexpr double tcp_base      = 0.050;  
    constexpr double tcp_per_route = 0.0005;  
    constexpr double tdp_per_fib   = 0.00030;  
    constexpr double tl_per_link   = 0.010;  

    
    constexpr int    hello_ivl  = 10; 
    constexpr int    dead_ivl   = 40;  
    constexpr double c_light    = 3.0e5; 
    
    constexpr double dijk_coeff = 1.5e-6;   


    constexpr double bw_link   = 100.0; 
    constexpr double queue_ms  = 0.5; 
}


struct Vec3 {
    double x{0}, y{0}, z{0};
    double norm()      const { return std::sqrt(x*x+y*y+z*z); }
    double dist(const Vec3& o) const {
        double dx=x-o.x,dy=y-o.y,dz=z-o.z;
        return std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    double lat_rad() const { double r=norm(); return r>0?std::asin(z/r):0.0; }
    double lon_deg() const { return std::atan2(y,x)*180.0/C::PI; }
    double lat_deg() const { return lat_rad()*180.0/C::PI; }
    std::string str() const {
        std::ostringstream ss;
        ss<<std::fixed<<std::setprecision(0)
          <<"("<<x<<","<<y<<","<<z<<")";
        return ss.str();
    }
};


// формула не из статьи, это стандартное преобразование географических координат
// (широта/долгота) в ECEF. Нужна для привязки наземных станций к ближайшему VID.
// x = r*cos(lat)*cos(lon)
// y = r*cos(lat)*sin(lon)
// z = r*sin(lat)
Vec3 llToECEF(double la_d,double lo_d,double r=C::Re){
    double la=la_d*C::PI/180.0, lo=lo_d*C::PI/180.0;
    return {r*std::cos(la)*std::cos(lo),r*std::cos(la)*std::sin(lo),r*std::sin(la)};
}

double Lv_km(){ return std::sqrt(2.0)*C::R_orb*std::sqrt(1.0-std::cos(2.0*C::PI/C::N)); }
// В статье Lh зависит от широты lat и фазового фактора F.
// Здесь F=0 для Iridium
double Lh_km(double lat){ return std::sqrt(2.0)*C::R_orb*std::sqrt((1.0-std::cos(C::PI/C::N))*std::cos(lat)); }
double maxCommTime_s(){
    double theta=C::theta_min;
    double g=std::acos((C::Re/C::R_orb)*std::cos(theta))-theta;
    return 2.0*C::R_orb*g/std::sqrt(3.986e5/C::R_orb);
}



namespace ISL {
//U_SNR считается как 1-exp(-Pt_eff*(Lv/L)^2)
//сохраняется зависимость от расстояния, но не воспроизводится полная модель мощности, шума и усиления антенн.
static double Lv_ref = 0.0;
double uSNR(double L){
    if(L<=0||Lv_ref<=0) return 0.0;
    return 1.0-std::exp(-C::Pt_eff*(Lv_ref/L)*(Lv_ref/L));
}
double uDur(bool intra,const Vec3& a,const Vec3& b){
    if(intra) return 1.0;
    double lA=a.lat_rad(),lB=b.lat_rad();
    if(std::abs(lA)>=C::polar_B||std::abs(lB)>=C::polar_B) return 0.0;
    double ml=std::max(std::abs(lA),std::abs(lB));
    return (C::polar_B-ml)/C::polar_B;
}
double uTotal(double L,bool intra,const Vec3& a,const Vec3& b){
    double us=uSNR(L),ul=uDur(intra,a,b);
    if(ul<=0.0) return 0.0;
    return std::pow(us,C::ws)*std::pow(ul,C::wl);
}
}

using SID=std::string;
using VID=std::string;

enum class RouterStatus{RUNNING,SUSPENDED,MIGRATING};
inline std::string rsStr(RouterStatus s){
    switch(s){
        case RouterStatus::RUNNING:   return "RUNNING";
        case RouterStatus::SUSPENDED: return "SUSPENDED";
        case RouterStatus::MIGRATING: return "MIGRATING";
    } return "?";
}

struct RouterConfig{
    int    helloIvl=C::hello_ivl;
    int    deadIvl =C::dead_ivl;
    double linkBW  =C::bw_link;
    int    area    =0;
    int    lsaSeq  =0;
};
struct RIBEntry{ VID nextHop; double metric; };
struct FIBEntry{ VID outIface; double cost;  };
struct RouterState{
    RouterStatus             status    =RouterStatus::RUNNING;
    std::map<VID,RIBEntry>  rib;
    std::map<VID,FIBEntry>  fib;
    RouterConfig             cfg;
    double                   downStart =-1.0;
    double                   totalDown = 0.0;

    void suspend(double now){
        if(status!=RouterStatus::SUSPENDED){status=RouterStatus::SUSPENDED;downStart=now;}
    }
    void restart(double now){
        if(downStart>=0.0){totalDown+=now-downStart;downStart=-1.0;}
        status=RouterStatus::RUNNING;
    }
    void updateFIBFromRIB(){
        fib.clear();
        for(auto&[d,e]:rib) fib[d]={e.nextHop,e.metric};
    }
};

struct Satellite{
    SID    sid;
    int    plane,idx;
    double raan,initAnom;
    Vec3   pos;
    void update(double t){
        double om=2.0*C::PI/C::T_orb,an=initAnom+om*t,in=C::incl,R=C::R_orb;
        pos.x=R*(std::cos(raan)*std::cos(an)-std::sin(raan)*std::sin(an)*std::cos(in));
        pos.y=R*(std::sin(raan)*std::cos(an)+std::cos(raan)*std::sin(an)*std::cos(in));
        pos.z=R*std::sin(an)*std::sin(in);
    }
    bool inPolar() const{return std::abs(pos.lat_rad())>=C::polar_B;}
};

struct VLink{ VID u,v; bool intra; };

struct VirtualNode{
    VID         vid;
    int         plane,idx;
    Vec3        ctr;
    bool        active   =false;
    SID         bound    ="";         // привязка после Complete
    SID         pendingBound="";      // после Response
    SID         prevBound="";
    RouterState router;
};

struct GS{
    std::string id;
    Vec3        pos;
    VID         nearVID;
    double      nearDist=1e18;
};


struct LSAEntry{
    VID  origin; int seq=0;
    std::vector<std::pair<VID,double>> links;
    bool sameLinks(const LSAEntry&o)const{return links==o.links;}
};
using LSDB=std::map<VID,LSAEntry>;
using RIBMap=std::map<VID,RIBEntry>;

RIBMap dijkstra(const LSDB& lsdb,const VID& src){
    constexpr double INF=std::numeric_limits<double>::infinity();
    std::map<VID,double> dist; std::map<VID,VID> prev;
    for(auto&[v,_]:lsdb){dist[v]=INF;prev[v]="";}
    dist[src]=0.0;
    using PQE=std::pair<double,VID>;
    std::priority_queue<PQE,std::vector<PQE>,std::greater<PQE>> pq;
    pq.push({0.0,src});
    while(!pq.empty()){
        auto[d,u]=pq.top();pq.pop();
        if(d>dist[u]) continue;
        auto it=lsdb.find(u); if(it==lsdb.end()) continue;
        for(auto&[v,w]:it->second.links){
            if(w<=0||w>=INF) continue;
            double nd=dist[u]+w;
            if(nd<dist[v]){dist[v]=nd;prev[v]=u;pq.push({nd,v});}
        }
    }
    RIBMap rib;
    for(auto&[dst,_]:lsdb){
        if(dst==src) continue;
        if(dist[dst]>=INF){rib[dst]={"UNREACHABLE",INF};continue;}
        VID cur=dst;
        while(!prev[cur].empty()&&prev[cur]!=src) cur=prev[cur];
        rib[dst]={cur,dist[dst]};
    }
    return rib;
}

std::vector<VID> getPath(const std::map<VID,RIBMap>& ribs,const VID& src,const VID& dst){
    std::vector<VID> path;
    if(!ribs.count(src)) return path;
    path.push_back(src);
    VID cur=src; std::set<VID> vis;
    while(cur!=dst&&!vis.count(cur)){
        vis.insert(cur);
        auto ri=ribs.find(cur); if(ri==ribs.end()) break;
        auto di=ri->second.find(dst);
        if(di==ri->second.end()||di->second.nextHop=="UNREACHABLE") break;
        cur=di->second.nextHop; path.push_back(cur);
    }
    return path;
}

class ServiceCubeModel {
public:
    // service cube задаётся как ближайший фиксированный центр VN в 3D Voronoi-разбиении
    VID findOwner(const Satellite& sat,
                  const std::vector<VirtualNode>& vns,
                  double /*t_unused*/) const
    {
        double minD = std::numeric_limits<double>::max();
        VID    owner;
        for(const auto& vn : vns){
            double d = sat.pos.dist(vn.ctr);
            if(d < minD){ minD = d; owner = vn.vid; }
        }
        return owner;
    }

    bool isNewEntry(const Satellite& sat,
                    const std::string& prevOwner,
                    const std::vector<VirtualNode>& vns,
                    double t) const
    { return findOwner(sat,vns,t) != prevOwner; }

    struct ServiceCell {
        double lat_ctr,  lon_ctr;      
        double along_km, along_deg;
        double cross_km, cross_deg;
        
        double lat_N, lat_S, lon_E, lon_W;
    };

    ServiceCell getServiceCell(const VirtualNode& vn) const {
        ServiceCell sc;
        sc.lat_ctr = vn.ctr.lat_deg();
        sc.lon_ctr = vn.ctr.lon_deg();

    
        double Lv        = Lv_km();
        double scaleOrb  = C::Re / C::R_orb;        
        sc.along_km  = (Lv / 2.0) * scaleOrb;         
        sc.along_deg = sc.along_km / C::Re * 180.0 / C::PI;

   
        double latR   = vn.ctr.lat_rad();
        double Lh     = Lh_km(latR);
        sc.cross_km   = (Lh / 2.0) * scaleOrb;
        double cosLat = std::max(0.01, std::cos(latR));
        sc.cross_deg  = sc.cross_km / (C::Re * cosLat) * 180.0 / C::PI;

        sc.lat_N = sc.lat_ctr + sc.along_deg;
        sc.lat_S = sc.lat_ctr - sc.along_deg;
        sc.lon_E = sc.lon_ctr + sc.cross_deg;
        sc.lon_W = sc.lon_ctr - sc.cross_deg;
        return sc;
    }

    
    double angDist_pub(const Satellite& sat, const VID& vid,
                       const std::vector<VirtualNode>& vns,
                       double /*t*/) const {
        
        for(const auto& vn : vns)
            if(vn.vid == vid) return sat.pos.dist(vn.ctr);
        return std::numeric_limits<double>::max();
    }

    static double normAnom(double a){
        a = std::fmod(a, 2.0*C::PI);
        return a < 0 ? a + 2.0*C::PI : a;
    }
    static double angDist(double a, double b){
        double d = std::abs(a-b);
        return std::min(d, 2.0*C::PI - d);
    }
};

struct OverheadModel {
    // Формула C(Q) не задана в статье явно, там O(t)=C(Q(t)).
    // C(Q) = βbase*Q + βcong*Q^2/Qsat + βchange*|deltaQ|
    //Линейный член моделирует базовую стоимость каждой SID/VID или link-привязки, 
    //квадратичный член - рост нагрузки управляющей плоскости при большом числе привязок, 
    //|deltaQ| - дополнительную стоимость резких изменений между временными слотами.
    double beta_base   = 1.0;    
    double beta_cong   = 0.5;      
    double beta_change = 2.0;      
    int    Q_sat       = C::M * C::N;  

    double compute(int Q, int deltaQ) const {
        return beta_base * Q
             + beta_cong * (double)Q * Q / Q_sat
             + beta_change * std::abs(deltaQ);
    }
};


struct OSPFStats{
    double routeCalc_ms=0, converge_ms=0;
    int    lsaPkts=0, helloPkts=0, changedLinks=0;
};

class OSPFSim {
public:
    LSDB                  lsdb;
    std::map<VID,int>     seqNum;
    std::map<VID,RIBMap>  ribs;
    int                   totalLSA=0;

    OSPFStats runEpoch(const std::vector<VirtualNode>& vns,
                       const std::vector<VLink>&       links,
                       const std::vector<Satellite>&   sats,
                       double)
    {
        OSPFStats st; LSDB newLSDB; int V=0;
        for(const auto& vn:vns){
            if(!vn.active) continue; V++;
            LSAEntry e; e.origin=vn.vid; e.seq=seqNum[vn.vid]++;
            const Satellite* sA=fSat(sats,vn.bound);
            for(const auto& lk:links){
                VID nb; bool intra=lk.intra;
                if(lk.u==vn.vid) nb=lk.v;
                else if(lk.v==vn.vid) nb=lk.u;
                else continue;
                const VirtualNode* vnB=fVN(vns,nb);
                if(!vnB||!vnB->active) continue;
                if(!intra&&isSeam(vn.plane,vnB->plane)) continue;
                const Satellite* sB=fSat(sats,vnB->bound);
                if(!sA||!sB) continue;
                double L=sA->pos.dist(sB->pos);
                double U=ISL::uTotal(L,intra,sA->pos,sB->pos);
                if(U<=0.0) continue;
                e.links.push_back({nb,1.0/U});
            }
            newLSDB[vn.vid]=e;
        }
        for(auto&[vid,ne]:newLSDB){
            auto it=lsdb.find(vid);
            if(it==lsdb.end()){
                st.changedLinks++;         
            } else {

                std::set<VID> oldNb, newNb;
                for(auto&[nb,_]:it->second.links) oldNb.insert(nb);
                for(auto&[nb,_]:ne.links)          newNb.insert(nb);
                if(oldNb!=newNb) st.changedLinks++;
            }
        }
        st.helloPkts=V*4;
        st.lsaPkts=st.changedLinks*V; totalLSA+=st.lsaPkts;
        int E=0; for(auto&[v,e]:newLSDB) E+=(int)e.links.size(); E/=2;
        int diam=C::M+C::N/2;
        double hopMs=2034.0/C::c_light*1000.0;
        st.converge_ms=diam*hopMs+(double)V*(V+E)*std::log2(V+1.0)*C::dijk_coeff*1000.0;
        st.routeCalc_ms=(double)V*(V+E)*std::log2(V+1.0)*C::dijk_coeff*1000.0;
        lsdb=newLSDB;
        ribs.clear();
        for(const auto& vn:vns) if(vn.active) ribs[vn.vid]=dijkstra(lsdb,vn.vid);
        return st;
    }
    std::pair<VID,double> route(const VID& s,const VID& d)const{
        auto it=ribs.find(s); if(it==ribs.end()) return{"NO_RIB",-1};
        auto jt=it->second.find(d); if(jt==it->second.end()) return{"NOT_FOUND",-1};
        return{jt->second.nextHop,jt->second.metric};
    }
private:
    const VirtualNode* fVN(const std::vector<VirtualNode>& v,const VID& id)const{
        for(auto& n:v) if(n.vid==id) return &n; return nullptr;}
    const Satellite* fSat(const std::vector<Satellite>& s,const SID& id)const{
        for(auto& n:s) if(n.sid==id) return &n; return nullptr;}
    bool isSeam(int a,int b)const{return(a==0&&b==C::M-1)||(a==C::M-1&&b==0);}
};


struct MigRes{
    SID    from,to; VID vid;
    double tcp,tdp,tl,Tm,downtime;
};
// Упрощение: Tm=tcp+tdp+tl как в статье но
// tcp/tdp/tl вычисляются по числу RIB/FIB-записей и линков,
// а не измеряются на Ryu/Open vSwitch/VMware-прототипе
MigRes performMigration(const SID& from,const SID& to,const VID& vid,
                         const RouterState& rs,int nLinks,double){
    int nR=rs.rib.size(), nF=std::max((int)rs.fib.size(),nR);
    MigRes m; m.from=from;m.to=to;m.vid=vid;
    m.tcp=C::tcp_base+nR*C::tcp_per_route;
    m.tdp=nF*C::tdp_per_fib;
    m.tl =nLinks*C::tl_per_link;
    m.Tm=m.tcp+m.tdp+m.tl;
    m.downtime=m.tcp/3.0;
    return m;
}


struct MappingEntry{int slot;SID sid;VID vid;bool active;Vec3 loc;double ts;};
struct InquireResult{bool found=false,active=false;VID vid;SID oldSID;};

class MapServer {
public:
    std::map<int,std::vector<MappingEntry>> history;
    std::map<VID,SID> curVIDtoSID;
    std::map<SID,VID> curSIDtoVID;
    std::map<VID,bool> actv;

    void bind(int slot,const SID& sid,const VID& vid,const Vec3& loc,double now){
        auto ov=curSIDtoVID.find(sid);
        if(ov!=curSIDtoVID.end()&&ov->second!=vid){
            curVIDtoSID.erase(ov->second); actv[ov->second]=false;
        }
        curSIDtoVID[sid]=vid; curVIDtoSID[vid]=sid; actv[vid]=true;
        history[slot].push_back({slot,sid,vid,true,loc,now});
    }
    void unbindSID(const SID& sid){
        auto it=curSIDtoVID.find(sid); if(it==curSIDtoVID.end()) return;
        VID v=it->second; curSIDtoVID.erase(it); curVIDtoSID.erase(v); actv[v]=false;
    }
    SID getCurrentSID(const VID& vid)const{
        auto it=curVIDtoSID.find(vid); return it==curVIDtoSID.end()?SID():it->second;
    }
    bool isActive(const VID& vid)const{
        auto it=actv.find(vid); return it!=actv.end()&&it->second;
    }
    InquireResult inquire(const SID& sid,const Vec3& loc,int slot,
                          const std::vector<VirtualNode>& vns,
                          const ServiceCubeModel& scm,double t)const{
        (void)sid;(void)slot;
        
        InquireResult r;
        double best=1e18;
        for(auto& vn:vns){
            double d=vn.ctr.dist(loc);
            if(d<best){best=d;r.vid=vn.vid;}
        }
        (void)scm;(void)t;
        if(r.vid.empty()){r.found=false;return r;}
        r.found=true; r.oldSID=getCurrentSID(r.vid); r.active=!r.oldSID.empty();
        return r;
    }
    void print(int mx=2)const{
        std::cout<<"\n== SID/VID Mapping Table (last "<<mx<<" slots) ==\n";
        std::cout<<std::left<<std::setw(6)<<"Slot"<<std::setw(16)<<"SID"
                 <<std::setw(16)<<"VID"<<"Active\n"<<std::string(46,'-')<<"\n";
        int c=0;
        for(auto rit=history.rbegin();rit!=history.rend()&&c<mx;++rit,++c)
            for(auto& e:rit->second)
                std::cout<<std::setw(6)<<e.slot<<std::setw(16)<<e.sid
                         <<std::setw(16)<<e.vid<<(e.active?"1":"0")<<"\n";
    }
};

enum class MsgType{
    Inquire_FromSID,Response_FromMapSrv,Shift_FromMapSrv,
    Request_FromSID,Confirm_FromSID,Migrate_FromSID,Complete_FromSID
};
inline std::string msgStr(MsgType t){
    switch(t){
        case MsgType::Inquire_FromSID:     return "Inquire";
        case MsgType::Response_FromMapSrv: return "Response";
        case MsgType::Shift_FromMapSrv:    return "Shift";
        case MsgType::Request_FromSID:     return "Request";
        case MsgType::Confirm_FromSID:     return "Confirm";
        case MsgType::Migrate_FromSID:     return "Migrate";
        case MsgType::Complete_FromSID:    return "Complete";
    } return "?";
}
struct Msg{
    MsgType type; std::string from,to;
    VID vid; SID oldSID,newSID;
    RouterState rs;
    double createdAt=0; bool ok=true; std::string note;
};
enum class Alg1Status{OK,IGNORED,FAILED_NO_VID,FAILED_NO_OLD,FAILED_NO_NEW,FAILED_NO_RS};
inline std::string a1Str(Alg1Status s){
    switch(s){
        case Alg1Status::OK:             return "OK";
        case Alg1Status::IGNORED:        return "IGNORED";
        case Alg1Status::FAILED_NO_VID:  return "FAIL_NO_VID";
        case Alg1Status::FAILED_NO_OLD:  return "FAIL_NO_OLD";
        case Alg1Status::FAILED_NO_NEW:  return "FAIL_NO_NEW";
        case Alg1Status::FAILED_NO_RS:   return "FAIL_NO_RS";
    } return "?";
}
struct Alg1Trace{
    int slot; double time; MsgType type;
    SID from,to,oldSID,newSID; VID vid;
    Alg1Status status; std::string note;
};
struct Agent{
    SID sid,curVID,lastCubeVID,targetVID;
    RouterState router;
};

class DDTMPlane {
public:
    MapServer              mapSrv;
    ServiceCubeModel       scm;
    std::map<SID,Agent>    agents;
    std::map<VID,RouterState> vidRS;
    std::map<VID,SID>      vidCurSID;
    std::vector<MigRes>    migrations;
    std::vector<Alg1Trace> trace;
    int fn_sum=0,fl_sum=0;

    
    int runSlot(int slot,double now,
                std::vector<VirtualNode>& vns,
                const std::vector<Satellite>& sats,
                const std::vector<VLink>& vlinks,
                double /*cubeR_unused*/)
    {
        fn_sum=fl_sum=0; migrations.clear(); trace.clear();

        std::map<VID,std::pair<SID,double>> bestForVN; 

        for(const auto& sat:sats){
            Agent& ag=agents[sat.sid]; ag.sid=sat.sid;
            VID owner=scm.findOwner(sat,vns,now);
            if(owner.empty()) continue;
            double d=scm.angDist_pub(sat,owner,vns,now);
            auto it=bestForVN.find(owner);
            if(it==bestForVN.end()||d<it->second.second)
                bestForVN[owner]={sat.sid,d};
        }

        
        std::vector<Msg> Q;
        for(auto&[vid,pair]:bestForVN){
            const SID& newSID=pair.first;
            Agent& ag=agents[newSID]; ag.sid=newSID;
            if(ag.lastCubeVID==vid) continue;  
            ag.lastCubeVID=vid; ag.targetVID=vid;
            SID oldSID=mapSrv.getCurrentSID(vid);
            Msg m; m.type=MsgType::Inquire_FromSID;
            m.from=newSID; m.to="MapSrv"; m.vid=vid;
            m.oldSID=oldSID; m.newSID=newSID; m.createdAt=now;
            m.note=oldSID.empty()?"free":"occupied by "+oldSID;
            Q.push_back(m);
        }

        for(size_t i=0;i<Q.size();i++){
            const Msg msg=Q[i];
            bindingProc(msg,Q,slot,now,vns,sats);
        }

        for(auto& vn:vns){
            vn.prevBound=vn.bound;

            if(!vn.bound.empty()||!vn.pendingBound.empty())
                vn.active=true;
        }

        vidCurSID.clear();
        for(auto& vn:vns)
            if(!vn.bound.empty()) vidCurSID[vn.vid]=vn.bound;

        for(auto& vn:vns) if(vn.active&&!vn.bound.empty()) fn_sum++;

        std::set<std::pair<VID,VID>> seen;
        for(auto& lk:vlinks){
            auto key=std::make_pair(std::min(lk.u,lk.v),std::max(lk.u,lk.v));
            if(seen.count(key)) continue;
            const VirtualNode* a=fVN(vns,lk.u); const VirtualNode* b=fVN(vns,lk.v);
            if(!a||!b||a->bound.empty()||b->bound.empty()) continue;
            const Satellite* sA=fSat(sats,a->bound); const Satellite* sB=fSat(sats,b->bound);
            if(!sA||!sB) continue;
            if(!lk.intra&&(sA->inPolar()||sB->inPolar())) continue;
            if(!lk.intra&&isSeam(a->plane,b->plane)) continue;
            seen.insert(key); fl_sum++;
        }
        return fn_sum+fl_sum;
    }

    double avgTm()const{
        if(migrations.empty()) return 0.0;
        double s=0; for(auto& m:migrations) s+=m.Tm;
        return s/migrations.size();
    }

    void printTrace(const VID& vid)const{
        std::cout<<"\n[Alg1 trace for "<<vid<<"]\n";
        for(auto& e:trace){
            if(e.vid!=vid) continue;
            std::cout<<"  t="<<std::fixed<<std::setprecision(3)<<e.time
                     <<"  "<<msgStr(e.type)
                     <<"  "<<e.from<<"->"<<e.to
                     <<"  old="<<(e.oldSID.empty() ? "none" : e.oldSID)
                     <<"  new="<<(e.newSID.empty() ? "none" : e.newSID)
                     <<"  "<<a1Str(e.status);
            if(!e.note.empty()) std::cout<<"  "<<e.note;
            std::cout<<"\n";
        }
    }

private:
    void bindingProc(const Msg& msg,std::vector<Msg>& Q,
                     int slot,double now,
                     std::vector<VirtualNode>& vns,
                     const std::vector<Satellite>& sats)
    {
        trace.push_back({slot,now,msg.type,msg.from,msg.to,
                         msg.oldSID,msg.newSID,msg.vid,Alg1Status::OK,msg.note});

        switch(msg.type){

        case MsgType::Inquire_FromSID:{
            const SID nSID=msg.newSID; const SID oSID=msg.oldSID; const VID vid=msg.vid;

            Msg rsp; rsp.type=MsgType::Response_FromMapSrv;
            rsp.from="MapSrv";rsp.to=nSID;rsp.vid=vid;
            rsp.oldSID=oSID;rsp.newSID=nSID;rsp.createdAt=now;
            rsp.note="MapSrv response";
            Q.push_back(rsp);
            if(!oSID.empty()&&oSID!=nSID){
                Msg sh; sh.type=MsgType::Shift_FromMapSrv;
                sh.from="MapSrv";sh.to=oSID;sh.vid=vid;
                sh.oldSID=oSID;sh.newSID=nSID;sh.createdAt=now;
                sh.note="migrate to "+nSID;
                Q.push_back(sh);
            }
            break;
        }

        case MsgType::Response_FromMapSrv:{
            //В статье: для свободного VID нет старого SID, поэтому цепочка Shift→Request→Confirm→Migrate→Complete не запускается
            //Здесь я считаю такую первичную привязку сразу финальной
            const SID nSID=msg.newSID.empty()?msg.to:msg.newSID;
            const VID vid=msg.vid;
            Agent& ag=agents[nSID]; ag.sid=nSID; ag.curVID=vid;
            const Satellite* sat=fSat(sats,nSID);
            Vec3 loc=sat?sat->pos:Vec3{};
            mapSrv.bind(slot,nSID,vid,loc,now);
            VirtualNode* vn=fVNm(vns,vid);
            if(vn){
                vn->pendingBound=nSID;         
                if(msg.oldSID.empty()){        
                    vn->bound  = nSID;
                    vn->active = true;
                }
            }
            break;
        }

        case MsgType::Shift_FromMapSrv:{
            const SID oSID=msg.oldSID; const SID nSID=msg.newSID; const VID vid=msg.vid;
            if(oSID.empty()){addT(slot,now,msg,Alg1Status::FAILED_NO_OLD,"no oldSID");break;}
            if(nSID.empty()){addT(slot,now,msg,Alg1Status::FAILED_NO_NEW,"no newSID");break;}
            RouterState rs=vidRS.count(vid)?vidRS[vid]:
                           (agents.count(oSID)?agents[oSID].router:RouterState{});
    
            Msg req; req.type=MsgType::Request_FromSID;
            req.from=oSID;req.to=nSID;req.vid=vid;
            req.oldSID=oSID;req.newSID=nSID;req.rs=rs;req.createdAt=now;
            req.note="request suspend+confirm";
            Q.push_back(req);
            break;
        }

        case MsgType::Request_FromSID:{
            const SID oSID=msg.oldSID; const SID nSID=msg.newSID; const VID vid=msg.vid;
            if(nSID.empty()){addT(slot,now,msg,Alg1Status::FAILED_NO_NEW,"no newSID");break;}

            agents[nSID].router.suspend(now);
            Msg conf; conf.type=MsgType::Confirm_FromSID;
            conf.from=nSID;conf.to=oSID;conf.vid=vid;
            conf.oldSID=oSID;conf.newSID=nSID;conf.createdAt=now;
            conf.note="new SID suspended";
            Q.push_back(conf);
            break;
        }

        case MsgType::Confirm_FromSID:{
            const SID oSID=msg.oldSID; const SID nSID=msg.newSID; const VID vid=msg.vid;
            RouterState rs=vidRS.count(vid)?vidRS[vid]:
                           (agents.count(oSID)?agents[oSID].router:RouterState{});
         
            MigRes mr=performMigration(oSID,nSID,vid,rs,4,now);
            migrations.push_back(mr);
            Msg mig; mig.type=MsgType::Migrate_FromSID;
            mig.from=oSID;mig.to=nSID;mig.vid=vid;
            mig.oldSID=oSID;mig.newSID=nSID;mig.rs=rs;mig.createdAt=now+mr.downtime;
            mig.note="tcp="+std::to_string(mr.tcp).substr(0,5)+"s";
            Q.push_back(mig);
            break;
        }

        case MsgType::Migrate_FromSID:{
        
            const SID oSID=msg.oldSID; const SID nSID=msg.newSID; const VID vid=msg.vid;
            Agent& ag=agents[nSID]; ag.sid=nSID; ag.curVID=vid;
            ag.router=msg.rs; ag.router.updateFIBFromRIB(); ag.router.restart(now);
            vidRS[vid]=ag.router;
            Msg done; done.type=MsgType::Complete_FromSID;
            done.from=nSID;done.to=oSID;done.vid=vid;
            done.oldSID=oSID;done.newSID=nSID;done.createdAt=now;
            done.note="migration complete";
            Q.push_back(done);
            break;
        }

        case MsgType::Complete_FromSID:{
            const SID oSID=msg.oldSID; const SID nSID=msg.newSID; const VID vid=msg.vid;
            VirtualNode* vn=fVNm(vns,vid);
            if(vn){
                vn->bound=nSID;        
                vn->active=true;
            }

            if(agents.count(oSID)&&agents[oSID].curVID==vid)
                agents[oSID].curVID.clear();
            if(mapSrv.curSIDtoVID.count(oSID)&&mapSrv.curSIDtoVID[oSID]==vid)
                mapSrv.curSIDtoVID.erase(oSID);
            break;
        }

        } 
    }

    void addT(int sl,double t,const Msg& m,Alg1Status s,const std::string& n){
        trace.push_back({sl,t,m.type,m.from,m.to,m.oldSID,m.newSID,m.vid,s,n});
    }
    const VirtualNode* fVN(const std::vector<VirtualNode>& v,const VID& id)const{
        for(auto& n:v) if(n.vid==id) return &n; return nullptr;}
    VirtualNode* fVNm(std::vector<VirtualNode>& v,const VID& id)const{
        for(auto& n:v) if(n.vid==id) return &n; return nullptr;}
    const Satellite* fSat(const std::vector<Satellite>& s,const SID& id)const{
        for(auto& n:s) if(n.sid==id) return &n; return nullptr;}
    bool isSeam(int a,int b)const{return(a==0&&b==C::M-1)||(a==C::M-1&&b==0);}
};


struct UnitTestResult {
    bool passed;
    std::string details;
};

UnitTestResult testAlg1MsgSequence() {

    const VID testVID = "VID_TEST";
    const SID oldSID  = "SID_OLD";
    const SID newSID  = "SID_NEW";

    std::vector<VirtualNode> vns;
    VirtualNode vn0;
    vn0.vid=testVID; vn0.plane=0; vn0.idx=0;
    vn0.ctr={C::R_orb,0,0};  
    vns.push_back(vn0);

    std::vector<Satellite> sats;

    Satellite sNew;
    sNew.sid="SID_NEW"; sNew.plane=0; sNew.idx=10;
    sNew.raan=0; sNew.initAnom=0.0;  
    sNew.update(0.0); sats.push_back(sNew);
    Satellite sOld;
    sOld.sid="SID_OLD"; sOld.plane=0; sOld.idx=1;
    sOld.raan=0; sOld.initAnom=2*C::PI/C::N;
    sOld.update(0.0); sats.push_back(sOld);

    DDTMPlane ddtm;

    ddtm.mapSrv.bind(0, oldSID, testVID, sOld.pos, 0.0);
    Agent& agOld=ddtm.agents[oldSID];
    agOld.sid=oldSID; agOld.curVID=testVID; agOld.lastCubeVID=testVID;
    for(int i=0;i<5;i++) agOld.router.rib["VID_"+std::to_string(i)]={"next",1.0};
    ddtm.vidRS[testVID]=agOld.router;

    Agent& agNew=ddtm.agents[newSID];
    agNew.sid=newSID; agNew.lastCubeVID="";

    bool initOldBound=(ddtm.mapSrv.getCurrentSID(testVID)==oldSID);
    bool initNewFree =(agNew.curVID!=testVID);

    std::vector<VLink> vlinks;
    ddtm.runSlot(1, 600.0, vns, sats, vlinks, 2200.0);


    const std::vector<MsgType> expected = {
        MsgType::Inquire_FromSID,
        MsgType::Response_FromMapSrv,
        MsgType::Shift_FromMapSrv,
        MsgType::Request_FromSID,
        MsgType::Confirm_FromSID,
        MsgType::Migrate_FromSID,
        MsgType::Complete_FromSID
    };


    std::vector<MsgType> actual;
    for(auto& e:ddtm.trace)
        if(e.vid==testVID) actual.push_back(e.type);

    bool seqOK = (actual.size() >= expected.size());
    if(seqOK) for(size_t i=0;i<expected.size();i++) seqOK&=(actual[i]==expected[i]);

    bool finalBoundOK=false, oldFreedOK=false;
    for(auto& vn:vns){
        if(vn.vid==testVID){
            finalBoundOK = (vn.bound==newSID);
            break;
        }
    }
    oldFreedOK = (ddtm.agents[oldSID].curVID!=testVID);


    std::ostringstream rpt;
    rpt << "  Initial state: oldBound=" << initOldBound
        << " newFree="    << initNewFree << "\n";
    rpt << "  Actual msg sequence: ";
    for(auto t:actual) rpt<<msgStr(t)<<" -> ";
    rpt << "\n";
    rpt << "  Expected sequence match: " << (seqOK?"PASS":"FAIL") << "\n";
    rpt << "  Final binding (after Complete) correct: "
        << (finalBoundOK?"PASS":"FAIL") << "\n";
    rpt << "  Old SID freed after Complete: "
        << (oldFreedOK?"PASS":"FAIL") << "\n";

    bool all = seqOK && finalBoundOK && oldFreedOK && initOldBound && initNewFree;
    rpt << "  OVERALL: " << (all?"PASS":"FAIL") << "\n";
    return {all, rpt.str()};
}


// В статье задается внешним STK-сценарием, я задаю круговой орбитой
std::vector<Satellite> buildIridium(){
    std::vector<Satellite> s; s.reserve(C::M*C::N);
    for(int m=0;m<C::M;m++){
        double raan=2*C::PI*m/C::M;
        for(int n=0;n<C::N;n++){
            double ph=C::PI*m/(C::M*C::N),an=2*C::PI*n/C::N+ph;
            Satellite sat; sat.sid="SID_"+std::to_string(m)+"_"+std::to_string(n);
            sat.plane=m;sat.idx=n;sat.raan=raan;sat.initAnom=an;sat.update(0.0);
            s.push_back(sat);
        }
    } return s;
}
std::vector<VirtualNode> buildVNs(const std::vector<Satellite>& sats){
    std::vector<VirtualNode> v; v.reserve(sats.size());
    for(auto& s:sats){
        VirtualNode vn; vn.vid="VID_"+std::to_string(s.plane)+"_"+std::to_string(s.idx);
        vn.plane=s.plane;vn.idx=s.idx;vn.ctr=s.pos;
        v.push_back(vn);
    } return v;
}
std::vector<VLink> buildVLinks(int M,int N){
    std::vector<VLink> lks;
    for(int m=0;m<M;m++) for(int n=0;n<N;n++){
        VID cur="VID_"+std::to_string(m)+"_"+std::to_string(n);
        VID nxt="VID_"+std::to_string(m)+"_"+std::to_string((n+1)%N);
        lks.push_back({cur,nxt,true});
        if(m<M-1){
            VID nxp="VID_"+std::to_string(m+1)+"_"+std::to_string(n);
            lks.push_back({cur,nxp,false});
        }
    } return lks;
}


//Использую nearVID(gs) = argmin_v ||pos(gs) - ctr(v)||
//pos(gs) - координаты наземной станции,
//ctr(v)  - фиксированный центр виртуального узла,
//||*||  - евклидово расстояние в 3D.
void attachGS(std::vector<GS>& gss,const std::vector<VirtualNode>& vns){
    for(auto& gs:gss){gs.nearVID="";gs.nearDist=1e18;
        for(auto& vn:vns){
            if(vn.bound.empty()) continue;
            double d=gs.pos.dist(vn.ctr);
            if(d<gs.nearDist){gs.nearDist=d;gs.nearVID=vn.vid;}
        }
    }
}


struct TopoMetrics{
    // Q(t) из eq.(5), но O(t)=C(Q(t))
    int    fn=0,fl=0,Q=0,Qprev=0;
    double O=0,Oprev=0,dO=0;
    OverheadModel ovm;

    void update(int fn_,int fl_){
        fn=fn_;fl=fl_;Qprev=Q;Q=fn+fl;
        Oprev=O;
        O=ovm.compute(Q, Q-Qprev);
        dO=std::abs(O-Oprev);
    }

    void print(int t)const{
        std::cout<<std::fixed<<std::setprecision(1)
                 <<"fn="<<fn<<" fl="<<fl
                 <<" Q("<<t<<")="<<Q
                 <<" C(Q)="<<O          
                 <<" delta0("<<t<<")="<<dO<<"\n";
    }
};

struct SlotRecord {
    int slot;
    double time_s;

    int activeVNs;
    int fn;
    int fl;
    int Q;

    double C_Q;
    double dO;

    int migrations;
    double avgTm;
    double avgDowntime;

    double routeCalc_ms;
    double converge_ms;
    int lsaPkts;
    int helloPkts;
    int changedLinks;

    static constexpr double NaN_SENTINEL = -1.0;
};

class CSVWriter {
public:
    explicit CSVWriter(const std::string& fname)
        : ofs_(fname) {
        if (!ofs_) {
            std::cerr << "[CSV] Cannot open " << fname << "\n";
        }
    }

    void writeHeader() {
        ofs_ << "slot,time_s,activeVNs,fn,fl,Q,C_Q,dO,"
                "migrations,avgTm_s,avgDown_s,"
                "routeCalc_ms,converge_ms,lsaPkts,helloPkts,changedLinks\n";
    }

    std::string fmt(double v, int prec = 3,
                    double sentinel = SlotRecord::NaN_SENTINEL) const {
        if (v <= sentinel + 1e-9) return "NaN";
        std::ostringstream s;
        s << std::fixed << std::setprecision(prec) << v;
        return s.str();
    }

    void write(const SlotRecord& r) {
        ofs_ << r.slot << ","
             << r.time_s << ","
             << r.activeVNs << ","
             << r.fn << ","
             << r.fl << ","
             << r.Q << ","
             << fmt(r.C_Q) << ","
             << fmt(r.dO) << ","
             << r.migrations << ","
             << fmt(r.avgTm, 4) << ","
             << fmt(r.avgDowntime, 4) << ","
             << fmt(r.routeCalc_ms, 3) << ","
             << fmt(r.converge_ms, 3) << ","
             << r.lsaPkts << ","
             << r.helloPkts << ","
             << r.changedLinks << "\n";
    }

private:
    std::ofstream ofs_;
};

struct SummaryStats {
    double mean;
    double stddev;
    double minv;
    double maxv;
    int n_valid;

    static SummaryStats of(const std::vector<double>& v,
                           double sentinel = SlotRecord::NaN_SENTINEL) {
        std::vector<double> valid;
        for (double x : v) {
            if (x > sentinel + 1e-9) valid.push_back(x);
        }

        if (valid.empty()) return {0, 0, 0, 0, 0};

        double sum = 0;
        for (double x : valid) sum += x;

        double mean = sum / valid.size();
        double ss = 0;
        double mn = valid[0];
        double mx = valid[0];

        for (double x : valid) {
            ss += (x - mean) * (x - mean);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
        }

        return {
            mean,
            std::sqrt(ss / valid.size()),
            mn,
            mx,
            static_cast<int>(valid.size())
        };
    }

    std::string str() const {
        std::ostringstream o;
        o << std::fixed << std::setprecision(3)
          << "mean=" << mean
          << " std=" << stddev
          << " min=" << minv
          << " max=" << maxv
          << " n=" << n_valid;
        return o.str();
    }
};

void writePlotScript(const std::string& csvFile) {
    const std::string pyFile = "plot_dyna_stn.py";
    std::ofstream py(pyFile);

    if (!py) {
        std::cerr << "[Plot] Cannot write " << pyFile << "\n";
        return;
    }

    auto L = [&](const std::string& s) {
        py << s << "\n";
    };

    L("#!/usr/bin/env python3");
    L("import csv");
    L("import math");
    L("import sys");
    L("");
    L("def to_float(x):");
    L("    try:");
    L("        v = float(x)");
    L("        return None if math.isnan(v) else v");
    L("    except:");
    L("        return None");
    L("");
    L("def col(rows, name):");
    L("    return [to_float(r[name]) for r in rows]");
    L("");
    L("def clean(values):");
    L("    return [v for v in values if v is not None]");
    L("");
    L("with open('" + csvFile + "', newline='') as f:");
    L("    rows = list(csv.DictReader(f))");
    L("");
    L("slots = clean(col(rows, 'slot'))");
    L("");
    L("print(f'Slots recorded: {len(rows)}')");
    L("for k in ['Q', 'C_Q', 'dO', 'converge_ms', 'routeCalc_ms', 'lsaPkts', 'avgTm_s']:");
    L("    v = clean(col(rows, k))");
    L("    if v:");
    L("        print(f'{k:18s}: mean={sum(v)/len(v):.3f} min={min(v):.3f} max={max(v):.3f}')");
    L("");
    L("try:");
    L("    import matplotlib.pyplot as plt");
    L("except ImportError:");
    L("    print('matplotlib not available')");
    L("    sys.exit(0)");
    L("");
    L("fig, axes = plt.subplots(3, 2, figsize=(13, 10))");
    L("fig.suptitle('Dyna-STN implementation metrics', fontsize=14)");
    L("");
    L("ax = axes[0, 0]");
    L("ax.plot(slots, clean(col(rows, 'Q')), marker='o', label='Q(t) = fn + fl')");
    L("ax.plot(slots, clean(col(rows, 'C_Q')), marker='s', label='C(Q) overhead')");
    L("ax.set_title('Dynamic bindings and overhead')");
    L("ax.set_xlabel('Slot')");
    L("ax.legend()");
    L("");
    L("ax = axes[0, 1]");
    L("ax.plot(slots, clean(col(rows, 'fn')), marker='o', label='Node bindings fn')");
    L("ax.plot(slots, clean(col(rows, 'fl')), marker='s', label='Link bindings fl')");
    L("ax.plot(slots, clean(col(rows, 'activeVNs')), marker='^', label='Active VNs')");
    L("ax.set_title('DDTM binding state')");
    L("ax.set_xlabel('Slot')");
    L("ax.legend()");
    L("");
    L("ax = axes[1, 0]");
    L("ax.plot(slots, clean(col(rows, 'routeCalc_ms')), marker='o', label='Route calculation')");
    L("ax.plot(slots, clean(col(rows, 'converge_ms')), marker='s', label='Convergence')");
    L("ax.set_title('Routing management plane time')");
    L("ax.set_xlabel('Slot')");
    L("ax.set_ylabel('ms')");
    L("ax.legend()");
    L("");
    L("ax = axes[1, 1]");
    L("ax.plot(slots, clean(col(rows, 'lsaPkts')), marker='o', label='LSA packets')");
    L("ax.plot(slots, clean(col(rows, 'helloPkts')), marker='s', label='Hello packets')");
    L("ax.set_title('Routing control messages')");
    L("ax.set_xlabel('Slot')");
    L("ax.legend()");
    L("");
    L("ax = axes[2, 0]");
    L("ax.plot(slots, clean(col(rows, 'changedLinks')), marker='o', label='Changed links')");
    L("ax.plot(slots, clean(col(rows, 'dO')), marker='s', label='Delta overhead')");
    L("ax.set_title('Topology dynamics')");
    L("ax.set_xlabel('Slot')");
    L("ax.legend()");
    L("");
    L("ax = axes[2, 1]");
    L("tm = col(rows, 'avgTm_s')");
    L("down = col(rows, 'avgDown_s')");
    L("ax.plot(slots, tm, marker='o', label='Avg migration time')");
    L("ax.plot(slots, down, marker='s', label='Avg downtime')");
    L("ax.set_title('Service migration')");
    L("ax.set_xlabel('Slot')");
    L("ax.set_ylabel('seconds')");
    L("ax.legend()");
    L("");
    L("plt.tight_layout()");
    L("out = 'dyna_stn_results.png'");
    L("plt.savefig(out, dpi=120, bbox_inches='tight')");
    L("print(f'Plot saved to {out}')");
    L("plt.show()");
    L("");

    std::cout << "[Plot script written to " << pyFile << "]\n";
}

int main() {
    
    std::cout << "=== Unit Test: Algorithm 1 Message Sequence ===\n";
    UnitTestResult utr = testAlg1MsgSequence();
    std::cout << utr.details;
    std::cout << (utr.passed ? "[UNIT TEST PASSED]\n" : "[UNIT TEST FAILED]\n");
    std::cout << std::string(55,'-') << "\n\n";

    const int    T_SLOTS = 12;   
    const double DT      = 600.0;  

    auto sats   = buildIridium();
    auto vns    = buildVNs(sats);
    auto vlinks = buildVLinks(C::M, C::N);

    ISL::Lv_ref = Lv_km();
    double Lh   = Lh_km(0.0);

    ServiceCubeModel scm;
    std::cout << "=== Constellation ===\n"
              << C::M<<"x"<<C::N<<" = "<<C::M*C::N<<" sats, "
              << vns.size()<<" VNs, "<<vlinks.size()<<" virtual links\n"
              << "Lv(intra)="<<std::fixed<<std::setprecision(0)<<ISL::Lv_ref
              <<" km  Lh(equator)="<<Lh<<" km\n"
              << "T_comm="<<std::setprecision(1)<<maxCommTime_s()<<" s\n"
              << "U_SNR(Lv)="<<std::setprecision(4)<<ISL::uSNR(ISL::Lv_ref)<<"  "
              << "U_SNR(Lh)="<<ISL::uSNR(Lh)<<"\n"
              << "ServiceCubeModel: angular Voronoi (no sphere radius)\n\n";


    std::cout << "=== Sample Service Cells  ===\n";
    for(int i=0;i<3;i++){
        auto sc=scm.getServiceCell(vns[i]);
        std::cout<<"  "<<vns[i].vid
                 <<"  subsat=("<<std::setprecision(1)<<sc.lat_ctr<<"deg,"<<sc.lon_ctr<<"deg)"
                 <<"  along=+-"<<std::setprecision(0)<<sc.along_km<<"km(+-"
                 <<std::setprecision(1)<<sc.along_deg<<"deg)"
                 <<"  cross=+-"<<std::setprecision(0)<<sc.cross_km<<"km(+-"
                 <<std::setprecision(1)<<sc.cross_deg<<"deg)"
                 <<"  bbox=["<<std::setprecision(1)<<sc.lat_S<<"deg,"<<sc.lat_N<<"deg]x["
                 <<sc.lon_W<<"deg,"<<sc.lon_E<<"deg]\n";
    }
    std::cout<<"\n";

    OverheadModel ovm;
    std::cout<<"=== Overhead Model C(Q) ===\n"
             <<"  C(Q) = "<<ovm.beta_base<<"*Q + "<<ovm.beta_cong
             <<"*Q^2/"<<ovm.Q_sat<<" + "<<ovm.beta_change<<"*|deltaQ|\n"
             <<"  Example: Q=66 -> C(Q)="
             <<std::setprecision(2)<<ovm.compute(66,10)
             <<" (vs simple Q=66)\n\n";

    DDTMPlane   ddtm;
    OSPFSim     ospf;
    TopoMetrics topo; topo.ovm=ovm;

    const std::string csvFile = "dyna_stn_slots.csv";
    CSVWriter csv(csvFile);
    csv.writeHeader();

    int    totalMig=0;
    double totalTm=0, totalDown=0;
    std::vector<double> hist_Q;
    std::vector<double> hist_CQ;
    std::vector<double> hist_conv;
    std::vector<double> hist_Tm;

    for(int t=1;t<=T_SLOTS;t++){
        double now=t*DT;
        std::cout<<"\n|- Slot t="<<t<<"  (t="<<std::fixed<<std::setprecision(0)
                 <<now<<" s) ------------------------------------\n";

        for(auto& s:sats) s.update(now);

        ddtm.runSlot(t,now,vns,sats,vlinks,0/*unused*/);
        topo.update(ddtm.fn_sum,ddtm.fl_sum);

        int actVN=0;
        for(auto& vn:vns) if(vn.active) actVN++;

        std::cout<<"| [DDTM] "; topo.print(t);
        std::cout<<"|        Active VNs: "<<actVN<<"/"<<vns.size()<<"\n";

        double slotTm = -1.0;
        double slotDown = -1.0;
        if(!ddtm.migrations.empty()){
            totalMig+=(int)ddtm.migrations.size();
            double sumTm=0, sumD=0;
            for(auto& m:ddtm.migrations){totalTm+=m.Tm;totalDown+=m.downtime;
                                          sumTm+=m.Tm;sumD+=m.downtime;}
            slotTm  = sumTm  / ddtm.migrations.size();
            slotDown= sumD   / ddtm.migrations.size();
            auto& m0=ddtm.migrations[0];
            std::cout<<"| [Mig] count="<<ddtm.migrations.size()
                     <<" avgTm="<<std::setprecision(4)<<ddtm.avgTm()<<"s"
                     <<"  ("<<m0.from<<"->"<<m0.to<<" Tm="<<m0.Tm<<"s)\n";
        } else {
            std::cout<<"| [Mig] none this slot\n";
        }

        auto ospfSt=ospf.runEpoch(vns,vlinks,sats,now);
        for(auto& vn:vns){
            if(!vn.active) continue;
            auto it=ospf.ribs.find(vn.vid);
            if(it!=ospf.ribs.end()){vn.router.rib=it->second;vn.router.updateFIBFromRIB();}
        }
        std::cout<<"| [OSPF] calc="<<std::setprecision(2)<<ospfSt.routeCalc_ms
                 <<" ms  conv="<<ospfSt.converge_ms
                 <<" ms  LSA="<<ospfSt.lsaPkts
                 <<"  deltaLinks="<<ospfSt.changedLinks<<"\n";
        SlotRecord rec{};
        rec.slot = t;
        rec.time_s = now;

        rec.activeVNs = actVN;
        rec.fn = topo.fn;
        rec.fl = topo.fl;
        rec.Q = topo.Q;
        rec.C_Q = topo.O;
        rec.dO = topo.dO;

        rec.migrations = static_cast<int>(ddtm.migrations.size());
        rec.avgTm = slotTm;
        rec.avgDowntime = slotDown;

        rec.routeCalc_ms = ospfSt.routeCalc_ms;
        rec.converge_ms = ospfSt.converge_ms;
        rec.lsaPkts = ospfSt.lsaPkts;
        rec.helloPkts = ospfSt.helloPkts;
        rec.changedLinks = ospfSt.changedLinks;

        csv.write(rec);

        hist_Q.push_back(topo.Q); hist_CQ.push_back(topo.O);
        hist_conv.push_back(ospfSt.converge_ms);
        if(slotTm>0) hist_Tm.push_back(slotTm);

        std::cout<<"| [Bindings] ";
        int sh=0; for(auto& vn:vns){if(sh++>=5){std::cout<<"...";break;}
            std::cout<<vn.vid<<":"<<(vn.bound.empty()?"none":vn.bound)<<" ";}
        std::cout<<"\n-------------------------------------------------\n";
    }

    std::cout << "\n=== Final summary ===\n";
    std::cout << "  Total migrations: " << totalMig << "\n";
    std::cout << "  Avg Tm: " 
            << (totalMig > 0 ? totalTm / totalMig : 0) << " s\n";
    std::cout << "  Avg downtime: "
            << (totalMig > 0 ? totalDown / totalMig : 0) << " s\n";
    std::cout << "  Cumulative LSA overhead: "
            << ospf.totalLSA << " packets\n";

    writePlotScript(csvFile);
    std::cout << "  CSV written to " << csvFile << "\n";
    std::cout << "  Run: python plot_dyna_stn.py\n";

    ddtm.mapSrv.print(2);

    for(auto& vn:vns){ if(!vn.bound.empty()){ddtm.printTrace(vn.vid);break;} }

    std::cout<<"\nDone.\n";
    return 0;
}

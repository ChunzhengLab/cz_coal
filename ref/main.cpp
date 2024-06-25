#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <valarray>
#include <cmath>
#include <limits>
#include <ctime>
#include <cstdlib>
#include <random>
#include <numeric>
#include <utility>
#include <algorithm>
#include <map>
#include <unordered_map>

using namespace std;

bool readZPCData(ifstream&);
void outputHadronData(int, ofstream&);
void randomQuarkList();
void doCoalescence();
void wignerCoalescence();
void updateHadronList();
bool checkIndex(const vector<int>&);
void eraseIndex(const vector<int>&);
inline double squareVector(const valarray<double>&);
inline double wigner2(int, int, const double&, const double&);
inline double wigner3(int, int, int, const double&, const double&);
inline double wigner4(int, int, int, int, const double&, const double&);
inline double wignerN(vector<int>, const double&, const double&);
double mesonCoalescence(vector<int>&);
double antimesonCoalescence(vector<int>&);
double baryonCoalescence(vector<int>&);
double antibaryonCoalescence(vector<int>&);
double tetraquarkCoalescence(vector<int>&);
inline void lorentz(const double&, const valarray<double>&,
                          double&,       valarray<double>&, const valarray<double>&);
double calculateCoordinateDistance(vector<double>&);
double calculateMomentDistance(vector<double>&);

namespace parameters {
    const double pi = acos(-1.0);
    int NSG, MUL, MesonNum, BaryonNum, Quark_num, AntiQuark_num;
    int NELP, NINP, NELT, NINT;
    double phiRP, bimp;

    int opt_f0 = 0; // 0 for close, 1 for open
}

namespace space {
    vector<int> types, orders;
    vector<double> pxs, pys, pzs, ms, es, xs, ys, zs, ts;
    vector<vector<int>> hadrons;
}

namespace flags {
    valarray<bool> flags_all;
    int quarks_num, antiquarks_num;
    int ds_num, us_num, ss_num, cs_num, bs_num;
    int dbs_num, ubs_num, sbs_num, cbs_num, bbs_num;
    vector<int> quarks_index, antiquarks_index;
    vector<int> ds_index, us_index, ss_index, cs_index, bs_index;
    vector<int> dbs_index, ubs_index, sbs_index, cbs_index, bbs_index;

    vector<vector<vector<pair<int,int>>>> bins;
}

namespace wigner {
    const double _hbar2 = 1.0/(0.19732697*0.19732697);
    const double _small = numeric_limits<double>::lowest();

    const double g_meson = 1.0/36.0;   // pion
    const double g_baryon = 2.0/216.0; // proton
    const double g_tetra = 1.0/1236.0; // tetraquark // f0(980)

    const double rms_pion2 = 0.61; // Phys. Lett. B 51, 402 (1974).
    const double rms_K2    = 0.34; // Phys. Lett. B 178, 435 (1986).
    const double rms_proton2 = 0.7071; // PDG: 2022
    // Phys Rev C. 101. 024908 from sigma = 2.4 fm, rms_4 = 2.93 fm
    const double rms_tetra2 = 1.0*1.0;

    const double _ln8 = log(8.0*g_meson);
    const double _ln64 = log(64.0*g_baryon);
    const double _ln512 = log(512.0*g_tetra);

    multimap<double,vector<int>,greater<double>> wigners;
}

namespace test {
    int _t_2=0, _t_3=0, _t_4=0, _t_t=0;
}

int main() {

    cout << "Wigner function mode." << endl;

    time_t timer1, timer2;
    time(&timer1);

    ifstream zpc("./input/zpc.dat");   // zpc data input
    ofstream out("./ana/init_combination.dat");  // hadron after zpc

    int event = 0;
    while(true) {

        cout << "------------------------ event:   " << event << endl;

        test::_t_2=0; test::_t_3=0; test::_t_4=0; test::_t_t=0;

        if(!readZPCData(zpc)) break;

        randomQuarkList(); // test

        doCoalescence();

        if(test::_t_2+test::_t_3+test::_t_4==0) {
            cerr << "Something is wrong, particle number is 0. " << endl;
            exit(0);
        }

        outputHadronData(event, out);
        vector<double> dist_c, dist_p;
        double dist_c_final = calculateCoordinateDistance(dist_c);
        double dist_m_final = calculateMomentDistance(dist_p);

        cout << "_t_2:     " << test::_t_2 << endl;
        cout << "_t_3:     " << test::_t_3 << endl;
        cout << "_t_4:     " << test::_t_4 << endl;
        cout << "_t_t:     " << test::_t_t << endl;
        cout << "Final Coordinate distance:   " << dist_c_final << endl;
        cout << "Final  Momentum  distance:   " << dist_m_final << endl;
        cout << "Meson:   " << dist_c[0] << ",  " << dist_p[0] << endl;
        cout << "Baryon:  " << dist_c[1] << ",  " << dist_p[1] << endl;
        cout << "Tetra:   " << dist_c[2] << ",  " << dist_p[2] << endl;

        ++event;
    }

    time(&timer2);
    cout << "Duration time: " << float(difftime(timer2,timer1))/60.0 << " minutes "<< endl;

    zpc.close();
    out.close();
    return 0;
}

void randomQuarkList() {
    using namespace parameters;
    using space::hadrons;
    using space::types;

    hadrons.clear();
    vector<int> quarks, aquarks;
    for(int i=0; i<MUL; ++i) {
        int type = types[i];
        if(type > 0) quarks.push_back(i);
        if(type < 0) aquarks.push_back(i);
    }

    shuffle(quarks.begin(), quarks.end(), default_random_engine(1));
    shuffle(aquarks.begin(), aquarks.end(), default_random_engine(2));
    for(int i=0; i<MesonNum; ++i) {
        hadrons.push_back(vector<int>());
        hadrons.back().push_back(quarks[i]);
        hadrons.back().push_back(aquarks[i]);
    }
    for(int i=MesonNum; i<Quark_num; i+=3) {
        hadrons.push_back(vector<int>());
        for(int j=0; j<3; ++j) {
            hadrons.back().push_back(quarks[i+j]);
        }
    }
    for(int i=MesonNum; i<AntiQuark_num; i+=3) {
        hadrons.push_back(vector<int>());
        for(int j=0; j<3; ++j) {
            hadrons.back().push_back(aquarks[i+j]);
        }
    }
    cout << "random hadron:     " << hadrons.size() << endl;
    vector<double> dist_c, dist_p;
    double dist_c_final = calculateCoordinateDistance(dist_c);
    double dist_m_final = calculateMomentDistance(dist_p);
    cout << "Random Coordinate distance:  " << dist_c_final << endl;
    cout << "Random  Momentum  distance:  " << dist_m_final << endl;
    cout << "Meson:   " << dist_c[0] << ",  " << dist_p[0] << endl;
    cout << "Baryon:  " << dist_c[1] << ",  " << dist_p[1] << endl;
    cout << "Tetra:   " << dist_c[2] << ",  " << dist_p[2] << endl;
}

void doCoalescence() {
    using namespace parameters;
    using namespace flags;
    using space::hadrons;
    using wigner::wigners;

    hadrons.clear();
    flags_all = false;
    // main loop
    while(true) {

        wigners.clear();

        wignerCoalescence();
        updateHadronList();

        cout << "hadron:   " << hadrons.size() << endl;
        cout << "q_rest:   " << quarks_num << endl;
        cout << "a_rest:   " << antiquarks_num << endl;
        if(quarks_num==0 && antiquarks_num==0) break;
        if((quarks_num<=2 && antiquarks_num==0) ||
           (quarks_num==0 && antiquarks_num<=2)) {
            cerr << "Error, rest quark cannot form hadrons. " << endl;
            cerr << "q_rest:   " << quarks_num << endl;
            cerr << "a_rest:   " << antiquarks_num << endl;
            exit(0);
        }
    }

    if(test::_t_2*2+test::_t_3*3+test::_t_4*4 != Quark_num+AntiQuark_num) {
        cerr << "Error, quark number is not conserved." << endl;
        exit(0);
    }
}

bool checkIndex(const vector<int>& qs) {
    using namespace flags;

    for(const int& qi : qs) {
        if(flags_all[qi]) return false;
    }
    return true;
};

void eraseIndex(const vector<int>& qs) {
    using namespace flags;
    using space::types;

    int type;

    for(const int& qi : qs) {
        type = types[qi];
        flags_all[qi] = true;

        switch(type) {
        case 1:
            --ds_num; break;
        case 2:
            --us_num; break;
        case 3:
            --ss_num; break;
        case 4:
            --cs_num; break;
        case 5:
            --bs_num; break;
        case -1:
            --dbs_num; break;
        case -2:
            --ubs_num; break;
        case -3:
            --sbs_num; break;
        case -4:
            --cbs_num; break;
        case -5:
            --bbs_num; break;
        default:
            cerr << "quark type cannot be " << type << endl;
            exit(0);
        }

        if(type > 0) {
            --quarks_num;
        } else if(type < 0) {
            --antiquarks_num;
        }
    }
};

void wignerCoalescence() {
    using namespace flags;
    using namespace wigner;
    using parameters::opt_f0;

    double ws;

    valarray<bool> temps1(flags_all);
/// meson
    if(quarks_num!=0 && antiquarks_num!=0) {
        for(const int& qi : quarks_index) { // quark first
            if(temps1[qi]) continue;
            vector<int> qs = {qi};
            ws = mesonCoalescence(qs);
            wigners.insert({ws, qs});

            temps1[qs[1]] = true;
        }

        cout << "size of wigner-mesons:       " << wigners.size() << endl;

        for(const int& qi : antiquarks_index) { // anti-quark first
            if(temps1[qi]) continue;
            vector<int> qs = {qi};
            ws = antimesonCoalescence(qs);
            wigners.insert({ws, qs});
        }

        cout << "size of wigner-anti-mesons:  " << wigners.size() << endl;
    }

/// baryon
    valarray<bool> temps2(flags_all);
    if(quarks_num>=3) {
        for(const int& qi : quarks_index) {
            if(temps2[qi]) continue;
            vector<int> qs = {qi};
            ws = baryonCoalescence(qs);
            wigners.insert({ws, qs});

            temps2[qs[1]] = true;
            temps2[qs[2]] = true;
        }

        cout << "size of wigner-baryons:      " << wigners.size() << endl;
    }

/// antibaryon
    valarray<bool> temps3(flags_all);
    if(antiquarks_num>=3) {
        for(const int& qi : antiquarks_index) {
            if(temps3[qi]) continue;
            vector<int> qs = {qi};
            ws = antibaryonCoalescence(qs);
            wigners.insert({ws, qs});

            temps3[qs[1]] = true;
            temps3[qs[2]] = true;
        }

        cout << "size of wigner-anti-baryons: " << wigners.size() << endl;
    }

    if(opt_f0==0) return;
/// tetraquark
    if(ss_num!=0 && sbs_num!=0 && us_num!=0 && ubs_num!=0) {
        for(const int& qi : ss_index) { // s-sbar-u-ubar
            if(flags_all[qi]) continue;
            vector<int> qs = {qi};
            ws = tetraquarkCoalescence(qs);
            wigners.insert({ws, qs});
        }

        cout << "size of wigner-tetraquarks:  " << wigners.size() << endl;
    }
}

void updateHadronList() {
    using namespace test;
    using namespace wigner;
    using space::hadrons;

    for(const auto& qs_pair : wigners) {
        if(checkIndex(qs_pair.second)) {
            hadrons.push_back(qs_pair.second);
            eraseIndex(qs_pair.second);

            switch(qs_pair.second.size()) {
            case 2:
                ++_t_2; break;
            case 3:
                ++_t_3; break;
            case 4:
                ++_t_4; break;
            default:
                cerr << "quark number of Hadron cannot be " <<
                        qs_pair.second.size() << endl;
                exit(0);
            }
        }
    }

}

double mesonCoalescence(vector<int>& qs) { // notes: 2023.7.4
    using namespace wigner;
    using namespace flags;

    int q1 = qs[0], q2;
    double _ws, _temp_w;
    _ws = _small;
    for(const int& _q2 : antiquarks_index) { // check distance between q1 and q2?
        if(flags_all[_q2]) continue;
        _temp_w = wigner2(q1, _q2, rms_pion2, _ln8);
        if(_temp_w > _ws) {
            q2 = _q2;
            _ws = _temp_w;
        }
    }
    qs.push_back(q2);
    return _ws;
}

double antimesonCoalescence(vector<int>& qs) {
    using namespace wigner;
    using namespace flags;

    int q1 = qs[0], q2;
    double _ws, _temp_w;
    _ws = _small;
    for(const int& _q2 : quarks_index) {
        if(flags_all[_q2]) continue;
        _temp_w = wigner2(q1, _q2, rms_pion2, _ln8);
        if(_temp_w > _ws) {
            q2 = _q2;
            _ws = _temp_w;
        }
    }
    qs.push_back(q2);
    swap(qs[0], qs[1]); // quark in first position
    return _ws;
}

double baryonCoalescence(vector<int>& qs) {
    using namespace wigner;
    using namespace flags;

    int q1 = qs[0], q2, q3;
    double _ws, _temp_w;
    _ws = _small;
    flags_all[q1] = true;
    for(const int& _q2 : quarks_index) {
        if(flags_all[_q2]) continue;
        _temp_w = wigner2(q1, _q2, rms_proton2, _ln64); // zwhtest
        if(_temp_w > _ws) {
            q2 = _q2;
            _ws = _temp_w;
        }
    }
    qs.push_back(q2);

    _ws = _small;
    flags_all[q2] = true;
    for(const int& _q3 : quarks_index) {
        if(flags_all[_q3]) continue;
        _temp_w = wigner3(q1, q2, _q3, rms_proton2, _ln64);
        if(_temp_w > _ws) {
            q3 = _q3;
            _ws = _temp_w;
        }
    }
    qs.push_back(q3);
    flags_all[q1] = false;
    flags_all[q2] = false;
    return _ws;
}

double antibaryonCoalescence(vector<int>& qs) {
    using namespace wigner;
    using namespace flags;

    int q1 = qs[0], q2, q3;
    double _ws, _temp_w;
    _ws = _small;
    flags_all[q1] = true;
    for(const int& _q2 : antiquarks_index) {
        if(flags_all[_q2]) continue;
        _temp_w = wigner2(q1, _q2, rms_proton2, _ln64); // zwhtest
        if(_temp_w > _ws) {
            q2 = _q2;
            _ws = _temp_w;
        }
    }
    qs.push_back(q2);

    _ws = _small;
    flags_all[q2] = true;
    for(const int& _q3 : antiquarks_index) {
        if(flags_all[_q3]) continue;
        _temp_w = wigner3(q1, q2, _q3, rms_proton2, _ln64);
        if(_temp_w > _ws) {
            q3 = _q3;
            _ws = _temp_w;
        }
    }
    qs.push_back(q3);
    flags_all[q1] = false;
    flags_all[q2] = false;
    return _ws;
}

double tetraquarkCoalescence(vector<int>& qs) { // s-sbar-u-ubar
    using namespace wigner;
    using namespace flags;

    int q1 = qs[0], q2, q3, q4;
    double _ws, _temp_w;
    _ws = _small;
    for(const int& _q2 : sbs_index) {
        if(flags_all[_q2]) continue;
        _temp_w = wigner2(q1, _q2, rms_tetra2, _ln512); // zwhtest
        if(_temp_w > _ws) {
            q2 = _q2;
            _ws = _temp_w;
        }
    }
    qs.push_back(q2);

    _ws = _small;
    for(const int& _q3 : us_index) {
        if(flags_all[_q3]) continue;
        _temp_w = wigner3(q1, q2, _q3, rms_tetra2, _ln512);
        if(_temp_w > _ws) {
            q3 = _q3;
            _ws = _temp_w;
        }
    }
    qs.push_back(q3);

    _ws = _small;
    for(const int& _q4 : ubs_index) {
        if(flags_all[_q4]) continue;
        _temp_w = wigner4(q1, q2, q3, _q4, rms_tetra2, _ln512);
        if(_temp_w > _ws) {
            q4 = _q4;
            _ws = _temp_w;
        }
    }
    qs.push_back(q4);
    return _ws;
}

inline double squareVector(const valarray<double>& _array) {
    return _array[0]*_array[0] +
           _array[1]*_array[1] +
           _array[2]*_array[2];
};

inline double wigner2(int q1, int q2,
                      const double& _rms, const double& _lnN) {
    using namespace space;
    using wigner::_hbar2;

    double es1 = es[q1], es2 = es[q2], es1cm, es2cm;
    double ts1 = ts[q1], ts2 = ts[q2], ts1cm, ts2cm;
    double ms1 = ms[q1], ms2 = ms[q2];
    valarray<double> ps1(3), ps2(3), ps1cm(3), ps2cm(3);
    valarray<double> rs1(3), rs2(3), rs1cm(3), rs2cm(3);
    ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
    ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
    rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
    rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];

    double _m2 = ms1 + ms2;
    double _t_mu_1 = 1.0/2.0*(1.0/ms1 + 1.0/ms2);
    double _omega  = 3.0/(4.0*_rms)*(1.0/ms1 + 1.0/ms2 - 2.0/_m2);
    double _sigma1 = _t_mu_1/_omega;

    valarray<double> _beta = (ps1 + ps2)/(es1 + es2);

    lorentz(es1, ps1, es1cm, ps1cm, _beta);
    lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
    lorentz(es2, ps2, es2cm, ps2cm, _beta);
    lorentz(ts2, rs2, ts2cm, rs2cm, _beta);

    double _t_max = max(ts1cm, ts2cm);
    rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
    rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);

    valarray<double> _rho_r = rs1cm - rs2cm;
    valarray<double> _rho_p = (ms2*ps1cm - ms1*ps2cm)/_m2;

    double _r2n = squareVector(_rho_r);
    double _p2n = squareVector(_rho_p);

    return _lnN -(1.0/2.0)*_r2n/_sigma1 -(2.0/1.0)*_p2n*_sigma1*_hbar2; // ln(Wigner)
}

inline double wigner3(int q1, int q2, int q3,
                      const double& _rms, const double& _lnN) {
    using namespace space;
    using wigner::_hbar2;

    double es1 = es[q1], es2 = es[q2], es3 = es[q3], es1cm, es2cm, es3cm;
    double ts1 = ts[q1], ts2 = ts[q2], ts3 = ts[q3], ts1cm, ts2cm, ts3cm;
    double ms1 = ms[q1], ms2 = ms[q2], ms3 = ms[q3];
    valarray<double> ps1(3), ps2(3), ps3(3), ps1cm(3), ps2cm(3), ps3cm(3);
    valarray<double> rs1(3), rs2(3), rs3(3), rs1cm(3), rs2cm(3), rs3cm(3);
    ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
    ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
    ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];
    rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
    rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];
    rs3[0] = xs[q3];  rs3[1] = ys[q3];  rs3[2] = zs[q3];

    double _m2 = ms1 + ms2;
    double _m3 = _m2 + ms3;
    double _t_mu_1 = 1.0/2.0*(1.0/ms1 + 1.0/ms2);
    double _t_mu_2 = 2.0/3.0*(1.0/_m2 + 1.0/ms3);
    double _omega  = 3.0/(6.0*_rms)*(1.0/ms1 + 1.0/ms2 + 1.0/ms3 - 3.0/_m3);
    double _sigma1 = _t_mu_1/_omega;
    double _sigma2 = _t_mu_2/_omega;

    valarray<double> _beta = (ps1 + ps2 + ps3)/(es1 + es2 + es3);

    lorentz(es1, ps1, es1cm, ps1cm, _beta);
    lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
    lorentz(es2, ps2, es2cm, ps2cm, _beta);
    lorentz(ts2, rs2, ts2cm, rs2cm, _beta);
    lorentz(es3, ps3, es3cm, ps3cm, _beta);
    lorentz(ts3, rs3, ts3cm, rs3cm, _beta);

    double _t_max = max(max(ts1cm, ts2cm),ts3cm);
    rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
    rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);
    rs3cm = rs3cm + ps3cm/es3cm*(_t_max-ts3cm);

    valarray<double> _rho_r    = rs1cm - rs2cm;
    valarray<double> _lambda_r = (ms1*rs1cm + ms2*rs2cm)/_m2 - rs3cm;
    valarray<double> _rho_p    = (ms2*ps1cm - ms1*ps2cm)/_m2;
    valarray<double> _lambda_p = (ms3*(ps1cm+ps2cm) - _m2*ps3cm)/_m3;

    double _r2n = squareVector(_rho_r);
    double _p2n = squareVector(_rho_p);
    double _r3n = squareVector(_lambda_r);
    double _p3n = squareVector(_lambda_p);

    return _lnN -(1.0/2.0)*_r2n/_sigma1 -(2.0/1.0)*_p2n*_sigma1*_hbar2
                -(2.0/3.0)*_r3n/_sigma2 -(3.0/2.0)*_p3n*_sigma2*_hbar2;  // ln(Wigner)
}

inline double wigner4(int q1, int q2, int q3, int q4,
                      const double& _rms, const double& _lnN) {
    using namespace space;
    using wigner::_hbar2;

    double es1 = es[q1], es2 = es[q2], es3 = es[q3], es4 = es[q4], es1cm, es2cm, es3cm, es4cm;
    double ts1 = ts[q1], ts2 = ts[q2], ts3 = ts[q3], ts4 = ts[q4], ts1cm, ts2cm, ts3cm, ts4cm;
    double ms1 = ms[q1], ms2 = ms[q2], ms3 = ms[q3], ms4 = ms[q4];
    valarray<double> ps1(3), ps2(3), ps3(3), ps4(3), ps1cm(3), ps2cm(3), ps3cm(3), ps4cm(3);
    valarray<double> rs1(3), rs2(3), rs3(3), rs4(3), rs1cm(3), rs2cm(3), rs3cm(3), rs4cm(3);
    ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
    ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
    ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];
    ps4[0] = pxs[q4]; ps4[1] = pys[q4]; ps4[2] = pzs[q4];
    rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
    rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];
    rs3[0] = xs[q3];  rs3[1] = ys[q3];  rs3[2] = zs[q3];
    rs4[0] = xs[q4];  rs4[1] = ys[q4];  rs4[2] = zs[q4];

    double _m2 = ms1 + ms2;
    double _m3 = _m2 + ms3;
    double _m4 = _m3 + ms4;
    double _t_mu_1 = 1.0/2.0*(1.0/ms1 + 1.0/ms2);
    double _t_mu_2 = 2.0/3.0*(1.0/_m2 + 1.0/ms3);
    double _t_mu_3 = 3.0/4.0*(1.0/_m3 + 1.0/ms4);
    double _omega  = 3.0/(8.0*_rms)*(1.0/ms1 + 1.0/ms2 + 1.0/ms3 + 1.0/ms4 - 4.0/_m4);
    double _sigma1 = _t_mu_1/_omega;
    double _sigma2 = _t_mu_2/_omega;
    double _sigma3 = _t_mu_3/_omega;

    valarray<double> _beta = (ps1 + ps2 + ps3 + ps4)/(es1 + es2 + es3 + es4);

    lorentz(es1, ps1, es1cm, ps1cm, _beta);
    lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
    lorentz(es2, ps2, es2cm, ps2cm, _beta);
    lorentz(ts2, rs2, ts2cm, rs2cm, _beta);
    lorentz(es3, ps3, es3cm, ps3cm, _beta);
    lorentz(ts3, rs3, ts3cm, rs3cm, _beta);
    lorentz(es4, ps4, es4cm, ps4cm, _beta);
    lorentz(ts4, rs4, ts4cm, rs4cm, _beta);

    double _t_max = max(max(max(ts1cm, ts2cm),ts3cm),ts4cm);
    rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
    rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);
    rs3cm = rs3cm + ps3cm/es3cm*(_t_max-ts3cm);
    rs4cm = rs4cm + ps4cm/es4cm*(_t_max-ts4cm);

    valarray<double> _rho_r    = rs1cm - rs2cm;
    valarray<double> _lambda_r = (ms1*rs1cm + ms2*rs2cm)/_m2 - rs3cm;
    valarray<double> _theta_r  = (ms1*rs1cm + ms2*rs2cm + ms3*rs3cm)/_m3 - rs4cm;
    valarray<double> _rho_p    = (ms2*ps1cm - ms1*ps2cm)/_m2;
    valarray<double> _lambda_p = (ms3*(ps1cm+ps2cm) - _m2*ps3cm)/_m3;
    valarray<double> _theta_p  = (ms4*(ps1cm+ps2cm+ps3cm) - _m3*ps4cm)/_m4;

    double _r2n = squareVector(_rho_r);
    double _p2n = squareVector(_rho_p);
    double _r3n = squareVector(_lambda_r);
    double _p3n = squareVector(_lambda_p);
    double _r4n = squareVector(_theta_r);
    double _p4n = squareVector(_theta_p);

    return _lnN -(1.0/2.0)*_r2n/_sigma1 -(2.0/1.0)*_p2n*_sigma1*_hbar2
                -(2.0/3.0)*_r3n/_sigma2 -(3.0/2.0)*_p3n*_sigma2*_hbar2
                -(3.0/4.0)*_r4n/_sigma3 -(4.0/3.0)*_p4n*_sigma3*_hbar2;   // ln(Wigner)
}

inline double wignerN(vector<int> qs, const double& _rms, const double& _lnN) {
// zwhtest haven't tested yet.
    using namespace space;
    using wigner::_hbar2;

    int qs_n = qs.size();
    double *ess = new double[qs_n];
    double *esscm = new double[qs_n];
    double *tss = new double[qs_n];
    double *tsscm = new double[qs_n];
    double *mss = new double[qs_n];
    valarray<double> *pss = new valarray<double>[qs_n];
    valarray<double> *psscm = new valarray<double>[qs_n];
    valarray<double> *rss = new valarray<double>[qs_n];
    valarray<double> *rsscm = new valarray<double>[qs_n];

    int qs_i;
    for(int i=0; i<qs_n; ++i) {
        qs_i = qs[i];
        pss[i].resize(3);
        rss[i].resize(3);

        mss[i]    = ms[qs_i];
        ess[i]    = es[qs_i];
        tss[i]    = ts[qs_i];
        pss[i][0] = pxs[qs_i];
        pss[i][1] = pys[qs_i];
        pss[i][2] = pzs[qs_i];
        rss[i][0] = xs[qs_i];
        rss[i][1] = ys[qs_i];
        rss[i][2] = zs[qs_i];
    }

    double *_mN = new double[qs_n];
    _mN[0] = mss[0];
    for(int i=1; i<qs_n; ++i) {
        _mN[i] = _mN[i-1] + mss[i];
    }

    double _omega = 0.0;
    for(int i=0; i<qs_n; ++i) {
        _omega = _omega + 1.0/mss[i];
    }
    _omega = (_omega - double(qs_n)/_mN[qs_n-1])*3.0/(2.0*qs_n*_rms);
    double *_t_mu = new double[qs_n];
    double *_sigma = new double[qs_n];
    for(int i=1; i<qs_n; ++i) {
        _t_mu[i] = double(i)/double(i+1)*(1.0/_mN[i-1] + 1.0/mss[i]);
        _sigma[i] = _t_mu[i]/_omega;
    }

    valarray<double> _beta(0.0, 3);
    for(int i=0; i<qs_n; ++i) {
        _beta = _beta + pss[i];
    }
    _beta = _beta/accumulate(ess,ess+qs_n,0.0);

    for(int i=0; i<qs_n; ++i) {
        psscm[i].resize(3);
        lorentz(ess[i], pss[i], esscm[i], psscm[i], _beta);
        rsscm[i].resize(3);
        lorentz(tss[i], rss[i], tsscm[i], rsscm[i], _beta);
    }

    double _t_max = *max_element(tsscm,tsscm+qs_n);
    for(int i=0; i<qs_n; ++i) {
        rsscm[i] = rsscm[i] + psscm[i]/esscm[i]*(_t_max-tsscm[i]);
    }

    valarray<double> *_rrn = new valarray<double>[qs_n];
    for(int i=0; i<qs_n; ++i) _rrn[i].resize(3);
    _rrn[0] = 0.0;
    for(int i=1; i<qs_n; ++i) {
        _rrn[i] = _rrn[i-1] + mss[i-1]*rsscm[i-1];
    }
    for(int i=1; i<qs_n; ++i) {
        _rrn[i] = (_rrn[i]/_mN[i-1] - rsscm[i])*sqrt(double(i)/double(i+1));
    }

    valarray<double> *_ppn = new valarray<double>[qs_n];
    for(int i=0; i<qs_n; ++i) _ppn[i].resize(3);
    _ppn[0] = 0.0;
    for(int i=1; i<qs_n; ++i) {
        _ppn[i] = _ppn[i-1] + psscm[i-1];
    }
    for(int i=1; i<qs_n; ++i) {
        _ppn[i] = (_ppn[i]*mss[i] - psscm[i]*_mN[i-1])/_mN[i]*sqrt(double(i+1)/double(i));
    }

    double *_rNn = new double[qs_n];
    double *_pNn = new double[qs_n];
    for(int i=1; i<qs_n; ++i) {
        _rNn[i] = 0.0;
        _pNn[i] = 0.0;
        for(int j=0; j<3; ++j) {
            _rNn[i] = _rNn[i] + _rrn[i][j]*_rrn[i][j];
            _pNn[i] = _pNn[i] + _ppn[i][j]*_ppn[i][j];
        }
    }

    double _w = _lnN;
    for(int i=1; i<qs_n; ++i) {
        _w = _w - _rNn[i]/_sigma[i] - _pNn[i]*_sigma[i]*_hbar2;
    }

    return _w;  // ln(Wigner)
}

inline void lorentz(const double& es1, const valarray<double>& ps1,
                          double& es2,       valarray<double>& ps2,
                                       const valarray<double>& beta) {
    double beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    if(beta2 < 1.0e-9) {
        ps2 = ps1;
        es2 = es1;
        return;
    } else if(beta2 > 0.999999999999999) {
        beta2 = 0.999999999999999;
    }

    double ga = 1.0/sqrt(1.0 - beta2);
    double _dot = beta[0]*ps1[0] + beta[1]*ps1[1] + beta[2]*ps1[2];

    es2 = ga*(es1 - _dot);
    ps2 = ps1 - beta*(ga*es1 - (ga - 1.0)/beta2*_dot);
}

double calculateCoordinateDistance(vector<double>& dists) {
    using namespace space;

    valarray<double> ps1(3), ps2(3), ps3(3), ps4(3);
    valarray<double> ps1cm(3), ps2cm(3), ps3cm(3), ps4cm(3);
    valarray<double> rs1(3), rs2(3), rs3(3), rs4(3);
    valarray<double> rs1cm(3), rs2cm(3), rs3cm(3), rs4cm(3);
    double es1, es2, es3, es4, es1cm, es2cm, es3cm, es4cm;
    double ts1, ts2, ts3, ts4, ts1cm, ts2cm, ts3cm, ts4cm;
    double _t_max;
    valarray<double> _beta(3);

    int count_d = 0;
    int count_d2 = 0, count_d3 = 0, count_d4 = 0;
    int q1, q2, q3, q4;
    double dist = 0.0, distt;
    double dist2 = 0.0, dist3 = 0.0, dist4 = 0.0;
    for(auto hadron:hadrons) {
        switch(hadron.size()) {
        case 2:
            q1 = hadron[0];
            q2 = hadron[1];

            es1 = es[q1], es2 = es[q2];
            ts1 = ts[q1], ts2 = ts[q2];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
            rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
            rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];

            _beta = (ps1 + ps2)/(es1 + es2);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);
            lorentz(ts2, rs2, ts2cm, rs2cm, _beta);

            _t_max = max(ts1cm, ts2cm);
            rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
            rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);

            distt = sqrt((rs1cm[0]-rs2cm[0])*(rs1cm[0]-rs2cm[0])
                       + (rs1cm[1]-rs2cm[1])*(rs1cm[1]-rs2cm[1])
                       + (rs1cm[2]-rs2cm[2])*(rs1cm[2]-rs2cm[2]));
            dist += distt;
            dist2 += distt;

            count_d += 1;
            count_d2 += 1;
            break;
        case 3:
            q1 = hadron[0];
            q2 = hadron[1];
            q3 = hadron[2];

            es1 = es[q1], es2 = es[q2], es3 = es[q3];
            ts1 = ts[q1], ts2 = ts[q2], ts3 = ts[q3];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
            ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];
            rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
            rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];
            rs3[0] = xs[q3];  rs3[1] = ys[q3];  rs3[2] = zs[q3];

            _beta = (ps1 + ps2 + ps3)/(es1 + es2 + es3);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);
            lorentz(ts2, rs2, ts2cm, rs2cm, _beta);
            lorentz(es3, ps3, es3cm, ps3cm, _beta);
            lorentz(ts3, rs3, ts3cm, rs3cm, _beta);

            _t_max = max(max(ts1cm, ts2cm),ts3cm);
            rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
            rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);
            rs3cm = rs3cm + ps3cm/es3cm*(_t_max-ts3cm);

            distt = sqrt((rs1cm[0]-rs2cm[0])*(rs1cm[0]-rs2cm[0])
                       + (rs1cm[1]-rs2cm[1])*(rs1cm[1]-rs2cm[1])
                       + (rs1cm[2]-rs2cm[2])*(rs1cm[2]-rs2cm[2]))
                  + sqrt((rs1cm[0]-rs3cm[0])*(rs1cm[0]-rs3cm[0])
                       + (rs1cm[1]-rs3cm[1])*(rs1cm[1]-rs3cm[1])
                       + (rs1cm[2]-rs3cm[2])*(rs1cm[2]-rs3cm[2]))
                  + sqrt((rs2cm[0]-rs3cm[0])*(rs2cm[0]-rs3cm[0])
                       + (rs2cm[1]-rs3cm[1])*(rs2cm[1]-rs3cm[1])
                       + (rs2cm[2]-rs3cm[2])*(rs2cm[2]-rs3cm[2]));
            dist += distt;
            dist3 += distt;

            count_d += 3;
            count_d3 += 3;
            break;
        case 4:
            // tetraquark
            q1 = hadron[0];
            q2 = hadron[1];
            q3 = hadron[2];
            q4 = hadron[3];

            es1 = es[q1], es2 = es[q2], es3 = es[q3], es4 = es[q4];
            ts1 = ts[q1], ts2 = ts[q2], ts3 = ts[q3], ts4 = ts[q4];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
            ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];
            ps4[0] = pxs[q4]; ps4[1] = pys[q4]; ps4[2] = pzs[q4];
            rs1[0] = xs[q1];  rs1[1] = ys[q1];  rs1[2] = zs[q1];
            rs2[0] = xs[q2];  rs2[1] = ys[q2];  rs2[2] = zs[q2];
            rs3[0] = xs[q3];  rs3[1] = ys[q3];  rs3[2] = zs[q3];
            rs4[0] = xs[q4];  rs4[1] = ys[q4];  rs4[2] = zs[q4];

            _beta = (ps1 + ps2 + ps3 + ps4)/(es1 + es2 + es3 + es4);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(ts1, rs1, ts1cm, rs1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);
            lorentz(ts2, rs2, ts2cm, rs2cm, _beta);
            lorentz(es3, ps3, es3cm, ps3cm, _beta);
            lorentz(ts3, rs3, ts3cm, rs3cm, _beta);
            lorentz(es4, ps4, es4cm, ps4cm, _beta);
            lorentz(ts4, rs4, ts4cm, rs4cm, _beta);

            _t_max = max(max(max(ts1cm, ts2cm),ts3cm),ts4cm);
            rs1cm = rs1cm + ps1cm/es1cm*(_t_max-ts1cm);
            rs2cm = rs2cm + ps2cm/es2cm*(_t_max-ts2cm);
            rs3cm = rs3cm + ps3cm/es3cm*(_t_max-ts3cm);
            rs4cm = rs4cm + ps4cm/es4cm*(_t_max-ts4cm);

            distt = sqrt((rs1cm[0]-rs2cm[0])*(rs1cm[0]-rs2cm[0])
                       + (rs1cm[1]-rs2cm[1])*(rs1cm[1]-rs2cm[1])
                       + (rs1cm[2]-rs2cm[2])*(rs1cm[2]-rs2cm[2]))
                  + sqrt((rs1cm[0]-rs3cm[0])*(rs1cm[0]-rs3cm[0])
                       + (rs1cm[1]-rs3cm[1])*(rs1cm[1]-rs3cm[1])
                       + (rs1cm[2]-rs3cm[2])*(rs1cm[2]-rs3cm[2]))
                  + sqrt((rs1cm[0]-rs4cm[0])*(rs1cm[0]-rs4cm[0])
                       + (rs1cm[1]-rs4cm[1])*(rs1cm[1]-rs4cm[1])
                       + (rs1cm[2]-rs4cm[2])*(rs1cm[2]-rs4cm[2]))
                  + sqrt((rs2cm[0]-rs3cm[0])*(rs2cm[0]-rs3cm[0])
                       + (rs2cm[1]-rs3cm[1])*(rs2cm[1]-rs3cm[1])
                       + (rs2cm[2]-rs3cm[2])*(rs2cm[2]-rs3cm[2]))
                  + sqrt((rs2cm[0]-rs4cm[0])*(rs2cm[0]-rs4cm[0])
                       + (rs2cm[1]-rs4cm[1])*(rs2cm[1]-rs4cm[1])
                       + (rs2cm[2]-rs4cm[2])*(rs2cm[2]-rs4cm[2]))
                  + sqrt((rs3cm[0]-rs4cm[0])*(rs3cm[0]-rs4cm[0])
                       + (rs3cm[1]-rs4cm[1])*(rs3cm[1]-rs4cm[1])
                       + (rs3cm[2]-rs4cm[2])*(rs3cm[2]-rs4cm[2]));
            dist += distt;
            dist4 += distt;

            count_d += 6;
            count_d4 += 6;
            break;
        default:
            cerr << "hadron have more than 4 quarks!" << endl;
            exit(0);
        }
    }
    dists.push_back(dist2/count_d2);
    dists.push_back(dist3/count_d3);
    dists.push_back(dist4/count_d4);
    return dist/count_d;
}

double calculateMomentDistance(vector<double>& dists) {
    using namespace space;

    valarray<double> ps1(3), ps2(3), ps3(3), ps4(3);
    valarray<double> ps1cm(3), ps2cm(3), ps3cm(3), ps4cm(3);
    double es1, es2, es3, es4, es1cm, es2cm, es3cm, es4cm;
    valarray<double> _beta(3);

    int count_d = 0;
    int count_d2 = 0, count_d3 = 0, count_d4 = 0;
    int q1, q2, q3, q4;
    double dist = 0.0, distt;
    double dist2 = 0.0, dist3 = 0.0, dist4 = 0.0;
    for(auto hadron:hadrons) {
        switch(hadron.size()) {
        case 2:
            q1 = hadron[0];
            q2 = hadron[1];

            es1 = es[q1], es2 = es[q2];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];

            _beta = (ps1 + ps2)/(es1 + es2);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);

            distt = sqrt((ps1cm[0]-ps2cm[0])*(ps1cm[0]-ps2cm[0])
                       + (ps1cm[1]-ps2cm[1])*(ps1cm[1]-ps2cm[1])
                       + (ps1cm[2]-ps2cm[2])*(ps1cm[2]-ps2cm[2]));
            dist += distt;
            dist2 += distt;

            count_d += 1;
            count_d2 += 1;
            break;
        case 3:
            q1 = hadron[0];
            q2 = hadron[1];
            q3 = hadron[2];

            es1 = es[q1], es2 = es[q2], es3 = es[q3];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
            ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];

            _beta = (ps1 + ps2 + ps3)/(es1 + es2 + es3);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);
            lorentz(es3, ps3, es3cm, ps3cm, _beta);

            distt = sqrt((ps1cm[0]-ps2cm[0])*(ps1cm[0]-ps2cm[0])
                       + (ps1cm[1]-ps2cm[1])*(ps1cm[1]-ps2cm[1])
                       + (ps1cm[2]-ps2cm[2])*(ps1cm[2]-ps2cm[2]))
                  + sqrt((ps1cm[0]-ps3cm[0])*(ps1cm[0]-ps3cm[0])
                       + (ps1cm[1]-ps3cm[1])*(ps1cm[1]-ps3cm[1])
                       + (ps1cm[2]-ps3cm[2])*(ps1cm[2]-ps3cm[2]))
                  + sqrt((ps2cm[0]-ps3cm[0])*(ps2cm[0]-ps3cm[0])
                       + (ps2cm[1]-ps3cm[1])*(ps2cm[1]-ps3cm[1])
                       + (ps2cm[2]-ps3cm[2])*(ps2cm[2]-ps3cm[2]));
            dist += distt;
            dist3 += distt;

            count_d += 3;
            count_d3 += 3;
            break;
        case 4:
            // tetraquark
            q1 = hadron[0];
            q2 = hadron[1];
            q3 = hadron[2];
            q4 = hadron[3];

            es1 = es[q1], es2 = es[q2], es3 = es[q3], es4 = es[q4];
            ps1[0] = pxs[q1]; ps1[1] = pys[q1]; ps1[2] = pzs[q1];
            ps2[0] = pxs[q2]; ps2[1] = pys[q2]; ps2[2] = pzs[q2];
            ps3[0] = pxs[q3]; ps3[1] = pys[q3]; ps3[2] = pzs[q3];
            ps4[0] = pxs[q4]; ps4[1] = pys[q4]; ps4[2] = pzs[q4];

            _beta = (ps1 + ps2 + ps3 + ps4)/(es1 + es2 + es3 + es4);

            lorentz(es1, ps1, es1cm, ps1cm, _beta);
            lorentz(es2, ps2, es2cm, ps2cm, _beta);
            lorentz(es3, ps3, es3cm, ps3cm, _beta);
            lorentz(es4, ps4, es4cm, ps4cm, _beta);

            distt = sqrt((ps1cm[0]-ps2cm[0])*(ps1cm[0]-ps2cm[0])
                       + (ps1cm[1]-ps2cm[1])*(ps1cm[1]-ps2cm[1])
                       + (ps1cm[2]-ps2cm[2])*(ps1cm[2]-ps2cm[2]))
                  + sqrt((ps1cm[0]-ps3cm[0])*(ps1cm[0]-ps3cm[0])
                       + (ps1cm[1]-ps3cm[1])*(ps1cm[1]-ps3cm[1])
                       + (ps1cm[2]-ps3cm[2])*(ps1cm[2]-ps3cm[2]))
                  + sqrt((ps1cm[0]-ps4cm[0])*(ps1cm[0]-ps4cm[0])
                       + (ps1cm[1]-ps4cm[1])*(ps1cm[1]-ps4cm[1])
                       + (ps1cm[2]-ps4cm[2])*(ps1cm[2]-ps4cm[2]))
                  + sqrt((ps2cm[0]-ps3cm[0])*(ps2cm[0]-ps3cm[0])
                       + (ps2cm[1]-ps3cm[1])*(ps2cm[1]-ps3cm[1])
                       + (ps2cm[2]-ps3cm[2])*(ps2cm[2]-ps3cm[2]))
                  + sqrt((ps2cm[0]-ps4cm[0])*(ps2cm[0]-ps4cm[0])
                       + (ps2cm[1]-ps4cm[1])*(ps2cm[1]-ps4cm[1])
                       + (ps2cm[2]-ps4cm[2])*(ps2cm[2]-ps4cm[2]))
                  + sqrt((ps3cm[0]-ps4cm[0])*(ps3cm[0]-ps4cm[0])
                       + (ps3cm[1]-ps4cm[1])*(ps3cm[1]-ps4cm[1])
                       + (ps3cm[2]-ps4cm[2])*(ps3cm[2]-ps4cm[2]));
            dist += distt;
            dist4 += distt;

            count_d += 6;
            count_d4 += 6;
            break;
        default:
            cerr << "hadron have more than 4 quarks!" << endl;
            exit(0);
        }
    }
    dists.push_back(dist2/count_d2);
    dists.push_back(dist3/count_d3);
    dists.push_back(dist4/count_d4);
    return dist/count_d;
}

bool readZPCData(ifstream& zpc) {
    using namespace parameters;
    using namespace space;
    using namespace flags;

    int temp_int;
    zpc >> temp_int >> temp_int >> MUL
        >> bimp >> NELP >> NINP >> NELT >> NINT
        >> phiRP >> NSG >> BaryonNum >> MesonNum;

    if(zpc.eof()) return false;

/*
        0   1   2   3   4  5  6  7  8  9  10
        ib, px, py, pz, m, x, y, z, t, i, j
*/
    Quark_num = 0; AntiQuark_num = 0;

    int type;
    double px, py, pz, m, x, y, z, t;

    quarks_num = 0; antiquarks_num = 0;
    ds_num = 0; us_num = 0; ss_num = 0; cs_num = 0; bs_num = 0;
    dbs_num = 0; ubs_num = 0; sbs_num = 0; cbs_num = 0; bbs_num = 0;
    types.resize(MUL); es.resize(MUL); orders.resize(MUL); flags_all.resize(MUL);
    pxs.resize(MUL); pys.resize(MUL); pzs.resize(MUL); ms.resize(MUL);
    xs.resize(MUL);   ys.resize(MUL);  zs.resize(MUL); ts.resize(MUL);

    quarks_index.clear();
    antiquarks_index.clear();
    ds_index.clear(); dbs_index.clear();
    us_index.clear(); ubs_index.clear();
    ss_index.clear(); sbs_index.clear();
    cs_index.clear(); cbs_index.clear();
    bs_index.clear(); bbs_index.clear();
    for(int i=0;i<MUL;++i) {
        zpc >> type >> px >> py >> pz >> m >> x  >> y  >> z  >> t  // data
            >> temp_int >> temp_int;

        types[i] = type; orders[i] = i;
        pxs[i] = px; pys[i] = py; pzs[i] = pz; ms[i] = m;
        xs[i]  = x;  ys[i]  = y;  zs[i]  = z;  ts[i] = t;
        es[i] = sqrt(m*m + px*px + py*py + pz*pz);
        if(type > 0) {
            Quark_num++;
            quarks_num++;
            quarks_index.push_back(i);
        }
        if(type < 0) {
            AntiQuark_num++;
            antiquarks_num++;
            antiquarks_index.push_back(i);
        }

        switch(type) {
        case 1:
            ds_index.push_back(i); ds_num++; break;
        case 2:
            us_index.push_back(i); us_num++; break;
        case 3:
            ss_index.push_back(i); ss_num++; break;
        case 4:
            cs_index.push_back(i); cs_num++; break;
        case 5:
            bs_index.push_back(i); bs_num++; break;
        case -1:
            dbs_index.push_back(i); dbs_num++; break;
        case -2:
            ubs_index.push_back(i); ubs_num++; break;
        case -3:
            sbs_index.push_back(i); sbs_num++; break;
        case -4:
            cbs_index.push_back(i); cbs_num++; break;
        case -5:
            bbs_index.push_back(i); bbs_num++; break;
        default:
            cerr << "quark flavour number cannot be " << type << "!" << endl;
            exit(0);
        }
    }

    cout << "particles:         " << MUL << endl;
    cout << "meson number:      " << MesonNum << endl;
    cout << "baryon number:     " << BaryonNum << endl;
    cout << "quark number:      " << Quark_num << endl;
    cout << "anti-quark number: " << AntiQuark_num << endl;

    return true;
}

void outputHadronData(int event, ofstream& out) {
    using space::hadrons;
    using space::orders;

    out << event << "  " << hadrons.size() << endl;
    for(unsigned int i=0;i<hadrons.size();++i) {
        out << setw(7) << i+1 << setw(7) << hadrons[i].size();
        for(unsigned int j=0;j<hadrons[i].size();++j) {
            out << setw(7) << orders[hadrons[i][j]]+1; // initial order
        }
        out << endl;
    }
}
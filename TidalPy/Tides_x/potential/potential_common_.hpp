#pragma once

#include <cstdint>
#include <vector>

#include "intmap_.hpp"
#include "keys_.hpp"
#include "numerics_.hpp"
#include "constants_.hpp"

struct c_FrequencyStorage
{
    size_t num_instances;
    double frequency;

    c_FrequencyStorage(double freq) :
        num_instances(1),
        frequency(freq)
    {
    }
};

struct c_ModeStorage
{
    double mode;
    double mode_strength;
    int n_coeff;
    int o_coeff;

    c_ModeStorage(double mode_, double mode_strength_, int n_coeff_, int o_coeff_) :
        mode(mode_),
        mode_strength(mode_strength_),
        n_coeff(n_coeff_),
        o_coeff(o_coeff_)
    {
    }

    c_ModeStorage(int n_coeff_, int o_coeff_) :
        n_coeff(n_coeff_),
        o_coeff(o_coeff_)
    {
    }

};

typedef c_IntMap<c_Key4, c_ModeStorage> c_ModeMap;
typedef c_IntMap<c_Key4, size_t> c_UniqueFreqIndexMap;
typedef std::vector<c_FrequencyStorage> c_UniqueFreqMap;


bool record_unique_frequencies(
        c_Key4& lmpq_key,
        double frequency,
        c_UniqueFreqIndexMap& frequency_index_map,
        c_UniqueFreqMap& frequency_map)
{
    bool nonzero_freq = not c_isclose(frequency, 0.0, 1.0e-9, tidalpy_config_ptr->d_MIN_FREQUENCY);

    // TODO: Do we want to keep zero frequencies? I don't think so...
    if (nonzero_freq)
    {    
        // Checks if this frequency has been recorded yet. If not then it adds it to the frequency map. 
        bool found = false;
        c_FrequencyStorage* frequency_storage_ptr = nullptr;
        for (size_t i = 0; i < frequency_map.size(); i++)
        {
            frequency_storage_ptr = &frequency_map[i];
            if (c_isclose(frequency, frequency_storage_ptr->frequency, 1.0e-9, 0))
            {   
                // Increment the number of times this frequency has shown up.
                frequency_storage_ptr->num_instances += 1;
                // Record that this mode corresponds to this unique frequency.
                frequency_index_map.set(lmpq_key, i);
                found = true;
                // Break out of the frequency loop.
                break;
            }
        }

        if (!found)
        {
            // Did not find the frequency, so this is a unique one!
            frequency_index_map.set(lmpq_key, frequency_map.size());
            frequency_map.emplace_back(frequency);
        }
    }

    return nonzero_freq;
}

inline c_IntMap<c_Key2, double>& c_get_lm_coeff_map() {
    static c_IntMap<c_Key2, double> lm_coeff_map;
    static bool initialized = false;

    if (!initialized) {
        lm_coeff_map.reserve(512);
        lm_coeff_map.set(c_Key2(1, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(1, 1), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(2, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(2, 1), 3.333333333333333e-1);
        lm_coeff_map.set(c_Key2(2, 2), 8.333333333333333e-2);
        lm_coeff_map.set(c_Key2(3, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(3, 1), 1.666666666666667e-1);
        lm_coeff_map.set(c_Key2(3, 2), 1.666666666666667e-2);
        lm_coeff_map.set(c_Key2(3, 3), 2.777777777777778e-3);
        lm_coeff_map.set(c_Key2(4, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(4, 1), 1.000000000000000e-1);
        lm_coeff_map.set(c_Key2(4, 2), 5.555555555555556e-3);
        lm_coeff_map.set(c_Key2(4, 3), 3.968253968253968e-4);
        lm_coeff_map.set(c_Key2(4, 4), 4.960317460317460e-5);
        lm_coeff_map.set(c_Key2(5, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(5, 1), 6.666666666666667e-2);
        lm_coeff_map.set(c_Key2(5, 2), 2.380952380952381e-3);
        lm_coeff_map.set(c_Key2(5, 3), 9.920634920634921e-5);
        lm_coeff_map.set(c_Key2(5, 4), 5.511463844797178e-6);
        lm_coeff_map.set(c_Key2(5, 5), 5.511463844797178e-7);
        lm_coeff_map.set(c_Key2(6, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(6, 1), 4.761904761904762e-2);
        lm_coeff_map.set(c_Key2(6, 2), 1.190476190476190e-3);
        lm_coeff_map.set(c_Key2(6, 3), 3.306878306878307e-5);
        lm_coeff_map.set(c_Key2(6, 4), 1.102292768959436e-6);
        lm_coeff_map.set(c_Key2(6, 5), 5.010421677088344e-8);
        lm_coeff_map.set(c_Key2(6, 6), 4.175351397573620e-9);
        lm_coeff_map.set(c_Key2(7, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(7, 1), 3.571428571428571e-2);
        lm_coeff_map.set(c_Key2(7, 2), 6.613756613756614e-4);
        lm_coeff_map.set(c_Key2(7, 3), 1.322751322751323e-5);
        lm_coeff_map.set(c_Key2(7, 4), 3.006253006253006e-7);
        lm_coeff_map.set(c_Key2(7, 5), 8.350702795147240e-9);
        lm_coeff_map.set(c_Key2(7, 6), 3.211808767364323e-10);
        lm_coeff_map.set(c_Key2(7, 7), 2.294149119545945e-11);
        lm_coeff_map.set(c_Key2(8, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(8, 1), 2.777777777777778e-2);
        lm_coeff_map.set(c_Key2(8, 2), 3.968253968253968e-4);
        lm_coeff_map.set(c_Key2(8, 3), 6.012506012506013e-6);
        lm_coeff_map.set(c_Key2(8, 4), 1.002084335417669e-7);
        lm_coeff_map.set(c_Key2(8, 5), 1.927085260418594e-9);
        lm_coeff_map.set(c_Key2(8, 6), 4.588298239091890e-11);
        lm_coeff_map.set(c_Key2(8, 7), 1.529432746363963e-12);
        lm_coeff_map.set(c_Key2(8, 8), 9.558954664774771e-14);
        lm_coeff_map.set(c_Key2(9, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(9, 1), 2.222222222222222e-2);
        lm_coeff_map.set(c_Key2(9, 2), 2.525252525252525e-4);
        lm_coeff_map.set(c_Key2(9, 3), 3.006253006253006e-6);
        lm_coeff_map.set(c_Key2(9, 4), 3.854170520837188e-8);
        lm_coeff_map.set(c_Key2(9, 5), 5.505957886910268e-10);
        lm_coeff_map.set(c_Key2(9, 6), 9.176596478183780e-12);
        lm_coeff_map.set(c_Key2(9, 7), 1.911790932954954e-13);
        lm_coeff_map.set(c_Key2(9, 8), 5.622914508691042e-15);
        lm_coeff_map.set(c_Key2(9, 9), 3.123841393717245e-16);
        lm_coeff_map.set(c_Key2(10, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(10, 1), 1.818181818181818e-2);
        lm_coeff_map.set(c_Key2(10, 2), 1.683501683501684e-4);
        lm_coeff_map.set(c_Key2(10, 3), 1.618751618751619e-6);
        lm_coeff_map.set(c_Key2(10, 4), 1.651787366073080e-8);
        lm_coeff_map.set(c_Key2(10, 5), 1.835319295636756e-10);
        lm_coeff_map.set(c_Key2(10, 6), 2.294149119545945e-12);
        lm_coeff_map.set(c_Key2(10, 7), 3.373748705214625e-14);
        lm_coeff_map.set(c_Key2(10, 8), 6.247682787434491e-16);
        lm_coeff_map.set(c_Key2(10, 9), 1.644127049324866e-17);
        lm_coeff_map.set(c_Key2(10, 10), 8.220635246624330e-19);
        lm_coeff_map.set(c_Key2(11, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(11, 1), 1.515151515151515e-2);
        lm_coeff_map.set(c_Key2(11, 2), 1.165501165501166e-4);
        lm_coeff_map.set(c_Key2(11, 3), 9.250009250009250e-7);
        lm_coeff_map.set(c_Key2(11, 4), 7.708341041674375e-9);
        lm_coeff_map.set(c_Key2(11, 5), 6.882447358637835e-11);
        lm_coeff_map.set(c_Key2(11, 6), 6.747497410429250e-13);
        lm_coeff_map.set(c_Key2(11, 7), 7.497219344921389e-15);
        lm_coeff_map.set(c_Key2(11, 8), 9.864762295949196e-17);
        lm_coeff_map.set(c_Key2(11, 9), 1.644127049324866e-18);
        lm_coeff_map.set(c_Key2(11, 10), 3.914588212678252e-20);
        lm_coeff_map.set(c_Key2(11, 11), 1.779358278490115e-21);
        lm_coeff_map.set(c_Key2(12, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(12, 1), 1.282051282051282e-2);
        lm_coeff_map.set(c_Key2(12, 2), 8.325008325008325e-5);
        lm_coeff_map.set(c_Key2(12, 3), 5.550005550005550e-7);
        lm_coeff_map.set(c_Key2(12, 4), 3.854170520837188e-9);
        lm_coeff_map.set(c_Key2(12, 5), 2.833948912380285e-11);
        lm_coeff_map.set(c_Key2(12, 6), 2.249165803476417e-13);
        lm_coeff_map.set(c_Key2(12, 7), 1.972952459189839e-15);
        lm_coeff_map.set(c_Key2(12, 8), 1.972952459189839e-17);
        lm_coeff_map.set(c_Key2(12, 9), 2.348752927606951e-19);
        lm_coeff_map.set(c_Key2(12, 10), 3.558716556980229e-21);
        lm_coeff_map.set(c_Key2(12, 11), 7.736340341261368e-23);
        lm_coeff_map.set(c_Key2(12, 12), 3.223475142192237e-24);
        lm_coeff_map.set(c_Key2(13, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(13, 1), 1.098901098901099e-2);
        lm_coeff_map.set(c_Key2(13, 2), 6.105006105006105e-5);
        lm_coeff_map.set(c_Key2(13, 3), 3.468753468753469e-7);
        lm_coeff_map.set(c_Key2(13, 4), 2.040443216913805e-9);
        lm_coeff_map.set(c_Key2(13, 5), 1.259532849946793e-11);
        lm_coeff_map.set(c_Key2(13, 6), 8.286400328597324e-14);
        lm_coeff_map.set(c_Key2(13, 7), 5.918857377569517e-16);
        lm_coeff_map.set(c_Key2(13, 8), 4.697505855213903e-18);
        lm_coeff_map.set(c_Key2(13, 9), 4.270459868376275e-20);
        lm_coeff_map.set(c_Key2(13, 10), 4.641804204756821e-22);
        lm_coeff_map.set(c_Key2(13, 11), 6.446950284384473e-24);
        lm_coeff_map.set(c_Key2(13, 12), 1.289390056876895e-25);
        lm_coeff_map.set(c_Key2(13, 13), 4.959192526449595e-27);
        lm_coeff_map.set(c_Key2(14, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(14, 1), 9.523809523809524e-3);
        lm_coeff_map.set(c_Key2(14, 2), 4.578754578754579e-5);
        lm_coeff_map.set(c_Key2(14, 3), 2.244487538605186e-7);
        lm_coeff_map.set(c_Key2(14, 4), 1.133579564952114e-9);
        lm_coeff_map.set(c_Key2(14, 5), 5.966208236590074e-12);
        lm_coeff_map.set(c_Key2(14, 6), 3.314560131438930e-14);
        lm_coeff_map.set(c_Key2(14, 7), 1.972952459189839e-16);
        lm_coeff_map.set(c_Key2(14, 8), 1.281137960512883e-18);
        lm_coeff_map.set(c_Key2(14, 9), 9.283608409513642e-21);
        lm_coeff_map.set(c_Key2(14, 10), 7.736340341261368e-23);
        lm_coeff_map.set(c_Key2(14, 11), 7.736340341261368e-25);
        lm_coeff_map.set(c_Key2(14, 12), 9.918385052899190e-27);
        lm_coeff_map.set(c_Key2(14, 13), 1.836737972759109e-28);
        lm_coeff_map.set(c_Key2(14, 14), 6.559778474139676e-30);
        lm_coeff_map.set(c_Key2(15, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(15, 1), 8.333333333333333e-3);
        lm_coeff_map.set(c_Key2(15, 2), 3.501400560224090e-5);
        lm_coeff_map.set(c_Key2(15, 3), 1.496325025736790e-7);
        lm_coeff_map.set(c_Key2(15, 4), 6.562829060249081e-10);
        lm_coeff_map.set(c_Key2(15, 5), 2.983104118295037e-12);
        lm_coeff_map.set(c_Key2(15, 6), 1.420525770616684e-14);
        lm_coeff_map.set(c_Key2(15, 7), 7.174372578872142e-17);
        lm_coeff_map.set(c_Key2(15, 8), 3.899115531995730e-19);
        lm_coeff_map.set(c_Key2(15, 9), 2.320902102378410e-21);
        lm_coeff_map.set(c_Key2(15, 10), 1.547268068252274e-23);
        lm_coeff_map.set(c_Key2(15, 11), 1.190206206347903e-25);
        lm_coeff_map.set(c_Key2(15, 12), 1.102042783655466e-27);
        lm_coeff_map.set(c_Key2(15, 13), 1.311955694827935e-29);
        lm_coeff_map.set(c_Key2(15, 14), 2.261992577289543e-31);
        lm_coeff_map.set(c_Key2(15, 15), 7.539975257631811e-33);
        lm_coeff_map.set(c_Key2(16, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(16, 1), 7.352941176470588e-3);
        lm_coeff_map.set(c_Key2(16, 2), 2.723311546840959e-5);
        lm_coeff_map.set(c_Key2(16, 3), 1.023801333398857e-7);
        lm_coeff_map.set(c_Key2(16, 4), 3.937697436149449e-10);
        lm_coeff_map.set(c_Key2(16, 5), 1.562578347678353e-12);
        lm_coeff_map.set(c_Key2(16, 6), 6.456935320984928e-15);
        lm_coeff_map.set(c_Key2(16, 7), 2.807363183036925e-17);
        lm_coeff_map.set(c_Key2(16, 8), 1.299705177331910e-19);
        lm_coeff_map.set(c_Key2(16, 9), 6.498525886659549e-22);
        lm_coeff_map.set(c_Key2(16, 10), 3.570618619043708e-24);
        lm_coeff_map.set(c_Key2(16, 11), 2.204085567310931e-26);
        lm_coeff_map.set(c_Key2(16, 12), 1.574346833793522e-28);
        lm_coeff_map.set(c_Key2(16, 13), 1.357195546373726e-30);
        lm_coeff_map.set(c_Key2(16, 14), 1.507995051526362e-32);
        lm_coeff_map.set(c_Key2(16, 15), 2.432250083107036e-34);
        lm_coeff_map.set(c_Key2(16, 16), 7.600781509709487e-36);
        lm_coeff_map.set(c_Key2(17, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(17, 1), 6.535947712418301e-3);
        lm_coeff_map.set(c_Key2(17, 2), 2.149982800137599e-5);
        lm_coeff_map.set(c_Key2(17, 3), 7.166609333791996e-8);
        lm_coeff_map.set(c_Key2(17, 4), 2.437622222378230e-10);
        lm_coeff_map.set(c_Key2(17, 5), 8.523154623700105e-13);
        lm_coeff_map.set(c_Key2(17, 6), 3.088099501340618e-15);
        lm_coeff_map.set(c_Key2(17, 7), 1.169734659598719e-17);
        lm_coeff_map.set(c_Key2(17, 8), 4.678938638394875e-20);
        lm_coeff_map.set(c_Key2(17, 9), 1.999546426664477e-22);
        lm_coeff_map.set(c_Key2(17, 10), 9.257159382705911e-25);
        lm_coeff_map.set(c_Key2(17, 11), 4.723040501380567e-27);
        lm_coeff_map.set(c_Key2(17, 12), 2.714391092747452e-29);
        lm_coeff_map.set(c_Key2(17, 13), 1.809594061831635e-31);
        lm_coeff_map.set(c_Key2(17, 14), 1.459350049864222e-33);
        lm_coeff_map.set(c_Key2(17, 15), 1.520156301941897e-35);
        lm_coeff_map.set(c_Key2(17, 16), 2.303267124154390e-37);
        lm_coeff_map.set(c_Key2(17, 17), 6.774315071042324e-39);
        lm_coeff_map.set(c_Key2(18, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(18, 1), 5.847953216374269e-3);
        lm_coeff_map.set(c_Key2(18, 2), 1.719986240110079e-5);
        lm_coeff_map.set(c_Key2(18, 3), 5.119006666994283e-8);
        lm_coeff_map.set(c_Key2(18, 4), 1.551214141513419e-10);
        lm_coeff_map.set(c_Key2(18, 5), 4.817435222091364e-13);
        lm_coeff_map.set(c_Key2(18, 6), 1.544049750670309e-15);
        lm_coeff_map.set(c_Key2(18, 7), 5.146832502234363e-18);
        lm_coeff_map.set(c_Key2(18, 8), 1.799591783998029e-20);
        lm_coeff_map.set(c_Key2(18, 9), 6.665154755548256e-23);
        lm_coeff_map.set(c_Key2(18, 10), 2.644902680773117e-25);
        lm_coeff_map.set(c_Key2(18, 11), 1.140044258953930e-27);
        lm_coeff_map.set(c_Key2(18, 12), 5.428782185494904e-30);
        lm_coeff_map.set(c_Key2(18, 13), 2.918700099728443e-32);
        lm_coeff_map.set(c_Key2(18, 14), 1.824187562330277e-34);
        lm_coeff_map.set(c_Key2(18, 15), 1.381960274492634e-36);
        lm_coeff_map.set(c_Key2(18, 16), 1.354863014208465e-38);
        lm_coeff_map.set(c_Key2(18, 17), 1.935518591726378e-40);
        lm_coeff_map.set(c_Key2(18, 18), 5.376440532573273e-42);
        lm_coeff_map.set(c_Key2(19, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(19, 1), 5.263157894736842e-3);
        lm_coeff_map.set(c_Key2(19, 2), 1.392369813422445e-5);
        lm_coeff_map.set(c_Key2(19, 3), 3.722913939632206e-8);
        lm_coeff_map.set(c_Key2(19, 4), 1.011661396639186e-10);
        lm_coeff_map.set(c_Key2(19, 5), 2.810170546219962e-13);
        lm_coeff_map.set(c_Key2(19, 6), 8.029058703485606e-16);
        lm_coeff_map.set(c_Key2(19, 7), 2.375461154877398e-18);
        lm_coeff_map.set(c_Key2(19, 8), 7.331670231103081e-21);
        lm_coeff_map.set(c_Key2(19, 9), 2.380412412695806e-23);
        lm_coeff_map.set(c_Key2(19, 10), 8.208318664468295e-26);
        lm_coeff_map.set(c_Key2(19, 11), 3.040118023877146e-28);
        lm_coeff_map.set(c_Key2(19, 12), 1.225854041885946e-30);
        lm_coeff_map.set(c_Key2(19, 13), 5.472562686990831e-33);
        lm_coeff_map.set(c_Key2(19, 14), 2.763920548985268e-35);
        lm_coeff_map.set(c_Key2(19, 15), 1.625835617050158e-37);
        lm_coeff_map.set(c_Key2(19, 16), 1.161311155035827e-39);
        lm_coeff_map.set(c_Key2(19, 17), 1.075288106514655e-41);
        lm_coeff_map.set(c_Key2(19, 18), 1.453092035830614e-43);
        lm_coeff_map.set(c_Key2(19, 19), 3.823926410080564e-45);
        lm_coeff_map.set(c_Key2(20, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(20, 1), 4.761904761904762e-3);
        lm_coeff_map.set(c_Key2(20, 2), 1.139211665527455e-5);
        lm_coeff_map.set(c_Key2(20, 3), 2.751718998858587e-8);
        lm_coeff_map.set(c_Key2(20, 4), 6.744409310927909e-11);
        lm_coeff_map.set(c_Key2(20, 5), 1.686102327731977e-13);
        lm_coeff_map.set(c_Key2(20, 6), 4.323339301876865e-16);
        lm_coeff_map.set(c_Key2(20, 7), 1.143740556052081e-18);
        lm_coeff_map.set(c_Key2(20, 8), 3.142144384758463e-21);
        lm_coeff_map.set(c_Key2(20, 9), 9.029150530915125e-24);
        lm_coeff_map.set(c_Key2(20, 10), 2.736106221489432e-26);
        lm_coeff_map.set(c_Key2(20, 11), 8.826149101578812e-29);
        lm_coeff_map.set(c_Key2(20, 12), 3.064635104714865e-31);
        lm_coeff_map.set(c_Key2(20, 13), 1.160846630573813e-33);
        lm_coeff_map.set(c_Key2(20, 14), 4.877506851150473e-36);
        lm_coeff_map.set(c_Key2(20, 15), 2.322622310071654e-38);
        lm_coeff_map.set(c_Key2(20, 16), 1.290345727817585e-40);
        lm_coeff_map.set(c_Key2(20, 17), 8.718552214983686e-43);
        lm_coeff_map.set(c_Key2(20, 18), 7.647852820161128e-45);
        lm_coeff_map.set(c_Key2(20, 19), 9.804939513027087e-47);
        lm_coeff_map.set(c_Key2(20, 20), 2.451234878256772e-48);
        lm_coeff_map.set(c_Key2(21, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(21, 1), 4.329004329004329e-3);
        lm_coeff_map.set(c_Key2(21, 2), 9.410878976096367e-6);
        lm_coeff_map.set(c_Key2(21, 3), 2.063789249143940e-8);
        lm_coeff_map.set(c_Key2(21, 4), 4.586198331430978e-11);
        lm_coeff_map.set(c_Key2(21, 5), 1.037601432450448e-13);
        lm_coeff_map.set(c_Key2(21, 6), 2.401855167709369e-16);
        lm_coeff_map.set(c_Key2(21, 7), 5.718702780260403e-19);
        lm_coeff_map.set(c_Key2(21, 8), 1.408547482822759e-21);
        lm_coeff_map.set(c_Key2(21, 9), 3.611660212366050e-24);
        lm_coeff_map.set(c_Key2(21, 10), 9.708764011736693e-27);
        lm_coeff_map.set(c_Key2(21, 11), 2.758171594243379e-29);
        lm_coeff_map.set(c_Key2(21, 12), 8.358095740131451e-32);
        lm_coeff_map.set(c_Key2(21, 13), 2.731403836644265e-34);
        lm_coeff_map.set(c_Key2(21, 14), 9.755013702300946e-37);
        lm_coeff_map.set(c_Key2(21, 15), 3.871037183452756e-39);
        lm_coeff_map.set(c_Key2(21, 16), 1.743710442996737e-41);
        lm_coeff_map.set(c_Key2(21, 17), 9.177423384193353e-44);
        lm_coeff_map.set(c_Key2(21, 18), 5.882963707816252e-46);
        lm_coeff_map.set(c_Key2(21, 19), 4.902469756513543e-48);
        lm_coeff_map.set(c_Key2(21, 20), 5.978621654284809e-50);
        lm_coeff_map.set(c_Key2(21, 21), 1.423481346258288e-51);
        lm_coeff_map.set(c_Key2(22, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(22, 1), 3.952569169960474e-3);
        lm_coeff_map.set(c_Key2(22, 2), 7.842399146746973e-6);
        lm_coeff_map.set(c_Key2(22, 3), 1.568479829349395e-8);
        lm_coeff_map.set(c_Key2(22, 4), 3.175060383298370e-11);
        lm_coeff_map.set(c_Key2(22, 5), 6.533046056169485e-14);
        lm_coeff_map.set(c_Key2(22, 6), 1.372488667262497e-16);
        lm_coeff_map.set(c_Key2(22, 7), 2.957949713927795e-19);
        lm_coeff_map.set(c_Key2(22, 8), 6.573221586506211e-22);
        lm_coeff_map.set(c_Key2(22, 9), 1.514567185830924e-24);
        lm_coeff_map.set(c_Key2(22, 10), 3.640786504401260e-27);
        lm_coeff_map.set(c_Key2(22, 11), 9.193905314144596e-30);
        lm_coeff_map.set(c_Key2(22, 12), 2.458263452979838e-32);
        lm_coeff_map.set(c_Key2(22, 13), 7.023609865656681e-35);
        lm_coeff_map.set(c_Key2(22, 14), 2.167780822733544e-37);
        lm_coeff_map.set(c_Key2(22, 15), 7.323583860586296e-40);
        lm_coeff_map.set(c_Key2(22, 16), 2.753227015258006e-42);
        lm_coeff_map.set(c_Key2(22, 17), 1.176592741563250e-44);
        lm_coeff_map.set(c_Key2(22, 18), 5.882963707816252e-47);
        lm_coeff_map.set(c_Key2(22, 19), 3.587172992570885e-49);
        lm_coeff_map.set(c_Key2(22, 20), 2.846962692516576e-51);
        lm_coeff_map.set(c_Key2(22, 21), 3.310421735484390e-53);
        lm_coeff_map.set(c_Key2(22, 22), 7.523685762464524e-55);
        lm_coeff_map.set(c_Key2(23, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(23, 1), 3.623188405797101e-3);
        lm_coeff_map.set(c_Key2(23, 2), 6.587615283267457e-6);
        lm_coeff_map.set(c_Key2(23, 3), 1.206522945653380e-8);
        lm_coeff_map.set(c_Key2(23, 4), 2.234301751209964e-11);
        lm_coeff_map.set(c_Key2(23, 5), 4.199815321823240e-14);
        lm_coeff_map.set(c_Key2(23, 6), 8.045623221883602e-17);
        lm_coeff_map.set(c_Key2(23, 7), 1.577573180761491e-19);
        lm_coeff_map.set(c_Key2(23, 8), 3.180591090244941e-22);
        lm_coeff_map.set(c_Key2(23, 9), 6.626231438010293e-25);
        lm_coeff_map.set(c_Key2(23, 10), 1.434249229006557e-27);
        lm_coeff_map.set(c_Key2(23, 11), 3.244907757933387e-30);
        lm_coeff_map.set(c_Key2(23, 12), 7.725970852222349e-33);
        lm_coeff_map.set(c_Key2(23, 13), 1.951002740460189e-35);
        lm_coeff_map.set(c_Key2(23, 14), 5.272980379622133e-38);
        lm_coeff_map.set(c_Key2(23, 15), 1.541807128544483e-40);
        lm_coeff_map.set(c_Key2(23, 16), 4.941689514565652e-43);
        lm_coeff_map.set(c_Key2(23, 17), 1.764889112344876e-45);
        lm_coeff_map.set(c_Key2(23, 18), 7.174345985141771e-48);
        lm_coeff_map.set(c_Key2(23, 19), 3.416355231019891e-50);
        lm_coeff_map.set(c_Key2(23, 20), 1.986253041290634e-52);
        lm_coeff_map.set(c_Key2(23, 21), 1.504737152492905e-54);
        lm_coeff_map.set(c_Key2(23, 22), 1.671930169436561e-56);
        lm_coeff_map.set(c_Key2(23, 23), 3.634630803122958e-58);
        lm_coeff_map.set(c_Key2(24, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(24, 1), 3.333333333333333e-3);
        lm_coeff_map.set(c_Key2(24, 2), 5.574136008918618e-6);
        lm_coeff_map.set(c_Key2(24, 3), 9.384067355081848e-9);
        lm_coeff_map.set(c_Key2(24, 4), 1.595929822292831e-11);
        lm_coeff_map.set(c_Key2(24, 5), 2.751603141884192e-14);
        lm_coeff_map.set(c_Key2(24, 6), 4.827373933130161e-17);
        lm_coeff_map.set(c_Key2(24, 7), 8.651207765466239e-20);
        lm_coeff_map.set(c_Key2(24, 8), 1.590295545122470e-22);
        lm_coeff_map.set(c_Key2(24, 9), 3.011923380913770e-25);
        lm_coeff_map.set(c_Key2(24, 10), 5.905732119438764e-28);
        lm_coeff_map.set(c_Key2(24, 11), 1.205251452946686e-30);
        lm_coeff_map.set(c_Key2(24, 12), 2.575323617407450e-33);
        lm_coeff_map.set(c_Key2(24, 13), 5.800278417584346e-36);
        lm_coeff_map.set(c_Key2(24, 14), 1.387626415690035e-38);
        lm_coeff_map.set(c_Key2(24, 15), 3.558016450487269e-41);
        lm_coeff_map.set(c_Key2(24, 16), 9.883379029131303e-44);
        lm_coeff_map.set(c_Key2(24, 17), 3.013225313759544e-46);
        lm_coeff_map.set(c_Key2(24, 18), 1.024906569305967e-48);
        lm_coeff_map.set(c_Key2(24, 19), 3.972506082581268e-51);
        lm_coeff_map.set(c_Key2(24, 20), 1.805684582991486e-53);
        lm_coeff_map.set(c_Key2(24, 21), 1.003158101661936e-55);
        lm_coeff_map.set(c_Key2(24, 22), 7.269261606245917e-58);
        lm_coeff_map.set(c_Key2(24, 23), 7.733257027921188e-60);
        lm_coeff_map.set(c_Key2(24, 24), 1.611095214150247e-61);
        lm_coeff_map.set(c_Key2(25, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(25, 1), 3.076923076923077e-3);
        lm_coeff_map.set(c_Key2(25, 2), 4.748338081671415e-6);
        lm_coeff_map.set(c_Key2(25, 3), 7.373195778992880e-9);
        lm_coeff_map.set(c_Key2(25, 4), 1.155673319591361e-11);
        lm_coeff_map.set(c_Key2(25, 5), 1.834402094589461e-14);
        lm_coeff_map.set(c_Key2(25, 6), 2.958713055789454e-17);
        lm_coeff_map.set(c_Key2(25, 7), 4.866304368074759e-20);
        lm_coeff_map.set(c_Key2(25, 8), 8.192431596085453e-23);
        lm_coeff_map.set(c_Key2(25, 9), 1.417375708665303e-25);
        lm_coeff_map.set(c_Key2(25, 10), 2.531028051188042e-28);
        lm_coeff_map.set(c_Key2(25, 11), 4.687088983681559e-31);
        lm_coeff_map.set(c_Key2(25, 12), 9.048434331431580e-34);
        lm_coeff_map.set(c_Key2(25, 13), 1.831666868710846e-36);
        lm_coeff_map.set(c_Key2(25, 14), 3.913818095535996e-39);
        lm_coeff_map.set(c_Key2(25, 15), 8.895041126218173e-42);
        lm_coeff_map.set(c_Key2(25, 16), 2.169522225906871e-44);
        lm_coeff_map.set(c_Key2(25, 17), 5.739476788113417e-47);
        lm_coeff_map.set(c_Key2(25, 18), 1.668452554684133e-49);
        lm_coeff_map.set(c_Key2(25, 19), 5.417053748974457e-52);
        lm_coeff_map.set(c_Key2(25, 20), 2.006316203323873e-54);
        lm_coeff_map.set(c_Key2(25, 21), 8.723113927495100e-57);
        lm_coeff_map.set(c_Key2(25, 22), 4.639954216752713e-59);
        lm_coeff_map.set(c_Key2(25, 23), 3.222190428300495e-61);
        lm_coeff_map.set(c_Key2(25, 24), 3.287949416633158e-63);
        lm_coeff_map.set(c_Key2(25, 25), 6.575898833266316e-65);
        lm_coeff_map.set(c_Key2(26, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(26, 1), 2.849002849002849e-3);
        lm_coeff_map.set(c_Key2(26, 2), 4.070004070004070e-6);
        lm_coeff_map.set(c_Key2(26, 3), 5.847706997132284e-9);
        lm_coeff_map.set(c_Key2(26, 4), 8.474937677003311e-12);
        lm_coeff_map.set(c_Key2(26, 5), 1.242659483431571e-14);
        lm_coeff_map.set(c_Key2(26, 6), 1.849195659868408e-17);
        lm_coeff_map.set(c_Key2(26, 7), 2.801811605861225e-20);
        lm_coeff_map.set(c_Key2(26, 8), 4.337169668515828e-23);
        lm_coeff_map.set(c_Key2(26, 9), 6.884396299231473e-26);
        lm_coeff_map.set(c_Key2(26, 10), 1.124901356083574e-28);
        lm_coeff_map.set(c_Key2(26, 11), 1.900171209600632e-31);
        lm_coeff_map.set(c_Key2(26, 12), 3.333633701053740e-34);
        lm_coeff_map.set(c_Key2(26, 13), 6.105556229036154e-37);
        lm_coeff_map.set(c_Key2(26, 14), 1.174145428660799e-39);
        lm_coeff_map.set(c_Key2(26, 15), 2.386474448497559e-42);
        lm_coeff_map.set(c_Key2(26, 16), 5.165529109302075e-45);
        lm_coeff_map.set(c_Key2(26, 17), 1.201285839372576e-47);
        lm_coeff_map.set(c_Key2(26, 18), 3.033550099425696e-50);
        lm_coeff_map.set(c_Key2(26, 19), 8.426528053960266e-53);
        lm_coeff_map.set(c_Key2(26, 20), 2.616934178248530e-55);
        lm_coeff_map.set(c_Key2(26, 21), 9.279908433505425e-58);
        lm_coeff_map.set(c_Key2(26, 22), 3.866628513960594e-60);
        lm_coeff_map.set(c_Key2(26, 23), 1.972769649979895e-62);
        lm_coeff_map.set(c_Key2(26, 24), 1.315179766653263e-64);
        lm_coeff_map.set(c_Key2(26, 25), 1.289391928091435e-66);
        lm_coeff_map.set(c_Key2(26, 26), 2.479599861714297e-68);
        lm_coeff_map.set(c_Key2(27, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(27, 1), 2.645502645502646e-3);
        lm_coeff_map.set(c_Key2(27, 2), 3.508624198279371e-6);
        lm_coeff_map.set(c_Key2(27, 3), 4.678165597705828e-9);
        lm_coeff_map.set(c_Key2(27, 4), 6.287856986163747e-12);
        lm_coeff_map.set(c_Key2(27, 5), 8.543283948592047e-15);
        lm_coeff_map.set(c_Key2(27, 6), 1.176760874461714e-17);
        lm_coeff_map.set(c_Key2(27, 7), 1.648124474036015e-20);
        lm_coeff_map.set(c_Key2(27, 8), 2.354463534337164e-23);
        lm_coeff_map.set(c_Key2(27, 9), 3.442198149615737e-26);
        lm_coeff_map.set(c_Key2(27, 10), 5.168465690113719e-29);
        lm_coeff_map.set(c_Key2(27, 11), 8.000720882528976e-32);
        lm_coeff_map.set(c_Key2(27, 12), 1.282166808097592e-34);
        lm_coeff_map.set(c_Key2(27, 13), 2.136944680162654e-37);
        lm_coeff_map.set(c_Key2(27, 14), 3.722900139656191e-40);
        lm_coeff_map.set(c_Key2(27, 15), 6.818498424278739e-43);
        lm_coeff_map.set(c_Key2(27, 16), 1.321414423309833e-45);
        lm_coeff_map.set(c_Key2(27, 17), 2.730195089483126e-48);
        lm_coeff_map.set(c_Key2(27, 18), 6.067100198851392e-51);
        lm_coeff_map.set(c_Key2(27, 19), 1.465483139819177e-53);
        lm_coeff_map.set(c_Key2(27, 20), 3.897561542072279e-56);
        lm_coeff_map.set(c_Key2(27, 21), 1.159988554188178e-58);
        lm_coeff_map.set(c_Key2(27, 22), 3.945539299959790e-61);
        lm_coeff_map.set(c_Key2(27, 23), 1.578215719983916e-63);
        lm_coeff_map.set(c_Key2(27, 24), 7.736351568548607e-66);
        lm_coeff_map.set(c_Key2(27, 25), 4.959199723428594e-68);
        lm_coeff_map.set(c_Key2(27, 26), 4.678490305121315e-70);
        lm_coeff_map.set(c_Key2(27, 27), 8.663870935409843e-72);
        lm_coeff_map.set(c_Key2(28, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(28, 1), 2.463054187192118e-3);
        lm_coeff_map.set(c_Key2(28, 2), 3.040807638508788e-6);
        lm_coeff_map.set(c_Key2(28, 3), 3.772714191698248e-9);
        lm_coeff_map.set(c_Key2(28, 4), 4.715892739622810e-12);
        lm_coeff_map.set(c_Key2(28, 5), 5.954410024776275e-15);
        lm_coeff_map.set(c_Key2(28, 6), 7.614335070046388e-18);
        lm_coeff_map.set(c_Key2(28, 7), 9.888746844216088e-21);
        lm_coeff_map.set(c_Key2(28, 8), 1.308035296853980e-23);
        lm_coeff_map.set(c_Key2(28, 9), 1.767615266018892e-26);
        lm_coeff_map.set(c_Key2(28, 10), 2.448220590053867e-29);
        lm_coeff_map.set(c_Key2(28, 11), 3.487493718025451e-32);
        lm_coeff_map.set(c_Key2(28, 12), 5.128667232390369e-35);
        lm_coeff_map.set(c_Key2(28, 13), 7.818090293278002e-38);
        lm_coeff_map.set(c_Key2(28, 14), 1.240966713218730e-40);
        lm_coeff_map.set(c_Key2(28, 15), 2.061406500363340e-43);
        lm_coeff_map.set(c_Key2(28, 16), 3.603857518117727e-46);
        lm_coeff_map.set(c_Key2(28, 17), 6.673810218736531e-49);
        lm_coeff_map.set(c_Key2(28, 18), 1.318934825837259e-51);
        lm_coeff_map.set(c_Key2(28, 19), 2.806244310292041e-54);
        lm_coeff_map.set(c_Key2(28, 20), 6.495935903453798e-57);
        lm_coeff_map.set(c_Key2(28, 21), 1.657126505983112e-59);
        lm_coeff_map.set(c_Key2(28, 22), 4.734647159951748e-62);
        lm_coeff_map.set(c_Key2(28, 23), 1.547270313709721e-64);
        lm_coeff_map.set(c_Key2(28, 24), 5.951039668114313e-67);
        lm_coeff_map.set(c_Key2(28, 25), 2.807094183072789e-69);
        lm_coeff_map.set(c_Key2(28, 26), 1.732774187081969e-71);
        lm_coeff_map.set(c_Key2(28, 27), 1.575249260983608e-73);
        lm_coeff_map.set(c_Key2(28, 28), 2.812945108899300e-75);
        lm_coeff_map.set(c_Key2(29, 0), 1.000000000000000e+0);
        lm_coeff_map.set(c_Key2(29, 1), 2.298850574712644e-3);
        lm_coeff_map.set(c_Key2(29, 2), 2.648445362572170e-6);
        lm_coeff_map.set(c_Key2(29, 3), 3.065330280754827e-9);
        lm_coeff_map.set(c_Key2(29, 4), 3.572646014865765e-12);
        lm_coeff_map.set(c_Key2(29, 5), 4.203112958665606e-15);
        lm_coeff_map.set(c_Key2(29, 6), 5.003705903173341e-18);
        lm_coeff_map.set(c_Key2(29, 7), 6.043123071465387e-21);
        lm_coeff_map.set(c_Key2(29, 8), 7.423984117279345e-24);
        lm_coeff_map.set(c_Key2(29, 9), 9.303238242204694e-27);
        lm_coeff_map.set(c_Key2(29, 10), 1.192722851564704e-29);
        lm_coeff_map.set(c_Key2(29, 11), 1.569372173111453e-32);
        lm_coeff_map.set(c_Key2(29, 12), 2.126520559771617e-35);
        lm_coeff_map.set(c_Key2(29, 13), 2.978320111724953e-38);
        lm_coeff_map.set(c_Key2(29, 14), 4.328953650763013e-41);
        lm_coeff_map.set(c_Key2(29, 15), 6.559020682974263e-44);
        lm_coeff_map.set(c_Key2(29, 16), 1.041114394122899e-46);
        lm_coeff_map.set(c_Key2(29, 17), 1.740993970105182e-49);
        lm_coeff_map.set(c_Key2(29, 18), 3.086868741321245e-52);
        lm_coeff_map.set(c_Key2(29, 19), 5.846342313108418e-55);
        lm_coeff_map.set(c_Key2(29, 20), 1.193131084307840e-57);
        lm_coeff_map.set(c_Key2(29, 21), 2.651402409572979e-60);
        lm_coeff_map.set(c_Key2(29, 22), 6.498535317580830e-63);
        lm_coeff_map.set(c_Key2(29, 23), 1.785311900434294e-65);
        lm_coeff_map.set(c_Key2(29, 24), 5.614188366145579e-68);
        lm_coeff_map.set(c_Key2(29, 25), 2.079329024498362e-70);
        lm_coeff_map.set(c_Key2(29, 26), 9.451495565901647e-73);
        lm_coeff_map.set(c_Key2(29, 27), 5.625890217798600e-75);
        lm_coeff_map.set(c_Key2(29, 28), 4.934991419121579e-77);
        lm_coeff_map.set(c_Key2(29, 29), 8.508605895037205e-79);
    }
    
    return lm_coeff_map;
};
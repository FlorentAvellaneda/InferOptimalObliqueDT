

#ifndef MALIB__RAND__H
#define MALIB__RAND__H 

#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>

namespace MaLib {

    class MonRand {
        static std::mt19937& getGen() {
            static std::mt19937 gen = std::mt19937();
            return gen;
        }
    public:
        template< class Sseq >
        static void seed(const Sseq& seq) {
            MonRand::getGen().seed(seq);
        }

        static bool getBool(double proba=0.5) {
            assert(MonRand::getGen().min() == 0);
            auto v = MonRand::getGen().operator()();
            if( v > MonRand::getGen().max() / 2 ) {
                v--;
            }

            return v < MonRand::getGen().max() * proba;
        }

        static unsigned int get() {
            return MonRand::getGen().operator()();
        }

        template<class T>
        static const T& get(const std::set<T> &data) {
            auto it = data.begin();
            std::advance(it, get(0, data.size()-1));
            return *it;
        }

        template<class T>
        static const T& get(const std::vector<T> &data) {
            return data[ get(0, data.size()-1) ];
        }

        template<class T1, class T2>
        static T1 getIndex(const std::map<T1, T2> &data) {
            auto it = data.begin();
            std::advance(it, get(0, data.size()-1));
            return it->first;
        }

        template<class T>
        static unsigned int getIndex(const std::vector<T> &data) {
            return get(0, data.size()-1);
        }

        template<class T1, class T2>
        static const T2& get(const std::map<T1, T2> &data) {
            auto it = data.begin();
            std::advance(it, get(0, data.size()-1));
            return it->second;
        }

        static unsigned int get(unsigned int min, unsigned int max) {
            std::uniform_int_distribution<> distrib(min, max);
            return distrib(MonRand::getGen());
        }

        static double getDouble(double min, double max) {
            std::uniform_real_distribution<> distrib(min, max);
            return distrib(MonRand::getGen());
        }

        static double getDouble(double max) {
            std::uniform_real_distribution<> distrib(0.0, max);
            return distrib(MonRand::getGen());
        }

        template <class T>
        static void shuffle(T &v) {
            std::shuffle(v.begin(), v.end(), getGen());
        }
        template <class RandomIt>
        static void shuffle(RandomIt first, RandomIt second) {
            std::shuffle(first, second, getGen());
        }

        template <class T, class T2>
        static std::vector<T> sampling_with_replacement(const std::vector<T> &v, const std::vector<T2> &weights, unsigned int nb) {
            std::vector<T> res;
            std::discrete_distribution<> d(weights.begin(), weights.end());
            for(unsigned int i=0; i<nb; i++) {
                res.push_back(v[d(getGen())]);
            }
            return res;
        }
    };

}


#endif

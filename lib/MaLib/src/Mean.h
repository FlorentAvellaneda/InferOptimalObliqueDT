#ifndef MALIB_MEAN_H
#define MALIB_MEAN_H


#include <iostream>
#include <vector>
#include <cmath>

/*
  Example:
    Mean M;
    for(auto obj: A) {
        M.add( f(obj) );
    }
    M.print();
*/

namespace MaLib {
    // Welford's Online Algorithm
    class Mean
    {
        unsigned int nb = 0;
        double mean = 0.0;
        double M2 = 0.0;
        double val_t = 3.0; // Pour le % de certitude (1 => 68; 2 => 95% ; 3 => 99,7)

        std::string title;

    public:

        Mean(std::string title="") : title(title) {

        }

        void clear() {
            nb = 0;
            mean = 0.0;
            M2 = 0.0;
        }

        // "t" implique le % de certitude (1 => 68%; 2 => 95% ; 3 => 99.7%)
        void setT(double t) {
                val_t = t;
        }

        //
        void add(double val) {
            ++nb;
            double delta = val - mean;
            mean += delta/(double)nb;
            M2 += delta*(val - mean);
        }

        double getMeanMin() const {
            return getMean() - val_t*(getEcartType() / sqrt(nb));
        }

        double getMeanMax()  const  {
            return getMean() + val_t*(getEcartType() / sqrt(nb));
        }

        double getMean() const  {
            return mean;
        }

        double getVariance() const {
            return M2/((double)nb-1.0);
        }

        double getEcartType() const  {
            return sqrt(getVariance());
        }

        double size() {
            return nb;
        }

        void print(std::string unite="") {
            std::cout << "Average " << title << ": " << getMean() << unite << " +- " << (getMeanMax()-getMean()) << unite << std::endl;
        }

        ~Mean() {

        }
    };

    static std::ostream& operator<<(std::ostream& os, const Mean& m) {
        os << m.getMean() << " ± " << (m.getMeanMax()-m.getMean());
        return os;
    }
}




#endif



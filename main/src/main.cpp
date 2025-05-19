
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/format.h>
#include <optional>
#include <functional>

#include <limits>
#include <libsvm/svm.h>

#include <z3++.h>
#include <gmpxx.h>

#include "csv.hpp"

#include "Chrono.h"
#include "rand.h"
#include "xdot.h"
#include "Mean.h"

#include "cadical.hpp"

#include "MesData.h"

#include "EvalMaxSAT.h"

using namespace std;

int verbosity = 1;
string GLOBAL_externalSolver="";
bool GLOBAL_exact = false;
bool GLOBAL_no_show = false;


const double GLOBAL_MAX_ROUND_ERROR = 1e-11;


double orientation(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y) {
    return (p2y - p1y) * (p3x - p2x) - (p2x - p1x) * (p3y - p2y);
    double val = (p2y - p1y) * (p3x - p2x) - (p2x - p1x) * (p3y - p2y);
    if (val == 0) return 0;  // colinear
    return (val > 0) ? 1 : -1; // clock or counterclock wise
}


std::vector<double> findSeparatorbyEnumeration(const std::vector<std::vector<double>> &data1, const std::vector<std::vector<double>> &data2, double r = 1) {
    
    for(unsigned int i=0; i<data1.size(); i++) {
        mpq_class x1 = data1[i][0];
        mpq_class y1 = data1[i][1];
        for(unsigned int j=0; j<data2.size(); j++) {
            mpq_class x2 = data2[j][0];
            mpq_class y2 = data2[j][1];

            mpq_class x3 = (x1+x2)/2.0;
            mpq_class y3 = (y1+y2)/2.0;

            mpq_class a = (y1-y2);
            mpq_class b = (x2-x1);
            mpq_class c = x1*y2 - x2*y1;
            bool inverser = false;

            mpq_class A = a*cos(1.5) - b*sin(1.5);
            mpq_class B = b*cos(1.5) + a*sin(1.5);
            mpq_class C = -x3*( a * cos(1.5) - b * sin(1.5) ) - y3*( a * sin(1.5) + b * cos(1.5) ); 

            if( A*x1 + B*y1 + C < 0 ) {
                A = a*cos(-1.5) - b*sin(-1.5);
                B = b*cos(-1.5) + a*sin(-1.5);
                C = -x3*( a * cos(-1.5) - b * sin(-1.5) ) - y3*( a * sin(-1.5) + b * cos(-1.5) ); 
                inverser = true;
            }
            
            std::vector<double> result = {a.get_d(), b.get_d(), c.get_d()};

            bool ok = true;
            for(unsigned int k=0; k<data1.size(); k++) {
                mpq_class x = data1[k][0];
                mpq_class y = data1[k][1];
                if( x * result[0] + y * result[1] + result[2] < -GLOBAL_MAX_ROUND_ERROR ) {
                    ok = false;
                    break;
                }
            }
            if(ok == false) {
                continue;
            }
            for(unsigned int k=0; k<data2.size(); k++) {
                mpq_class x = data2[k][0];
                mpq_class y = data2[k][1];
                if( x * result[0] + y * result[1] + result[2] > GLOBAL_MAX_ROUND_ERROR ) {
                    ok = false;
                    break;
                }
            }
            if(ok == false) {
                continue;
            }
            if(ok == true) {
                unsigned int bestData1 = i;
                for(unsigned int k=0; k<data1.size(); k++) {
                    mpq_class x = data1[k][0];
                    mpq_class y = data1[k][1];
                    if( abs(a*x + b*y + c) <= GLOBAL_MAX_ROUND_ERROR ) {
                        if( (x-x2)*(x-x2) + (y-y2)*(y-y2) < (data1[bestData1][0]-x2)*(data1[bestData1][0]-x2) + (data1[bestData1][1]-y2)*(data1[bestData1][1]-y2) ) {
                            bestData1 = k;
                        }
                    }
                }
                unsigned int bestData2 = j;
                for(unsigned int k=0; k<data2.size(); k++) {
                    mpq_class x = data2[k][0];
                    mpq_class y = data2[k][1];
                    if( abs(a*x + b*y + c) <= GLOBAL_MAX_ROUND_ERROR ) {
                        if( (x-x1)*(x-x1) + (y-y1)*(y-y1) < (data2[bestData2][0]-x1)*(data2[bestData2][0]-x1) + (data2[bestData2][1]-y1)*(data2[bestData2][1]-y1) ) {
                            bestData2 = k;
                        }
                    }
                }

                x1 = data1[bestData1][0];
                y1 = data1[bestData1][1];
                x2 = data2[bestData2][0];
                y2 = data2[bestData2][1];

                x3 = (x1+x2)/2.0;
                y3 = (y1+y2)/2.0;

                a = (y1-y2);
                b = (x2-x1);
                c = x1*y2 - x2*y1;

                mpq_class d = 0.0000001;
                do {
                    d = d*2;
                    if(inverser) {
                        A = a*cos(-d.get_d()) - b*sin(-d.get_d());
                        B = b*cos(-d.get_d()) + a*sin(-d.get_d());
                        C = -x3*( a * cos(-d.get_d()) - b * sin(-d.get_d()) ) - y3*( a * sin(-d.get_d()) + b * cos(-d.get_d()) ); 
                    } else {
                        A = a*cos(d.get_d()) - b*sin(d.get_d());
                        B = b*cos(d.get_d()) + a*sin(d.get_d());
                        C = -x3*( a * cos(d.get_d()) - b * sin(d.get_d()) ) - y3*( a * sin(d.get_d()) + b * cos(d.get_d()) ); 
                    }
                } while( abs(A*x1 + B*y1 + C) < 10*GLOBAL_MAX_ROUND_ERROR );

                result[0] = A.get_d();
                result[1] = B.get_d();
                result[2] = C.get_d();
                
                return result;
            }
        }
    }

    std::cout << "BUG: no linearly separable solution found, but it should be (probably due to rounding errors)." << std::endl;

    return std::vector<double>();

}

std::vector<double> findSeparator(const std::vector<std::vector<double>> &data1, const std::vector<std::vector<double>> &data2) {
    int num_features;
    if(data1.size() > data2.size())
        num_features = data1[0].size();
    else
        num_features = data2[0].size();

    int num_data = data1.size() + data2.size();

    // Prepar data for LibSVM
    svm_problem prob;
    prob.l = num_data;
    std::vector<double> labels;
    for(int i = 0; i < data1.size(); i++)
        labels.push_back(1.0);
    for(int i = 0; i < data2.size(); i++)
        labels.push_back(-1.0);
    prob.y = labels.data();

    svm_node* x = new svm_node[num_data * (num_features + 1)];
    prob.x = new svm_node*[num_data];

    int node_index = 0;
    for (int i = 0; i < data1.size(); ++i) {
        prob.x[i] = &x[node_index];
        for (int j = 0; j < num_features; ++j) {
            x[node_index].index = j + 1;
            x[node_index].value = data1[i][j];
            node_index++;
        }
        x[node_index].index = -1;
        node_index++;
    }
    for (int i = 0; i < data2.size(); ++i) {
        prob.x[i + data1.size()] = &x[node_index];
        for (int j = 0; j < num_features; ++j) {
            x[node_index].index = j + 1;
            x[node_index].value = data2[i][j];
            node_index++;
        }
        x[node_index].index = -1;
        node_index++;
    }

    // Parameter for hard margin
    svm_parameter param;
    param.svm_type = C_SVC;
    param.kernel_type = LINEAR;
    param.C = std::numeric_limits<double>::infinity();
    param.eps = 1e-3;
    param.gamma = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    param.nu = 0.001;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;

    svm_model* model = svm_train(&prob, &param);
    std::vector<double> result(num_features, 0);

    // Verif
    for(auto &p: data1) {
        svm_node test_point[3];
        test_point[0].index = 1;
        test_point[0].value = p[0];
        test_point[1].index = 2;
        test_point[1].value = p[1];
        test_point[2].index = -1;

        double prediction = svm_predict(model, test_point);
        if(prediction != 1) {
            return findSeparatorbyEnumeration(data1, data2);
        }
    }
    for(auto &p: data2) {
        svm_node test_point[3];
        test_point[0].index = 1;
        test_point[0].value = p[0];
        test_point[1].index = 2;
        test_point[1].value = p[1];
        test_point[2].index = -1;

        double prediction = svm_predict(model, test_point);
        if(prediction != -1) {
            return findSeparatorbyEnumeration(data1, data2);
        }
    }

    assert(model->param.kernel_type == LINEAR);

    const int num_sv = model->l;

    for (int i = 0; i < num_sv; ++i) {
        double alpha = model->sv_coef[0][i];
        for (int j = 0; j < num_features; ++j) {
            result[j] += alpha * model->SV[i][j].value;
        }
    }
    result.push_back(-model->rho[0]);

    svm_free_and_destroy_model(&model);
    delete[] x;
    delete[] prob.x;

    return result;
}

double vectorialProduct(double ax, double ay, double bx, double by) {
    if((ax == std::numeric_limits<double>::infinity()) && (ay == std::numeric_limits<double>::infinity())) {
        ax = 1;
        ay = 1;
    }
    if((ax == -std::numeric_limits<double>::infinity()) && (ay == std::numeric_limits<double>::infinity())) {
        ax = -1;
        ay = 1;
    }
    if((ax == std::numeric_limits<double>::infinity()) && (ay == -std::numeric_limits<double>::infinity())) {
        ax = 1;
        ay = -1;
    }
    if((ax == -std::numeric_limits<double>::infinity()) && (ay == -std::numeric_limits<double>::infinity())) {
        ax = -1;
        ay = -1;
    }
    if((bx == std::numeric_limits<double>::infinity()) && (by == std::numeric_limits<double>::infinity())) {
        bx = 1;
        by = 1;
    }
    if((bx == -std::numeric_limits<double>::infinity()) && (by == std::numeric_limits<double>::infinity())) {
        bx = -1;
        by = 1;
    }
    if((bx == std::numeric_limits<double>::infinity()) && (by == -std::numeric_limits<double>::infinity())) {
        bx = 1;
        by = -1;
    }
    if((bx == -std::numeric_limits<double>::infinity()) && (by == -std::numeric_limits<double>::infinity())) {
        bx = -1;
        by = -1;
    }
    assert(ax != std::numeric_limits<double>::infinity());
    assert(ay != std::numeric_limits<double>::infinity());
    assert(bx != std::numeric_limits<double>::infinity());
    assert(by != std::numeric_limits<double>::infinity());
    assert(ax != -std::numeric_limits<double>::infinity());
    assert(ay != -std::numeric_limits<double>::infinity());
    assert(bx != -std::numeric_limits<double>::infinity());
    assert(by != -std::numeric_limits<double>::infinity());

    return ax*by - ay*bx;
}

bool isInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) {
    assert((cx != std::numeric_limits<double>::infinity()) && (cy != std::numeric_limits<double>::infinity()) && (cx != -std::numeric_limits<double>::infinity()) && (cy != -std::numeric_limits<double>::infinity()));
    assert((px != std::numeric_limits<double>::infinity()) && (py != std::numeric_limits<double>::infinity()) && (px != -std::numeric_limits<double>::infinity()) && (py != -std::numeric_limits<double>::infinity()));

    if(ax == -std::numeric_limits<double>::infinity() && ay == -std::numeric_limits<double>::infinity()) {
        double Bb = by - bx;
        double Bc = cy - cx;
        double Bp = py - px;

        bool basgauche;
        if(Bc + GLOBAL_MAX_ROUND_ERROR < Bb) {
            basgauche = orientation(bx, by, cx, cy, px, py) >= -GLOBAL_MAX_ROUND_ERROR;
        } else if(Bc - GLOBAL_MAX_ROUND_ERROR > Bb) {
            basgauche = orientation(bx, by, cx, cy, px, py) <= GLOBAL_MAX_ROUND_ERROR;
        } else {
            if( abs( Bp - Bc ) > GLOBAL_MAX_ROUND_ERROR ) {
                return false;
            }
            return (px - GLOBAL_MAX_ROUND_ERROR <= max(bx, cx)) && (py - GLOBAL_MAX_ROUND_ERROR <= max(by, cy)) && (px + GLOBAL_MAX_ROUND_ERROR >= min(bx, cx)) && (py + GLOBAL_MAX_ROUND_ERROR >= min(by, cy));
        }
        return (Bp + GLOBAL_MAX_ROUND_ERROR >= min(Bb, Bc)) && (Bp - GLOBAL_MAX_ROUND_ERROR <= max(Bb, Bc)) && basgauche;
    }

    if(ax == std::numeric_limits<double>::infinity() && ay == std::numeric_limits<double>::infinity()) {
        double Bb = by - bx;
        double Bc = cy - cx;
        double Bp = py - px;

        bool hautdroite;
        if(Bc + GLOBAL_MAX_ROUND_ERROR < Bb) {
            hautdroite = orientation(bx, by, cx, cy, px, py) <= GLOBAL_MAX_ROUND_ERROR;
        } else if(Bc - GLOBAL_MAX_ROUND_ERROR > Bb) {
            hautdroite = orientation(bx, by, cx, cy, px, py) >= -GLOBAL_MAX_ROUND_ERROR;
        } else {
            if( abs( Bp - Bc ) > GLOBAL_MAX_ROUND_ERROR ) {
                return false;
            }
            return (px + GLOBAL_MAX_ROUND_ERROR >= min(bx, cx)) && (py + GLOBAL_MAX_ROUND_ERROR >= min(by, cy)) && (px - GLOBAL_MAX_ROUND_ERROR <= max(bx, cx)) && (py - GLOBAL_MAX_ROUND_ERROR <= max(by, cy));
        }
        return (Bp + GLOBAL_MAX_ROUND_ERROR >= min(Bb, Bc)) && (Bp - GLOBAL_MAX_ROUND_ERROR <= max(Bb, Bc)) && hautdroite;
    }

    double v1 = vectorialProduct(bx-ax, by-ay, px-ax, py-ay);
    double v2 = vectorialProduct(cx-bx, cy-by, px-bx, py-by);
    double v3 = vectorialProduct(ax-cx, ay-cy, px-cx, py-cy);

    if(abs(v1) <= GLOBAL_MAX_ROUND_ERROR && abs(v2) <= GLOBAL_MAX_ROUND_ERROR && abs(v3) <= GLOBAL_MAX_ROUND_ERROR) {
        return (px - GLOBAL_MAX_ROUND_ERROR >= min(ax, min(bx, cx))) && (px + GLOBAL_MAX_ROUND_ERROR <= max(ax, max(bx, cx)));
    }

    return ((v1 - GLOBAL_MAX_ROUND_ERROR >= 0) && (v2 - GLOBAL_MAX_ROUND_ERROR >= 0) && (v3 - GLOBAL_MAX_ROUND_ERROR >= 0)) || ((v1 + GLOBAL_MAX_ROUND_ERROR <= 0) && (v2 + GLOBAL_MAX_ROUND_ERROR <= 0) && (v3 + GLOBAL_MAX_ROUND_ERROR <= 0));
}

struct VirtualDecisionNode {
    virtual ~VirtualDecisionNode() = default;
    virtual std::string toDot( const std::function< std::string(unsigned int) > & getNameFeature, const std::function< std::string(unsigned int, unsigned int) > & getValue ) const = 0;
    virtual unsigned int get(const std::vector< double > & feature_values) const = 0;
};

class LeafDecisionNode : public VirtualDecisionNode {
    unsigned int feature;
    unsigned int value;
    int nb_in_leaf;
    double acc_in_leaf;
public:
    LeafDecisionNode(unsigned int feature, unsigned int value, int nb_in_leaf=-1, double acc_in_leaf=-1) : feature(feature), value(value), nb_in_leaf(nb_in_leaf), acc_in_leaf(acc_in_leaf) {

    }

    virtual std::string toDot( const std::function< std::string(unsigned int) > & getNameFeature, const std::function< std::string(unsigned int, unsigned int) > & getValue ) const override {
        std::ostringstream oss;
        if(nb_in_leaf != -1) {
            oss << "Q" << this << " [shape=\"box\", label=\""<< getNameFeature(feature) << " = " << getValue(feature, value) << " (# = "<<nb_in_leaf<<", acc="<<acc_in_leaf<<")\"];" << endl;
        } else {
            oss << "Q" << this << " [shape=\"box\", label=\""<< getNameFeature(feature) << " = " << getValue(feature, value) << "\"];" << endl;
        }
        return oss.str();
    }

    virtual unsigned int get(const std::vector< double > & feature_values) const override {
        return value;
    }
};

class NumericDecisionNodeV2 : public VirtualDecisionNode {
    unsigned int feature;
    double realValue;
    std::unique_ptr<VirtualDecisionNode> left; // If the feature is less than the value
    std::unique_ptr<VirtualDecisionNode> right; // If the feature is greater or equal to the value
public:
    NumericDecisionNodeV2(unsigned int feature, double value, std::unique_ptr<VirtualDecisionNode> && gauche, std::unique_ptr<VirtualDecisionNode> && droite) : feature(feature), realValue(value), left(std::move(gauche)), right(std::move(droite)) {

    }

    virtual std::string toDot( const std::function< std::string(unsigned int) > & getNameFeature, const std::function< std::string(unsigned int, unsigned int) > & getValue ) const override {
        std::ostringstream oss;
        
        oss << "Q" << this << " [shape=\"box\", label=\"" << getNameFeature(feature) << " ≥ " << realValue << "\"];" << endl;
        if(right != nullptr) {
            oss << left->toDot(getNameFeature, getValue);
            oss << right->toDot(getNameFeature, getValue);
            oss << "Q" << this << " -> Q" << left << " [fontcolor=red, label=\"no\", color=red];" << std::endl;
            oss << "Q" << this << " -> Q" << right << " [fontcolor=green, label=\"yes\", color=green];" << std::endl;
        }

        return oss.str();
    }

    virtual unsigned int get(const std::vector< double > & feature_values) const override {
        if(!left.get()) {
            assert(false);
            return realValue;
        }

        if(feature_values[feature] >= realValue) {
            return right->get(feature_values);
        }
        return left->get(feature_values);
    }

    static std::unique_ptr<VirtualDecisionNode> random(const std::vector<Domain> & features, double proba_feuille, unsigned int indexOut, unsigned int max_depth) {
        if( (max_depth == 0) || MaLib::MonRand::getBool(proba_feuille) ) {
            unsigned int value = MaLib::MonRand::get(0, features[indexOut].size()-1);
            return std::make_unique<LeafDecisionNode>(indexOut , MaLib::MonRand::get(0, features[indexOut].size()-1));
        } else {
            unsigned int feature = MaLib::MonRand::getIndex(features);
            double value = 0;

            while( (feature == indexOut) || (features[feature].size() == 0)) {
                feature = MaLib::MonRand::getIndex(features);
            }
            assert( feature != indexOut );

            if(features[feature].size() > 1) {
                if( features[feature].is_numeric() ) {
                    value = features[feature].toDouble( MaLib::MonRand::get(1, features[feature].size()-1) );
                } else {
                    value = features[feature].toDouble( MaLib::MonRand::get(0, features[feature].size()-1) );
                }
            } else {
                value = 0;
            }
            
            assert( features[feature].is_numeric() == true );

            return std::make_unique<NumericDecisionNodeV2>(
                feature,
                value,
                NumericDecisionNodeV2::random( features, proba_feuille, indexOut, max_depth-1 ),
                NumericDecisionNodeV2::random( features, proba_feuille, indexOut, max_depth-1 )
            );
        }

        assert(false);
        return nullptr;
    }
};

class LinearDecisionNode : public VirtualDecisionNode {
    unsigned int feature1;
    double coef1;
    unsigned int feature2;
    double coef2;

    double threshold;

    std::unique_ptr<VirtualDecisionNode> left; // If the feature is less than the value
    std::unique_ptr<VirtualDecisionNode> right; // If the feature is greater or equal to the value
public:

    LinearDecisionNode(unsigned int feature1, double coef1, unsigned int feature2, double coef2, double seuil, std::unique_ptr<VirtualDecisionNode> && gauche, std::unique_ptr<VirtualDecisionNode> && droite) : feature1(feature1), coef1(coef1), feature2(feature2), coef2(coef2), threshold(seuil), left(std::move(gauche)), right(std::move(droite)) {

    }

    static std::unique_ptr<VirtualDecisionNode> random(const std::vector<Domain> & features, unsigned int indexOut, unsigned int depth) {
        if(depth == 0) {
            return std::make_unique<LeafDecisionNode>(indexOut , MaLib::MonRand::get(0, features[indexOut].size()-1));
        }

        unsigned int f1 = MaLib::MonRand::getIndex(features);
        while( (f1 == indexOut) || (features[f1].size() == 0)) {
            f1 = MaLib::MonRand::getIndex(features);
        }
        unsigned int f2 = MaLib::MonRand::getIndex(features);
        while( (f2 == indexOut) || (features[f2].size() == 0)) {
            f2 = MaLib::MonRand::getIndex(features);
        }

        return std::make_unique<LinearDecisionNode>(
            f1,
            MaLib::MonRand::getDouble(-10, 10),
            f2,
            MaLib::MonRand::getDouble(-10, 10),
            MaLib::MonRand::getDouble(-10, 10),
            LinearDecisionNode::random( features, indexOut, depth-1 ),
            LinearDecisionNode::random( features, indexOut, depth-1 )
        );
    }

    std::string toDot( const std::function< std::string(unsigned int) > & getNameFeature, const std::function< std::string(unsigned int, unsigned int) > & getValue ) const override {
        std::ostringstream oss;
        
        oss << "Q" << this << " [shape=\"box\", label=\"" << getNameFeature(feature1) << " * " << coef1 << " + " << getNameFeature(feature2) << " * " << coef2 << " ≤ " << threshold << "\"];" << endl;
        if(right != nullptr) {
            oss << left->toDot(getNameFeature, getValue);
            oss << right->toDot(getNameFeature, getValue);
            oss << "Q" << this << " -> Q" << left << " [fontcolor=red, label=\"no\", color=red];" << std::endl;
            oss << "Q" << this << " -> Q" << right << " [fontcolor=green, label=\"yes\", color=green];" << std::endl;
        }

        return oss.str();
    }

    unsigned int get(const std::vector< double > & feature_values) const override {
        if(!left.get()) {
            return MaLib::MonRand::getBool();
        }

        double v = feature_values[feature1] * coef1 + feature_values[feature2] * coef2;

        if(v <= threshold) {
            return right->get(feature_values);
        }
        return left->get(feature_values);
    }

};

class DecisionTree {
    std::unique_ptr<VirtualDecisionNode> root;
    std::vector<Domain> _vFeatures;
public:
    DecisionTree(std::unique_ptr<VirtualDecisionNode> && root, std::vector<Domain> vFeatures) : root(std::move(root)), _vFeatures(vFeatures) {

    }

    std::string toDot() {
        return fmt::format("digraph mon_graphe {{\n{}}}\n", root->toDot( 
            [this](unsigned int i) { return _vFeatures[i].getName(); },
            [this](unsigned int i, unsigned int j) { return _vFeatures[i].toString(j); }
        ));
    }

    double accuracy(const MesData &data) const {
        unsigned int acc = 0;
        for(unsigned int i=0; i<data.getData().size(); i++) {
            acc += (root->get(data.getDataDouble(i)) == data.getData()[i].back().value());
        }
        return 100.0*acc/(double)data.size();
    }

    unsigned int operator() (const std::vector< double > & feature_values) const {
        return root->get(feature_values);
    }

    static DecisionTree randomLinear(unsigned int nombre_features, unsigned int nombre_valeur_par_features, unsigned int nb_class, unsigned int depth) {

        std::vector<Domain> vFeatures;
        bool noFeature = true;
        for(unsigned int i=0; i<nombre_features; i++) {
            std::set<double> tmp;
            for(unsigned int j=0; j<nombre_valeur_par_features; j++) {
                tmp.insert(MaLib::MonRand::getDouble(100));
            }
            if(tmp.size() > 1) {
                noFeature = false;
            }
            vFeatures.push_back( Domain("Feature"+to_string(i), tmp) );
        }

        std::set<double> tmp;
        for(unsigned int c=0; c<nb_class; c++) {
            tmp.insert(c);
        }
        vFeatures.push_back( Domain("Out", tmp) );

        if(noFeature)
            depth = 0;

        return DecisionTree( LinearDecisionNode::random( vFeatures, vFeatures.size()-1, depth ), vFeatures );
    }

    static DecisionTree random(unsigned int nombre_features, unsigned int nombre_max_valeur_par_features=2, double proba_feuille=0.3, unsigned int nb_class=3, unsigned int max_depth=3) {
        assert(nb_class>=1);

        std::vector<Domain> vFeatures;
        bool noFeature = true;
        for(unsigned int i=0; i<nombre_features; i++) {
            std::set<double> tmp;
            unsigned int nombre_valeur_par_features = MaLib::MonRand::get(1, nombre_max_valeur_par_features); // TODO: changer 0 en 2
            for(unsigned int j=0; j<nombre_valeur_par_features; j++) {
                tmp.insert(MaLib::MonRand::getDouble(100));
            }
            if(tmp.size() > 1) {
                noFeature = false;
            }
            vFeatures.push_back( Domain("Feature"+to_string(i), tmp) );
        }

        std::set<double> tmp;
        for(unsigned int c=0; c<nb_class; c++) {
            tmp.insert(c);
        }
        vFeatures.push_back( Domain("Out", tmp) );

        if(noFeature)
            max_depth = 0;

        return DecisionTree( NumericDecisionNodeV2::random( vFeatures, proba_feuille, vFeatures.size()-1, max_depth ), vFeatures );
    }

    MesData generate_dataset(unsigned int size) {
        std::vector< std::vector< std::optional<unsigned int> > > data;
        data.reserve(size);

        for(unsigned int i=0; i<size; i++) {
            std::vector< std::optional<unsigned int> > exemple;
            exemple.reserve(_vFeatures.size());

            for(unsigned int j=0; j<_vFeatures.size()-1; j++) {
                if(_vFeatures[j].size()) {
                    exemple.push_back( MaLib::MonRand::get(0, _vFeatures[j].size()-1) );
                } else {
                    exemple.push_back( std::nullopt );
                }
            }

            std::vector<double> exemple_values;
            for(unsigned int i=0; i<exemple.size(); i++) {
                if(exemple[i].has_value()) {
                    exemple_values.push_back( _vFeatures[i].toDouble(exemple[i].value()) );
                } else {
                    exemple_values.push_back( 0 );
                }
            }

            exemple.push_back( root->get(exemple_values) );

            data.push_back( exemple );
        }

        return MesData(data, _vFeatures);
    }

};

struct AbstractSolveur {
    virtual void add(int lit) = 0;
    virtual void add(const std::vector<int> &clause) = 0;
    virtual void add_soft(unsigned int w) = 0;
    virtual bool solve() = 0;
    virtual bool getValue(int lit) = 0;
    virtual unsigned int getNbClause() {
        return 0;
    }
};

class SolverSMT : public AbstractSolveur {
    z3::context c;
    z3::model model;

    std::vector<z3::expr> constraints;
    std::vector<std::tuple<z3::expr, unsigned int>> constraints_soft;

    z3::expr tmp;

    public:

    SolverSMT() : c(), model(c), tmp(c.bool_val(false)) {

    }

    void add_side1(int lit1, int lit2, int idA, double x, int idB, double y, int idC) {

        z3::expr tmp = c.bool_val(false);
        if(lit1 > 0) {
            tmp = c.bool_const(std::to_string(lit1).c_str());
        } else {
            tmp = !c.bool_const(std::to_string(-lit1).c_str());
        }

        if(lit2 > 0) {
            tmp = tmp || c.bool_const(std::to_string(lit2).c_str());
        } else {
            tmp = tmp || !c.bool_const(std::to_string(-lit2).c_str());
        }
        
        constraints.push_back(  tmp || (
            (c.real_val(std::to_string(x).c_str()) * c.real_const(std::to_string(idA).c_str())) +
            (c.real_val(std::to_string(y).c_str()) * c.real_const(std::to_string(idB).c_str())) +
            (-c.real_const(std::to_string(idC).c_str())) < 0 // <= c.real_val("-0.001")
        ));
    }


    void add_side2(int lit1, int lit2, int idA, double x, int idB, double y, int idC) {

        z3::expr tmp = c.bool_val(false);
        if(lit1 > 0) {
            tmp = c.bool_const(std::to_string(lit1).c_str());
        } else {
            tmp = !c.bool_const(std::to_string(-lit1).c_str());
        }

        if(lit2 > 0) {
            tmp = tmp || c.bool_const(std::to_string(lit2).c_str());
        } else {
            tmp = tmp || !c.bool_const(std::to_string(-lit2).c_str());
        }
        
        constraints.push_back(  tmp || (
            (c.real_val(std::to_string(x).c_str()) * c.real_const(std::to_string(idA).c_str())) +
            (c.real_val(std::to_string(y).c_str()) * c.real_const(std::to_string(idB).c_str())) +
            (-c.real_const(std::to_string(idC).c_str())) > 0 // c.real_val("0.001")
        ));
    }

    void add_side1(int idA, double x) {
        tmp = tmp || ( 
            c.real_val(std::to_string(x).c_str()) < c.real_const(std::to_string(idA).c_str())
        );
    }

    void add_side1(int idA, double x, int idB, double y, int idC) {
        tmp = tmp || ( 
            (c.real_val(std::to_string(x).c_str()) * c.real_const(std::to_string(idA).c_str())) +
            (c.real_val(std::to_string(y).c_str()) * c.real_const(std::to_string(idB).c_str())) +
            (-c.real_const(std::to_string(idC).c_str())) < 0 // <= c.real_val("-0.001")
            );
    }

    void add_side2(int idA, double x) {
        tmp = tmp || ( 
            c.real_val(std::to_string(x).c_str()) >= c.real_const(std::to_string(idA).c_str())
        );
    }

    void add_side2(int idA, double x, int idB, double y, int idC) {
        tmp = tmp || ( 
            (c.real_val(std::to_string(x).c_str()) * c.real_const(std::to_string(idA).c_str())) +
            (c.real_val(std::to_string(y).c_str()) * c.real_const(std::to_string(idB).c_str())) +
            (-c.real_const(std::to_string(idC).c_str())) > 0 // c.real_val("0.001")
            );
    }

    void add_soft(unsigned int w) {
        constraints_soft.push_back({tmp, w});
        tmp = c.bool_val(false);
    }

    void add(int lit) {
        if(lit == 0) {
            constraints.push_back(tmp);
            tmp = c.bool_val(false);
        } else {
            if(lit > 0) {
                tmp = tmp || c.bool_const(std::to_string(lit).c_str());
            } else {
                tmp = tmp || !c.bool_const(std::to_string(-lit).c_str());
            }
        }
    }

    void add(const std::vector<int> &clause) {
        for(auto &l: clause) {
            add(l);
        }
        add(0);
    }

    bool getValue(int lit) {
        if(lit > 0) {
            return model.eval(c.bool_const(std::to_string(lit).c_str())).is_true();
        }
        return !model.eval(c.bool_const(std::to_string(-lit).c_str())).is_true();
    }

    double getRealValue(int lit) {
        return model.eval(c.real_const(std::to_string(lit).c_str())).as_double();
    }

    bool solve() {
        z3::optimize s(c);
        for(auto &e : constraints) {
            s.add(e);
        }
        for(auto &[e, w] : constraints_soft) {
            s.add_soft(e, w);
        }
        auto result = s.check();
        if(result == z3::sat) {
            model = s.get_model();
            return true;
        }
        return false;
    }
};

class ExternalSolver : public AbstractSolveur {
    std::string cmd;

    std::vector< std::tuple< unsigned int, std::vector<int> > > clauses;
    std::vector<int> tmpClause;

public:
    ExternalSolver(std::string cmd) : cmd(cmd) {
        this->cmd.append(" /tmp/formula3.wcnf");
    }

    virtual void add(int lit) {
        if(lit==0) {
            clauses.push_back({0, tmpClause});
            tmpClause.clear();
        } else {
            tmpClause.push_back(lit);
        }
    }

    virtual void add(const std::vector<int> &clause) {
        clauses.push_back({0, clause});
    }
    
    virtual void add_soft(unsigned int w) {
        clauses.push_back({w, tmpClause});
        tmpClause.clear();
    }


    void exportFormula(std::string vFile) {
        std::ofstream file( vFile );
        for(auto &line: clauses) {
            if(std::get<0>(line)==0) {
                file << "h ";
            } else {
                file << std::get<0>(line) << " ";
            }

            for(auto &l: std::get<1>(line)) {
                file << l << " ";
            }
            file << "0" << std::endl;
        }
    }

    t_weight cost=0;
    std::vector<bool> values;

    virtual bool solve() {

        unsigned int costSolver = 0;
        values.clear();

        exportFormula("/tmp/formula3.wcnf");

        std::string result = exec(cmd);
        std::istringstream issResult(result);

        std::string line;

        std::string res;
        while (std::getline(issResult, line))
        {
            std::istringstream iss(line);
            char action;
            if(!(iss >> action)) {
                continue;
            }
            switch (action) {

            case 'c':
            {
                break;
            }

            case 'o':
            {
                if(!(iss >> costSolver)) {
                    assert(false);
                }
                break;
            }

            case 's':
            {
                if(!(iss >> res)) {
                    assert(false);
                }
                if( (res.compare("OPTIMUM") != 0) && (res.compare("SATISFIABLE") != 0) ) {
                    return false;
                }
                break;
            }

            case 'v':
            {
                if(!(iss >> res)) {
                    assert(false);
                }
                values.push_back(0); // fake lit

                int res2;
                if(!(iss >> res2)) {
                    if(res.compare("-1") ==  0) {
                        res = "0"; // special case of a single variable
                    }

                    // New format
                    for(unsigned int i=0; i<res.size(); i++) {
                        values.push_back(res[i] == '1');
                    }
                } else {
                    // Old format
                    int lit = std::atoi(res.c_str());
                    values.push_back(lit > 0);
                    assert((values.size()-1) == abs(lit));

                    values.push_back(res2 > 0);

                    while(iss>>lit) {
                        values.push_back(lit > 0);
                        assert((values.size()-1) == abs(lit));
                    }
                }

                break;
            }

            default:
                assert(false);
            }

        }

        assert(costSolver != std::numeric_limits<unsigned int>::max());
        assert(values.size());

        cost += costSolver;

        return true;
    }

    bool getValue(int lit) {
        assert(values.size() > lit);
        return values[lit];
    }

    t_weight getCost() {
        return cost;
    }

    std::string exec(std::string cmd) {
        std::array<char, 128> buffer;
        std::string result;
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
        if (!pipe) {
            throw std::runtime_error("popen() failed!");
        }
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            result += buffer.data();
        }
        return result;
    }
};

class SolveurMaxSAT : public AbstractSolveur {
    EvalMaxSAT<Solver_cadical> solver;
    std::vector<int> clause;
    bool no_more_var = false;

    unsigned int nbClause = 0;

    public:

    SolveurMaxSAT() {
        solver.setTargetComputationTime(1);
        solver.unactivateMultiSolveStrategy();
        solver.unactivateUBStrategy();
    }

    void add(int lit) {
        while(solver.nVars() < abs(lit)) {
            assert(no_more_var == false);
            if(no_more_var) {
                fmt::print("Error: no more var\n");
                exit(-1);
            }
            solver.newVar();
        }
        if(lit==0) {
            solver.addClause(clause);
            nbClause++;
            clause.clear();
        } else {
            clause.push_back(lit);
        }
    }

    void add(const std::vector<int> &clause) {
        solver.addClause(clause);
        nbClause++;
    }

    bool solve() {
        return solver.solve();
    }

    bool getValue(int lit) {
        return solver.getValue(lit);
    }

    void add_soft(unsigned int w) {
        no_more_var=true;
        solver.addClause(clause, w);
        nbClause++;
        clause.clear();
    }

    unsigned int getNbClause() {
        return nbClause;
    }
};

// Algo Pouya
class InferDT_Pouya {
    std::vector< Domain > _vDomain;
    unsigned int _index_out;
    unsigned int _depth;

    AbstractSolveur *solver;

    struct Var {
        int id;
        static int nbVar;
        Var() : id(++nbVar) { }
    };

    std::vector< std::map<unsigned int, Var> > A; // A[n][j] = true if the feature j is chosen for the split at branching node n.
    std::vector< std::map<unsigned int, Var> > S; // S[i][n] = true if the exemple i is directed towards the left child, if it passes through branching node n.
    std::vector< std::map<unsigned int, Var> > Z; // Z[i][l] = true is the exemple i ends up at leaf node t.
    std::vector< std::map<unsigned int, Var> > C; // C[l][c] = true if the label c is assigned to leaf node l.

    std::vector< std::vector< std::optional<unsigned int> > > examples;
public:

    InferDT_Pouya(const std::vector< Domain > & _vDomain, unsigned int index_out, unsigned int depth, std::string externalSolver="") : _vDomain(_vDomain), _index_out(index_out), _depth(depth) {
        if(externalSolver.size() == 0) {
            solver = new SolveurMaxSAT();
        } else {
            solver = new ExternalSolver(externalSolver);
        }

        A.resize(std::pow(2, depth));
        C.resize(std::pow(2, depth));

        // exactly one feature is chosen at each branching node n
        for(unsigned int n = 1; n < A.size(); n++) {
            for(unsigned int j1 = 0; j1 < _vDomain.size(); j1++) {
                if(j1 == index_out) {
                    continue;
                }
                for(unsigned int j2 = j1+1; j2 < _vDomain.size(); j2++) {
                    if(j2 == index_out) {
                        continue;
                    }
                    solver->add(-A[n][j1].id);
                    solver->add(-A[n][j2].id);
                    solver->add(0); // (14)
                }
            }
        }

        for(unsigned int n = 1; n < A.size(); n++) {
            for(auto &[key, val]: A[n]) {
                if(key == index_out) {
                    assert(false);
                    continue;
                }
                solver->add(val.id);
            }
            solver->add(0); // (15)
        }

    }

    void use_all_examples( const std::vector< std::vector< std::optional<unsigned int> > > & examples ) {
        assert(this->examples.size() == 0);

        this->examples = examples;

        S.resize( examples.size() );
        Z.resize( examples.size() );

        auto get_order = [&](unsigned int feature_index){
            std::vector<int> order;
            order.reserve(examples.size());

            for(unsigned int i=0; i<examples.size(); i++) {
                if(examples[i][feature_index].has_value())
                    order.push_back(i);
            }

            std::sort(order.begin(), order.end(), [&](int a, int b){
                return examples[a][feature_index].value() < examples[b][feature_index].value();
            });

            return order;
        };

        for(unsigned int n = 1; n < A.size(); n++) {
            for(auto &[j, val]: A[n]) {
                assert(j != _index_out);

                auto order = get_order(j);
                for(unsigned a = 1; a<order.size(); a++) {
                    unsigned int i1 = order[a-1];
                    assert(examples[i1][j].has_value());
                    unsigned int i2 = order[a];
                    assert(examples[i2][j].has_value());

                    if(_vDomain[j].is_numeric()) {
                        solver->add(-A[n][j].id);
                        solver->add(S[i1][n].id);
                        solver->add(-S[i2][n].id);
                        solver->add(0); // (16)
                    }

                    if( examples[i1][j].value() == examples[i2][j].value() ) {
                        solver->add(-A[n][j].id);
                        solver->add(-S[i1][n].id);
                        solver->add(S[i2][n].id);
                        solver->add(0); // (17)

                        if(_vDomain[j].is_categorical()) {
                            solver->add(-A[n][j].id);
                            solver->add(S[i1][n].id);
                            solver->add(-S[i2][n].id);
                            solver->add(0); // (18)
                        }
                    }
                }

                if(order.size()) {
                    solver->add(-A[n][j].id);
                    solver->add(S[order.front()][n].id);
                    solver->add(0); // (23)
                
                    if(_vDomain[j].is_numeric() && (examples[order.back()][j].value() != examples[order.front()][j].value())) {
                        solver->add(-A[n][j].id);
                        solver->add(-S[order.back()][n].id);
                        solver->add(0); // (24)
                    }
                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                unsigned int n = C.size() + l;
                std::vector<int> clause = { Z[i][l].id };

                while(n > 1) {
                    bool right = n%2;
                    n = n / 2;
                    if(!right) {
                        solver->add(-Z[i][l].id );
                        solver->add( S[i][n].id );
                        solver->add(0); // (19)

                        clause.push_back( -S[i][n].id );
                    } else {
                        solver->add(-Z[i][l].id );
                        solver->add( -S[i][n].id );
                        solver->add(0); // (20)
                        
                        clause.push_back( S[i][n].id );
                    }
                }
                solver->add(clause); // (21)
            }

        }

        for(unsigned int l=0; l<C.size(); l++) {
            for(unsigned int c1=0; c1<_vDomain[_index_out].size(); c1++) {
                for(unsigned int c2=c1+1; c2<_vDomain[_index_out].size(); c2++) {
                    solver->add(-C[l][c1].id);
                    solver->add(-C[l][c2].id);
                    solver->add(0); // (22)
                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                if( examples[i][ _index_out ].has_value() ) {
                    solver->add( -Z[i][l].id );
                    solver->add( C[l][ examples[i][ _index_out ].value() ].id );
                    if(GLOBAL_exact) {
                        solver->add(0);
                    } else {
                        solver->add_soft(1);
                    }
                }
            }
        }
    }

    bool infer() {
        return solver->solve();
    }

    DecisionTree getModel() {
        std::vector< unsigned int > IDs( examples.size() );
        for(unsigned int i=0; i<examples.size(); i++) {
            IDs[i] = i;
        }
        return DecisionTree( getModel(IDs, 1), _vDomain );
    }
private:

    std::unique_ptr<VirtualDecisionNode> getModel(const std::vector< unsigned int > &IDs, unsigned int n) {
        if(n >= A.size()) {
            for(auto &[key, value]: C[n - A.size()]) {
                if( solver->getValue(value.id)) {
                    return std::make_unique<LeafDecisionNode>( _index_out, key );
                }
            }
            assert(false);
            return std::make_unique<LeafDecisionNode>(_index_out, 0 );
        }

        unsigned int index_feature=-1;
        assert(index_feature == -1);
        for(auto &[key, value]: A[n]) {
            if( solver->getValue(value.id)) {
                index_feature = key;
            }
        }
        assert(index_feature != -1);

        std::vector< unsigned int > IDs_left;
        std::vector< unsigned int > IDs_right;

        unsigned int maxValue_left = 0;
        unsigned int minValue_right = _vDomain[index_feature].size()-1;

        for(auto id: IDs) {
            if( solver->getValue(S[id][n].id) ) {
                IDs_left.push_back(id);
                if( examples[id][index_feature].has_value() ) {
                    if( examples[id][index_feature].value() > maxValue_left ) {
                        maxValue_left = examples[id][index_feature].value();
                    }
                }
            } else {
                IDs_right.push_back(id);
                if( examples[id][index_feature].has_value() ) {
                    if( examples[id][index_feature].value() < minValue_right ) {
                        minValue_right = examples[id][index_feature].value();
                    }
                }
            }
        }

        if(IDs_left.size() == 0) {
            return getModel(IDs_right, n*2+1);
        }
        if(IDs_right.size() == 0) {
            return getModel(IDs_left, n*2);
        }
        

        unsigned int index_value = minValue_right;

        return std::make_unique<NumericDecisionNodeV2>(index_feature, _vDomain[index_feature].toDouble(index_value), getModel(IDs_left, n*2), getModel(IDs_right, n*2+1) );
        
    }

};
int InferDT_Pouya::Var::nbVar = 0;

// Our encoding with Z3
class InferODT_Z3 {
    std::vector< Domain > _vDomain;
    unsigned int _index_out;
    unsigned int _depth;

    SolverSMT solver;

    struct Var {
        int id;
        static int nbVar;
        Var() : id(++nbVar) { }
    };

    std::vector< std::map<unsigned int, Var> > A; // A[n][j] = true if the feature j is chosen for the split at branching node n.
    std::vector< std::map<unsigned int, Var> > S; // S[i][n] = true if the exemple i is directed towards the left child, if it passes through branching node n.
    std::vector< std::map<unsigned int, Var> > Z; // Z[i][l] = true is the exemple i ends up at leaf node t.
    std::vector< std::map<unsigned int, Var> > C; // C[l][c] = true if the label c is assigned to leaf node l.

    std::vector< std::tuple<Var, Var, Var> > EQ; // EQ[n] = (a, b, c) if the equation on node n is : a * x + b * y > c

    std::vector< std::map< std::tuple<unsigned int, unsigned int>, Var> > AA; // AA[n][j1][j2] = true if the feature j1 and j2 are chosen for the split at branching node n.


    std::vector< std::vector< std::optional<unsigned int> > > examples;


    void init_constraints() {

        A.resize(std::pow(2, _depth));
        AA.resize(std::pow(2, _depth));
        C.resize(std::pow(2, _depth));
        EQ.resize(std::pow(2, _depth));

        // exactly one feature is chosen at each branching node n
        for(unsigned int n = 1; n < A.size(); n++) {
            for(unsigned int j1 = 0; j1 < _vDomain.size(); j1++) {
                if(j1 == _index_out) {
                    continue;
                }
                for(unsigned int j2 = j1+1; j2 < _vDomain.size(); j2++) {
                    if(j2 == _index_out) {
                        continue;
                    }
                    solver.add(-A[n][j1].id);
                    solver.add(-A[n][j2].id);
                    solver.add(0); // (14)
                }
            }
        }

        for(unsigned int n = 1; n < A.size(); n++) {
            for(unsigned int j1 = 0; j1 < _vDomain.size(); j1++) {
                if(j1 == _index_out) {
                    continue;
                }
                for(unsigned int j2 = 0; j2 < _vDomain.size(); j2++) {
                    if(j2 == _index_out) {
                        continue;
                    }
                    for(unsigned int j3 = j2+1; j3 < _vDomain.size(); j3++) {
                        if(j3 == _index_out) {
                            continue;
                        }

                        solver.add(-A[n][j1].id);
                        solver.add(-AA[n][{j2,j3}].id);
                        solver.add(0); // (14-bis)
                    }
                }

            }
        }

        for(unsigned int n = 1; n < A.size(); n++) {
            for(unsigned int j1 = 0; j1 < _vDomain.size(); j1++) {
                if(j1 == _index_out) {
                    continue;
                }
                for(unsigned int j2 = j1+1; j2 < _vDomain.size(); j2++) {
                    if(j2 == _index_out) {
                        continue;
                    }

                    for(unsigned int j3 = 0; j3 < _vDomain.size(); j3++) {
                        if(j3 == _index_out) {
                            continue;
                        }
                        for(unsigned int j4 = j3+1; j4 < _vDomain.size(); j4++) {
                            if(j4 == _index_out) {
                                continue;
                            }

                            if( std::make_tuple(j1, j2) >= std::make_tuple(j3, j4) )
                                continue;

                            solver.add(-AA[n][{j1,j2}].id);
                            solver.add(-AA[n][{j3,j4}].id);
                            solver.add(0); // (14-bis)
                        }
                    }
                }
            }
        }

        for(unsigned int n = 1; n < A.size(); n++) {
           for(auto &[key, val]: A[n]) {
                assert(key != _index_out);
                solver.add(val.id);
            }
            for(auto &[key, val]: AA[n]) {
                assert(std::get<0>(key) != _index_out);
                assert(std::get<1>(key) != _index_out);
                solver.add(val.id);
            }
            solver.add(0); // (15)
        }
    }

public:

    InferODT_Z3(const std::vector< Domain > & domain, unsigned int index_out, unsigned int depth) : _vDomain(domain), _index_out(index_out), _depth(depth) {
        init_constraints();
    }

    void use_all_examples( const std::vector< std::vector< std::optional<unsigned int> > > & examples ) {
        assert(this->examples.size() == 0);

        this->examples = examples;

        S.resize( examples.size() );
        Z.resize( examples.size() );

        auto get_order = [&](unsigned int feature_index){
            std::vector<int> order;
            order.reserve(examples.size());

            for(unsigned int i=0; i<examples.size(); i++) {
                if(examples[i][feature_index].has_value())
                    order.push_back(i);
            }

            std::sort(order.begin(), order.end(), [&](int a, int b){
                return examples[a][feature_index].value() < examples[b][feature_index].value();
            });

            return order;
        };

        for(unsigned int n = 1; n < A.size(); n++) {
            for(auto &[j, val]: A[n]) {
                assert(j != _index_out);

                auto order = get_order(j);
                for(unsigned a = 1; a<order.size(); a++) {
                    unsigned int i1 = order[a-1];
                    assert(examples[i1][j].has_value());
                    unsigned int i2 = order[a];
                    assert(examples[i2][j].has_value());

                    if(_vDomain[j].is_numeric()) {
                        solver.add(-A[n][j].id);
                        solver.add(S[i1][n].id);
                        solver.add(-S[i2][n].id);
                        solver.add(0); // (16)
                    }

                    if( examples[i1][j].value() == examples[i2][j].value() ) {
                        solver.add(-A[n][j].id);
                        solver.add(-S[i1][n].id);
                        solver.add(S[i2][n].id);
                        solver.add(0); // (17)

                        if(_vDomain[j].is_categorical()) {
                            solver.add(-A[n][j].id);
                            solver.add(S[i1][n].id);
                            solver.add(-S[i2][n].id);
                            solver.add(0); // (18)
                        }
                    }
                }

                if(order.size()) {
                    solver.add(-A[n][j].id);
                    solver.add(S[order.front()][n].id);
                    solver.add(0); // (23)

                    if(_vDomain[j].is_numeric() && (examples[order.back()][j].value() != examples[order.front()][j].value())) {
                        solver.add(-A[n][j].id);
                        solver.add(-S[order.back()][n].id);
                        solver.add(0); // (24)
                    }
                }
            }

            for(auto &[key, val]: AA[n]) {
                int j1 = std::get<0>(key);
                int j2 = std::get<1>(key);
                
                for(unsigned int e = 0; e < examples.size(); e++) {
                    assert(examples[e][j1].has_value());
                    assert(examples[e][j2].has_value());

                    solver.add(-AA[n][key].id);
                    solver.add(-S[e][n].id);
                    solver.add_side1( std::get<0>(EQ[n]).id, _vDomain[j1].toDouble( examples[e][j1].value() ), std::get<1>(EQ[n]).id, _vDomain[j2].toDouble( examples[e][j2].value() ), std::get<2>(EQ[n]).id );
                    solver.add(0);
                    
                    solver.add(-AA[n][key].id);
                    solver.add(S[e][n].id);
                    solver.add_side2( std::get<0>(EQ[n]).id, _vDomain[j1].toDouble( examples[e][j1].value() ), std::get<1>(EQ[n]).id, _vDomain[j2].toDouble( examples[e][j2].value() ), std::get<2>(EQ[n]).id );
                    solver.add(0);

                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                unsigned int n = C.size() + l;
                std::vector<int> clause = { Z[i][l].id };

                while(n > 1) {
                    bool right = n%2;
                    n = n / 2;
                    if(!right) {
                        solver.add(-Z[i][l].id );
                        solver.add( S[i][n].id );
                        solver.add(0); // (19)

                        clause.push_back( -S[i][n].id );
                    } else {
                        solver.add(-Z[i][l].id );
                        solver.add( -S[i][n].id );
                        solver.add(0); // (20)
                        
                        clause.push_back( S[i][n].id );
                    }
                }
                solver.add(clause); // (21)
            }
        }

        for(unsigned int l=0; l<C.size(); l++) {
            for(unsigned int c1=0; c1<_vDomain[_index_out].size(); c1++) {
                for(unsigned int c2=c1+1; c2<_vDomain[_index_out].size(); c2++) {
                    solver.add(-C[l][c1].id);
                    solver.add(-C[l][c2].id);
                    solver.add(0); // (22)
                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                if( examples[i][ _index_out ].has_value() ) {
                    solver.add( -Z[i][l].id );
                    solver.add( C[l][ examples[i][ _index_out ].value() ].id );
                    if(GLOBAL_exact) {
                        solver.add(0);
                    } else {
                        solver.add_soft(1);
                    }
                }
            }
        }
    }

    bool infer() {
        std::cout << "solve..." << std::endl;
        return solver.solve();
    }

    DecisionTree getModel() {
        std::vector< unsigned int > IDs;//( examples.size() );
        for(unsigned int i=0; i<examples.size(); i++) {
            bool mal_classe = false;
            for(unsigned int l=0; l<C.size(); l++) {
                assert( examples[i][ _index_out ].has_value() );
                if( (solver.getValue( Z[i][l].id ) == true) && (solver.getValue( C[l][ examples[i][ _index_out ].value() ].id ) == false) ) {
                    mal_classe = true;
                    break;
                }
            }
            if(mal_classe)
                continue;

            IDs.push_back(i);
        }
        return DecisionTree( getModel(IDs, 1), _vDomain );
    }
private:

    std::unique_ptr<VirtualDecisionNode> getModel(const std::vector< unsigned int > &IDs, unsigned int n) {
        if(n >= A.size()) {
            for(auto &[key, value]: C[n - A.size()]) {
                if( solver.getValue(value.id)) {
                    return std::make_unique<LeafDecisionNode>( _index_out, key );
                }
            }
            assert(false);
            return std::make_unique<LeafDecisionNode>(_index_out, 0 );
        }

        unsigned int index_feature=-1;
        assert(index_feature == -1);
        for(auto &[key, value]: A[n]) {
            if( solver.getValue(value.id)) {
                index_feature = key;
            }
        }
        if(index_feature != -1) {
            std::vector< unsigned int > IDs_left;
            std::vector< unsigned int > IDs_right;

            unsigned int maxValue_left = 0;
            unsigned int minValue_right = _vDomain[index_feature].size()-1;

            for(auto id: IDs) {
                if( solver.getValue(S[id][n].id) ) {
                    IDs_left.push_back(id);
                    if( examples[id][index_feature].has_value() ) {
                        if( examples[id][index_feature].value() > maxValue_left ) {
                            maxValue_left = examples[id][index_feature].value();
                        }
                    }
                } else {
                    IDs_right.push_back(id);
                    if( examples[id][index_feature].has_value() ) {
                        if( examples[id][index_feature].value() < minValue_right ) {
                            minValue_right = examples[id][index_feature].value();
                        }
                    }
                }
            }

            if(IDs_left.size() == 0) {
                return getModel(IDs_right, n*2+1);
            }
            if(IDs_right.size() == 0) {
                return getModel(IDs_left, n*2);
            }
            
            unsigned int index_value = minValue_right;

            return std::make_unique<NumericDecisionNodeV2>(index_feature, _vDomain[index_feature].toDouble(index_value), getModel(IDs_left, n*2), getModel(IDs_right, n*2+1) );
        }
        for(auto &[key, value]: AA[n]) {
            int j1 = std::get<0>(key);
            int j2 = std::get<1>(key);
            if( solver.getValue(value.id) ) {
            std::vector< unsigned int > IDs_left;
                std::vector<std::vector<double>> gauche;
            std::vector< unsigned int > IDs_right;
                std::vector<std::vector<double>> droite;
                for(auto id: IDs) {
                    if( solver.getValue(S[id][n].id) ) {
                        IDs_left.push_back(id);
                        gauche.push_back({_vDomain[j1].toDouble(examples[id][j1].value()), _vDomain[j2].toDouble(examples[id][j2].value())});
                    } else {
                        IDs_right.push_back(id);
                        droite.push_back({_vDomain[j1].toDouble(examples[id][j1].value()), _vDomain[j2].toDouble(examples[id][j2].value())});
                    }
                }
                if(gauche.size() == 0) {
                    return getModel(IDs, n*2+1);
                }
                if(droite.size() == 0) {
                    return getModel(IDs, n*2);
                }
                auto separator = findSeparator(gauche, droite);

                if(separator.size() != 3) {
                    return std::make_unique<LinearDecisionNode>(j1, -solver.getRealValue( std::get<0>(EQ[n]).id ), j2, -solver.getRealValue( std::get<1>(EQ[n]).id ), -solver.getRealValue( std::get<2>(EQ[n]).id ), getModel(IDs_left, n*2), getModel(IDs_right, n*2+1));
                }

                return std::make_unique<LinearDecisionNode>(j1, separator[0], j2, separator[1], -separator[2], getModel(IDs_left, n*2), getModel(IDs_right, n*2+1));
            }
        }
        assert(false);
        return nullptr;
    }
};
int InferODT_Z3::Var::nbVar = 0;

// Our encoding with SAT
class InferODT_SAT {
    std::vector< Domain > _vDomain;
    unsigned int _index_out;
    unsigned int _depth;

    AbstractSolveur *solver;

    struct Var {
        int id;
        static int nbVar;
        Var() : id(++nbVar) { }
    };

    std::vector< std::map<unsigned int, Var> > S; // S[i][n] = true if the exemple i is directed towards the left child, if it passes through branching node n.
    std::vector< std::map<unsigned int, Var> > Z; // Z[i][l] = true is the exemple i ends up at leaf node t.
    std::vector< std::map<unsigned int, Var> > C; // C[l][c] = true if the label c is assigned to leaf node l.

    std::vector< std::map< std::tuple<unsigned int, unsigned int>, Var> > AA; // AA[n][j1][j2] = true if the feature j1 and j2 are chosen for the split at branching node n.

    std::vector< std::vector< std::optional<unsigned int> > > examples;

    void init_constraints() {
        AA.resize(std::pow(2, _depth));
        C.resize(std::pow(2, _depth));

        for(unsigned int n = 1; n < AA.size(); n++) {
            for(unsigned int j1 = 0; j1 < _vDomain.size(); j1++) {
                if(j1 == _index_out) {
                    continue;
                }
                for(unsigned int j2 = j1+1; j2 < _vDomain.size(); j2++) {
                    if(j2 == _index_out) {
                        continue;
                    }
                    AA[n][{j1,j2}];

                    for(unsigned int j3 = 0; j3 < _vDomain.size(); j3++) {
                        if(j3 == _index_out) {
                            continue;
                        }
                        for(unsigned int j4 = j3+1; j4 < _vDomain.size(); j4++) {
                            if(j4 == _index_out) {
                                continue;
                            }

                            if( std::make_tuple(j1, j2) >= std::make_tuple(j3, j4) )
                                continue;

                            solver->add(-AA[n][{j1,j2}].id);
                            solver->add(-AA[n][{j3,j4}].id);
                            solver->add(0); // (14-bis)
                        }
                    }
                }
            }
        }

        for(unsigned int n = 1; n < AA.size(); n++) {
            for(auto &[key, val]: AA[n]) {
                assert(std::get<0>(key) != _index_out);
                assert(std::get<1>(key) != _index_out);
                solver->add(val.id);
            }
            solver->add(0); // (15)
        }
    }

public:

    virtual ~InferODT_SAT() {
        delete solver;
    }


    InferODT_SAT(const std::vector< Domain > & domain, unsigned int index_out, unsigned int depth) : _vDomain(domain), _index_out(index_out), _depth(depth) {
        if(GLOBAL_externalSolver.size() == 0) {
            solver = new SolveurMaxSAT();
        } else {
            solver = new ExternalSolver(GLOBAL_externalSolver);
        }
        init_constraints();
    }

    void use_all_examples( const std::vector< std::vector< std::optional<unsigned int> > > & examples ) {
        assert(this->examples.size() == 0);

        this->examples = examples;

        S.resize( examples.size() );
        Z.resize( examples.size() );

        auto getAllTriangleBasGauche = [&](unsigned int e1, unsigned int e2, unsigned int j1, unsigned int j2) {
            std::unordered_set< unsigned int > pointsInTriangle;
            for(unsigned int p=0; p<examples.size(); p++) {
                if(p==e1 || p==e2) {
                    continue;
                }

                if( isInTriangle(
                    -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 
                    _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ), 
                    _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                    _vDomain[j1].toDouble( examples[p][j1].value() ), _vDomain[j2].toDouble( examples[p][j2].value() )
                ) ) {
                    pointsInTriangle.insert(p);
                }
            }
            return pointsInTriangle;
        };

        auto getAllTriangleHautDroite = [&](unsigned int e1, unsigned int e2, unsigned int j1, unsigned int j2) {
            std::unordered_set< unsigned int > pointsInTriangle;
            for(unsigned int p=0; p<examples.size(); p++) {
                if(p==e2 || p==e1) {
                    continue;
                }

                if( isInTriangle(
                    std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 
                    _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ), 
                    _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                    _vDomain[j1].toDouble( examples[p][j1].value() ), _vDomain[j2].toDouble( examples[p][j2].value() )
                ) ) {
                    pointsInTriangle.insert(p);
                }
            }
            return pointsInTriangle;
        };

        unsigned int nbTriangle = 0;
        for(auto &[key, _]: AA[1]) {
            int j1 = std::get<0>(key);
            int j2 = std::get<1>(key);

            std::vector< unsigned int > firsts;
            std::map< std::tuple<double, double >, unsigned int > dejaVue;

            for(unsigned int e1=0; e1<examples.size(); e1++) {
                std::tuple<double, double> pos = { _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ) };
                if( dejaVue.count(pos) ) {
                    unsigned int e = dejaVue[pos];
                    for(unsigned int n = 1; n < AA.size(); n++) {
                        solver->add(-AA[n][key].id);
                        solver->add(-S[e1][n].id);
                        solver->add(S[e][n].id);
                        solver->add(0); // (16)

                        solver->add(-AA[n][key].id);
                        solver->add(S[e1][n].id);
                        solver->add(-S[e][n].id);
                        solver->add(0); // (16)
                    }
                } else {
                    dejaVue[pos] = e1;
                }
            }

            auto getAllTriangleBasGauche2 = [&](unsigned int e1, unsigned int e2, unsigned int j1, unsigned int j2) {
                std::vector< unsigned int > pointsInTriangle;
                for(auto [_, p]: dejaVue) {
                    if(p == e1 || p == e2)
                        continue;
                    if( isInTriangle(
                        -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 
                        _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ), 
                        _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                        _vDomain[j1].toDouble( examples[p][j1].value() ), _vDomain[j2].toDouble( examples[p][j2].value() )
                    ) ) {
                        pointsInTriangle.push_back(p);
                    }
                }
                return pointsInTriangle;
            };

            auto getAllTriangleHautDroite2 = [&](unsigned int e1, unsigned int e2, unsigned int j1, unsigned int j2) {
                std::vector< unsigned int > pointsInTriangle;
                for(auto [_, p]: dejaVue) {
                    if(p == e2 || p == e1)
                        continue;
                    if( isInTriangle(
                        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 
                        _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ), 
                        _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                        _vDomain[j1].toDouble( examples[p][j1].value() ), _vDomain[j2].toDouble( examples[p][j2].value() )
                    ) ) {
                        pointsInTriangle.push_back(p);
                    }
                }
                return pointsInTriangle;
            };

            auto getCoord = [&](unsigned int e, unsigned int j1, unsigned int j2) {
                return std::make_tuple(_vDomain[j1].toDouble( examples[e][j1].value() ), _vDomain[j2].toDouble( examples[e][j2].value() ));
            };

            for(auto it=dejaVue.begin(); it!=dejaVue.end(); ++it) {
                unsigned int e1 = it->second;
                auto it2 = it;
                for(++it2; it2!=dejaVue.end(); ++it2) {
                    unsigned int e2 = it2->second;

                    auto allIncludedBasGauche = getAllTriangleBasGauche2(e1, e2, j1, j2);
                    std::unordered_set<unsigned int> removed_exemples;
                    unsigned int nb_utilise = 0;
                    for(auto e3: allIncludedBasGauche) {                       
                        bool a_utiliser = true;
                        for(auto e4: allIncludedBasGauche) {
                            if(e3 == e4)
                                continue;
                            if(removed_exemples.count(e4))
                                continue;

                            if( isInTriangle(
                                -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 
                                _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ),  
                                _vDomain[j1].toDouble( examples[e4][j1].value() ), _vDomain[j2].toDouble( examples[e4][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e3][j1].value() ), _vDomain[j2].toDouble( examples[e3][j2].value() )
                            ) ) {
                                a_utiliser = false;
                                removed_exemples.insert(e3);
                                break;
                            }
                            
                            if( isInTriangle(
                                -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 
                                _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e4][j1].value() ), _vDomain[j2].toDouble( examples[e4][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e3][j1].value() ), _vDomain[j2].toDouble( examples[e3][j2].value() )
                            ) ) {
                                a_utiliser = false;
                                removed_exemples.insert(e3);
                                break;
                            }
                            
                        }
                        if(a_utiliser) {
                            for(unsigned int n = 1; n < AA.size(); n++) {
                                solver->add(-AA[n][key].id);
                                solver->add(-S[e1][n].id);
                                solver->add(-S[e2][n].id);
                                solver->add(S[e3][n].id);
                                solver->add(0); // (16)
                            }
                            nb_utilise++;
                        }
                    }
                    //fmt::println("nb_utilise: {}", nb_utilise);

                    auto allIncludedHautDroite = getAllTriangleHautDroite2(e1, e2, j1, j2);
                    removed_exemples.clear();
                    for(auto e3: allIncludedHautDroite) {
                        bool a_utiliser = true;
                        for(auto e4: allIncludedHautDroite) {
                            if(e3 == e4)
                                continue;
                            if(removed_exemples.count(e4))
                                continue;

                            if( isInTriangle(
                                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 
                                _vDomain[j1].toDouble( examples[e1][j1].value() ), _vDomain[j2].toDouble( examples[e1][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e4][j1].value() ), _vDomain[j2].toDouble( examples[e4][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e3][j1].value() ), _vDomain[j2].toDouble( examples[e3][j2].value() )
                            ) ) {
                                a_utiliser = false;
                                removed_exemples.insert(e3);
                                break;
                            }
                            if( isInTriangle(
                                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 
                                _vDomain[j1].toDouble( examples[e2][j1].value() ), _vDomain[j2].toDouble( examples[e2][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e4][j1].value() ), _vDomain[j2].toDouble( examples[e4][j2].value() ), 
                                _vDomain[j1].toDouble( examples[e3][j1].value() ), _vDomain[j2].toDouble( examples[e3][j2].value() )
                            ) ) {
                                a_utiliser = false;
                                removed_exemples.insert(e3);
                                break;
                            }
                            
                        }
                        if(a_utiliser) {
                            for(unsigned int n = 1; n < AA.size(); n++) {
                                solver->add(-AA[n][key].id);
                                solver->add(S[e1][n].id);
                                solver->add(S[e2][n].id);
                                solver->add(-S[e3][n].id);
                                solver->add(0); // (16)  
                            }
                        }
                    }
                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                unsigned int n = C.size() + l;
                std::vector<int> clause = { Z[i][l].id };

                while(n > 1) {
                    bool right = n%2;
                    n = n / 2;
                    if(!right) {
                        solver->add(-Z[i][l].id );
                        solver->add( S[i][n].id );
                        solver->add(0); // (19)

                        clause.push_back( -S[i][n].id );
                    } else {
                        solver->add(-Z[i][l].id );
                        solver->add( -S[i][n].id );
                        solver->add(0); // (20)
                        
                        clause.push_back( S[i][n].id );
                    }
                }
                solver->add(clause); // (21)
            }
        }

        for(unsigned int l=0; l<C.size(); l++) {
            for(unsigned int c1=0; c1<_vDomain[_index_out].size(); c1++) {
                for(unsigned int c2=c1+1; c2<_vDomain[_index_out].size(); c2++) {
                    solver->add(-C[l][c1].id);
                    solver->add(-C[l][c2].id);
                    solver->add(0); // (22)
                }
            }
        }

        for(unsigned int i = 0; i < examples.size(); i++) {
            for(unsigned int l=0; l<C.size(); l++) {
                if( examples[i][ _index_out ].has_value() ) {
                    solver->add( -Z[i][l].id );
                    solver->add( C[l][ examples[i][ _index_out ].value() ].id );
                    if(GLOBAL_exact) {
                        solver->add(0);
                    } else {
                        solver->add_soft(1);
                    }
                }
            }
        }
    }


    bool infer() {
        fmt::print("solver.solve()...\n");
        std::cout << "getNbClause(): " << solver->getNbClause() << std::endl;
        return solver->solve();
    }


    DecisionTree getModel() {

        Chrono Chro("Temps pour décodage");

        std::vector<int> nb_in_leaf(C.size(), 0);
        std::vector<int> err_in_leaf(C.size(), 0);

        std::vector< unsigned int > IDs;
        for(unsigned int i=0; i<examples.size(); i++) {
            bool mal_classe = false;
            for(unsigned int l=0; l<C.size(); l++) {
                assert( examples[i][ _index_out ].has_value() );
                if( (solver->getValue( Z[i][l].id ) == true) ) {
                    nb_in_leaf[l]++;
                }
                if( (solver->getValue( Z[i][l].id ) == true) && (solver->getValue( C[l][ examples[i][ _index_out ].value() ].id ) == false) ) {
                    mal_classe = true;
                    err_in_leaf[l]++;
                }
            }
            if(mal_classe)
                continue;

            IDs.push_back(i);
        }
        return DecisionTree( getModel(IDs, 1, nb_in_leaf, err_in_leaf), _vDomain );
    }
private:

    std::unique_ptr<VirtualDecisionNode> getModel(const std::vector< unsigned int > &IDs, unsigned int n, const std::vector<int> &nb_in_leaf, const std::vector<int> &err_in_leaf) {
        if(n >= AA.size()) {
            for(auto &[key, value]: C[n - AA.size()]) {
                if( solver->getValue(value.id)) {
                    return std::make_unique<LeafDecisionNode>( _index_out, key, nb_in_leaf[n - AA.size()], 1 - err_in_leaf[n - AA.size()] / (double)nb_in_leaf[n - AA.size()] );
                }
            }
            assert(false);
            return std::make_unique<LeafDecisionNode>(_index_out, 0, -1);
        }

        for(auto &[key, value]: AA[n]) {
            int j1 = std::get<0>(key);
            int j2 = std::get<1>(key);
            if( solver->getValue(value.id) ) {
                std::vector< unsigned int > IDs_left;
                std::vector<std::vector<double>> gauche;
                std::vector< unsigned int > IDs_right;
                std::vector<std::vector<double>> droite;
                for(auto id: IDs) {
                    if( solver->getValue(S[id][n].id) ) {
                        IDs_left.push_back(id);
                        gauche.push_back({_vDomain[j1].toDouble(examples[id][j1].value()), _vDomain[j2].toDouble(examples[id][j2].value())});
                    } else {
                        IDs_right.push_back(id);
                        droite.push_back({_vDomain[j1].toDouble(examples[id][j1].value()), _vDomain[j2].toDouble(examples[id][j2].value())});
                    }
                }
                if(gauche.size() == 0) {
                    return getModel(IDs, n*2+1, nb_in_leaf, err_in_leaf);
                }
                if(droite.size() == 0) {
                    return getModel(IDs, n*2, nb_in_leaf, err_in_leaf);
                }
                auto separator = findSeparator(gauche, droite);

                if(separator.size() != 3) {
                    assert(false);
                    std::cerr << "separator.size() != 3" << std::endl;
                }

                return std::make_unique<LinearDecisionNode>(j1, separator[0], j2, separator[1], -separator[2], getModel(IDs_left, n*2, nb_in_leaf, err_in_leaf), getModel(IDs_right, n*2+1, nb_in_leaf, err_in_leaf));
            }
        }
        assert(false);
        return nullptr;
    }

};
int InferODT_SAT::Var::nbVar = 0;

template<class Algo>
double test_synthetic(unsigned int numberExemples=64, unsigned int numberFeature = 4, unsigned int NumberValuePerFeature=100, unsigned int numberClass=2, unsigned int depth=3, unsigned int nbItt=10) {
    Mean moyenne;
    Chrono C;

    auto arbre = DecisionTree::randomLinear(numberFeature, NumberValuePerFeature, numberClass, depth);

    for(unsigned int it=0; it<nbItt; it++) {
        unsigned int seed = MaLib::MonRand::get();
        srand(seed);
        MaLib::MonRand::seed(seed);

        auto data = arbre.generate_dataset(numberExemples);

        C.tic();
        Algo infer(data.getDomain(), data.numberFeature()-1, depth);
        infer.use_all_examples(data.getData());
        auto res = infer.infer();
        moyenne.add( C.tacSec() );
        if(!res) {
            fmt::println("==================  seed = {} =====================", seed);
            fmt::println("No solution :(");

            fmt::println("data = {}", data);
            xdot::show(arbre.toDot());
            assert(false);
            exit(-1);
        }

    }
    return moyenne.getMean();
}


template<class Algo>
void test(unsigned int numberFeature_Max = 4, unsigned int NumberValuePerFeature_Max = 10, unsigned int numberClass_Max = 5, unsigned int maxDepth_Max = 4) {

    for(unsigned int i=0; i<10000; i++) {
        unsigned int seed = MaLib::MonRand::get();
        srand(seed);
        MaLib::MonRand::seed(seed);
        fmt::println("==================  seed = {} =====================", seed);

        unsigned int nombre_features = MaLib::MonRand::get(2, numberFeature_Max);
        unsigned int nombre_valeur_max_par_features = MaLib::MonRand::get(3, NumberValuePerFeature_Max);
        double proba_feuille = 0;
        unsigned int nb_class= MaLib::MonRand::get(2, numberClass_Max);
        unsigned int max_depth= MaLib::MonRand::get(2, maxDepth_Max);

        auto arbre = DecisionTree::random(nombre_features, nombre_valeur_max_par_features, proba_feuille, nb_class, max_depth);

        for(unsigned int nbData=2; ; nbData*=2) {
            auto data = arbre.generate_dataset(nbData);
            fmt::println("nbData = {}", nbData);
            fmt::println("accuracy oracle: {}%", arbre.accuracy(data));
            if( arbre.accuracy(data) != 100 ) {
                fmt::println("data = {}", data);
                xdot::show(arbre.toDot());
                assert(false);
                exit(-1);
            }

            Algo infer(data.getDomain(), data.numberFeature()-1, max_depth);

            infer.use_all_examples(data.getData());

            auto res = infer.infer();
            if(!res) {
                fmt::println("==================  seed = {} =====================", seed);
                fmt::println("No solution :(");

                fmt::println("data = {}", data);
                xdot::show(arbre.toDot());
                assert(false);
                exit(-1);
            }

            fmt::println("seed = {}", seed);
            auto dt = infer.getModel();

            fmt::println("Training accuracy inferDT : {}%", dt.accuracy(data));
            if( dt.accuracy(data) != 100 ) {
                fmt::println("data = {}", data);
                xdot::show(arbre.toDot(), false);
                xdot::show(dt.toDot());
                fmt::println("seed = {}", seed);
                assert(false);
                fmt::println("accuracy inferDT : {}%", dt.accuracy(data));
                exit(-1);
            }

            auto data2 = arbre.generate_dataset(1000);
            auto testAccuray = dt.accuracy(data2);
            fmt::println("accuracy test : {}%", testAccuray);
            if(testAccuray==100)
                break;
        }
        fmt::println("");
    }
}


template<class Algo>
void test2(unsigned int numberFeature_Max = 4, unsigned int NumberValuePerFeature_Max = 10, unsigned int numberClass_Max = 5, unsigned int maxDepth_Max = 4) {

    for(unsigned int i=0; i<10000; i++) {
        unsigned int seed = MaLib::MonRand::get();
        srand(seed);
        MaLib::MonRand::seed(seed);
        fmt::println("==================  seed = {} =====================", seed);

        unsigned int nombre_features = MaLib::MonRand::get(2, numberFeature_Max);
        unsigned int nombre_valeur_max_par_features = MaLib::MonRand::get(3, NumberValuePerFeature_Max);
        double proba_feuille = 0;
        unsigned int nb_class= 2;
        unsigned int max_depth= MaLib::MonRand::get(3, maxDepth_Max);

        auto arbre = DecisionTree::random(nombre_features, nombre_valeur_max_par_features, proba_feuille, nb_class, max_depth);

        unsigned int nbData=32;
        auto data = arbre.generate_dataset(nbData);
        fmt::println("nbData = {}", nbData);
        fmt::println("accuracy oracle: {}%", arbre.accuracy(data));
        if( arbre.accuracy(data) != 100 ) {
            fmt::println("data = {}", data);
            xdot::show(arbre.toDot());
            assert(false);
            exit(-1);
        }

        Algo infer(data.getDomain(), data.numberFeature()-1, 1);
        InferODT_Z3 inferCheck(data.getDomain(), data.numberFeature()-1, 1);

        infer.use_all_examples(data.getData());
        inferCheck.use_all_examples(data.getData());

        auto res = infer.infer();
        auto resCheck = inferCheck.infer();

        if(!res) {
            fmt::println("==================  seed = {} =====================", seed);
            fmt::println("No solution :(");

            fmt::println("data = {}", data);
            xdot::show(arbre.toDot());
            assert(false);
            exit(-1);
        }

        fmt::println("seed = {}", seed);
        auto dt = infer.getModel();
        auto dtCheck = inferCheck.getModel();

        fmt::println("Training accuracy inferDT : {}%", dt.accuracy(data));
        if( dt.accuracy(data) != dtCheck.accuracy(data) ) {
            fmt::println("data = {}", data);
            fmt::println("accuracy inferDT : {}%", dt.accuracy(data));
            fmt::println("accuracy inferDT Check: {}%", dtCheck.accuracy(data));
            xdot::show(dt.toDot());
            xdot::show(dtCheck.toDot());
            exit(-1);
        }
    }
}

template<class Algo>
int monMain(std::string file_path, unsigned int k=2, unsigned int numberSplit=10) {

    if(numberSplit == 0)
    {
        MesData data(file_path, false);
        fmt::println("data = {}", data);

        Algo infer(data.getDomain(), data.numberFeature()-1, k);

        infer.use_all_examples(data.getData());

        auto res = infer.infer();
        if(res) {
            auto arbre = infer.getModel();
            fmt::println("accuracy inferDT : {}%", arbre.accuracy(data));
            if(GLOBAL_no_show == false) {
                xdot::show(arbre.toDot());
            }
        } else {
            fmt::println("No Solution :(");
        }
        return 0;
    }

    if(numberSplit == -1) // 1 out cross validation
    {
        MaLib::Chrono C;
        
        MesData data(file_path, false);
        MaLib::Mean mean_train("Training Accuracy");
        MaLib::Mean mean_test("Testing Accuracy");
        
        for(unsigned int i=0; i<data.getData().size(); i++) {
            MesData train = data;
            MesData test({ data.getData()[i] }, data.getDomain());
            train.remove(i);

            Algo infer(train.getDomain(), train.numberFeature()-1, k);
            infer.use_all_examples(train.getData());

            auto res = infer.infer();
            if(res) {
                auto arbre = infer.getModel();
                fmt::println("train accuracy : {}%", arbre.accuracy(train));
                fmt::println("test inferDT : {}%", arbre.accuracy(test));

                mean_train.add( arbre.accuracy(train) );
                mean_test.add( arbre.accuracy(test) );
            } else {
                fmt::println("No solution :(");
            }
        }

        mean_train.print();
        mean_test.print();
        return 0;
    }

    // k cross validation
    {
        MesData data(file_path, false);
        MaLib::Mean mean("accuracy");

        assert(numberSplit <= data.getData().size());
        
        auto datas = data.split_for_cross_validation(numberSplit);
        for(auto &[train, test]: datas) {
            Algo infer(train.getDomain(), train.numberFeature()-1, k);


            for(auto &v: train.getData()) {
                if(v.back().value() == 0) {
                    std::cout << train.getDomain()[0].toDouble(v[0].value()) << "\t" << train.getDomain()[2].toDouble(v[2].value()) << std::endl;
                }
            }
            std::cout << "-----" << std::endl;
            for(auto &v: train.getData()) {
                if(v.back().value() == 1) {
                    std::cout << train.getDomain()[0].toDouble(v[0].value()) << "\t" << train.getDomain()[2].toDouble(v[2].value()) << std::endl;
                }
            }
            
            infer.use_all_examples(train.getData());

            auto res = infer.infer();
            if(res) {
                auto arbre = infer.getModel();
                fmt::println("train accuracy : {}%", arbre.accuracy(train));
                fmt::println("test inferDT : {}%", arbre.accuracy(test));

                //xdot::show(arbre.toDot());
                mean.add( arbre.accuracy(test) );
            } else {
                fmt::println("No Solution :(");
            }
        }

        mean.print();
        return 0;
    }


    return 0;
}


void exp1() {
    fmt::println("SAT, numberExemples");
    std::vector< std::tuple<unsigned int, double> > res;

    unsigned int numberExemples=64;
    unsigned int numberFeature = 4;
    unsigned int NumberValuePerFeature=100;
    unsigned int numberClass=2;
    unsigned int depth=3;
    unsigned int nbItt=10;
    for(unsigned int numberExemples=8; numberExemples<=512; numberExemples*=2) {
        auto r = test_synthetic<InferODT_SAT>(numberExemples, numberFeature, NumberValuePerFeature, numberClass, depth, nbItt);
        res.push_back({numberExemples, r});
        fmt::println("res = {}", res);
    }
}

void exp2() {
    fmt::println("Z3, numberExemples");
    std::vector< std::tuple<unsigned int, double> > res;

    unsigned int numberExemples=64;
    unsigned int numberFeature = 4;
    unsigned int NumberValuePerFeature=100;
    unsigned int numberClass=2;
    unsigned int depth=3;
    unsigned int nbItt=10;
    for(unsigned int numberExemples=8; numberExemples<=512; numberExemples*=2) {
        auto r = test_synthetic<InferODT_Z3>(numberExemples, numberFeature, NumberValuePerFeature, numberClass, depth, nbItt);
        res.push_back({numberExemples, r});
        fmt::println("res = {}", res);
    }
}

void exp3() {
    fmt::println("SAT, depth, feature=depth+1");
    std::vector< std::tuple<unsigned int, double> > res;

    unsigned int numberExemples=64;
    unsigned int numberFeature = 4;
    unsigned int NumberValuePerFeature=100;
    unsigned int numberClass=2;
    unsigned int depth=3;
    unsigned int nbItt=10;
    for(unsigned int depth=1; depth<=8; depth++) {
        auto r = test_synthetic<InferODT_SAT>(numberExemples, depth+1, NumberValuePerFeature, numberClass, depth, nbItt);
        res.push_back({numberExemples, r});
        fmt::println("res = {}", res);
    }
}

void exp4() {
    fmt::println("Z3, depth, feature=depth+1");
    std::vector< std::tuple<unsigned int, double> > res;

    unsigned int numberExemples=64;
    unsigned int numberFeature = 4;
    unsigned int NumberValuePerFeature=100;
    unsigned int numberClass=2;
    unsigned int depth=3;
    unsigned int nbItt=10;
    for(unsigned int depth=1; depth<=8; depth++) {
        auto r = test_synthetic<InferODT_Z3>(numberExemples, depth+1, NumberValuePerFeature, numberClass, depth, nbItt);
        res.push_back({numberExemples, r});
        fmt::println("res = {}", res);
    }
}




int main(int argc, char *argv[])
{
    assert([](){std::cout << "c Assertion activated. For better performance, compile the project with assertions disabled. (-DNDEBUG)" << std::endl; return true;}());
    MaLib::Chrono C("Total time");

    CLI::App app("InferOptimalObliqueDT");

    unsigned int exp=0;
    app.add_option("--exp", exp, "Experiment for Fig2 and Fig3 (1, 2, 3 or 4)");

    unsigned int seed = 0;
    app.add_option("--seed", seed, "seed (0 means random seed)");

    app.add_option("-v", verbosity, "vebosity (default: 1)");

    unsigned int depth = 2;
    app.add_option("-d", depth, "depth (default: 2)");

    GLOBAL_exact=false;
    app.add_flag("--exact", GLOBAL_exact, "Exact classification");

    bool Z3=false;
    app.add_flag("--Z3", Z3, "Use Z3 encodage");

    int kcross = 0;
    app.add_option("--kcross", kcross, "k cross validation (default: 0 (no cross validation), -1 for one out cross validation)");

    bool shati=false;
    app.add_flag("--shati", shati, "Use Shati et al. encodage");

    GLOBAL_no_show=false; // Ne pas afficher l'arbre
    app.add_flag("--noshow", GLOBAL_no_show, "Don't show the tree");

    GLOBAL_externalSolver = "";
    app.add_option("--solver", GLOBAL_externalSolver, "External solver cmd (internal solver by default)");

    std::string file = "";
    app.add_option("CSV_file", file, "CSV file")->check(CLI::ExistingFile);//->required();

    CLI11_PARSE(app, argc, argv);

    if(seed==0) {
        seed = time(NULL);
        fmt::println("seed = {}", seed);
    }

    srand(seed);
    MaLib::MonRand::seed(seed);

    if(exp) {
        if(exp==1) {
            exp1();
        } else if(exp==2) {
            exp2();
        } else if(exp==3) {
            exp3();
        } else if(exp==4) {
            exp4();
        }
        return 0;
    }

    if(file.size() == 0) {
        fmt::println("No csv file specified.");
        exit(-1);
    }

    if(Z3) {
        if(shati) {
            fmt::print("Z3 and Shati are mutually exclusive.\n");
            exit(-1);
        }
        return monMain<InferODT_Z3>(file, depth, kcross);
    }
    if(shati) {
        return monMain<InferDT_Pouya>(file, depth, kcross);
    }

    return monMain<InferODT_SAT>(file, depth, kcross);
}


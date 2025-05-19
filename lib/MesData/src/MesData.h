#ifndef MESDATA_06KS11IZE
#define MESDATA_06KS11IZE


#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/format.h>
#include <optional>


#include "csv.hpp"

#include "Domain.h"

class MesData {
    std::vector< std::vector< std::optional<unsigned int> > > _vvDATA;
    std::vector< Domain > _vDomain;


    void initialize(const std::vector<std::string> &feature_names, const std::vector<std::vector<std::string>> &data, const std::vector< std::set<std::string> > &domain) {
        _vDomain.reserve(domain.size());
        std::vector< std::map<std::string, int> > mCat2int;
        mCat2int.resize(domain.size());
        for(unsigned int j=0; j<domain.size(); j++) {
            _vDomain.emplace_back(feature_names[j], domain[j]);
        }

        _vvDATA.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            _vvDATA[i].reserve(data[i].size());
            for(unsigned int j=0; j<data[i].size(); j++) {
                if(data[i][j].size() == 0) {
                    _vvDATA[i].push_back(std::nullopt);
                } else {
                    _vvDATA[i].push_back( _vDomain[j].getIndex(data[i][j]) );
                }
            }
        }
    }

public:

    MesData(const std::vector< std::vector< std::optional<unsigned int> > > &data, const std::vector< Domain > &domain) : _vvDATA(data), _vDomain(domain) {

    }

    MesData(const std::string & file_path, bool one_hot_encoding = false) {
        assert(one_hot_encoding == false);
        std::vector<std::vector<std::string>> data;
        std::vector< std::set<std::string> > domain;
        std::vector<std::string> feature_names;
        {
            csv::CSVReader reader(file_path);

            feature_names = reader.get_col_names();
            domain.resize(feature_names.size());

            for(auto& row: reader) {
                if(row.size() == 0)
                    continue;

                std::vector<std::string> features;
                assert(row.size() == feature_names.size());
                for(unsigned int i=0; i<row.size(); i++) {
                    std::string v = row[i].get<>();
                    features.push_back(v);

                    if(v.size()) {
                        domain[i].insert(v);
                    }
                }

                data.push_back(features);
            }
        }
        initialize(feature_names, data, domain);
    }

    const std::vector< Domain > & getDomain() const {
        return _vDomain;
    }

    const unsigned int numberFeature() const {
        return _vDomain.size();
    }

    const std::vector< std::vector< std::optional<unsigned int> > > getData() const {
        return _vvDATA;
    }

    void remove(unsigned int i) {
        _vvDATA.erase(_vvDATA.begin() + i);
    }

    const std::vector< std::vector< double > > getDataDouble() const {
        std::vector< std::vector< double > > res;
        res.reserve(_vvDATA.size());
        for(auto & row: _vvDATA) {
            std::vector< double > row_double;
            row_double.reserve(row.size());
            for(auto & v: row) {
                if(v.has_value()) {
                    row_double.push_back( _vDomain[v.value()].toDouble(v.value()) );
                } else {
                    row_double.push_back( 0 );
                }
            }
            res.push_back(row_double);
        }
        return res;
    }

    const std::vector< double > getDataDouble(int i) const {
        std::vector< double > res;
        res.reserve(_vvDATA[i].size());
        for(unsigned int j=0; j<_vvDATA[i].size(); j++) {
            if(_vvDATA[i][j].has_value()) {
                res.push_back( _vDomain[j].toDouble(_vvDATA[i][j].value()) );
            } else {
                res.push_back( 0 );
            }
        }
        return res;
    }

    std::vector< std::tuple<MesData, MesData> > split_for_cross_validation(unsigned int k) {
        std::vector< std::tuple<MesData, MesData> > res;
        std::random_shuffle(_vvDATA.begin(), _vvDATA.end());

        for(unsigned int i=0; i<k; i++) {
            std::vector< std::vector< std::optional<unsigned int> > > data_train;
            std::vector< std::vector< std::optional<unsigned int> > > data_test;

            for(unsigned int j=0; j<_vvDATA.size(); j++) {
                if(j%k == i) {
                    data_test.push_back(_vvDATA[j]);
                } else {
                    data_train.push_back(_vvDATA[j]);
                }
            }

            res.push_back( std::make_tuple(MesData(data_train, _vDomain), MesData(data_test, _vDomain)) );
        }

        return res;
    }


    unsigned int size() const {
        return _vvDATA.size();
    }

    template <typename T>
    friend struct fmt::formatter;
};


template <> struct fmt::formatter<MesData> {
  constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

  template <typename FormatContext>
  auto format(const MesData &d, FormatContext &ctx) const {
    fmt::format_to(ctx.out(), "Features:\n");
    for(auto &dom: d._vDomain) {
        fmt::format_to(ctx.out(), "- {}\n", dom);
    }

    std::vector<unsigned int> max_size_text(d._vDomain.size(), 0);
    for(unsigned int i=0; i<d._vvDATA.size(); i++) {
        for(unsigned int j=0; j<d._vvDATA[i].size(); j++) {
            if(d._vvDATA[i][j].has_value()) {
                unsigned int size = d._vDomain[j].toString(d._vvDATA[i][j].value()).size()
                + std::to_string( d._vvDATA[i][j].value() ).size() + 5;
                if(max_size_text[j] < size) {
                    max_size_text[j] = size;
                }
            }
        }
    }

    unsigned int size_num = std::to_string(d._vvDATA.size()-1).size() + 5;

    fmt::format_to(ctx.out(), "Data:\n");

    fmt::format_to(ctx.out(), "  {:<{}}", "", size_num);
    for(unsigned int j=0; j<d._vDomain.size(); j++) {
        fmt::format_to(ctx.out(), "{:<{}} ", d._vDomain[j].getName(), max_size_text[j]);
    }
    fmt::format_to(ctx.out(), "\n");

    for(unsigned int i=0; i<d._vvDATA.size(); i++) {
        fmt::format_to(ctx.out(), "  {:<{}}", std::to_string(i) + ":", size_num);
        for(unsigned int j=0; j<d._vvDATA[i].size(); j++) {
            if(d._vvDATA[i][j].has_value()) {
                fmt::format_to(ctx.out(), "{:<{}}", fmt::format("{}", d._vvDATA[i][j].value()) + " (" + d._vDomain[j].toString(d._vvDATA[i][j].value() ) + ")", max_size_text[j]);
            } else {
                fmt::format_to(ctx.out(), "{:<{}}", "?", max_size_text[j]);
            }
            fmt::format_to(ctx.out(), " ");
        }
        fmt::format_to(ctx.out(), "\n");
    }

    return ctx.out();
  }
};

#endif

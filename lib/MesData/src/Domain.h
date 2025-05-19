#ifndef DOMAIN_8SKGI93
#define DOMAIN_8SKGI93

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <optional>

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/format.h>

class Domain {
    std::string name;
    bool isNumeric;

    std::vector<std::string> values;
    std::map<std::string, unsigned int> mValues2index;

    std::vector<double> values_double;
    std::map<double, unsigned int> mValuesDouble2index;

public:

    Domain(const std::string & name, const std::set<std::string> & sValues) : name(name), isNumeric(true) {
        assert(sValues.size() > 0);

        try {
            for(auto & v : sValues) {
                std::size_t pos;
                std::stod(v, &pos);  

                if(pos != v.length()) {
                    isNumeric = false;
                    break;
                }
            }
        } catch (const std::invalid_argument&) {
            isNumeric = false;
        } catch (const std::out_of_range&) {
            isNumeric = false;
        }

        if(isNumeric) {
            values_double.reserve(sValues.size());
            for(auto & v : sValues) {
                values_double.push_back(std::stod(v));
            }
            std::sort(values_double.begin(), values_double.end());
            values_double.erase(std::unique(values_double.begin(), values_double.end()), values_double.end());
            for(unsigned int i = 0; i<values_double.size(); i++) {
                mValuesDouble2index[ values_double[i] ] = i;
            }
        } else {
            values = std::vector<std::string>(sValues.begin(), sValues.end());
            values_double.reserve(values.size());
            for(unsigned int i=0; i<values.size(); i++) {
                values_double.push_back(i);
                mValues2index[values[i]] = i;
                mValuesDouble2index[i] = i;
            }
        }
    }

    Domain(const std::string & name, const std::set<double> & sValues) : name(name), isNumeric(true) {
        assert(sValues.size() > 0);
        values_double.reserve(sValues.size());
        for(auto & v : sValues) {
            values_double.push_back(v);
            mValuesDouble2index[v] = values_double.size()-1;
        }
    }

    std::string getName() const {
        return name;
    }

    unsigned int getIndex(const double & v) const {
        assert(mValuesDouble2index.find(v) != mValuesDouble2index.end());
        return mValuesDouble2index.at(v);
    }

    unsigned int getIndex(const std::string & v) const {
        if(isNumeric) {
            return getIndex(std::stod(v));
        }
        assert(mValues2index.find(v) != mValues2index.end());
        return mValues2index.at(v);
    }

    std::string toString(unsigned int index) const {
        if(isNumeric) {
            assert(index < values_double.size());
            return fmt::format("{}", values_double[index]);
        }
        assert(index < values.size());
        return values.at(index);
    }
    
    double toDouble(unsigned int index) const {
        assert(values_double.size() > index);
        if(isNumeric) {
            return values_double[index];
        }
        return values_double.at(index);
    }


    unsigned int size() const {
        return values_double.size();
    }

    bool is_categorical() const {
        return !isNumeric;
    }

    bool is_numeric() const {
        return isNumeric;
    }

    bool isMin(std::string v) const {
        assert(isNumeric);
        assert(values_double.size() > 0);
        return std::stod(v) == values_double.front();
    }

    bool isMax(std::string v) const {
        assert(isNumeric);
        assert(values_double.size() > 0);
        return std::stod(v) == values_double.back();
    }

    template <typename T>
    friend struct fmt::formatter;
};


template <> struct fmt::formatter<Domain> {
  constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

  template <typename FormatContext>
  auto format(const Domain &d, FormatContext &ctx) const {
    assert(d.size() > 0);
    if(d.is_numeric()) {
        return fmt::format_to(ctx.out(), "{} (numeric) = [{}, {}]", d.getName(), d.toDouble(0), d.toDouble(d.size()-1));
    }
    return fmt::format_to(ctx.out(), "{} (categorical) = {}", d.getName(), d.values);
  }
};


#endif

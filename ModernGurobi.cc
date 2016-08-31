#include "ModernGurobi.hh"

namespace ModernGurobi {

LinearExpr operator+(LinearExpr x, const LinearExpr &y)
{ return x += y; }

LinearExpr operator-(LinearExpr x, const LinearExpr &y)
{ return x -= y; }

LinearExpr operator*(const LinearExpr &x, double d)
{
    LinearExpr ret;
    LinearExpr::CoeffMap &cmap = ret.coeffmap();
    for(auto entry: x.coeffmap()) {
        cmap[entry.first] = d * entry.second;
    }
    return ret;
}

LinearExpr operator*(double d, const LinearExpr &x)
{ return x*d; }

LinearConstraint operator<=(const LinearExpr &expr, double rhs)
{ return LinearConstraint({expr, GRB_LESS_EQUAL, rhs}); }

LinearConstraint operator<=(double lhs, const LinearExpr &expr)
{ return expr >= lhs; }

LinearConstraint operator>=(const LinearExpr &expr, double rhs)
{ return LinearConstraint({expr, GRB_GREATER_EQUAL, rhs}); }

LinearConstraint operator>=(double lhs, const LinearExpr &expr)
{ return expr <= lhs; }

LinearConstraint operator==(const LinearExpr &expr, double rhs)
{ return LinearConstraint({expr, GRB_EQUAL, rhs}); }

LinearConstraint operator==(double lhs, const LinearExpr &expr)
{ return expr == lhs; }

void LinearConstraint::add_to_model(GRBmodel *model) const {
    std::vector<int> indices;
    std::vector<double> coeffs;
    const LinearExpr::CoeffMap &parts = expr_.coeffmap();
    for(auto &entry: parts) {
        indices.push_back(static_cast<int>(entry.first->idx()));
        coeffs.push_back(entry.second);
    }
    int error = GRBaddconstr(model,
                             static_cast<int>(parts.size()),
                             &indices[0],
                             &coeffs[0],
                             sense_, rhs_, nullptr);
    if (error) {
        throw GurobiException("LinearConstraint::add_to_model(): GRBaddconstr() failed", error);
    }
}

} // namespace

#include "kernels.h"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
    try {
        require(pc_abi_version()==1,"ABI version");
        // Two cells: volume=.3, face=.45. H(12)-H(8)=36,
        // H(16)-H(14)=12; exchange=.45*(16-12)=1.8.
        std::array<double,2> old_t{8,14}, t{12,16}, q{0,0}, r{}, d{};
        std::array<double,1> e{};
        auto assemble=[&](int side) { return pc_assemble(2,6,24,10,14,1.5,1,.3,1,
            0,0,0,0,side,old_t.data(),t.data(),q.data(),r.data(),d.data(),e.data()); };
        require(assemble(0)==0,"assembly status");
        require(std::abs(r[0]-9)<1e-13 && std::abs(r[1]-5.4)<1e-13,"two cell energy");
        require(std::abs(d[0]-4.05)<1e-13 && std::abs(d[1]-2.25)<1e-13
            && std::abs(e[0]+.45)<1e-13,"hand Jacobian");
        t={10,14}; require(assemble(0)==0,"right endpoint");
        require(std::abs(d[0]-4.05)<1e-13 && std::abs(d[1]-2.25)<1e-13,"right tangent");
        require(assemble(1)==0,"left endpoint");
        require(std::abs(d[0]-2.25)<1e-13 && std::abs(d[1]-4.05)<1e-13,"left tangent");
        // Positive definite matrix [[4,-1,0],[-1,4,-1],[0,-1,3]], x=[1,2,3].
        std::array<double,3> diag{4,4,3}, rhs{2,4,7}, x{};
        std::array<double,2> off{-1,-1};
        require(pc_solve(3,diag.data(),off.data(),rhs.data(),x.data())==0,"LDL status");
        for (std::size_t i=0;i<3;++i) require(std::abs(x[i]-static_cast<double>(i+1))<1e-14,"known solve");
        double scalar_d=2, scalar_rhs=6, scalar_x=0;
        require(pc_solve(1,&scalar_d,nullptr,&scalar_rhs,&scalar_x)==0 && scalar_x==3,"one cell solve");
        x.fill(99); diag[1]=-1;
        require(pc_solve(3,diag.data(),off.data(),rhs.data(),x.data())==3,"indefinite rejected");
        require(x[0]==99 && x[1]==99 && x[2]==99,"failed solve leaves outputs untouched");
        r.fill(99); t[1]=std::numeric_limits<double>::quiet_NaN();
        require(assemble(0)==1 && r[0]==99 && r[1]==99,"invalid trial atomic rejection");
        require(pc_solve(0,nullptr,nullptr,nullptr,nullptr)==1,"empty solve rejected");
        std::cout << "Independent energy, endpoint, LDL and failure invariants passed.\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}

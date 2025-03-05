#include "HeadVMD.hpp"

template<d_types SignalDT>
VMDDecomposition VMD(std::vector<SignalDT> f, double alpha_var, double tau, int k, bool DC, int init, double tol) {
    const F<SignalDT> f_hat_plus(f);
    int T = f_hat_plus.getT();
    Fregs fregs(T);
    NumbIter N;
    alpha Alpha(k, alpha_var);

    Omega omega_plus(N(), k, init, f_hat_plus.getFs(), DC);
    Lamb lambda_hat(N(), fregs());
    Epsilon uDiff(tol);
    int n = 0;
    Sum_UK<SignalDT> sum_uk(T);
    Ud<SignalDT> u_hat_plus(N(), k, fregs());

    while (uDiff() > tol && n < N() - 1) {
        sum_uk.update_sum_first(T, k, n, u_hat_plus);
        u_hat_plus.update_uhat(T, 0, n, f_hat_plus.getF(), sum_uk, lambda_hat, Alpha, fregs, omega_plus);

        if (!DC) {
            omega_plus.update_omega_plus(T, 0, n, u_hat_plus, fregs);
        }

        for (int i = 1; i < k; ++i) {
            sum_uk.update_sum_uk(T, i, n, u_hat_plus);
            u_hat_plus.update_uhat(T, i, n, f_hat_plus.getF(), sum_uk, lambda_hat, Alpha, fregs, omega_plus);
            omega_plus.update_omega_plus(T, i, n, u_hat_plus, fregs);
        }
        if (!u_hat_plus.isFinite() || !omega_plus.isFinite())
        {
            break;
        }

        lambda_hat.update_lamb(T, n, tau, u_hat_plus, f_hat_plus.getF());
        n = n + 1;
        uDiff.update_udiff(T, k, n, u_hat_plus);
    }

    //Niter = std::min(Niter, n);
    N.SetMinToN(n);

    omega_plus.count_omega(N(), k);

    u_hat_plus.u(N(), T, k);

    return new VMDDecomposition<SignalDT>(u_hat_plus.getU_hat(), u_hat_plus.getU(), omega_plus());
}

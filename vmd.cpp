#include "HeadVMD.hpp"

template<typename SignalDT>
VMDDecomposition VMD(std::vector<SignalDT> f, double alpha_var, double tau, int k, bool DC, int init, double tol) {
    F f_hat_plus<SignalDT>(f);
    int T = f_hat_plus.getT();
    Fregs fregs(T);
    NumbIter N;
    alpha Alpha(k, alpha_var);
    Omega omega_plus(N(), k, init, f_hat_plus.getFs(), DC);
    Lamb lambda_hat(N(), fregs());
    Epsilon uDiff(tol);

    int n = 0;
    Sum_UK sum_uk(T);
    Ud u_hat_plus(N(), k, fregs());

    while (uDiff() > tol && n < N() - 1) {
        sum_uk.update_sum_first(T, k, n, u_hat_plus);
        u_hat_plus.update_uhat(T, 0, n, f_hat_plus.getF(), sum_uk, lambda_hat, Alpha, fregs, omega_plus);

        if (!DC) {
            omega_plus.update_omega_plus(T, 0, n, u_hat_plus, fregs);
        }

        for (int i = 1; i < k; i++) {
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

    return new VMDDecomposition(u_hat_plus.getU_hat(), u_hat_plus.getU(), omega_plus());
}

# Copyright 2022, Michael Schneider, All rights reserved.
#!/usr/bin/env python3
"""
Test Suite for Variational Inference Module
"""
import unittest
import numpy as np
import scAbsolute
import sys
import os


class TestModule(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestModule, self).__init__(*args, **kwargs)
        self.dir = os.path.dirname(os.path.abspath(__file__))

    def build_dataset(self, N, K, D, prob, zero_prob, dist):
        """Build a mixture dataset for testing purposes

           Args:
                N (int): number of elements
                K (int): number of mixture components
                D (int): dimensionality of each component
                prob (array-like), shape
                zero_prob (float): dropout probability
                dist: distribution to draw samples from

           Returns:
                X: The return value. True for success, False otherwise.
                label:
        """
        label = np.random.choice(range(K), size=N, replace=True, p=prob)
        assert(prob.shape[0] == K)

        X = np.zeros((N, D))
        for i in range(N):
            new_element = dist()
            assert(new_element.shape == (D, K))
            X[i, :] = new_element[:, label[i]]

        # dropout events
        if zero_prob > 0:
            idx_zero = np.random.choice(range(N), size=int(zero_prob*N),
                                        replace=False)
            for i in idx_zero:
                X[i, :] = np.zeros((D, 1))

        return X, label

    def test_dpgmm(self):
        """Test traditional DPGMM model"""
        np.random.seed(11)

        N = 50000  # Number of data points
        K = 5  # Number of components
        D = 1  # dimensionality of data
        prob = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        mu = np.array([2, 3, 4, 6, 9])
        sigma = np.array([0.1, 0.05, 0.05, 0.1, 0.1])
        f_gen = lambda: np.random.normal(size=(D, K), loc=mu, scale=sigma)
        zero_prob = 0

        M = 3
        data = np.zeros((M, N, D))
        L = np.zeros((M, N))
        for j in range(M):
            X, label = self.build_dataset(N, K, D, prob, zero_prob, f_gen)
            data[j, :, :] = X
            L[j, :] = label

        data[1, :, :] = np.copy(data[0, :, :])
        assert(np.all(data[1, :, :] == data[0, :, :]))
        T = 64
        model = scAbsolute.dpgmm(T = T, means_prior = np.arange(1,T+1)[:, np.newaxis],
                            verbose=0, verbose_interval=50)

        vlb, inf_params = model.fit(data, batch_size=4, n_iterations=1, n_steps=20, n_cores=2)
        vlbs, delta_vlbs = vlb

        peaks = np.array([model.means[0][i] for i in range(model.means.shape[1]) if model.weights[0][i] > 1e-2])
        dif = np.min(np.abs(mu - peaks), axis=1)
        np.testing.assert_array_less(np.sort(dif), 0.1 * np.ones(peaks.shape[0]))

    #  def test_realdata(self):
        #  np.random.seed(11)

        #  print("DEBUG")
        #  data = np.loadtxt(os.path.expandvars("$HOME/scAbsolute/example-data/debug_unnormalized.txt"))
        #  data = data[:, :, np.newaxis]
        #  print(data.shape)
        #  print(np.min(data, axis=1))
        #  print(np.max(data, axis=1))
        #  #  plt.hist(data[0, :, :], bins=100)
        #  #  plt.hist(Y[0], bins=1000)
        #  #  plt.show()
        #  #  sys.exit(2)
        #  T = 128
        #  model = scAbsolute.dpgmm(T = T, means_prior = np.arange(1,T+1)[:, np.newaxis],
                                 #  verbose=0, verbose_interval=50)
        #  vlb, inf_params = model.fit(data, batch_size=4, n_iterations=1, n_steps=3, n_cores=4)
        #  vlbs, delta_vlbs = vlb
        #  result, solution = scAbsolute.run_scAbsolute(data, T=T, batch_size=4, n_steps=3, n_init=1)

        #  return model, vlbs, delta_vlbs, data, result, solution


    def test_scrun(self):
        np.random.seed(11)

        N = 20000  # Number of data points
        K = 5  # Number of components
        D = 1  # dimensionality of data
        prob = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        mu = np.array([1.05, 2.1, 2.95, 5.03, 7.95])
        sigma = np.array([0.05, 0.15, 0.2, 0.1, 0.3])
        scale = 2.79
        f_gen = lambda: np.random.normal(size=(D, K), loc=mu*scale, scale=sigma)
        zero_prob = 0.0

        M = 3
        data = np.zeros((M, N, D))
        L = np.zeros((M, N))
        for j in range(M):
            X, label = self.build_dataset(N, K, D, prob, zero_prob, f_gen)
            data[j, :, :] = X
            L[j, :] = label

        data[1, :, :] = np.copy(data[0, :, :])
        assert(np.all(data[1, :, :] == data[0, :, :]))

        T = 128
        data = data[0, :, :]
        data = data[np.newaxis, :, :]
        #  data = data + np.random.normal(0, 0.3, data.shape)
        #  print(data.shape)
        result, solution = scAbsolute.run_scAbsolute(data, T=T, batch_size=4, n_steps=20, verbose=2, n_init=2)

        means = solution["means"]
        weights = solution["weights"]
        prescale = solution["prescale"]
        peaks = np.array([means[0][i] for i in range(means.shape[1]) if weights[0][i] > 0.1])
        mu = mu *  scale * prescale

        dif = np.abs(mu - peaks)
        np.testing.assert_array_less(np.sort(dif), 2*np.ones(peaks.shape[0]))

    def test_elbo(self):
        # modified from https://github.com/mcusi/tf_dpgmm
        import scipy.stats
        from scipy.special import digamma, gamma
        verbose = True

        np.random.seed(2018)
        D = 1
        K = 4

        nu = np.random.randn(K)
        zeta = np.array([0.2, 0.2, 0.1, 0.5])

        a = np.ones(K)
        b = 1.5*np.ones(K)

        lambda1 = np.ones(K)
        lambda2 = 2.*np.ones(K)

        alpha = 3.
        p_phi = scipy.stats.beta(1, alpha)
        p_mu = scipy.stats.norm
        p_Lambda = scipy.stats.gamma(a=1., scale=1.)
        def log_p_z(phi, z):
            p = np.concatenate([[1], np.cumprod(1-phi[:-1])]) * phi
            return np.log(p[z])
        def p_x(z, mu, Lambda): return scipy.stats.norm(loc=mu[z], scale=np.sqrt(1./Lambda[z]))

        q_phi = scipy.stats.beta(lambda1, lambda2)
        q_mu = scipy.stats.norm(loc=nu)
        q_Lambda = scipy.stats.gamma(a=a, scale=1./b) #Lambda is precision!
        q_z = scipy.stats.rv_discrete(values=(range(K), zeta))

        x = 3.
        N = 40001

        ### ### ### ### ### ### ### ### ### ###
        # test derivations of variational lower bound ##

        # V ##
        if verbose:
            print("V")
        # Analytical
        lgamma = lambda x: np.log(gamma(x))
        bound = K*(lgamma(1. + alpha) - lgamma(alpha)) + \
                    sum(((alpha - 1.)*(digamma(l2_k) - digamma(l1_k + l2_k))
                    - lgamma(l1_k + l2_k) + lgamma(l1_k) + lgamma(l2_k)
                    - (l1_k - 1.)*(digamma(l1_k) - digamma(l1_k + l2_k))
                    - (l2_k - 1.)*(digamma(l2_k) - digamma(l1_k + l2_k)))
                    for (l1_k, l2_k) in zip(lambda1, lambda2))

        # Monte Carlo
        bounds = []
        for i in range(N):
            phi = q_phi.rvs()

            bounds.append(sum(p_phi.logpdf(phi) - q_phi.logpdf(phi))) #Sum over K for MC estimate
            if verbose and i%5000 == 0: print(i, "\t", np.mean(bounds))

        if verbose:
            print("V: - Analytical term:\t %.5f \n   - MC-estimate:\t %.5f" % (bound, np.mean(bounds)))
        np.testing.assert_array_almost_equal(np.mean(bounds), bound, decimal=1, err_msg="V term")

        # mu ##
        if verbose:
            print("mu")
        # Analytical
        # we test with m_{0,k} = 0
        bound = sum(-0.5*nu_k**2 for nu_k in nu)

        # Monte Carlo
        bounds = []
        for i in range(N):
            mu = q_mu.rvs()

            bounds.append(sum(p_mu.logpdf(mu) - q_mu.logpdf(mu))) #Sum over K for MC estimate
            if verbose and i%5000 == 0: print(i, "\t", np.mean(bounds))

        if verbose:
            print("mu: - Analytical term:\t %.5f \n    - MC-estimate:\t %.5f" % (bound, np.mean(bounds)))
        np.testing.assert_array_almost_equal(np.mean(bounds), bound, decimal=1, err_msg="mu term")

        # Lambda ##
        if verbose:
            print("Lambda")
        # Analytical
        bound = sum(lgamma(a_k) - (a_k-1.)*digamma(a_k) - np.log(b_k) + a_k - np.divide(a_k,b_k)
                    for (a_k, b_k) in zip(a, b))

        # Monte Carlo
        bounds = []
        for i in range(N):
            Lambda = q_Lambda.rvs()

            bounds.append(sum(p_Lambda.logpdf(Lambda) - q_Lambda.logpdf(Lambda))) #Sum over K for MC estimate
            if verbose and i%5000 == 0: print(i, "\t", np.mean(bounds))

        if verbose:
            print("Lambda: - Analytical term:\t %.5f \n        - MC-estimate:\t %.5f" % (bound, np.mean(bounds)))
        np.testing.assert_array_almost_equal(np.mean(bounds), bound, decimal=1, err_msg="Lambda term")

        # Z ##
        if verbose:
            print("Z")
        # Analytical
        bound = sum(zeta_k*(
                        - np.log(zeta_k)
                        + digamma(l1_k) - digamma(l1_k+l2_k)
                        + sum(digamma(lambda2[j]) - digamma(lambda1[j]+lambda2[j]) for j in range(k)))
                    for (l1_k, l2_k, zeta_k, k) in zip(lambda1, lambda2, zeta, range(K)))

        # Monte Carlo
        bounds = []
        for i in range(N):
            phi = q_phi.rvs()
            z = q_z.rvs()

            bounds.append(log_p_z(phi, z) - q_z.logpmf(z)) #There's only a single datapoint, so no need for sum
            if verbose and i%5000 == 0: print(i, "\t", np.mean(bounds))

        if verbose:
            print("Z: - Analytical term:\t %.5f \n   - MC-estimate:\t %.5f" % (bound, np.mean(bounds)))
        np.testing.assert_array_almost_equal(np.mean(bounds), bound, decimal=1, err_msg="Z term")

        # X ##
        if verbose:
            print("X")
        # Analytical
        bound = sum(zeta_k * (
                    -1/2. * (np.log(2 * np.pi) - digamma(ak) + np.log(bk))
                    - (82519129./65840739) * (ak/float(bk))
                     * (2*np.pi)**(-1/2.) * ((nu_k - x)**2 + 1)
                ) for (ak, bk, nu_k, zeta_k) in zip(a, b, nu, zeta))
        #  print("DEBUG", bound)
        bound = sum(zeta_k * (
                    -1/2. * (np.log(2 * np.pi) - digamma(ak) + np.log(bk)
                    + (ak/float(bk))
                     * ((nu_k - x)**2 + 1))
                ) for (ak, bk, nu_k, zeta_k) in zip(a, b, nu, zeta))
        #  print("DEBUG", bound)

        # Monte Carlo
        bounds = []
        for i in range(N):
            mu = q_mu.rvs()
            Lambda = q_Lambda.rvs()
            z = q_z.rvs()

            bounds.append(p_x(z, mu, Lambda).logpdf(x)) #There's only a single datapoint, so no need for sum
            if verbose and i%5000 == 0: print(i, "\t", np.mean(bounds))

        if verbose:
            print("X: - Analytical term:\t %.5f \n   - MC-estimate:\t %.5f" % (bound, np.mean(bounds)))
        np.testing.assert_array_almost_equal(np.mean(bounds), bound, decimal=1, err_msg="X term")


if __name__ == "__main__":
    unittest.main()

    #  test = TestModule()
    #  test.test_dpgmm()
    #  test.test_scrun()

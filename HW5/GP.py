import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import cholesky, det, lstsq, inv
from scipy.optimize import minimize

beta=5
beta_inv=(1/beta)

####Rational quadratic
def kernel(xn, xm, l, sigma, a):
    xn_minus_xm_aquare = (xn**2).reshape(-1, 1) + np.sum(xm**2, 1) - 2 * np.dot(xn, xm.T)
    return sigma**2*(1+xn_minus_xm_aquare /(2*a*l**2))**(a*-1)

####posterior 
def posterior(x_star, x_train, y_train, l, sigma, a):
    C = kernel(x_train, x_train, l, sigma, a) + beta_inv * np.eye(len(x_train))
    K_star = kernel(x_train, x_star, l, sigma, a)
    K_starstar = kernel(x_star, x_star, l, sigma, a) + beta_inv * np.eye(len(x_star))
    C_inv = inv(C)
    
    mu_s = K_star.T.dot(C_inv).dot(y_train)
    cov_s = K_starstar - K_star.T.dot(C_inv).dot(K_star)
    
    return mu_s, cov_s

####minus log likelihood
def mll(x_train, y_train):
    y_train = y_train.ravel()
    
    def mll_result(theta):
        print(theta)
        K = kernel(x_train, x_train, l=theta[0], sigma=theta[1],a=theta[2]) + beta_inv * np.eye(len(x_train))
        return 0.5*np.log(det(K)) + 0.5*y_train.dot(inv(K).dot(y_train)) + 0.5*len(x_train)*np.log(2*np.pi)

    return mll_result

####plot function
def plot(mu, cov, x_axi, x_train, y_train):
    x_axi = x_axi.ravel()
    mu = mu.ravel()
    #### 95% confidence
    uncertainty = 1.96 * np.sqrt(np.diag(cov))

    plt.fill_between(x_axi, mu + uncertainty, mu - uncertainty, alpha=0.1)
    plt.plot(x_axi, mu, label='Mean')
    plt.plot(x_train, y_train, 'rx')
    plt.legend()
    plt.show()

####Load data
x_data = []
y_data = []
file = open("input.data")
for line in file:
    tmp = line.split(' ')
    x_data.append(float(tmp[0]))
    y_data.append(float(tmp[1].split('\n')[0]))
x_train = np.array(x_data).reshape(-1, 1)
y_train = np.array(y_data).reshape(-1, 1)

####x-axi
x_axi = np.arange(-50, 50, 0.2).reshape(-1, 1)

####before optimize
l_init=1
sigma_init=1
a_init=1
mu, cov = posterior(x_axi, x_train, y_train, l_init, sigma_init, a_init)
plot(mu, cov, x_axi, x_train, y_train)

####optimize
minimal = minimize(mll(x_train, y_train), [l_init, sigma_init, a_init], bounds=((0.00001, None), (0.00001, None), (0.00001, None)), method='L-BFGS-B')
l_opt, sigma_opt, a_opt = minimal.x

####after optimize
mu, cov = posterior(x_axi, x_train, y_train, l_opt, sigma_opt, a_opt)
plot(mu, cov, x_axi, x_train, y_train)

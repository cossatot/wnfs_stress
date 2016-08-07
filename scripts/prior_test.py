import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def make_pdf(vals, probs, n_interp=1000):

    val_min = np.min(vals)
    val_max = np.max(vals)
    
    # if the PDF is just a point (no uncertainty)
    if val_min == val_max: 
        return val_min, 1.
    
    # if not...
    else:
        pdf_range = np.linspace(val_min, val_max, n_interp)

        pmf = interp1d(vals, probs, bounds_error=False, fill_value=0.)
        pmf_samples = pmf(pdf_range)
        pdf_probs = pmf_samples / np.sum(pmf_samples) # normalize

    return pdf_range, pdf_probs


def make_cdf(pdf_range, pdf_probs):
    return (pdf_range, np.cumsum(pdf_probs))


def inverse_transform_sample(vals, probs, n_samps, n_interp=1000):
    
    pdf_range, pdf_probs = make_pdf(vals, probs, n_interp)
    cdf_range, cdf_probs = make_cdf(pdf_range, pdf_probs)

    if len(cdf_probs) == 1:
        return np.ones(n_samps) * pdf_range
    
    else:
        cdf_interp = interp1d(cdf_probs, cdf_range, bounds_error=False,
                              fill_value=0.)
        samps = np.random.rand(n_samps)

        return cdf_interp(samps)


#s3_xs = np.array([-1, 0, 1])
#s3_ys = np.array([1, 0, 1])

n_samps = int(5e5)

s3_xs = np.linspace(0.000001, 0.9999, 10000)
s3_ys = np.log(1/(1-s3_xs**500))
#s3_ys = 1 / np.exp(-s3_xs)



s1s = np.random.uniform(size=n_samps)
s3ps = inverse_transform_sample(s3_xs, s3_ys, n_samps)

#s3ps = np.random.uniform(size=n_samps)

s3s = s3ps * s1s

#s3s[::2] *= -1

plt.figure(1)
plt.plot(s3_xs, s3_ys)

plt.figure(2)
plt.hist(s3s, bins=100, normed=True)

s1z = np.random.uniform(size=n_samps)
s3z = np.random.uniform(size=n_samps) * -1

s1s = np.hstack((s1s, s1z))
s3s = np.hstack((s3s, s3z))

plt.figure(3)

plt.subplot(131)
plt.hist(s1s, bins=100, normed=True)
plt.xlabel('s1')

plt.subplot(132)
plt.hist(s3s, bins=100, normed=True)
plt.xlabel('s3')

plt.subplot(133)
plt.plot(s1s, s3s, ',')

plt.show()


import numpy as np
import matplotlib.pyplot as plt

# Read the csvs
linearResult_linear = np.loadtxt("linear_uniform.csv", delimiter=",")[:, 1]

plt.plot(1/linearResult_linear, label="Linear Uniform")

plt.xlabel("Number of Samples for the MC estimator")
plt.ylabel("Empirical variance of the MC estimator")
plt.title("Monte Carlo Integration Results")
plt.legend()
plt.grid(True)
plt.show()
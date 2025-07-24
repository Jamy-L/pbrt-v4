#include <cmath>
#include <iostream>
#include <functional>
#include <random>
#include <fstream>
#include <vector>

class RNG {
public:
    virtual float sample() = 0;
    virtual ~RNG() = default;
    virtual double pdf(double x) = 0;
};

class UniformRNG : public RNG {
    float xmin, xmax;
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist;
public:
    UniformRNG(float xmin, float xmax)
      : xmin(xmin), xmax(xmax), gen(std::random_device{}()), dist(xmin, xmax) {}
    float sample() override { return dist(gen); }
    double pdf(double x) override {
        return (x >= xmin && x <= xmax) ? inv_range : 0.0f;
    }
    private:
    double inv_range = 1.0f / (xmax - xmin);
};

class NormalRNG : public RNG {
    std::mt19937 gen;
    std::normal_distribution<float> dist;
public:
    NormalRNG(float mean, float stddev)
      : gen(std::random_device{}()), dist(mean, stddev) {}
    float sample() override { return dist(gen); }
    double pdf(double x) override {
        double exponent = -0.5 * std::pow((x - dist.mean()) * inv_stddev, 2);
        return coeff * std::exp(exponent);
    }
private:
    double coeff = 1.0 / (std::sqrt(2 * M_PI) * dist.stddev());
    double inv_stddev = 1.0 / dist.stddev();
};


std::vector<double> monteCarloIntegrator(std::function<float(float)> func,
                            RNG& rng, float xmin, float xmax, int n_samples = 100'000,
                        bool verbose = false) {
    double sum = 0.0f;
        std::vector<double> estimates;
        float x;
    for (int i = 0; i < n_samples; ++i) {
        x = rng.sample();
        sum += ((xmin < x) && (x < xmax) ? func(x) : 0.0f) / (rng.pdf(x) * n_samples);
        estimates.push_back(sum * n_samples / (i+1)); // Division is made here to avoid a stupid float overflow

        if (verbose && i % (n_samples / 10) == 0) {  // every 10%
            std::cout << "Progress: " << (100 * i / n_samples) << "%\r" << std::flush;
        }
    }
    if (verbose) std::cout << "Progress: 100%\n";
    return estimates;
}

//----- A few basc functions to test the Monte Carlo integration -----
double linearFunction(double x) {
    return x;
}

float quadraticFunction(float x) {
    return x * x;
}

float cubicFunction(float x) {
    return x * x * x;
}

#include <fstream>
#include <iomanip>
#include <limits>

void saveEstimates(const std::vector<double>& estimates, const std::string& filename) {
    std::ofstream out(filename);
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;
    for (size_t i = 0; i < estimates.size(); ++i) {
        out << i << "," << estimates[i] << "\n";
    }
}

std::vector<double> estimateVariance(std::function<std::vector<double>()> func, double gt_value, uint n_samples, bool verbose = true) {
    std::vector<double> estimates;
    // initialize with the first sample
    auto sample = func();
    for (size_t i = 0; i < sample.size(); ++i) {
        double error = (sample[i] - gt_value) * (sample[i] - gt_value);
        estimates.push_back(error);
    }

    for (uint iter = 1; iter < n_samples; ++iter) {
        sample = func();
        for (size_t i = 0; i < sample.size(); ++i) {
            double error = (sample[i] - gt_value) * (sample[i] - gt_value);
            estimates[i] += error;
        }

        if (verbose && iter % (n_samples / 10) == 0) {  // every 10%
            std::cout << "Progress: " << (100 * iter / n_samples) << "%\r" << std::flush;
        }
    }
    
    if (verbose) std::cout << "Progress: 100%\n";


    // Normalize the estimates
    for (size_t i = 0; i < estimates.size(); ++i) {
        estimates[i] /= n_samples - 1;
    }
    return estimates;

}
int main() {
    auto uniformRNG = UniformRNG(0.0f, 1.0f);
    auto normalRNG = NormalRNG(1.0f, 0.7f);
    constexpr int N_trials_var = 10'000;
    constexpr float linearResult_gt = 0.5f;

    // Estimate variance for linear function
    auto acc_linearResult_linear = estimateVariance(
        [&]() { return monteCarloIntegrator(linearFunction, uniformRNG, 0.0f, 1.0f); },
        linearResult_gt,
        N_trials_var);
        
    saveEstimates(acc_linearResult_linear, "linear_uniform.csv");

    std::cout << "Monte Carlo integration completed and results saved." << std::endl;
    return 0;
}

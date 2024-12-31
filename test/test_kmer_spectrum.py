import pytest
from dedup.kmer_spectrum import gauss, bimodal, find_minimum_between_peaks, fit_kmer_spectrum

def test_gauss():
    mu = 0
    sigma = 1
    A = 1
    x = 0
    expected = 1  # Calculated value
    assert pytest.approx(gauss(x, mu, sigma, A), 0.0001) == expected

def test_bimodal():
    mu1 = 1
    sigma1 = 1
    A1 = 1
    sigma2 = 1
    A2 = 1
    x = 2
    expected = 1.6065306597126334  # Calculated value
    assert pytest.approx(bimodal(x, mu1, sigma1, A1, sigma2, A2), 0.0001) == expected

def test_find_minimum_between_peaks():
    mu1 = 1
    sigma1 = 1
    A1 = 1
    sigma2 = 1
    A2 = 1
    min_x = find_minimum_between_peaks(mu1, sigma1, A1, sigma2, A2)
    assert 1.0 < min_x < 2.0  # The minimum should be between the two peaks

def test_fit_kmer_spectrum():
    # Sample histogram data with two peaks
    data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 20, 50, 100, 50, 20, 10, 0, 0, 0, 0, 0, 10, 20, 50, 100, 50, 20, 10, 0, 0, 0, 0, 0, 0]
    min_kmer_depth = 5
    max_kmer_depth = 50
    lower_bound, upper_bound = fit_kmer_spectrum(data, min_kmer_depth, max_kmer_depth)
    print(lower_bound, upper_bound)
    assert 15 <= lower_bound <= 25
    assert 25 <= upper_bound <= 35

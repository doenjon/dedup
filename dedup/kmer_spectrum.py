import os
import sys
import subprocess
import logging
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar, differential_evolution


logger = logging.getLogger("dedup_logger")

def get_homozygous_kmer_range(kmer_db, tmp_dir, min_kmer_depth, max_kmer_depth):
    '''
    Get the range of kmer frequencies that are homozygous

    Args:
        kmer_db (str): path to kmc kmer database

    Returns:
        tuple: (min, max) kmer frequency
    '''

    # get kmer histogram
    kmer_histo_data = get_kmer_histogram_data(kmer_db, tmp_dir)

    # Fit model to kmer spectrum
    lower_bound, upper_bound = fit_kmer_spectrum(kmer_histo_data, min_kmer_depth, max_kmer_depth)

    logger.info(f"Set homozygous kmer range to ({lower_bound}, {upper_bound})")
    return (lower_bound, upper_bound)


def get_kmer_histogram_data(kmer_db, tmp_dir):
    """
    Retrieves the k-mer histogram data from a given k-mer database.

    Args:
        kmer_db (str): The path to the kmc k-mer database.

    Returns:
        list: The k-mer histogram data.

    Raises:
        FileNotFoundError: If the k-mer histogram file does not exist.

    """
    histo_file = os.path.join(tmp_dir, 'kmer_counts.histo')
    # cmd = f"kmc_tools transform {kmer_db} -ci{min_kmer_depth} -cx{max_kmer_depth} histogram {histo_file}"
    cmd = f"kmc_tools transform {kmer_db} histogram {histo_file}"
    logger.info(cmd)
        
    if not os.path.exists(histo_file):
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()
        if retval:
            logger.critical(f"make_kmer_db ret: {retval}")
            sys.exit(retval)
    else:
        logger.debug(f"Skipping because results already exist")

    data = []
    with open(histo_file, 'r') as f:
        for line in f:
            x, y = line.strip().split()
            data.append(int(y))

    return data


# Gaussian function
def gauss(x, mu, sigma, A):
    """
    Implements gaussian function

    Parameters:
    x (float): The input value.
    mu (float): The mean of the Gaussian distribution.
    sigma (float): The standard deviation of the Gaussian distribution.
    A (float): The amplitude of the Gaussian function.

    Returns:
    float: The value of the Gaussian function at the given point.
    """
    return A / np.sqrt(sigma) * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))
    

# Mixture of Gaussians
def bimodal(x, mu1, sigma1, A1, sigma2, A2):
    """
    Implements mixutre of two gaussian functions with constraint that second peak is exactly twice the first.
    Allows seperate amplitude and standard deviation for each peak.

    Parameters:
    x (float): The input value.
    mu1 (float): Mean of the first peak.
    sigma1 (float): Standard deviation of the first peak.
    A1 (float): Amplitude of the first peak.
    sigma2 (float): Standard deviation of the second peak.
    A2 (float): Amplitude of the second peak.

    Returns:
    float: The value of the bimodal distribution function at x.
    """
    mu2 = 2*mu1 # second peak is exactly twice the first
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

def find_minimum_between_peaks(mu1, sigma1, A1, sigma2, A2):
    """
    Finds the minimum point between two Gaussian peaks.
    
    Parameters:
        mu1, sigma1, A1: Parameters for the first Gaussian peak.
        sigma2, A2: Parameters for the second Gaussian peak, with mu2 = 2*mu1.

    Returns:
        The x value where the bimodal function reaches its minimum between the two peaks.
    """
    mu2 = 2 * mu1  # second peak is exactly twice the first

    def objective(x):
        return bimodal(x, mu1, sigma1, A1, sigma2, A2)

    # Assuming the minimum is between mu1 and mu2
    result = minimize_scalar(objective, bounds=(mu1, mu2), method='bounded')
    min_x = result.x
    
    return min_x

def loss_function(params, *args):
    """
    Objective function to compute the sum of squared residuals
    between the model predictions and the observed data.

    Parameters:
        params (tuple): The parameters for the bimodal function. Values are (mu1, sigma1, A1, sigma2, A2).
        args (tuple): The x and y data to fit the curve to.

    Returns:
        float: The sum of squared residuals between the model predictions and the observed data.
    """
    x, y = args
    residuals = y - bimodal(x, *params)
    return np.sum(residuals**2)

def fit_kmer_spectrum(data, min_kmer_depth, max_kmer_depth):
    """
    Fits a mixture of two Gaussian curves to the given data and returns the 
    mean and standard deviation of the homogeneous peak.

    Parameters:
        data (array-like): The histogram data to fit the curve to.
        init_params (tuple): The initial parameters for the fitting. Values are (mu1, sigma1, A1, sigma2, A2).
        min_kmer_depth (int): The minimum depth of kmers to consider. Values below are set to zero.
        max_kmer_depth (int): The maximum depth of kmers to consider. Values above are removed.
    
    Returns:
        tuple: bounds for homozygous peak

    """
    # Normalize and remove data above and below bounds
    data = data[:max_kmer_depth]
    
    # Zero out error kmers below min_kmer_depth
    for i in range(min_kmer_depth):
        data[i] = 0 
    
    # Normalize the data to make fitting easier
    data = np.array(data) / np.max(data)

    # add x data
    x = np.arange(len(data))
    y = data

    bounds = [(min_kmer_depth, max_kmer_depth),  # mu1
              (1e-6, 100),  # sigma1
              (1e-6, 10),  # A1
              (1e-6, 100),  # sigma2
              (1e-6, 10)]  # A2


    # Perform global optimization
    result = differential_evolution(loss_function, bounds, args=(x, y))

    if not result.success:
        logger.error(f"Optimization failed: {result.message}")
        logger.error(f"Consider providing homozygous_lower_bound and homozygous_upper_bound manually.")
        raise RuntimeError(f"Optimization failed: {result.message}")

    warning_threshold = 1  
    if result.fun > warning_threshold:
        logger.warning(f"Optimizer may not have found a good model of kmer-spectrum -- suggest manually checking curve fit (kmer_spectrum_fit.png)")
        logger.warning(f"Consider providing homozygous_lower_bound and homozygous_upper_bound manually.")

    params = result.x
    
    # Graph the data and fit to check quality
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 6))
    sns.barplot(x=np.arange(len(data)), y=data, color="skyblue")

    plt.title('Kmer histogram')
    plt.xlabel('Kmer depth')
    plt.ylabel('Relative Frequency')

    # Add the fitted Gaussian curve
    x_vals = np.linspace(0, len(data), 1000)
    fitted_curve = bimodal(x_vals, params[0], params[1], params[2], params[3], params[4])
    plt.plot(x_vals, fitted_curve, color='red', label='Fitted Gaussian Curve')

    # Set the x label to show only 1 in 10 labels
    plt.xticks(np.arange(0, len(data), 10))

    plt.legend()

    output_file = 'kmer_spectrum_fit.png'  
    plt.savefig(output_file)
    
    mean_het, std_het, _, std_homo, _ = params
    mean_homo = 2*mean_het # enforced in fitting

    # Set bounds at minimum between peaks and 2 standard deviations above homozygous mean
    homo_left_bound = round(find_minimum_between_peaks(*params))
    homo_right_bound = round(mean_homo + 2*std_homo)

    return homo_left_bound, homo_right_bound

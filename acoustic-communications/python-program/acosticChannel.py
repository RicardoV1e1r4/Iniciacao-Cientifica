# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:06:29 2026

@author: Ricardo Alexandre
"""

"""
Abstract:
    This function generates the impulse response of an underwater acoustic
    (UWA) channel, taking into account the Doppler effect.

    This function uses the UWA channel model presented in [1] and the
    methodology presented in [2], which is based on the UWA channel models
    described in [1].

References:
    [1] S. Zhou and Z. Wang, "OFDM for Underwater Acoustic Communications."
        Chichester: John Wiley & Sons, May 2014.

    [2] R. S. Chaves, "Modeling and Simulation of Underwater Acoustic
        Communication Systems (in Portuguese)", Graduation Thesis,
        Federal University of Rio de Janeiro, Rio de Janeiro, 2016.
"""

import numpy as np
from scipy.stats import expon, rayleigh
from fractions import Fraction


def acoustic_channel(setup, gain, delay, doppler):
    """
    Generate the impulse response of an underwater acoustic (UWA) channel.

    Parameters
    ----------
    setup : dict
        Dictionary containing the channel parameters.

        setup["Ts"]
            Sampling period of the UWA channel, in seconds.

        setup["paths"]
            Number of multipath components. Must be a positive integer.

        setup["delayspread"]
            Delay spread of the UWA channel.

    gain : dict
        Dictionary containing the gain parameters.

        gain["attenuation"]
            Attenuation occurring during the delay spread, in dB.

    delay : dict
        Dictionary containing the delay parameters.

        delay["mean"]
            Mean value of the exponential probability density function
            used to generate the delay intervals.

    doppler : dict
        Dictionary containing the Doppler parameters.

        doppler["type"]
            Type of Doppler effect. Possible values are:

                "none"
                    No Doppler scaling factor (DSF).

                "uniform"
                    Uniform Doppler scaling factor.

                "non-uniform"
                    Non-uniform Doppler scaling factor.

        doppler["velocity"]
            Relative velocity between transmitter and receiver.

            For "uniform", this must be a scalar.

            For "non-uniform", this must be an array containing one
            velocity value for each multipath component.

    Returns
    -------
    h : numpy.ndarray
        Discrete impulse response of the UWA channel.

    delay_bar : numpy.ndarray
        Delay associated with each multipath component after Doppler
        correction.

    gain_tap : numpy.ndarray
        Gain associated with each multipath component.

    Q : int
        Downsampling factor associated with the Doppler scaling factor.

    M : int
        Upsampling factor associated with the Doppler scaling factor.
    """

    # ==============================================================
    # 1. GENERATING THE DISTRIBUTION OF DELAYS
    # ==============================================================

    # Create an exponential probability distribution for the delay
    # intervals (Delta_tau).
    #
    # MATLAB:
    # distribution_delay = makedist('exponential','mu',delay.mean);

    distribution_delay = expon(scale=delay["mean"])

    # Generate exponentially distributed delay intervals.
    #
    # MATLAB:
    # Dtau = random(distribution_delay, setup.paths, 1);

    Dtau = distribution_delay.rvs(size=setup["paths"])

    # Convert the delay intervals from continuous time to discrete-time
    # indexes.
    #
    # MATLAB:
    # Dtau_index = ceil(Dtau/setup.Ts);

    Dtau_index = np.ceil(Dtau / setup["Ts"]).astype(int)

    # Calculate the discrete-time indexes of the multipath components
    # using the cumulative sum of the delay intervals.
    #
    # MATLAB:
    # delay_index = cumsum(Dtau_index);

    delay_index = np.cumsum(Dtau_index)

    # Convert the discrete indexes back to time.
    #
    # MATLAB:
    # delay = (delay_index - 1)*setup.Ts;

    delay_values = (delay_index - 1) * setup["Ts"]

    # ==============================================================
    # 2. GENERATING THE DISTRIBUTION OF GAINS
    # ==============================================================

    # Calculate the exponential attenuation power factor.
    #
    # The attenuation is specified in dB.
    #
    # MATLAB:
    # alpha = log(10^(gain.attenuation/10))/setup.delayspread;

    alpha = (
        np.log(10 ** (gain["attenuation"] / 10))
        / setup["delayspread"]
    )

    # Calculate the gain variance associated with each multipath
    # component.
    #
    # MATLAB:
    # gain_variance = exp(-alpha*delay);

    gain_variance = np.exp(-alpha * delay_values)

    # Generate Rayleigh-distributed gains.
    #
    # MATLAB:
    # gain_tap = raylrnd(
    #     sqrt(gain_variance*2/(4-pi))
    # );

    rayleigh_scale = np.sqrt(
        gain_variance * 2 / (4 - np.pi)
    )

    # scipy.stats.rayleigh uses the parameter "scale".
    gain_tap = rayleigh.rvs(scale=rayleigh_scale)

    # Normalize the gains so that the total energy of the channel
    # taps is equal to one.
    #
    # MATLAB:
    # gain_tap = gain_tap/norm(gain_tap);

    gain_tap = gain_tap / np.linalg.norm(gain_tap)

    # ==============================================================
    # 3. DOPPLER EFFECT
    # ==============================================================

    # Approximate speed of sound in water, in meters per second.
    c = 1500.0

    # Get the Doppler type.
    doppler_type = doppler["type"]

    # --------------------------------------------------------------
    # 3.1 No Doppler
    # --------------------------------------------------------------

    if doppler_type == "none":

        # No Doppler scaling factor.
        a_max = 0.0

        # No upsampling/downsampling is required.
        Q = 1
        M = 1

    # --------------------------------------------------------------
    # 3.2 Uniform Doppler
    # --------------------------------------------------------------

    elif doppler_type == "uniform":

        # Calculate the Doppler scaling factor:
        #
        # a = v/c

        a_max = doppler["velocity"] / c

        # The MATLAB function rat() is used to approximate
        # (1 + a_max) by a rational number Q/M.
        #
        # MATLAB:
        # [Q,M] = rat(1+a_max);

        fraction = Fraction(1 + a_max).limit_denominator()

        Q = fraction.numerator
        M = fraction.denominator

    # --------------------------------------------------------------
    # 3.3 Non-uniform Doppler
    # --------------------------------------------------------------

    elif doppler_type == "non-uniform":

        # Calculate the Doppler scaling factor for each multipath.
        #
        # a = v/c

        a = np.asarray(doppler["velocity"]) / c

        # Use the maximum Doppler scaling factor.
        a_max = np.max(a)

        # Approximate (1 + a_max) by Q/M.
        fraction = Fraction(1 + a_max).limit_denominator()

        Q = fraction.numerator
        M = fraction.denominator

    # --------------------------------------------------------------
    # Invalid Doppler type
    # --------------------------------------------------------------

    else:
        raise ValueError(
            "Invalid Doppler type. "
            "Use 'none', 'uniform' or 'non-uniform'."
        )

    # ==============================================================
    # 4. CORRECTING THE DELAYS DUE TO DOPPLER
    # ==============================================================

    # Calculate the new discrete-time delay indexes after applying
    # the Doppler scaling factor.
    #
    # MATLAB:
    # delay_index_bar = ceil(delay_index./(1 + a_max));

    delay_index_bar = np.ceil(
        delay_index / (1 + a_max)
    ).astype(int)

    # Convert the corrected discrete indexes back to time.
    #
    # MATLAB:
    # delay_bar = (delay_index_bar - 1)*setup.Ts;

    delay_bar = (
        delay_index_bar - 1
    ) * setup["Ts"]

    # ==============================================================
    # 5. GENERATING THE CHANNEL IMPULSE RESPONSE
    # ==============================================================

    # Create the discrete impulse response.
    #
    # MATLAB automatically expands the vector when assigning values
    # to positions that do not yet exist. In Python, we need to
    # allocate the array explicitly.

    h = np.zeros(np.max(delay_index_bar))

    # Place each multipath gain at its corresponding delay position.
    #
    # MATLAB:
    # h(delay_index_bar) = gain_tap;
    #
    # MATLAB indexes start at 1, while Python indexes start at 0.
    # Therefore, delay_index_bar - 1 is used here.

    h[delay_index_bar - 1] = gain_tap

    # Return the channel impulse response, corrected delays,
    # multipath gains, and Doppler resampling factors.
    return h, delay_bar, gain_tap, Q, M

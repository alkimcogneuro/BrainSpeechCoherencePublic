function x = SimulateColoredNoise(L, fs, PSDfun)
    % ===============================================================================================================
    % Real-valued Gaussian noise of length L whose EXPECTED one-sided power spectral density is PSDfun(f).
    % ===============================================================================================================
    % Used by SpectrumVsPowerLaw.m to build its Monte-Carlo null (PSDfun = a fitted power law), and by
    % testscript_SpectrumVsPowerLaw_synthetic.m to make test signals with known spectra (power law, power law + peak,
    % knee, ...).
    %
    %   L      : number of samples (any positive integer, even or odd)
    %   fs     : sampling rate, Hz
    %   PSDfun : function handle, f (column vector of frequencies in Hz, > 0) -> one-sided PSD in units^2/Hz.
    %            Must be vectorized. Example: @(f) 0.1 * f.^-1.2
    %
    % Normalization is the inverse of CrossSpectralDensity.m's: that function returns
    %   PSD_one_sided(f_k) = 2 |X_k|^2 / (fs * sum(w.^2))   (0 < k < L/2; the DC and Nyquist bins are not doubled),
    % and for white noise of variance s2, E|X_k|^2 = s2 * sum(w.^2), so the one-sided PSD is 2*s2/fs. Here the DFT
    % coefficients are drawn directly with E|X_k|^2 = L * fs * S1(f_k) / 2 (0 < k < L/2), and L * fs * S1(f_N) at the
    % Nyquist bin (real), with X_0 = 0 (no DC), then inverse-transformed. Sanity check (done in the synthetic test
    % script): Welch PSDs of the output average back to PSDfun(f).
    %
    % The noise is circular (periodic) over the L samples; Welch windows much shorter than L see no wrap-around effect.
    % ===============================================================================================================
    n_bins = floor(L/2) + 1;
    f = (0:n_bins-1)' * (fs / L);
    S = zeros(n_bins, 1);
    S(2:end) = PSDfun(f(2:end));
    if any(S < 0) || any(~isfinite(S))
        error('SimulateColoredNoise:BadPSD', 'PSDfun must return finite, non-negative values for f > 0.');
    end
    Z = (randn(n_bins, 1) + 1i*randn(n_bins, 1)) .* sqrt(L * fs * S / 4);
    Z(1) = 0;
    if mod(L, 2) == 0
        Z(end) = sqrt(L * fs * S(end)) * randn;                 % Nyquist bin is real
        full_spectrum = [Z; conj(Z(end-1:-1:2))];
    else
        full_spectrum = [Z; conj(Z(end:-1:2))];
    end
    x = real(ifft(full_spectrum));
end

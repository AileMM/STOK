 %越接近于0，表明两个任务越相似，这里需要定义一个阈值
function KL_div = kl_divergence_gaussian(mu_P, Sigma_P, mu_Q, Sigma_Q)
    % mu_P, mu_Q: Mean vectors (1xD)
    % Sigma_P, Sigma_Q: Covariance matrices (DxD)

    d = length(mu_P);  % Dimensionality of the decision space

    % Ensure matrices are symmetric
    Sigma_P = (Sigma_P + Sigma_P') / 2;
    Sigma_Q = (Sigma_Q + Sigma_Q') / 2;

    % Compute log of determinants for numerical stability
    log_det_Sigma_P = log(det(Sigma_P));
    log_det_Sigma_Q = log(det(Sigma_Q));

    % Compute inverse of covariance matrix Q
    inv_Sigma_Q = inv(Sigma_Q);

    % Calculate each term in the KL divergence formula
    term1 = log_det_Sigma_Q - log_det_Sigma_P;
    term2 = trace(inv_Sigma_Q * Sigma_P);
    term3 = (mu_Q - mu_P) * inv_Sigma_Q * (mu_Q - mu_P)';

    % Final KL divergence computation
    KL_div = 0.5 * (term1 - d + term2 + term3);
end
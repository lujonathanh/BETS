Function: g_n, alpha_k, beta_k, n, k_k, h {
alpha_k + beta_k * n**h / (k_k**h + n**h)
}

Function: q_n, alpha_m, beta_m, n, k_m, p {
alpha_s + beta_m/ (1 + (n/k_m)**p)
}

Function: r_nm, delta, m, lambda_k {
delta * m + lambda_k
}

Function: s_nm, lambda_m {
lambda_m
}

Gen_n:
    $pool > n
    g_n

Deg_n:
    n > $pool
    r_nm

Gen_m:
    $pool > m
    q_n

Deg_m:
    m > $pool
    s_nm


# Parameters
alpha_k = 0.03
alpha_m = 0.00015
beta_k = 0.015
beta_m = 0.005
k_k = 10
k_m = 34
lambda_k = 0.0001
lambda_s = 0.0001
h = 2
p = 2
delta = 0.00002


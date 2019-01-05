Gen_z:
    $pool > Z
    alpha_z + beta_z * Z**p_z / (C_z**p_z + Z**p_z)

Deg_z:
    Z > $pool
    (delta_z * 1.0/(1 + Z* 1.0/gamma_z + A * 1.0/gamma_a) + lambda_z ) * Z


Gen_a:
    $pool > A
    alpha_a + beta_a * 1/(1 + (Z/C_a)**p_a)

Deg_a:
    A > $pool
    (delta_a * 1.0/(1 + Z* 1.0/gamma_z + A * 1.0/gamma_a) + lambda_a ) * A


Gen_x:
    $pool > X
    alpha_x + beta_x * Z**p_x /(C_x ** p_x + Z ** p_x)

Deg_x:
    X > $pool
    lambda_x * X


Gen_y:
    $pool > Y
    alpha_y + beta_y * Z ** p_y/(C_y ** p_y + Z ** p_y)

Deg_y:
    Y > $pool
    lambda_y * Y

# Variable Species
A = 0
Z = 0
X = 0
Y = 0

# Parameters
alpha_a = 0
alpha_z = 0.9 / 3600
alpha_x = 0
alpha_y = 0
beta_z = 0.03
beta_a = 0.003
beta_x = 0.003
beta_y = 0.003
C_z = 20
C_a = 3.3
C_x = 20
C_y = 10
p_a = 5
p_x = 5
p_y = 5
p_z = 2
lambda_a = 0.0001
lambda_z = 0.0001
lambda_x = 0.0001
lambda_y = 0.0001
delta_z = 0.001
delta_a = 0.001
gamma_z = 100
gamma_a = 1

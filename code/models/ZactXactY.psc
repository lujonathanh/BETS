Gen_z:
    $pool > Z
    alpha_z

Deg_Z:
    Z > $pool
    lambda_z * Z


Gen_x:
    $pool > X
    alpha_x + beta_x * Z**p_x /(C_x ** p_x + Z ** p_x)

Deg_x:
    X > $pool
    lambda_x * X


Gen_y:
    $pool > Y
    alpha_y + beta_y * X ** p_y/(C_y ** p_y + X ** p_y)

Deg_y:
    Y > $pool
    lambda_y * Y

# Variable Species
Z = 0
X = 0
Y = 0

# Parameters
alpha_z = 0.01
alpha_x = 0.01
alpha_y = 0.01
beta_x = 0.1
beta_y = 0.1
C_x = 50
C_y = 50
p_x = 2
p_y = 2
lambda_z = 0.001
lambda_x = 0.001
lambda_y = 0.001

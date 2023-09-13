n_iterations <- 2100
lr <- list(theta = .2, Sigma = .02, iterations_decrease_learning_rate = c(100, 200, 400), number_iterations_average = c(20, 40, 80), reduce_gradient_factor = .5)

# short chains to check code:
#n_iterations <- 30 
#lr <- list(theta = .2, Sigma = .02, iterations_decrease_learning_rate = c(6, 12, 18), number_iterations_average = c(2, 4, 8), reduce_gradient_factor = .5)

print("estimation 1/8")
ph2.stand <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.stand, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:7, random = list(8:46)), 
                 n_iterations, 
                 theta_start = model.stand$theta, Sigma_start = 0, 
                 fixed_Sigma = FALSE, fixed_theta = FALSE,
                 learning_rates = lr, model.stand$dinvv, min_Sigma = 0.0001)

print("estimation 2/8")
ph2.fullt <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.fullt, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:8, random = list(9:47)), 
                 n_iterations, 
                 theta_start = model.fullt$theta, Sigma_start = 0, 
                 fixed_Sigma = FALSE, fixed_theta = FALSE,
                 learning_rates = lr, model.fullt$dinvv, min_Sigma = 0.0001)

print("estimation 3/8")
ph2.notrt <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.notrt, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:6, random = list(7:45)), 
                 n_iterations, 
                 theta_start = model.notrt$theta, Sigma_start = 0, 
                 fixed_Sigma = FALSE, fixed_theta = FALSE,
                 learning_rates = lr, model.notrt$dinvv, min_Sigma = 0.0001)

print("estimation 4/8")
ph2.nosts <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.nosts, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:4, random = list(5:43)), 
                 n_iterations, 
                 theta_start = model.nosts$theta, Sigma_start = 0, 
                 fixed_Sigma = FALSE, fixed_theta = FALSE,
                 learning_rates = lr, model.nosts$dinvv, min_Sigma = 0.0001)



# control chains: no random effects

print("estimation 5/8")
p2c.stand <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.stand, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:7, random = list(8:46)), 
                 n_iterations, 
                 theta_start = model.stand$theta, Sigma_start = 0, 
                 fixed_Sigma = TRUE, fixed_theta = FALSE,
                 learning_rates = lr, model.stand$dinvv, min_Sigma = 0.0001)

print("estimation 6/8")
p2c.fullt <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.fullt, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:8, random = list(9:47)), 
                 n_iterations, 
                 theta_start = model.fullt$theta, Sigma_start = 0, 
                 fixed_Sigma = TRUE, fixed_theta = FALSE,
                 learning_rates = lr, model.fullt$dinvv, min_Sigma = 0.0001)

print("estimation 7/8")
p2c.notrt <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.notrt, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:6, random = list(7:45)), 
                 n_iterations, 
                 theta_start = model.notrt$theta, Sigma_start = 0, 
                 fixed_Sigma = TRUE, fixed_theta = FALSE,
                 learning_rates = lr, model.notrt$dinvv, min_Sigma = 0.0001)

print("estimation 8/8")
p2c.nosts <- RM2(input_forward_simulation = list(n_samples = 1, effects = effects.ranef.nosts, dataset = data), 
                 which_parameters = list(lambda = 1, theta = 1:4, random = list(5:43)), 
                 n_iterations, 
                 theta_start = model.nosts$theta, Sigma_start = 0, 
                 fixed_Sigma = TRUE, fixed_theta = FALSE,
                 learning_rates = lr, model.nosts$dinvv, min_Sigma = 0.0001)




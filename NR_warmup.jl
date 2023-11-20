using LinearAlgebra
using Printf
using Plots

X0 = [1.0, 0.5]  # initial guess
ϵ = 0.005  # tolerance
X_ans = X0[:]  # we will update x values which is our answers
iterations = 0
mismatch = zeros(10,1)

function deltaX(x1, x2)
    F = [exp(x1 * x2) - 1.2; cos(x1 + x2) - 0.5]
    return F
end

function compute_jacobian(x1, x2)
    Jac = [x2 * exp(x1 * x2) x1 * exp(x1 * x2); -sin(x1 + x2) -sin(x1 + x2)]
    adj_jac = [-sin(x1 + x2) -x1 * exp(x1 * x2); sin(x1 + x2) x2 * exp(x1 * x2)]
    det_jac = det(Jac)
    invJacobian = adj_jac / det_jac
    return invJacobian
end

function stopping_criteria(x_prev, x_ans, f_xp,i)
    max_diff_x = maximum(abs.(x_ans - x_prev) ./ abs.(x_prev))
    max_abs_f_xp = maximum(abs.(f_xp))
    mismatch[i] = max_abs_f_xp
    return max(max_diff_x, max_abs_f_xp), mismatch
end
# %% ------
function convergence(X_ans, iterations)
    while true
        X_prev = X_ans  
        J = compute_jacobian(X_ans[1], X_ans[2])

        if det(J) < eps() #Machine Epsilon for Float64
            println("Jacobian is non-invertible.")
            break
        end

        F_xp = deltaX(X_ans[1], X_ans[2])
        X_ans = X_ans - J * F_xp

        iterations= iterations+1
        tolerance, mismatch = stopping_criteria(X_prev, X_ans, F_xp,iterations)
        if tolerance < ϵ
            @printf("Converged to the solution in %d iterations.",iterations)
            @printf("\nx1 is: %0.5f ",X_ans[1])
            @printf("\nx2 is: %0.5f",X_ans[2])
            
            break
        end

    end
    return mismatch, iterations
end
mismatch, iterations = convergence(X_ans, iterations)
plot(1:iterations, mismatch[1:iterations],xlabel = "Iteration", ylabel = "Mismatch", legend = false)

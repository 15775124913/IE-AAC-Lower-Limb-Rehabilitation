function state = execution_dynamics_step(action, state, params, varargin)
%EXECUTION_DYNAMICS_STEP Optional coupled, nonideal joint execution surrogate.
% This is a lumped joint-space dynamic validation model, NOT an identified
% robot/patient model and NOT a mathematical stability or passivity proof.
% Governing equation (gravity nominally compensated):
%   M*qdd + D*qd = tau_act + tau_hri_signed,
%   tau_act_dot = (sat(tau_PD)-tau_act)/tau_act_const.
% This corrected stress-test implementation uses (i) a causal position-
% consistent desired velocity for the PD servo (NOT the unsynchronized raw
% outer-loop velocity when its position correction saturates), and
% (ii) an implicit dissipative load/velocity update to avoid numerical
% sample-to-sample sign-flip oscillations at high velocity-opposing load.
% A non-diagonal SPD M imposes cross-joint inertia coupling.
% The surrogate signed load is not an identified human/robot torque law.
% There is no inertial feedforward: position + derived velocity PD only.
%
% state = execution_dynamics_step('init', [], params, q_ref, dq_ref)
% state = execution_dynamics_step('step', state, params, index, q_des, dq_des, ddq_des, tau_mag)
% state = execution_dynamics_step('finish', state, params, q_terminal, dq_terminal, ddq_terminal)

    switch lower(char(action))
        case 'init'
            q_ref = varargin{1}; dq_ref = varargin{2};
            N = size(q_ref,1);
            cfg = params.execution_dynamics;
            I = cfg.inertia_diag(:);
            M = diag(I);
            % Symmetric, constant, explicitly positive-definite effective inertia.
            M(1,2) = cfg.coupling_fraction*sqrt(I(1)*I(2)); M(2,1) = M(1,2);
            M(2,3) = cfg.coupling_fraction*sqrt(I(2)*I(3)); M(3,2) = M(2,3);
            assert(all(eig(M)>0), 'Execution inertia M must be positive definite');
            assert(cfg.actuator_time_constant > 0 && cfg.delay_samples >= 0 && ...
                cfg.delay_samples == floor(cfg.delay_samples), 'Invalid actuator/delay');
            assert(all(cfg.torque_limit>0) && all(cfg.kp>0) && all(cfg.kd>0));
            assert(cfg.velocity_sign_eps>0 && all(cfg.viscous_damping>=0));
            state = struct();
            state.M = M;
            state.q = zeros(N,3); state.dq = zeros(N,3);
            state.ddq = zeros(N,3); state.tau_act = zeros(1,3);
            state.tau_act_hist = zeros(N,3);
            state.tau_target_hist = zeros(N,3);
            state.tau_hri_hist = zeros(N,3);
            state.torque_saturation = false(N,3);
            state.des_q = repmat(q_ref(1,:),N,1);
            state.des_dq = repmat(dq_ref(1,:),N,1);
            state.des_ddq = zeros(N,3);
            % Audit controller-level derivatives versus derivatives of the
            % bounded position sent to this execution plant.
            state.raw_des_dq = zeros(N,3);
            state.raw_des_ddq = zeros(N,3);
            state.velocity_solver_iterations = zeros(N,1);
            state.velocity_solver_residual_Nm = zeros(N,1);
            state.q(1,:) = q_ref(1,:);
            state.dq(1,:) = dq_ref(1,:);
            state.cfg = cfg;
        case 'step'
            i = varargin{1};
            q_des = reshape(varargin{2},1,3);
            dq_des = reshape(varargin{3},1,3);
            ddq_des = reshape(varargin{4},1,3);
            tau_mag = reshape(varargin{5},1,3);
            dt = params.dt;
            cfg = state.cfg;
            state.des_q(i,:) = q_des;
            state.raw_des_dq(i,:) = dq_des;
            state.raw_des_ddq(i,:) = ddq_des;
            % Independently derive the inner-loop velocity target from
            % the delivered bounded position sequence. This remains a
            % useful causal servo-interface check even after the V3
            % outer-loop correction q/v/a states have been unified.
            if i > 1
                state.des_dq(i,:) = (state.des_q(i,:) - state.des_q(i-1,:)) / dt;
            else
                state.des_dq(i,:) = dq_des;
            end
            % PD-only inner loop: no acceleration feed-forward, avoiding
            % impulsive feed-forward from position-limit corners.
            state.des_ddq(i,:) = 0;
            k_delayed = max(1,i-cfg.delay_samples);
            q_delayed = state.des_q(k_delayed,:);
            dq_delayed = state.des_dq(k_delayed,:);
            cmd = cfg.kp(:)'.*(q_delayed-state.q(i,:)) + ...
                  cfg.kd(:)'.*(dq_delayed-state.dq(i,:));
            lim = cfg.torque_limit(:)';
            saturated = abs(cmd) > lim;
            target = max(-lim,min(lim,cmd));
            % Backward Euler actuator-lag update (stable for any positive dt/tau).
            a = dt/(cfg.actuator_time_constant+dt);
            state.tau_act = state.tau_act + a*(target-state.tau_act);
            % IMPLICIT signed-load discretization: solve for the new
            % velocity before evaluating velocity-opposing interaction.
            % Treat the magnitude of interaction as frozen over this one
            % sample; it is recomputed from executed state at the next one.
            [vnew, tau_signed, nIter, residual] = solve_implicit_velocity( ...
                state.M, state.dq(i,:)', state.tau_act(:), abs(tau_mag(:)), ...
                cfg.viscous_damping(:), cfg.velocity_sign_eps, dt);
            state.dq(i+1,:) = vnew';
            state.ddq(i,:) = ((vnew-state.dq(i,:)')/dt)';
            state.q(i+1,:) = state.q(i,:) + dt*state.dq(i+1,:);
            state.velocity_solver_iterations(i) = nIter;
            state.velocity_solver_residual_Nm(i) = residual;
            state.tau_act_hist(i,:) = state.tau_act;
            state.tau_target_hist(i,:) = target;
            state.tau_hri_hist(i,:) = tau_signed;
            state.torque_saturation(i,:) = saturated;
            if any(~isfinite(state.q(i+1,:))) || any(~isfinite(state.dq(i+1,:)))
                error('Nonfinite execution-plant state at sample %d', i);
            end
        case 'finish'
            N = size(state.q,1);
            % N-1 integration transitions update the executed state through
            % sample N.  The terminal reference is a LOGGED snapshot only:
            % it is not sent to the actuator and does not add a time step.
            assert(N>=2 && numel(varargin)==3, ...
                'Finish requires terminal command position, velocity and acceleration.');
            q_terminal  = reshape(varargin{1},1,3);
            dq_terminal = reshape(varargin{2},1,3);
            ddq_terminal = reshape(varargin{3},1,3);
            assert(all(isfinite([q_terminal,dq_terminal,ddq_terminal])), ...
                'Nonfinite terminal command.');
            state.des_q(N,:) = q_terminal;
            state.raw_des_dq(N,:) = dq_terminal;
            state.raw_des_ddq(N,:) = ddq_terminal;
            % As for samples 2..N-1, the inner-loop desired velocity is
            % defined by the causal position difference (not raw dq_cmd).
            state.des_dq(N,:) = (state.des_q(N,:)-state.des_q(N-1,:)) / params.dt;
            state.des_ddq(N,:) = 0;  % PD-only actuator, no accel feedforward
            % Last-sample acceleration, torque and solver fields below are
            % endpoint display copies, not an additional dynamics update.
            state.ddq(N,:) = state.ddq(N-1,:);
            state.tau_act_hist(N,:) = state.tau_act_hist(N-1,:);
            state.tau_target_hist(N,:) = state.tau_target_hist(N-1,:);
            state.tau_hri_hist(N,:) = state.tau_hri_hist(N-1,:);
            state.torque_saturation(N,:) = state.torque_saturation(N-1,:);
            state.velocity_solver_iterations(N) = state.velocity_solver_iterations(N-1);
            state.velocity_solver_residual_Nm(N) = state.velocity_solver_residual_Nm(N-1);
        otherwise
            error('Unknown execution model action: %s',char(action));
    end
end

function [v, tau_signed, it, residual] = solve_implicit_velocity(M, vold, tau_act, mag, D, epsv, dt)
% Solve strongly monotone implicit velocity balance with damped Newton.
%   M*(v-vold)/dt + D.*v + mag.*tanh(v/epsv) = tau_act.
% Jacobian: M/dt + diag(D + mag/epsv .* sech(v/epsv)^2) is SPD.
% High load can violate the explicit-Euler local stability condition;
% the implicit solution prevents an artificial ± acceleration cycle.
    v = vold;
    tol = 1e-10*(1 + max(abs(tau_act)) + max(mag));
    max_iter = 40;
    for it = 1:max_iter
        tv = tanh(v/epsv);
        F = M*(v-vold)/dt + D.*v + mag.*tv - tau_act;
        residual = max(abs(F));
        if residual <= tol
            break;
        end
        J = M/dt + diag(D + (mag/epsv).*(1-tv.^2));
        dv = J\F;
        % Monotonicity of the underlying nonlinear force balance does not
        % guarantee an undamped Newton step reduces the residual norm.
        % Backtracking prevents a wrong branch or large Newton overshoot.
        alpha = 1;
        accepted = false;
        for ls = 1:25
            vc = v-alpha*dv;
            Fc = M*(vc-vold)/dt + D.*vc + mag.*tanh(vc/epsv) - tau_act;
            if max(abs(Fc)) < residual
                v = vc;
                accepted = true;
                break;
            end
            alpha = alpha/2;
        end
        if ~accepted
            error('Implicit execution-velocity solver line search failed at iter %d',it);
        end
    end
    tv = tanh(v/epsv);
    F = M*(v-vold)/dt + D.*v + mag.*tv - tau_act;
    residual = max(abs(F));
    if residual > tol
        error('Implicit execution-velocity solver did not converge: residual %.3e Nm',residual);
    end
    tau_signed = -mag.*tv;
end

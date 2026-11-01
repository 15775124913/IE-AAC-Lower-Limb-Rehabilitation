function [q_next, v_next, a_next] = bounded_admittance_step(q, v, a_prev, a_des, dt, q_lim, v_lim, a_lim, j_lim)
%BOUNDED_ADMITTANCE_STEP Shared scalar discrete outer-loop motion supervisor.
% One sample transition; used identically by all four controllers.
%
% State convention (at sample k): q[k], v[k], a[k].  The selected a[k+1]
% is the actual finite difference of v[k+1], not a pre-projection request.
%   v[k+1] = v[k] + dt*a[k+1]
%   q[k+1] = q[k] + dt*v[k+1]     (semi-implicit Euler)
%   jerk[k+1] = (a[k+1]-a[k])/dt.
%
% This function never clips q or v *after* integrating, which otherwise
% violates kinematic identities and creates artificial jerk spikes.
% The hard one-step reachable interval is intersected with a predictive
% jerk-limited stopping envelope near +/-q_lim.  The latter uses a discrete
% maximal-braking rollout, not a stability/passivity proof.
%
% If constraints become mutually infeasible, raise an explicit error rather
% than quietly violate a published bound.
    if ~(isscalar(dt) && dt>0 && q_lim>0 && v_lim>0 && a_lim>0 && j_lim>0)
        error('MotionSupervisor:InvalidLimits','Expected positive scalar dt and motion limits.');
    end
    if any(~isfinite([q,v,a_prev,a_des,dt,q_lim,v_lim,a_lim,j_lim]))
        error('MotionSupervisor:Nonfinite','Nonfinite motion state, target, or limit.');
    end

    % Feasible acceleration interval enforcing acceleration, jerk, velocity,
    % and next-sample position bounds simultaneously.
    lo = max([-a_lim, a_prev-j_lim*dt, (-v_lim-v)/dt, (-q_lim-q-v*dt)/(dt^2)]);
    hi = min([ a_lim, a_prev+j_lim*dt, ( v_lim-v)/dt, ( q_lim-q-v*dt)/(dt^2)]);
    if lo > hi+1e-9
        error('MotionSupervisor:Unreachable', ...
            'Constraint conflict at q=%.7f, v=%.7f, a=%.7f (lo %.7f > hi %.7f).', ...
            q,v,a_prev,lo,hi);
    end
    if hi<lo
        hi=lo; % only roundoff-scale interval overlap
    end
    candidate = min(max(a_des,lo),hi);

    % Predict whether the next state can still be brought to rest before
    % either position limit with maximally admissible jerk/acceleration.
    % Monotonicity w.r.t. candidate acceleration allows bisection.
    if ~stoppable_before_limit(q,v,candidate,dt,q_lim,v_lim,a_lim,j_lim,1)
        if ~stoppable_before_limit(q,v,lo,dt,q_lim,v_lim,a_lim,j_lim,1)
            error('MotionSupervisor:UpperBoundaryUnreachable', ...
                'Cannot stop before upper position bound with the available jerk/acceleration.');
        end
        left=lo; right=candidate;
        for iter=1:30
            middle=0.5*(left+right);
            if stoppable_before_limit(q,v,middle,dt,q_lim,v_lim,a_lim,j_lim,1)
                left=middle;
            else
                right=middle;
            end
        end
        candidate=left;
    end
    if ~stoppable_before_limit(q,v,candidate,dt,q_lim,v_lim,a_lim,j_lim,-1)
        if ~stoppable_before_limit(q,v,hi,dt,q_lim,v_lim,a_lim,j_lim,-1)
            error('MotionSupervisor:LowerBoundaryUnreachable', ...
                'Cannot stop before lower position bound with the available jerk/acceleration.');
        end
        left=candidate; right=hi;
        for iter=1:30
            middle=0.5*(left+right);
            if stoppable_before_limit(q,v,middle,dt,q_lim,v_lim,a_lim,j_lim,-1)
                right=middle;
            else
                left=middle;
            end
        end
        candidate=right;
    end

    a_next = candidate;
    v_next = v + a_next*dt;
    q_next = q + v_next*dt;
    epsq=2e-9; epsv=2e-9; epsa=2e-9; epsj=2e-9;
    if abs(q_next)>q_lim+epsq || abs(v_next)>v_lim+epsv || ...
            abs(a_next)>a_lim+epsa || abs(a_next-a_prev)>j_lim*dt+epsj || ...
            ~stoppable_before_limit(q,v,a_next,dt,q_lim,v_lim,a_lim,j_lim,1) || ...
            ~stoppable_before_limit(q,v,a_next,dt,q_lim,v_lim,a_lim,j_lim,-1)
        error('MotionSupervisor:InternalValidation', ...
            'Discrete motion transition did not satisfy all constraints.');
    end
end

function feasible=stoppable_before_limit(q,v,a,dt,q_lim,v_lim,a_lim,j_lim,sgn)
% Evaluate the maximum future position if braking starts *next* sample.
% sgn=+1 for the upper limit; sgn=-1 for the lower limit.
    x=sgn*(q+(v+a*dt)*dt);
    speed=sgn*(v+a*dt);
    accel=sgn*a;
    feasible=(x<=q_lim-1e-6);
    if ~feasible, return; end
    % Bound is generous; reachable v and a must be within their limits.
    nmax=ceil((2*a_lim/j_lim + 2*v_lim/a_lim + 2*dt)/dt)+5;
    for k=1:nmax
        if speed<=0 && accel<=0
            return;
        end
        accel=max(-a_lim,accel-j_lim*dt);
        speed=speed+accel*dt;
        % Record the forward-most sample (negative step cannot increase it).
        x=x+max(speed,0)*dt;
        if x>q_lim-1e-6
            feasible=false;
            return;
        end
    end
    error('MotionSupervisor:BrakingRollout','Stopping preview exceeded its bounded horizon.');
end

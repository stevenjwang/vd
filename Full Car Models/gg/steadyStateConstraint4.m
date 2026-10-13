function [c,ceq] = steadyStateConstraint4(state,P,latAccelTarget)
%STEADYSTATECONSTRAINT4 Constraint4 from a cached vehicle-equation state.

redline = stateValue(state,'redline',NaN);
if isscalar(redline) && isfinite(redline)
    rpmLimit = state.engineRpm-redline;
else
    rpmLimit = 1e6;
end
betaLimit = abs(state.beta)-20;
loadLimits = -state.Fzvirtual(1:4);
if ~isfinite(rpmLimit), rpmLimit = 1e6; end
if ~isfinite(betaLimit), betaLimit = 1e6; end
loadLimits(~isfinite(loadLimits)) = 1e6;
c = [rpmLimit,betaLimit,loadLimits];

if isfield(state,'aeroConverged')
    aeroConverged = logical(state.aeroConverged);
    outsideMap = logical(stateValue(state,'aeroOutsideMap',false));
    coverageValid = logical(stateValue(state,'aeroCoverageValid',true));
    aeroResidual = stateValue(state,'aeroResidualIn',0);
    aeroTolerance = stateValue(state,'aeroResidualToleranceIn',1e-9);
    invalidAero = ~isfinite(double(aeroConverged)) || ~aeroConverged || ...
        outsideMap || ~coverageValid || ~isfinite(double(aeroResidual)) || ...
        aeroResidual > aeroTolerance;
    c(end+1:end+3) = double([~aeroConverged, ...
        outsideMap || ~coverageValid,invalidAero]);
end

ceq = [P(3)*P(5)-latAccelTarget,state.latAccel,state.yawAccel, ...
    state.wheelAccel(1:4)];
ceq(~isfinite(ceq)) = 1e6;
end

function value = stateValue(state,name,defaultValue)
if isfield(state,name) && ~isempty(state.(name))
    value = state.(name);
else
    value = defaultValue;
end
end

function [stencil_F2C,BC_F2C,HcU,HfU] = coeff_orders2to2_C2F
stencil_F2C = [1.0./4.0,1.0./2.0,1.0./4.0];
if nargout > 1
    BC_F2C = [1.0./2.0,1.0./2.0];
end
if nargout > 2
    HcU = 1.0./2.0;
end
if nargout > 3
    HfU = 1.0./2.0;
end

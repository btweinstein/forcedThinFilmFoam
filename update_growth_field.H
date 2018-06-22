
forAll(mesh.C(), cell_i){
    dimensionedScalar cur_h ("cur_h", dimLength, h[cell_i]);
    if(cur_h < epsilon) growth_factor[cell_i] = 0; // No growth if less than epsilon
    else if(cur_h < hmax) growth_factor[cell_i] = (a*(1 - cur_h/hmax)).value();
    else growth_factor[cell_i] = 0; // No growth or shrinking above hmax
};

growth = h*growth_factor; // Effectively a*h(1-h)
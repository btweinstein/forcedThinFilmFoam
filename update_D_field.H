
forAll(mesh.C(), cell_i){
    dimensionedScalar cur_h ("cur_h", dimLength, h[cell_i]);
    if(cur_h < hmax) D[cell_i] = (Dbasal*pow(1 - cur_h/hmax, blocking_exponent)).value();
    else D[cell_i] = 0; // No diffusing when too dense
};
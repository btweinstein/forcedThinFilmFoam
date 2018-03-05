forAll(mesh.C(), cell_i){
    dimensionedScalar cur_h ("cur_h", dimensionSet(0,1,0,0,0,0,0), h[cell_i]);

    dimensionedScalar zero_h ("zero_h", dimensionSet(0,1,0,0,0,0,0), SMALL);

    if(cur_h < zero_h) growth[cell_i] = 0; // Zero check
    else if(cur_h < d) growth[cell_i] = (a*cur_h).value();
    else if (cur_h < (hmax - d)) growth[cell_i] = (a*d).value();
    else if (cur_h < hmax) growth[cell_i] = (a*(hmax - cur_h)).value(); // Slow growth as you approach hmax
    else growth[cell_i] = 0; // No growth past hmax
};
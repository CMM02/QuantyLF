function get_slater_integrals(ion, oxy)

    -- Slater integral values

    if ion == 29 then
        if oxy == 1 then
            nd = 10
        elseif oxy == 2 then
            nd = 9
        elseif oxy == 3 then
            nd = 8
        elseif oxy == 4 then
            nd = 7
        end
        zeta_3d = 0.102;
        F2dd = 12.854;
        F4dd = 7.980 --- initial state parameters
        zeta_2p = 13.498;
        F2pd = 8.177;
        G1pd = 6.169;
        G3pd = 3.510 ---  final  state parameters
        Xzeta_3d = 0.124;
        XF2dd = 13.611;
        XF4dd = 8.457 ---  final  state parameters

    elseif ion == 28 then
        if oxy == 2 then
            nd = 8
        elseif oxy == 3 then
            nd = 7
        elseif oxy == 4 then
            nd = 6
        end
        zeta_3d = 0.083;
        F2dd = 12.233;
        F4dd = 7.597
        zeta_2p = 11.507;
        F2pd = 7.720;
        G1pd = 5.783;
        G3pd = 3.290
        Xzeta_3d = 0.102;
        XF2dd = 13.005;
        XF4dd = 8.084

    elseif ion == 27 then
        if oxy == 2 then
            nd = 7
        elseif oxy == 3 then
            nd = 6
        elseif oxy == 4 then
            nd = 5
        end
        zeta_3d = 0.066;
        F2dd = 11.604;
        F4dd = 7.209
        zeta_2p = 9.748;
        F2pd = 7.259;
        G1pd = 5.394;
        G3pd = 3.068
        Xzeta_3d = 0.083;
        XF2dd = 12.395;
        XF4dd = 7.707

    elseif ion == 26 then
        if oxy == 2 then
            nd = 6
            zeta_3d = 0.052;
            F2dd = 10.965;
            F4dd = 6.815
            zeta_2p = 8.200;
            F2pd = 6.792;
            G1pd = 5.0;
            G3pd = 2.843
            Xzeta_3d = 0.067;
            XF2dd = 11.778;
            XF4dd = 7.327
        elseif oxy == 3 then
            nd = 5
            zeta_3d = 0.059;
            F2dd = 12.043;
            F4dd = 7.535
            zeta_2p = 8.199;
            F2pd = 7.446;
            G1pd = 5.566;
            G3pd = 3.166
            Xzeta_3d = 0.074;
            XF2dd = 12.818;
            XF4dd = 8.023
        elseif oxy == 4 then
            nd = 4
        end

    elseif ion == 25 then
        if oxy == 2 then
            nd = 5
            zeta_3d = 0.040;
            F2dd = 10.315;
            F4dd = 6.413
            zeta_2p = 6.846;
            F2pd = 6.320;
            G1pd = 4.603;
            G3pd = 2.617
            Xzeta_3d = 0.053;
            XF2dd = 11.154;
            XF4dd = 6.942
        elseif oxy == 3 then
            nd = 4
            zeta_3d = 0.046;
            F2dd = 11.415;
            F4dd = 7.148;
            zeta_2p = 6.845;
            F2pd = 6.988;
            G1pd = 5.179;
            G3pd = 2.945;
            Xzeta_3d = 0.059;
            XF2dd = 12.210;
            XF4dd = 7.649
        elseif oxy == 4 then
            nd = 3
            zeta_3d = 0.052;
            F2dd = 12.416;
            F4dd = 7.820;
            zeta_2p = 6.845;
            F2pd = 7.658;
            G1pd = 5.776;
            G3pd = 3.288
            Xzeta_3d = 0.066;
            XF2dd = 13.177;
            XF4dd = 8.299;
        elseif oxy == 7 then
            nd = 2
        end

    elseif ion == 24 then
        if oxy == 2 then
            nd = 4
            zeta_3d = 0.030;
            F2dd = 9.469;
            F4dd = 6.002
            zeta_2p = 5.668;
            F2pd = 5.841;
            G1pd = 4.204;
            G3pd = 2.388
            Xzeta_3d = 0.041;
            XF2dd = 10.521;
            XF4dd = 6.522
        elseif oxy == 3 then
            nd = 3
            zeta_3d = 0.035;
            F2dd = 10.777;
            F4dd = 6.755
            zeta_2p = 5.667;
            F2pd = 6.526;
            G1pd = 4.788;
            G3pd = 2.722
            Xzeta_3d = 0.047;
            XF2dd = 11.596;
            XF4dd = 7.270
        elseif oxy == 4 then
            nd = 2
        end

    elseif ion == 23 then
        if oxy == 2 then
            nd = 3
        elseif oxy == 3 then
            nd = 2
        elseif oxy == 4 then
            nd = 1
        end
        zeta_3d = 0.022;
        F2dd = 8.961;
        F4dd = 5.576
        zeta_2p = 4.650;
        F2pd = 5.351;
        G1pd = 3.792;
        G3pd = 2.154
        Xzeta_3d = 0.031;
        XF2dd = 9.875;
        XF4dd = 6.152

    elseif ion == 22 then
        if oxy == 2 then
            nd = 2
        elseif oxy == 3 then
            nd = 1
        elseif oxy == 4 then
            nd = 0
        end
        zeta_3d = 0.016;
        F2dd = 8.243;
        F4dd = 5.132
        zeta_2p = 3.776;
        F2pd = 4.849;
        G1pd = 3.376;
        G3pd = 1.917
        Xzeta_3d = 0.023;
        XF2dd = 9.213;
        XF4dd = 5.744

    elseif ion == 21 then
        if oxy == 2 then
            nd = 1
        elseif oxy == 3 then
            nd = 2
        end
        zeta_3d = 0.010;
        F2dd = 0;
        F4dd = 0
        zeta_2p = 3.032;
        F2pd = 4.332;
        G1pd = 2.950;
        G3pd = 1.674
        Xzeta_3d = 0.017;
        XF2dd = 8.530;
        XF4dd = 5.321

    else
        print("Could not recognize the ion name...")
        os.exit()
    end
    return nd, zeta_3d, F2dd, F4dd, zeta_2p, F2pd, G1pd, G3pd, Xzeta_3d, XF2dd, XF4dd
end

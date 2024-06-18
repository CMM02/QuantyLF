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
            zeta_3d = 0.083;
            F2dd = 12.234; --- Thesis de Groot
            F4dd = 7.598 --- Thesis de Groot
            zeta_2p = 11.507; --- Thesis de Groot
            F2pd = 7.721; --- Thesis de Groot
            G1pd = 5.787; --- Thesis de Groot
            G3pd = 3.291 --- Thesis de Groot
            Xzeta_3d = 0.102; --- Thesis de Groot
            XF2dd = 13.005;
            XF4dd = 8.084
        elseif oxy == 3 then
            nd = 7
            zeta_3d = 0.083;
            F2dd = 13.277; --- Thesis de Groot
            F4dd = 8.295; --- Thesis de Groot
            zeta_2p = 11.506; --- Thesis de Groot
            F2pd = 8.350; --- Thesis de Groot
            G1pd = 6.332; --- Thesis de Groot
            G3pd = 3.603 --- Thesis de Groot
            Xzeta_3d = 0.112; --- Thesis de Groot
            XF2dd = 14.022; --- Thesis de Groot
            XF4dd = 8.764 --- Thesis de Groot
        elseif oxy == 4 then
            nd = 6
            print("No data available for this ion and valence configuration...")
            os.exit()
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
            zeta_3d = 0.066;
            F2dd = 11.605; --- Thesis de Groot
            F4dd = 7.209 --- Thesis de Groot
            zeta_2p = 9.746; --- Thesis de Groot
            F2pd = 7.260; --- Thesis de Groot
            G1pd = 5.397; --- Thesis de Groot
            G3pd = 3.069 --- Thesis de Groot
            Xzeta_3d = 0.092; --- Thesis de Groot
            XF2dd = 12.396; --- Thesis de Groot
            XF4dd = 7.708 --- Thesis de Groot
        elseif oxy == 3 then
            nd = 6
            zeta_3d = 0.066;
            F2dd = 12.663; --- Thesis de Groot
            F4dd = 7.917 --- Thesis de Groot
            zeta_2p = 9.748; --- Thesis de Groot
            F2pd = 7.900; --- Thesis de Groot
            G1pd = 5.961; --- Thesis de Groot
            G3pd = 3.386 --- Thesis de Groot
            Xzeta_3d = 0.082; --- Thesis de Groot
            XF2dd = 13.422; --- Thesis de Groot
            XF4dd = 78.395 --- Thesis de Groot
        elseif oxy == 4 then
            nd = 5
            print("No data available for this ion and valence configuration...")
            os.exit()
        end

    elseif ion == 26 then
        if oxy == 2 then
            nd = 6
            zeta_3d = 0.052;
            F2dd = 10.966; --- Thesis de Groot
            F4dd = 6.815 --- Thesis de Groot
            zeta_2p = 8.200; --- Thesis de Groot
            F2pd = 6.793; --- Thesis de Groot
            G1pd = 5.004; --- Thesis de Groot
            G3pd = 2.844 --- Thesis de Groot
            Xzeta_3d = 0.067; --- Thesis de Groot
            XF2dd = 11.779; --- Thesis de Groot
            XF4dd = 7.327 --- Thesis de Groot
        elseif oxy == 3 then
            nd = 5
            zeta_3d = 0.059;
            F2dd = 12.043; --- Thesis de Groot
            F4dd = 7.535 --- Thesis de Groot
            zeta_2p = 8.199; --- Thesis de Groot
            F2pd = 7.446; --- Thesis de Groot
            G1pd = 5.566; --- Thesis de Groot
            G3pd = 3.166 --- Thesis de Groot
            Xzeta_3d = 0.074; --- Thesis de Groot
            XF2dd = 12.818; --- Thesis de Groot
            XF4dd = 8.023 --- Thesis de Groot
        elseif oxy == 4 then
            nd = 4
            print("No data available for this ion and valence configuration...")
            os.exit()
        end

    elseif ion == 25 then
        if oxy == 2 then
            nd = 5
            zeta_3d = 0.040;
            F2dd = 10.316; --- Thesis de Groot
            F4dd = 6.414 --- Thesis de Groot
            zeta_2p = 6.846; --- Thesis de Groot
            F2pd = 6.321; --- Thesis de Groot
            G1pd = 4.606; --- Thesis de Groot
            G3pd = 2.618 --- Thesis de Groot
            Xzeta_3d = 0.053; --- Thesis de Groot
            XF2dd = 11.155; --- Thesis de Groot
            XF4dd = 6.943 --- Thesis de Groot
        elseif oxy == 3 then
            nd = 4
            zeta_3d = 0.046;
            F2dd = 11.415; --- Thesis de Groot
            F4dd = 7.148; --- Thesis de Groot
            zeta_2p = 6.845; --- Thesis de Groot
            F2pd = 6.988; --- Thesis de Groot
            G1pd = 5.179; --- Thesis de Groot
            G3pd = 2.945; --- Thesis de Groot
            Xzeta_3d = 0.059; --- Thesis de Groot
            XF2dd = 12.210; --- Thesis de Groot
            XF4dd = 7.649 --- Thesis de Groot
        elseif oxy == 4 then
            nd = 3
            zeta_3d = 0.052;
            F2dd = 12.416; --- Thesis de Groot
            F4dd = 7.820; --- Thesis de Groot
            zeta_2p = 6.845; --- Thesis de Groot
            F2pd = 7.658; --- Thesis de Groot
            G1pd = 5.776; --- Thesis de Groot
            G3pd = 3.288 --- Thesis de Groot
            Xzeta_3d = 0.066; --- Thesis de Groot
            XF2dd = 13.177; --- Thesis de Groot
            XF4dd = 8.299; --- Thesis de Groot
        elseif oxy == 7 then
            nd = 2
            print("No data available for this ion and valence configuration...")
            os.exit()
        end

    elseif ion == 24 then
        if oxy == 2 then
            nd = 4
            zeta_3d = 0.030;
            F2dd = 9.649; --- Thesis de Groot
            F4dd = 6.002 --- Thesis de Groot
            zeta_2p = 5.668; --- Thesis de Groot
            F2pd = 5.841; --- Thesis de Groot
            G1pd = 4.204; --- Thesis de Groot
            G3pd = 2.388 --- Thesis de Groot
            Xzeta_3d = 0.041; --- Thesis de Groot
            XF2dd = 10.522; --- Thesis de Groot
            XF4dd = 6.552 --- Thesis de Groot
        elseif oxy == 3 then
            nd = 3
            zeta_3d = 0.035;
            F2dd = 10.777; --- Thesis de Groot
            F4dd = 6.755 --- Thesis de Groot
            zeta_2p = 5.667; --- Thesis de Groot
            F2pd = 6.526; --- Thesis de Groot
            G1pd = 4.788; --- Thesis de Groot
            G3pd = 2.722 --- Thesis de Groot
            Xzeta_3d = 0.047; --- Thesis de Groot
            XF2dd = 11.596; --- Thesis de Groot
            XF4dd = 7.270 --- Thesis de Groot
        elseif oxy == 4 then
            nd = 2
            print("No data available for this ion and valence configuration...")
            os.exit()
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

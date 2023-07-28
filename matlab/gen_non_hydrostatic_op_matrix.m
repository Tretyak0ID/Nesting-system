function D = gen_vertical_SBP_matrix(N, format)

    DV1 = [-2 3 -1 0 0];
    DV2 = [-1 1 0 0 0];
    DV3 = [1/24 -9/8 9/8 -1/24 0];
    DV4 = [-1/71 6/71 -83/71 81/71 -3/71];
    DVend1 = [0 0 1 -3 2];
    DVend2 = [0 0 0 -1 1];
    DVend3 = [0 1/24 -9/8 9/8 -1/24];
    DVend4 = [3/71 -81/71 83/71 -6/71 1/71];

    DP1 = [-79/78 27/26 -1/26 1/78 0];
    DP2 = [2/21 -9/7 9/7 -2/21 0];
    DP3 = [1/75 0 -27/25 83/75 -1/25];
    DPend1 = [0 -1/78 1/26 -27/26 79/78];
    DPend2 = [0 2/21 -9/7 9/7 -2/21];
    DPend3 = [1/25 -83/75 27/25 0 -1/75];

    e = ones(N, 1);
    DP = full(spdiags([1/24*e -9/8*e 9/8*e -1/24*e], -2:1, N, N));
    DP(1,1:5) = DP1;
    DP(2,1:5) = DP2;
    DP(3,1:5) = DP3;
    DP(end,end-4:end) = DPend1;
    DP(end-1,end-4:end) = DPend2;
    DP(end-2,end-4:end) = DPend3;

    DV = full(spdiags([1/24*e -9/8*e 9/8*e -1/24*e], -2:1, N, N));
    DV(1,1:5) = DV1;
    DV(2,1:5) = DV2;
    DV(3,1:5) = DV3;
    DV(4,1:5) = DV4;
    DV(end,end-4:end) = DVend1;
    DV(end-1,end-4:end) = DVend2;
    DV(end-2,end-4:end) = DVend3;
    DV(end-3,end-4:end) = DVend4;

    AP = eye(N);
    AV = eye(N);

    AP(1,1) = 7/18;
    AP(2,2) = 9/8;
    AP(3,3) = 1;
    AP(4,4) = 71/72;
    AP(end,end) = 7/18;
    AP(end-1,end-1) = 9/8;
    AP(end-2,end-2) = 1;
    AP(end-3,end-3) = 71/72;

    AV(1,1) = 13/12;
    AV(2,2) = 7/8;
    AV(3,3) = 25/24;
    AV(end,end) = 13/12;
    AV(end-1,end-1) = 7/8;
    AV(end-2,end-2) = 25/24;
    
    D = eye(N) - AV * DV * AP * DP;

    if format == "full"
        D = D;
    else if format == "spars"
            D = sparse(D)
        end
    end
    
end
function [stable_idx] = fixStableOnsetByHand(stable_idx)
% correct onset of stable period

stable_idx(2,98) = NaN;

stable_idx(7,324) = 544;
stable_idx(7,370) = 727;

stable_idx(8,210) = 453;

stable_idx(10,360) = 575;

stable_idx(12,372) = 910;

stable_idx(13,84) = 850;
stable_idx(13,223) = 908;

stable_idx(18,24) = 882;

stable_idx(19,383) = 1061;

stable_idx(22,228) = 779;

stable_idx(23,104) = 788;

stable_idx(25,11) = 640;
stable_idx(25,66) = 787;

stable_idx(27,359) = 998;
stable_idx(27,24) = 831;

stable_idx(28,303) = 962;
stable_idx(28,383) = 763;
stable_idx(28,326) = 911;

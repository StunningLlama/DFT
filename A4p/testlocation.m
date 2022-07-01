setup([2 8 14; 8 8 8], 1, 1, [48; 48; 48], diag([16 16 16]), [2; 2; 2], true, false);
global gbl_Vdual;
global gbl_r;
%gbl_Vdual = gbl_r(:,1).^2;
visualize2(real(-gbl_Vdual));
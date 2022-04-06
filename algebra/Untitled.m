syms rdd rd r thetadd thetad theta phidd phid phi er etheta ephi erd ethetad ephid accel
erd = thetad*etheta + sin(theta)*phid*ephi;
ethetad = -thetad*er + cos(theta)*phid*ephi;
ephid = -phid*sin(theta)*er - phid*cos(theta)*etheta;
accel = rdd*er + rd*erd + (rd*thetad + r*thetadd)*etheta + r*thetad*ethetad + (rd*sin(theta)*phid + r*cos(theta)*thetad*phid + r*sin(theta)*phidd)*ephi + r*sin(theta)*phid*ephid;
finaleq = simplify(expand(accel))
latexcode = latex(finaleq);
latexcode = strrep(latexcode, '\mathrm{rdd}', '\ddot{r}');
latexcode = strrep(latexcode, '\mathrm{rd}', '\dot{r}');
latexcode = strrep(latexcode, '\mathrm{thetadd}', '\ddot{\theta}');
latexcode = strrep(latexcode, '\mathrm{thetad}', '\dot{\theta}');
latexcode = strrep(latexcode, '\mathrm{theta}', '\theta ');
latexcode = strrep(latexcode, '\mathrm{phidd}', '\ddot{\varphi}');
latexcode = strrep(latexcode, '\mathrm{phid}', '\dot{\varphi}');
latexcode = strrep(latexcode, '\mathrm{phi}', '\varphi ');
latexcode = strrep(latexcode, '\mathrm{er}', 'e_r ');
latexcode = strrep(latexcode, '\mathrm{etheta}', 'e_\theta ');
latexcode = strrep(latexcode, '\mathrm{ephi}', 'e_\phi ');
latexcode = strrep(latexcode, '\,', '');
latexcode
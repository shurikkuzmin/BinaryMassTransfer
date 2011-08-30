size(15cm,5cm);

pair A=(-7.5,0), B=(-7.5,2), C=(7.5,2), D=(7.5,0);

draw (A--B--C--D--cycle);
draw("$C_{inlet}$",(A+B)*0.5,W);
draw("$\partial_y C = 0$",(C+B)*0.5,N);
draw("$\partial_y C = 0$",(A+D)*0.5,S);

pair bubble_vel_top=B+(C-B)*0.65;
pair bubble_vel_bottom=A+(D-A)*0.65;

draw("$U_{\mathrm{bubble}}$",bubble_vel_top,S);
draw("$U_{\mathrm{bubble}}$",bubble_vel_bottom,N);
draw(bubble_vel_top+(-0.5,-0.1)--bubble_vel_top+(0.75,-0.1),linewidth(0.4mm),Arrow(1mm));
draw(bubble_vel_bottom+(-0.5,0.1)--bubble_vel_bottom+(0.75,0.1),linewidth(0.4mm),Arrow(1mm));


draw("$\frac{\partial C}{\partial x}=0$",(C+D)*0.5,E);

pair AxisCenter=(-8,-0.5);

draw(AxisCenter--AxisCenter+(0.75,0),linewidth(0.4mm),Arrow(1mm));
draw(AxisCenter--AxisCenter+(0,0.75),linewidth(0.4mm),Arrow(1mm));

draw("$x$",AxisCenter+(0.5,0),S);
draw("$y$",AxisCenter+(0,0.5),W);
//draw(circle(0,r));
real rad=0.8;
real bubble_length=4;
pair right_sphere=A+5.0/10.0*(D-A)+(0.0,1.0);
pair left_sphere=right_sphere-(bubble_length,0);
path p1=arc(right_sphere,rad,-90,90);
path p2=arc(left_sphere,rad,90,270);

//path P=arc(0,R,step,0);
draw(p1);
draw(p2);
draw(right_sphere+(0,rad)--left_sphere+(0,rad));
draw(right_sphere+(0,-rad)--left_sphere+(0,-rad));
draw("$C_s$",(right_sphere+left_sphere)*0.5+(0,-rad),N);

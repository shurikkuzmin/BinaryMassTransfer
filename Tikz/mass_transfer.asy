size(15cm,5cm);

pair A=(-7.5,0), B=(-7.5,2), C=(7.5,2), D=(7.5,0);

pair left=0.5*A+(2,0);
pair right=0.5*D-(2,0);
draw (A--D);
draw (B--C);
pair O=(0,0);
real len=2;

draw(left--left+(0,2));
draw(right--right+(0,2));
draw(left+(0,-0.5)--right+(0,-0.5),Arrows(HookHead));
draw(0.5*(left+right)+(0,-0.3),"$\Delta x$",O);

draw("$U$",left+(0,1.5),W);
draw(left+(-len,1)--left+(0,1),Arrow(1mm));
draw("$C(x)$",left+(0,0.5),W);
draw(right+(0,1)--right+(len,1),Arrow(1mm));
draw("$U$",right+(0,1.5),E);
draw("$C(x+\Delta x)$",right+(0,0.5),E);

pair AxisCenter=(-8,-0.5);
draw(AxisCenter--AxisCenter+(0.75,0),linewidth(0.4mm),Arrow(1mm));
draw(AxisCenter--AxisCenter+(0,0.75),linewidth(0.4mm),Arrow(1mm));
draw("$x$",AxisCenter+(0.5,0),S);
draw("$y$",AxisCenter+(0,0.5),W);


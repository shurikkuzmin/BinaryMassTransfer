size(5cm,5cm);

pair A=(-3,0), B=(0,0), C=(0,-3);

pair left=A+(1,0);
pair right=B+(0,-2);
draw (A--B);
draw (B--C);
//dot(left);
//dot(right);


void cross(pair A) {
  real len=0.1;
  pair left_bottom=A-(len,len);
  pair right_top=A+(len,len);
  pair left_top=A+(-len,len);
  pair right_bottom=A+(len,-len);
  draw(left_bottom--right_top);
  draw(left_top--right_bottom);
}

cross(B);
cross(left);
cross(right);
dot(B+(2,2),black+5bp);
draw(B--B+(1,1),linewidth(0.5mm),Arrow(1mm));
draw(rotate(45)*"$\bar{n}$",0.5*(B+(1,1)),N);
//pair O=(0,0);
//real len=2;

//draw(left--left+(0,2));
//draw(right--right+(0,2));
//draw(left+(0,-0.5)--right+(0,-0.5),Arrows(HookHead));
//draw(0.5*(left+right)+(0,-0.3),"$\Delta x$",O);

//draw("$U$",left+(0,1.5),W);
//draw(left+(-len,1)--left+(0,1),Arrow(1mm));
//draw("$C(x)$",left+(0,0.5),W);
//draw(right+(0,1)--right+(len,1),Arrow(1mm));
//draw("$U$",right+(0,1.5),E);
//draw("$C(x+\Delta x)$",right+(0,0.5),E);

//pair AxisCenter=(-8,-0.5);
//draw(AxisCenter--AxisCenter+(0.75,0),linewidth(0.4mm),Arrow(1mm));
//draw(AxisCenter--AxisCenter+(0,0.75),linewidth(0.4mm),Arrow(1mm));
//draw("$x$",AxisCenter+(0.5,0),S);
//draw("$y$",AxisCenter+(0,0.5),W);

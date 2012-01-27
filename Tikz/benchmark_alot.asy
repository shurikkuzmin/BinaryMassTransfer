size(15cm,7.5cm);

pair AxisCenter=(-8,-0.5);
draw(AxisCenter--AxisCenter+(0.75,0),linewidth(0.4mm),Arrow(1mm));
draw(AxisCenter--AxisCenter+(0,0.75),linewidth(0.4mm),Arrow(1mm));
draw("$x$",AxisCenter+(0.5,0),S);
draw("$y$",AxisCenter+(0,0.5),W);


void one_bubble(real shift)
{
		pair A=(-7.5+shift,0), B=(-7.5+shift,2), C=(7.5+shift,2), D=(7.5+shift,0);

		draw (A--B--C--D--cycle);

		pair bubble_vel_top=B+(C-B)*0.8;
		pair bubble_vel_bottom=A+(D-A)*0.8;

		draw("$U_{\mathrm{bubble}}$",bubble_vel_top,S);
		draw("$U_{\mathrm{bubble}}$",bubble_vel_bottom,N);
		draw(bubble_vel_top+(-0.5,-0.1)--bubble_vel_top+(0.75,-0.1),linewidth(0.4mm),Arrow(1mm));
		draw(bubble_vel_bottom+(-0.5,0.1)--bubble_vel_bottom+(0.75,0.1),linewidth(0.4mm),Arrow(1mm));


		real rad=0.8;
		real bubble_length=4;
		pair right_sphere=A+5.0/10.0*(D-A)+bubble_length/2.0+(0.0,1.0);
		pair left_sphere=right_sphere-(bubble_length,0);
		path p1=arc(right_sphere,rad,-90,90);
		path p2=arc(left_sphere,rad,90,270);

		draw(p1);
		draw(p2);
		draw(right_sphere+(0,rad)--left_sphere+(0,rad));
		draw(right_sphere+(0,-rad)--left_sphere+(0,-rad));
		draw("$C^*$",(right_sphere+left_sphere)*0.5+(0,-rad),N);

}
draw((7.5,0)--(7.5,2));
one_bubble(0);
one_bubble(15);

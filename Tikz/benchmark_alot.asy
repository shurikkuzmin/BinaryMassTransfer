size(15cm,7.5cm);

pair AxisCenter=(-8,-0.5);
draw(AxisCenter--AxisCenter+(0.75,0),linewidth(0.4mm),Arrow(1mm));
draw(AxisCenter--AxisCenter+(0,0.75),linewidth(0.4mm),Arrow(1mm));
draw("$x$",AxisCenter+(0.5,0),S);
draw("$y$",AxisCenter+(0,0.5),W);


void one_bubble(real shift, bool left,bool right)
{
		pair A=(-7.5+shift,0), B=(-7.5+shift,2), C=(7.5+shift,2), D=(7.5+shift,0);

		//draw (A--B--C--D--cycle);
        if (left) 
            draw(A--B);
       
        draw(B--C);
        draw(A--D);
        
        if (right)
            draw(C--D);
      
		

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

import roundedpath;
void make_conjunction(real shift)
{
	pair bubble_vel_top=(7.5+shift,0);
	pair bubble_vel_bottom=(7.5+shift,2);
	pair middle_vel_top=(7.95+shift,1);
	pair middle_vel_bottom=(7.45+shift,1.5);
    draw(bubble_vel_top..middle_vel_top..middle_vel_bottom..bubble_vel_bottom,linetype("3 2 3 2"));
	//draw(roundedpath(bubble_vel_top--middle_vel_top--middle_vel_bottom--bubble_vel_bottom,0.75),dashed);
}
make_conjunction(-0.25);
make_conjunction(0.25);

pair bubble_vel_top=(7.5,0);
pair bubble_vel_bottom=(7.5,2);

draw("$U_{\mathrm{bubble}}$",bubble_vel_top,S);
draw("$U_{\mathrm{bubble}}$",bubble_vel_bottom,N);
draw(bubble_vel_top+(-0.5,-0.1)--bubble_vel_top+(0.75,-0.1),linewidth(0.4mm),Arrow(1mm));
draw(bubble_vel_bottom+(-0.5,0.1)--bubble_vel_bottom+(0.75,0.1),linewidth(0.4mm),Arrow(1mm));
draw("$\frac{\partial C}{\partial x}=0$",(-7.5,1),W);
draw("$\frac{\partial C}{\partial x}=0$",(22.5,1),E);


one_bubble(0,true,false);
one_bubble(15,false,true);
